#include "H2GCROC_Common.hxx"

// ============================================================================
// LOGGING SETUP
// ============================================================================

void set_easylogger(){
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.setGlobally(el::ConfigurationType::Format, "%datetime{%H:%m:%s}[%levshort] (%fbase) %msg");
    defaultConf.set(el::Level::Info, el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;34m%levshort\033[0m] (%fbase) %msg");
    defaultConf.set(el::Level::Warning, el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;33m%levshort\033[0m] (%fbase) %msg");
    defaultConf.set(el::Level::Error, el::ConfigurationType::Format, 
        "%datetime{%H:%m:%s}[\033[1;31m%levshort\033[0m] (%fbase) %msg");
    el::Loggers::reconfigureLogger("default", defaultConf);
}

// ============================================================================
// ARGUMENT PARSING
// ============================================================================

ScriptOptions parse_arguments_single_root(int argc, char **argv, const std::string& version) {
    ScriptOptions opts;

    opts.script_version = version;
    opts.script_name = argv[0];
    opts.script_name = opts.script_name.substr(opts.script_name.find_last_of("/\\") + 1);
    opts.script_name = opts.script_name.substr(0, opts.script_name.find_last_of("."));

    START_EASYLOGGINGPP(argc, argv);
    set_easylogger();

    argparse::ArgumentParser program(opts.script_name, opts.script_version);

    program.add_argument("-f", "--file").help("Input .root file").required();
    program.add_argument("-o", "--output").help("Output .root file").required();
    program.add_argument("-e", "--events").help("Number of events to process").default_value(std::string("-1"));
    program.add_argument("-v", "--verbose").help("Verbose mode").default_value(false).implicit_value(true);
    program.add_argument("--focal").help("Use FoCal mapping").default_value(false).implicit_value(true);

    try {
        program.parse_args(argc, argv);
        opts.input_file = program.get<std::string>("--file");
        opts.output_file = program.get<std::string>("--output");
        auto n_events_str = program.get<std::string>("--events");
        opts.n_events = std::stoi(n_events_str);
        opts.verbose = program.get<bool>("--verbose");
        opts.focal = program.get<bool>("--focal");
    } catch (const std::runtime_error& err) {
        LOG(ERROR) << err.what();
        LOG(INFO) << program;
        exit(1);
    }

    if (access(opts.input_file.c_str(), F_OK) == -1) {
        LOG(ERROR) << "Input file " << opts.input_file << " does not exist!";
        exit(1);
    }

    if (opts.input_file.substr(opts.input_file.find_last_of(".") + 1) != "root") {
        LOG(ERROR) << "Input file " << opts.input_file << " should end with .root!";
        exit(1);
    }

    if (opts.output_file.substr(opts.output_file.find_last_of(".") + 1) != "root") {
        LOG(ERROR) << "Output file " << opts.output_file << " should end with .root!";
        exit(1);
    }

    opts.output_folder = opts.output_file.substr(0, opts.output_file.find_last_of("/\\"));
    if (opts.output_folder.empty()) {
        opts.output_folder = "./dump/" + opts.script_name;
    }

    if (access(opts.output_folder.c_str(), F_OK) == -1) {
        LOG(INFO) << "Creating output folder " << opts.output_folder;
        if (mkdir(opts.output_folder.c_str(), 0777) == -1) {
            LOG(ERROR) << "Failed to create output folder " << opts.output_folder;
            exit(1);
        }
    }

    if (access(opts.output_file.c_str(), F_OK) != -1) {
        LOG(WARNING) << "Output file " << opts.output_file << " already exists!";
    }

    LOG(INFO) << "Script name: " << opts.script_name;
    LOG(INFO) << "Input file: " << opts.input_file;
    LOG(INFO) << "Output file: " << opts.output_file << " in " << opts.output_folder;
    LOG(INFO) << "Number of events: " << opts.n_events;

    return opts;
}

ScriptOptions parse_arguments_single_json(int argc, char **argv, const std::string& version) {
    ScriptOptions opts;

    opts.script_version = version;
    opts.script_name = argv[0];
    opts.script_name = opts.script_name.substr(opts.script_name.find_last_of("/\\") + 1);
    opts.script_name = opts.script_name.substr(0, opts.script_name.find_last_of("."));

    START_EASYLOGGINGPP(argc, argv);
    set_easylogger();

    argparse::ArgumentParser program(opts.script_name, opts.script_version);

    program.add_argument("-f", "--file").help("Input .json file").required();
    program.add_argument("-o", "--output").help("Output .root file").required();
    program.add_argument("-e", "--events").help("Number of events to process").default_value(std::string("-1"));
    program.add_argument("--focal").help("Use FoCal mapping").default_value(false).implicit_value(true);
    program.add_argument("-p", "--pedestal").help("Pedestal .root file").default_value(std::string(""));
    program.add_argument("-c", "--csv").help("Input .csv file").default_value(std::string(""));
    program.add_argument("-v", "--verbose").help("Verbose mode").default_value(false).implicit_value(true);

    try {
        program.parse_args(argc, argv);
        opts.input_file = program.get<std::string>("--file");
        opts.output_file = program.get<std::string>("--output");
        auto n_events_str = program.get<std::string>("--events");
        opts.n_events = std::stoi(n_events_str);
        opts.focal = program.get<bool>("--focal");
        opts.csv_file = program.get<std::string>("--csv");
        opts.verbose = program.get<bool>("--verbose");
        opts.pedestal_file = program.get<std::string>("--pedestal");
    } catch (const std::runtime_error& err) {
        LOG(ERROR) << err.what();
        LOG(INFO) << program;
        exit(1);
    }

    if (access(opts.input_file.c_str(), F_OK) == -1) {
        LOG(ERROR) << "Input file " << opts.input_file << " does not exist!";
        exit(1);
    }

    if (opts.input_file.substr(opts.input_file.find_last_of(".") + 1) != "json") {
        LOG(ERROR) << "Input file " << opts.input_file << " should end with .json!";
        exit(1);
    }

    if (opts.output_file.substr(opts.output_file.find_last_of(".") + 1) != "root") {
        LOG(ERROR) << "Output file " << opts.output_file << " should end with .root!";
        exit(1);
    }

    opts.output_folder = opts.output_file.substr(0, opts.output_file.find_last_of("/\\"));
    if (opts.output_folder.empty()) {
        opts.output_folder = "./dump/" + opts.script_name;
    }

    if (access(opts.output_folder.c_str(), F_OK) == -1) {
        LOG(INFO) << "Creating output folder " << opts.output_folder;
        if (mkdir(opts.output_folder.c_str(), 0777) == -1) {
            LOG(ERROR) << "Failed to create output folder " << opts.output_folder;
            exit(1);
        }
    }

    if (access(opts.output_file.c_str(), F_OK) != -1) {
        LOG(WARNING) << "Output file " << opts.output_file << " already exists!";
    }

    LOG(INFO) << "Script name: " << opts.script_name;
    LOG(INFO) << "Input file: " << opts.input_file;
    LOG(INFO) << "Output file: " << opts.output_file << " in " << opts.output_folder;
    LOG(INFO) << "Number of events: " << opts.n_events;

    return opts;
}

ScriptOptions parse_arguments_single_root_single_csv(int argc, char **argv, const std::string& version){
    ScriptOptions opts;

    opts.script_version = version;
    opts.script_name = argv[0];
    opts.script_name = opts.script_name.substr(opts.script_name.find_last_of("/\\") + 1);
    opts.script_name = opts.script_name.substr(0, opts.script_name.find_last_of("."));

    START_EASYLOGGINGPP(argc, argv);
    set_easylogger();

    argparse::ArgumentParser program(opts.script_name, opts.script_version);

    program.add_argument("-f", "--file").help("Input .root file").required();
    program.add_argument("-o", "--output").help("Output .root file").required();
    program.add_argument("-c", "--csv").help("Input .csv file").required();
    program.add_argument("--focal").help("Use FoCal mapping").default_value(false).implicit_value(true);
    program.add_argument("-e", "--events").help("Number of events to process").default_value(std::string("-1"));
    program.add_argument("-v", "--verbose").help("Verbose mode").default_value(false).implicit_value(true);
    program.add_argument("--timewalk").help("Enable timewalk correction").default_value(false).implicit_value(true);
    program.add_argument("-t", "--timewalk_file").help("Timewalk correction .json file").default_value(std::string(""));

    try {
        program.parse_args(argc, argv);
        opts.input_file = program.get<std::string>("--file");
        opts.output_file = program.get<std::string>("--output");
        auto n_events_str = program.get<std::string>("--events");
        opts.n_events = std::stoi(n_events_str);
        opts.csv_file = program.get<std::string>("--csv");
        opts.verbose = program.get<bool>("--verbose");
        opts.focal = program.get<bool>("--focal");
        opts.timewalk = program.get<bool>("--timewalk");
        opts.timewalk_file = program.get<std::string>("--timewalk_file");
    } catch (const std::runtime_error& err) {
        LOG(ERROR) << err.what();
        LOG(INFO) << program;
        exit(1);
    }

    if (access(opts.input_file.c_str(), F_OK) == -1) {
        LOG(ERROR) << "Input file " << opts.input_file << " does not exist!";
        exit(1);
    }

    if (opts.input_file.substr(opts.input_file.find_last_of(".") + 1) != "root") {
        LOG(ERROR) << "Input file " << opts.input_file << " should end with .root!";
        exit(1);
    }

    if (opts.output_file.substr(opts.output_file.find_last_of(".") + 1) != "root") {
        LOG(ERROR) << "Output file " << opts.output_file << " should end with .root!";
        exit(1);
    }

    if (opts.csv_file.substr(opts.csv_file.find_last_of(".") + 1) != "csv") {
        LOG(ERROR) << "CSV file " << opts.csv_file << " should end with .csv!";
        exit(1);
    }

    opts.output_folder = opts.output_file.substr(0, opts.output_file.find_last_of("/\\"));
    if (opts.output_folder.empty()) {
        opts.output_folder = "./dump/" + opts.script_name;
    }

    if (access(opts.output_folder.c_str(), F_OK) == -1) {
        LOG(INFO) << "Creating output folder " << opts.output_folder;
        if (mkdir(opts.output_folder.c_str(), 0777) == -1) {
            LOG(ERROR) << "Failed to create output folder " << opts.output_folder;
            exit(1);
        }
    }

    if (access(opts.output_file.c_str(), F_OK) != -1) {
        LOG(WARNING) << "Output file " << opts.output_file << " already exists!";
    }

    LOG(INFO) << "Script name: " << opts.script_name;
    LOG(INFO) << "Input file: " << opts.input_file;
    LOG(INFO) << "Output file: " << opts.output_file << " in " << opts.output_folder;
    LOG(INFO) << "Number of events: " << opts.n_events;

    return opts;
}

// ============================================================================
// CHANNEL UTILITIES
// ============================================================================

int get_valid_fpga_channel(int fpga_channel){
    if (fpga_channel < 0 || fpga_channel >= FPGA_CHANNEL_NUMBER) {
        LOG(ERROR) << "Invalid FPGA channel " << fpga_channel;
        return -1;
    }
    std::vector<int> CM_channels = {0, 38, 76, 114};
    for (auto CM_channel : CM_channels) {
        if (fpga_channel == CM_channel) {
            return -1;
        }
    }
    std::vector<int> Calib_channels = {19, 57, 95, 133};
    for (auto Calib_channel : Calib_channels) {
        if (fpga_channel == Calib_channel) {
            return -1;
        }
    }
    if (fpga_channel < 19) {
        return fpga_channel - 1;
    }
    if (fpga_channel < 38) {
        return fpga_channel - 2;
    }
    if (fpga_channel < 57) {
        return fpga_channel - 3;
    }
    if (fpga_channel < 76) {
        return fpga_channel - 4;
    }
    if (fpga_channel < 95) {
        return fpga_channel - 5;
    }
    if (fpga_channel < 114) {
        return fpga_channel - 6;
    }
    if (fpga_channel < 133) {
        return fpga_channel - 7;
    }
    if (fpga_channel < 152) {
        return fpga_channel - 8;
    }
    return -1;
}

int get_total_fpga_channel(int fpga_channel){
    if (fpga_channel < 0 || fpga_channel >= FPGA_CHANNEL_NUMBER_VALID) {
        LOG(ERROR) << "Invalid FPGA channel " << fpga_channel;
        return -1;
    }
    if (fpga_channel <= 17) {
        return fpga_channel + 1;
    }
    if (fpga_channel <= 35) {
        return fpga_channel + 2;
    }
    if (fpga_channel <= 53) {
        return fpga_channel + 3;
    }
    if (fpga_channel <= 71) {
        return fpga_channel + 4;
    }
    if (fpga_channel <= 89) {
        return fpga_channel + 5;
    }
    if (fpga_channel <= 107) {
        return fpga_channel + 6;
    }
    if (fpga_channel <= 125) {
        return fpga_channel + 7;
    }
    if (fpga_channel <= 143) {
        return fpga_channel + 8;
    }
    return -1;
}

// ============================================================================
// GLOBAL CHANNEL PAINTER IMPLEMENTATION
// ============================================================================

GlobalChannelPainter::GlobalChannelPainter(const std::string& mapping_file) {
    std::ifstream ifs(mapping_file);
    if (!ifs.is_open()) {
        LOG(ERROR) << "Failed to open mapping file " << mapping_file;
        exit(1);
    }
    if (this->mapping_json.empty()) {
        ifs >> this->mapping_json;
    } else {
        LOG(WARNING) << "Mapping file " << mapping_file << " overwritten!";
        mapping_json.clear();
        ifs >> this->mapping_json;
    }
    ifs.close();

    this->module_fpga_list = this->mapping_json["module_fpga"].get<std::vector<int>>();
    this->module_asic_list = this->mapping_json["module_asic"].get<std::vector<int>>();
    this->module_connector_list = this->mapping_json["module_connector"].get<std::vector<int>>();
    for (int _connector_index = 0; _connector_index < 4; _connector_index++) {
        this->connector_list_list.push_back(this->mapping_json["connector_" + std::to_string(_connector_index)].get<std::vector<int>>());
    }

    // check if the module fpga is based on 208
    for (int _module_index = 0; _module_index < 25; _module_index++) {
        if (this->module_fpga_list[_module_index] >= 208) {
            this->module_fpga_list[_module_index] -= 208;
        }
    }
    // check if the vector length is the same for all connectors
    if (this->connector_list_list.size() != 4) {
        LOG(ERROR) << "Connector vector length mismatch!";
        exit(1);
    }
    for (int _connector_index = 0; _connector_index < 4; _connector_index++) {
        if (this->connector_list_list[_connector_index].size() != 16) {
            LOG(ERROR) << "Connector vector length mismatch!";
            exit(1);
        }
    }
    // check if the vector length is the same for all asics
    if ((this->module_asic_list.size() != this->module_fpga_list.size()) || (this->module_asic_list.size() != this->module_connector_list.size())) {
        LOG(ERROR) << "ASIC vector length mismatch!";
        exit(1);
    }
    this->is_EEEMCal_mapping = true;
    this->is_FoCal_mapping = false;
}

GlobalChannelPainter::GlobalChannelPainter(const std::string& mapping_file, const std::string& channel_mapping_file){
    std::ifstream ifs(mapping_file);
    if (!ifs.is_open()) {
        LOG(ERROR) << "Failed to open mapping file " << mapping_file;
        exit(1);
    }
    if (this->focal_mapping_json.empty()) {
        ifs >> this->focal_mapping_json;
    } else {
        LOG(WARNING) << "Mapping file " << mapping_file << " overwritten!";
        focal_mapping_json.clear();
        ifs >> this->focal_mapping_json;
    }
    ifs.close();

    std::ifstream channel_ifs(channel_mapping_file);
    if (!channel_ifs.is_open()) {
        LOG(ERROR) << "Failed to open channel mapping file " << channel_mapping_file;
        exit(1);
    }
    if (this->focal_channel_mapping_json.empty()) {
        channel_ifs >> this->focal_channel_mapping_json;
    } else {
        LOG(WARNING) << "Channel mapping file " << channel_mapping_file << " overwritten!";
        focal_channel_mapping_json.clear();
        channel_ifs >> this->focal_channel_mapping_json;
    }
    channel_ifs.close();

    for(int _module_index = 0; _module_index < 9; _module_index++){
        this->focal_module_board_list.push_back(this->focal_mapping_json["module" + std::to_string(_module_index + 1) + "_board"].get<std::vector<int>>());
        this->focal_module_channel_list.push_back(this->focal_mapping_json["module" + std::to_string(_module_index + 1) + "_chn"].get<std::vector<int>>());
    }

    for (auto& [key, value] : this->focal_channel_mapping_json["CAEN_Ports"].items()) {
        int from = std::stoi(key);
        std::string to_str = value;
        int to = std::stoi(to_str);
        this->focal_channel_map[from] = to;
    }

    for (auto& [key, value] : this->focal_channel_mapping_json["CAEN_Boards"].items()) {
        int from = std::stoi(key);
        std::string to_str = value;
        int to = std::stoi(to_str);
        this->focal_fpga_map[to] = from;
    }

    this->is_EEEMCal_mapping = false;
    this->is_FoCal_mapping = true;
}

GlobalChannelPainter::~GlobalChannelPainter() {
    this->clear_canvas();
    if (this->painter_canvas != nullptr) {
        delete this->painter_canvas;
    }
}

void GlobalChannelPainter::clear_canvas() {
    if (!this->sub_canvas_list.empty()) {
        for (auto& sub_canvas : this->sub_canvas_list) {
            if (sub_canvas != nullptr) {
                delete sub_canvas;
                sub_canvas = nullptr;
            }
        }
        this->sub_canvas_list.clear();
    }
}

void GlobalChannelPainter::draw_global_channel_hists2D(std::vector <TH2D*> hists, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title){
    this->clear_canvas();
    this->painter_canvas = new TCanvas(hist_name.c_str(), hist_title.c_str(), 800, 600);

    if (is_EEEMCal_mapping) {
        LOG(INFO) << "Drawing global channel hists for EEEMCal";
        this->painter_canvas->Divide(5, 5, 0.0, 0.0);
        for (int _sub_canvas_index = 1; _sub_canvas_index <= 25; _sub_canvas_index++) {
            auto _pad = (TPad*) this->painter_canvas->GetPad(_sub_canvas_index);
            _pad->SetLeftMargin(0);
            _pad->SetRightMargin(0);
            _pad->SetTopMargin(0);
            _pad->SetBottomMargin(0);

            auto _pad_module_fpga = this->module_fpga_list[_sub_canvas_index - 1];
            auto _pad_module_asic = this->module_asic_list[_sub_canvas_index - 1];
            auto _pad_module_connector = this->module_connector_list[_sub_canvas_index - 1];
            auto _pad_connector = this->connector_list_list[_pad_module_connector-1];

            auto _pad_canvas = new TCanvas(("sub_canvas_" + std::to_string(_sub_canvas_index)).c_str(), ("Sub Canvas " + std::to_string(_sub_canvas_index)).c_str(), 800, 600);
            this->sub_canvas_list.push_back(_pad_canvas);
            _pad_canvas->Divide(4, 4, 0.0, 0.0);
            for (int _sub_sub_canvas_index = 1; _sub_sub_canvas_index <= 16; _sub_sub_canvas_index++) {
                auto _sub_pad = (TPad*) _pad_canvas->GetPad(_sub_sub_canvas_index);
                _sub_pad->SetLeftMargin(0);
                _sub_pad->SetRightMargin(0);
                _sub_pad->SetTopMargin(0);
                _sub_pad->SetBottomMargin(0);

                auto _unified_channel = _pad_connector[_sub_sub_canvas_index - 1] + _pad_module_fpga * FPGA_CHANNEL_NUMBER_VALID + _pad_module_asic * int(FPGA_CHANNEL_NUMBER_VALID / 2);
                auto _channel_index = hists_channel_map[_unified_channel];
                auto _hist = hists[_channel_index];

                if (_hist->GetEntries() == 0) {
                    _hist->SetBinContent(1, 0.0);
                }
                _sub_pad->cd();
                _hist->Draw("colz");
            }
            _pad->cd();
            _pad_canvas->DrawClonePad();
        }
    } else {
        if (is_FoCal_mapping) {
            LOG(INFO) << "Drawing FoCal mapping";
            this->painter_canvas->Divide(3, 3, 0.0, 0.0);
            for (int _sub_canvas_index = 1; _sub_canvas_index <= 9; _sub_canvas_index++) {
                auto _pad = (TPad*) this->painter_canvas->GetPad(_sub_canvas_index);
                _pad->SetLeftMargin(0);
                _pad->SetRightMargin(0);
                _pad->SetTopMargin(0);
                _pad->SetBottomMargin(0);

                auto _pad_module_board = this->focal_module_board_list[_sub_canvas_index - 1];
                auto _pad_module_channel = this->focal_module_channel_list[_sub_canvas_index - 1];

                auto _pad_canvas = new TCanvas(("sub_canvas_" + std::to_string(_sub_canvas_index)).c_str(), ("Sub Canvas " + std::to_string(_sub_canvas_index)).c_str(), 800, 600);
                this->sub_canvas_list.push_back(_pad_canvas);
                int _pad_divide_x_y;
                if (_sub_canvas_index == 5) {
                    _pad_divide_x_y = 7;
                } else {
                    _pad_divide_x_y = 5;
                }
                _pad_canvas->Divide(_pad_divide_x_y, _pad_divide_x_y, 0.0, 0.0);
                for (int _sub_sub_canvas_index = 1; _sub_sub_canvas_index <= _pad_divide_x_y*_pad_divide_x_y; _sub_sub_canvas_index++) {
                    auto _sub_pad = (TPad*) _pad_canvas->GetPad(_sub_sub_canvas_index);
                    _sub_pad->SetLeftMargin(0);
                    _sub_pad->SetRightMargin(0);
                    _sub_pad->SetTopMargin(0);
                    _sub_pad->SetBottomMargin(0);

                    auto _board_number = _pad_module_board[_sub_sub_canvas_index - 1];
                    auto _channel_number = _pad_module_channel[_sub_sub_canvas_index - 1];
                    if (this->focal_channel_map.find(_channel_number) == this->focal_channel_map.end()) {
                        LOG(ERROR) << "Channel number " << _channel_number << " not found in the channel map!";
                        exit(1);
                    }
                    auto _h2gcroc_channel_number = this->focal_channel_map[_channel_number];
                    if (this->focal_fpga_map.find(_board_number) == this->focal_fpga_map.end()) {
                        LOG(ERROR) << "Board number " << _board_number << " not found in the fpga map!";
                        exit(1);
                    }
                    auto _h2gcroc_board_number = this->focal_fpga_map[_board_number];

                    auto _unified_channel = _h2gcroc_board_number * (FPGA_CHANNEL_NUMBER_VALID / 2)  + _h2gcroc_channel_number;
                    auto _channel_index = hists_channel_map[_unified_channel];
                    auto _hist = hists[_channel_index];

                    if (_hist->GetEntries() == 0) {
                        _hist->SetBinContent(1, 0.0);
                    }
                    _sub_pad->cd();
                    _hist->Draw("colz");
                }
                _pad->cd();
                _pad_canvas->DrawClonePad();
            }
        }
    }
    
}

void GlobalChannelPainter::draw_global_channel_hists1D(std::vector <TH1D*> hists, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title) {
    this->clear_canvas();
    this->painter_canvas = new TCanvas(hist_name.c_str(), hist_title.c_str(), 800, 600);
    if (is_EEEMCal_mapping) {
        this->painter_canvas->Divide(5, 5, 0.0, 0.0);
        for (int _sub_canvas_index = 1; _sub_canvas_index <= 25; _sub_canvas_index++) {
            auto _pad = (TPad*) this->painter_canvas->GetPad(_sub_canvas_index);
            _pad->SetLeftMargin(0);
            _pad->SetRightMargin(0);
            _pad->SetTopMargin(0);
            _pad->SetBottomMargin(0);

            auto _pad_module_fpga = this->module_fpga_list[_sub_canvas_index - 1];
            auto _pad_module_asic = this->module_asic_list[_sub_canvas_index - 1];
            auto _pad_module_connector = this->module_connector_list[_sub_canvas_index - 1];
            auto _pad_connector = this->connector_list_list[_pad_module_connector-1];

            auto _pad_canvas = new TCanvas(("sub_canvas_" + std::to_string(_sub_canvas_index)).c_str(), ("Sub Canvas " + std::to_string(_sub_canvas_index)).c_str(), 800, 600);
            this->sub_canvas_list.push_back(_pad_canvas);
            _pad_canvas->Divide(4, 4, 0.0, 0.0);
            for (int _sub_sub_canvas_index = 1; _sub_sub_canvas_index <= 16; _sub_sub_canvas_index++) {
                auto _sub_pad = (TPad*) _pad_canvas->GetPad(_sub_sub_canvas_index);
                _sub_pad->SetLeftMargin(0);
                _sub_pad->SetRightMargin(0);
                _sub_pad->SetTopMargin(0);
                _sub_pad->SetBottomMargin(0);

                auto _unified_channel = _pad_connector[_sub_sub_canvas_index - 1] + _pad_module_fpga * FPGA_CHANNEL_NUMBER_VALID + _pad_module_asic * int(FPGA_CHANNEL_NUMBER_VALID / 2);
                auto _channel_index = hists_channel_map[_unified_channel];
                auto _hist = hists[_channel_index];

                if (_hist->GetEntries() == 0) {
                    _hist->SetBinContent(1, 0.0);
                }
                _sub_pad->cd();
                _hist->Draw("hist");
            }
            _pad->cd();
            _pad_canvas->DrawClonePad();
        }
    } else {
        if (is_FoCal_mapping) {
            this->painter_canvas->Divide(3, 3, 0.0, 0.0);
            for (int _sub_canvas_index = 1; _sub_canvas_index <= 9; _sub_canvas_index++) {
                auto _pad = (TPad*) this->painter_canvas->GetPad(_sub_canvas_index);
                _pad->SetLeftMargin(0);
                _pad->SetRightMargin(0);
                _pad->SetTopMargin(0);
                _pad->SetBottomMargin(0);

                auto _pad_module_board = this->focal_module_board_list[_sub_canvas_index - 1];
                auto _pad_module_channel = this->focal_module_channel_list[_sub_canvas_index - 1];

                auto _pad_canvas = new TCanvas(("sub_canvas_" + std::to_string(_sub_canvas_index)).c_str(), ("Sub Canvas " + std::to_string(_sub_canvas_index)).c_str(), 800, 600);
                this->sub_canvas_list.push_back(_pad_canvas);
                int _pad_divide_x_y;
                if (_sub_canvas_index == 5) {
                    _pad_divide_x_y = 7;
                } else {
                    _pad_divide_x_y = 5;
                }
                _pad_canvas->Divide(_pad_divide_x_y, _pad_divide_x_y, 0.0, 0.0);
                for (int _sub_sub_canvas_index = 1; _sub_sub_canvas_index <= _pad_divide_x_y*_pad_divide_x_y; _sub_sub_canvas_index++) {
                    auto _sub_pad = (TPad*) _pad_canvas->GetPad(_sub_sub_canvas_index);
                    _sub_pad->SetLeftMargin(0);
                    _sub_pad->SetRightMargin(0);
                    _sub_pad->SetTopMargin(0);
                    _sub_pad->SetBottomMargin(0);

                    auto _board_number = _pad_module_board[_sub_sub_canvas_index - 1];
                    auto _channel_number = _pad_module_channel[_sub_sub_canvas_index - 1];
                    if (this->focal_channel_map.find(_channel_number) == this->focal_channel_map.end()) {
                        LOG(ERROR) << "Channel number " << _channel_number << " not found in the channel map!";
                        exit(1);
                    }
                    auto _h2gcroc_channel_number = this->focal_channel_map[_channel_number];
                    if (this->focal_fpga_map.find(_board_number) == this->focal_fpga_map.end()) {
                        LOG(ERROR) << "Board number " << _board_number << " not found in the fpga map!";
                        exit(1);
                    }
                    auto _h2gcroc_board_number = this->focal_fpga_map[_board_number];

                    auto _unified_channel = _h2gcroc_board_number * (FPGA_CHANNEL_NUMBER_VALID / 2)  + _h2gcroc_channel_number;
                    auto _channel_index = hists_channel_map[_unified_channel];
                    auto _hist = hists[_channel_index];

                    if (_hist->GetEntries() == 0) {
                        _hist->SetBinContent(1, 0.0);
                    }

                    _sub_pad->cd();
                    _hist->Draw("hist");
                }
                _pad->cd();
                _pad_canvas->DrawClonePad();
            }
        }
    }
}

void GlobalChannelPainter::draw_global_channel_hists1D_run_group(std::vector <std::vector <TH1D*>> hists_list, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, std::vector <EColor> colors, std::vector <std::string> legend_labels){
    this->clear_canvas();
    this->painter_canvas = new TCanvas(hist_name.c_str(), hist_title.c_str(), 800, 600);
    if (is_EEEMCal_mapping) {
        this->painter_canvas->Divide(5, 5, 0.0, 0.0);
        for (int _sub_canvas_index = 1; _sub_canvas_index <= 25; _sub_canvas_index++) {
            auto _pad = (TPad*) this->painter_canvas->GetPad(_sub_canvas_index);
            _pad->SetLeftMargin(0);
            _pad->SetRightMargin(0);
            _pad->SetTopMargin(0);
            _pad->SetBottomMargin(0);

            auto _pad_module_fpga = this->module_fpga_list[_sub_canvas_index - 1];
            auto _pad_module_asic = this->module_asic_list[_sub_canvas_index - 1];
            auto _pad_module_connector = this->module_connector_list[_sub_canvas_index - 1];
            auto _pad_connector = this->connector_list_list[_pad_module_connector-1];

            auto _pad_canvas = new TCanvas(("sub_canvas_" + std::to_string(_sub_canvas_index)).c_str(), ("Sub Canvas " + std::to_string(_sub_canvas_index)).c_str(), 800, 600);
            this->sub_canvas_list.push_back(_pad_canvas);
            _pad_canvas->Divide(4, 4, 0.0, 0.0);
            for (int _sub_sub_canvas_index = 1; _sub_sub_canvas_index <= 16; _sub_sub_canvas_index++) {
                auto _sub_pad = (TPad*) _pad_canvas->GetPad(_sub_sub_canvas_index);
                _sub_pad->SetLeftMargin(0);
                _sub_pad->SetRightMargin(0);
                _sub_pad->SetTopMargin(0);
                _sub_pad->SetBottomMargin(0);

                auto _unified_channel = _pad_connector[_sub_sub_canvas_index - 1] + _pad_module_fpga * FPGA_CHANNEL_NUMBER_VALID + _pad_module_asic * int(FPGA_CHANNEL_NUMBER_VALID / 2);
                auto _channel_index = hists_channel_map[_unified_channel];
                // auto _hist_list = hists_list[_channel_index];

                TLegend* _legend = new TLegend(0.7, 0.7, 0.89, 0.89);
                _legend->SetFillColor(0);
                _legend->SetBorderSize(0);
                _legend->SetLineColorAlpha(0, 0);

                for (int _hist_index = 0; _hist_index < hists_list.size(); _hist_index++) {
                    auto _hist = hists_list[_hist_index][_channel_index];
                    _hist->SetLineColor(colors[_hist_index % colors.size()]);
                    if (_hist->GetEntries() == 0) {
                        _hist->SetBinContent(1, 0.0);
                    }
                    _hist->SetStats(0);
                    _sub_pad->cd();
                    if (_hist_index == 0) {
                        _hist->Draw();
                    } else {
                        _hist->Draw("same");
                    }
                    _legend->AddEntry(_hist, legend_labels[_hist_index].c_str(), "l");
                }
                _legend->Draw();
            }
            _pad->cd();
            _pad_canvas->DrawClonePad();
        }
    } else {
        if (is_FoCal_mapping) {
            this->painter_canvas->Divide(3, 3, 0.0, 0.0);
            for (int _sub_canvas_index = 1; _sub_canvas_index <= 9; _sub_canvas_index++) {
                auto _pad = (TPad*) this->painter_canvas->GetPad(_sub_canvas_index);
                _pad->SetLeftMargin(0);
                _pad->SetRightMargin(0);
                _pad->SetTopMargin(0);
                _pad->SetBottomMargin(0);

                auto _pad_module_board = this->focal_module_board_list[_sub_canvas_index - 1];
                auto _pad_module_channel = this->focal_module_channel_list[_sub_canvas_index - 1];

                auto _pad_canvas = new TCanvas(("sub_canvas_" + std::to_string(_sub_canvas_index)).c_str(), ("Sub Canvas " + std::to_string(_sub_canvas_index)).c_str(), 800, 600);
                this->sub_canvas_list.push_back(_pad_canvas);
                int _pad_divide_x_y;
                if (_sub_canvas_index == 5) {
                    _pad_divide_x_y = 7;
                } else {
                    _pad_divide_x_y = 5;
                }
                _pad_canvas->Divide(_pad_divide_x_y, _pad_divide_x_y, 0.0, 0.0);
                for (int _sub_sub_canvas_index = 1; _sub_sub_canvas_index <= _pad_divide_x_y*_pad_divide_x_y; _sub_sub_canvas_index++) {
                    auto _sub_pad = (TPad*) _pad_canvas->GetPad(_sub_sub_canvas_index);
                    _sub_pad->SetLeftMargin(0);
                    _sub_pad->SetRightMargin(0);
                    _sub_pad->SetTopMargin(0);
                    _sub_pad->SetBottomMargin(0);

                    auto _board_number = _pad_module_board[_sub_sub_canvas_index - 1];
                    auto _channel_number = _pad_module_channel[_sub_sub_canvas_index - 1];
                    if (this->focal_channel_map.find(_channel_number) == this->focal_channel_map.end()) {
                        LOG(ERROR) << "Channel number " << _channel_number << " not found in the channel map!";
                        exit(1);
                    }

                    auto _h2gcroc_channel_number = this->focal_channel_map[_channel_number];
                    if (this->focal_fpga_map.find(_board_number) == this->focal_fpga_map.end()) {
                        LOG(ERROR) << "Board number " << _board_number << " not found in the fpga map!";
                        exit(1);
                    }
                    auto _h2gcroc_board_number = this->focal_fpga_map[_board_number];

                    auto _unified_channel = _h2gcroc_board_number * (FPGA_CHANNEL_NUMBER_VALID / 2)  + _h2gcroc_channel_number;
                    auto _channel_index = hists_channel_map[_unified_channel];

                    TLegend* _legend = new TLegend(0.7, 0.7, 0.89, 0.89);
                    _legend->SetFillColor(0);
                    _legend->SetBorderSize(0);
                    _legend->SetLineColorAlpha(0, 0);

                    for (int _hist_index = 0; _hist_index < hists_list.size(); _hist_index++) {
                        auto _hist = hists_list[_hist_index][_channel_index];
                        _hist->SetLineColor(colors[_hist_index % colors.size()]);
                        if (_hist->GetEntries() == 0) {
                            _hist->SetBinContent(1, 0.0);
                        }
                        _hist->SetStats(0);
                        _sub_pad->cd();
                        if (_hist_index == 0) {
                            _hist->Draw();
                        } else {
                            _hist->Draw("same");
                        }
                        _legend->AddEntry(_hist, legend_labels[_hist_index].c_str(), "l");
                    }
                    _legend->Draw();

                }
                _pad->cd();
                _pad_canvas->DrawClonePad();
            }
        }
    }
}

void GlobalChannelPainter::draw_global_channel_hists1D_group(std::vector <std::vector <TH1D*>> hists_list, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, std::vector <EColor> colors, std::vector <std::string> legend_labels){
    this->clear_canvas();
    this->painter_canvas = new TCanvas(hist_name.c_str(), hist_title.c_str(), 800, 600);
    if (is_EEEMCal_mapping) {
        this->painter_canvas->Divide(5, 5, 0.0, 0.0);
        for (int _sub_canvas_index = 1; _sub_canvas_index <= 25; _sub_canvas_index++) {
            auto _pad = (TPad*) this->painter_canvas->GetPad(_sub_canvas_index);
            _pad->SetLeftMargin(0);
            _pad->SetRightMargin(0);
            _pad->SetTopMargin(0);
            _pad->SetBottomMargin(0);

            auto _pad_module_fpga = this->module_fpga_list[_sub_canvas_index - 1];
            auto _pad_module_asic = this->module_asic_list[_sub_canvas_index - 1];
            auto _pad_module_connector = this->module_connector_list[_sub_canvas_index - 1];
            auto _pad_connector = this->connector_list_list[_pad_module_connector-1];

            auto _pad_canvas = new TCanvas(("sub_canvas_" + std::to_string(_sub_canvas_index)).c_str(), ("Sub Canvas " + std::to_string(_sub_canvas_index)).c_str(), 800, 600);
            this->sub_canvas_list.push_back(_pad_canvas);
            _pad_canvas->Divide(4, 4, 0.0, 0.0);
            for (int _sub_sub_canvas_index = 1; _sub_sub_canvas_index <= 16; _sub_sub_canvas_index++) {
                auto _sub_pad = (TPad*) _pad_canvas->GetPad(_sub_sub_canvas_index);
                _sub_pad->SetLeftMargin(0);
                _sub_pad->SetRightMargin(0);
                _sub_pad->SetTopMargin(0);
                _sub_pad->SetBottomMargin(0);

                auto _unified_channel = _pad_connector[_sub_sub_canvas_index - 1] + _pad_module_fpga * FPGA_CHANNEL_NUMBER_VALID + _pad_module_asic * int(FPGA_CHANNEL_NUMBER_VALID / 2);
                auto _channel_index = hists_channel_map[_unified_channel];
                auto _hist_list = hists_list[_channel_index];

                TLegend* _legend = new TLegend(0.7, 0.7, 0.89, 0.89);
                _legend->SetFillColor(0);
                _legend->SetBorderSize(0);
                _legend->SetLineColorAlpha(0, 0);

                for (int _hist_index = 0; _hist_index < _hist_list.size(); _hist_index++) {
                    auto _hist = _hist_list[_hist_index];
                    _hist->SetLineColor(colors[_hist_index % colors.size()]);
                    if (_hist->GetEntries() == 0) {
                        _hist->SetBinContent(1, 0.0);
                    }
                    _hist->SetStats(0);
                    _sub_pad->cd();
                    if (_hist_index == 0) {
                        _hist->Draw();
                    } else {
                        _hist->Draw("same");
                    }
                    _legend->AddEntry(_hist, legend_labels[_hist_index].c_str(), "l");
                }
                _legend->Draw();
            }
            _pad->cd();
            _pad_canvas->DrawClonePad();
        }
    } else {
        if (is_FoCal_mapping) {
            this->painter_canvas->Divide(3, 3, 0.0, 0.0);
            for (int _sub_canvas_index = 1; _sub_canvas_index <= 9; _sub_canvas_index++) {
                auto _pad = (TPad*) this->painter_canvas->GetPad(_sub_canvas_index);
                _pad->SetLeftMargin(0);
                _pad->SetRightMargin(0);
                _pad->SetTopMargin(0);
                _pad->SetBottomMargin(0);

                auto _pad_module_board = this->focal_module_board_list[_sub_canvas_index - 1];
                auto _pad_module_channel = this->focal_module_channel_list[_sub_canvas_index - 1];

                auto _pad_canvas = new TCanvas(("sub_canvas_" + std::to_string(_sub_canvas_index)).c_str(), ("Sub Canvas " + std::to_string(_sub_canvas_index)).c_str(), 800, 600);
                this->sub_canvas_list.push_back(_pad_canvas);
                int _pad_divide_x_y;
                if (_sub_canvas_index == 5) {
                    _pad_divide_x_y = 7;
                } else {
                    _pad_divide_x_y = 5;
                }
                _pad_canvas->Divide(_pad_divide_x_y, _pad_divide_x_y, 0.0, 0.0);
                for (int _sub_sub_canvas_index = 1; _sub_sub_canvas_index <= _pad_divide_x_y*_pad_divide_x_y; _sub_sub_canvas_index++) {
                    auto _sub_pad = (TPad*) _pad_canvas->GetPad(_sub_sub_canvas_index);
                    _sub_pad->SetLeftMargin(0);
                    _sub_pad->SetRightMargin(0);
                    _sub_pad->SetTopMargin(0);
                    _sub_pad->SetBottomMargin(0);

                    auto _board_number = _pad_module_board[_sub_sub_canvas_index - 1];
                    auto _channel_number = _pad_module_channel[_sub_sub_canvas_index - 1];
                    if (this->focal_channel_map.find(_channel_number) == this->focal_channel_map.end()) {
                        LOG(ERROR) << "Channel number " << _channel_number << " not found in the channel map!";
                        exit(1);
                    }

                    auto _h2gcroc_channel_number = this->focal_channel_map[_channel_number];
                    if (this->focal_fpga_map.find(_board_number) == this->focal_fpga_map.end()) {
                        LOG(ERROR) << "Board number " << _board_number << " not found in the fpga map!";
                        exit(1);
                    }
                    auto _h2gcroc_board_number = this->focal_fpga_map[_board_number];

                    auto _unified_channel = _h2gcroc_board_number * (FPGA_CHANNEL_NUMBER_VALID / 2)  + _h2gcroc_channel_number;
                    auto _channel_index = hists_channel_map[_unified_channel];
                    auto _hist_list = hists_list[_channel_index];

                    TLegend* _legend = new TLegend(0.7, 0.7, 0.89, 0.89);
                    _legend->SetFillColor(0);
                    _legend->SetBorderSize(0);
                    _legend->SetLineColorAlpha(0, 0);

                    for (int _hist_index = 0; _hist_index < _hist_list.size(); _hist_index++) {
                        auto _hist = _hist_list[_hist_index];
                        _hist->SetLineColor(colors[_hist_index % colors.size()]);
                        if (_hist->GetEntries() == 0) {
                            _hist->SetBinContent(1, 0.0);
                        }
                        _hist->SetStats(0);
                        _sub_pad->cd();
                        if (_hist_index == 0) {
                            _hist->Draw();
                        } else {
                            _hist->Draw("same");
                        }
                        _legend->AddEntry(_hist, legend_labels[_hist_index].c_str(), "l");
                    }
                    _legend->Draw();

                }
                _pad->cd();
                _pad_canvas->DrawClonePad();
            }
        }
    }
}

void GlobalChannelPainter::draw_canvas_components(TCanvas* source_canvas, TPad* target_pad) {
    target_pad->cd();
    TList* primitives = source_canvas->GetListOfPrimitives();
    bool is_first = true;
    for (int i = 0; i < primitives->GetSize(); ++i) {
        TObject* component;
        try {
            component = primitives->At(i);
        } catch (const std::exception& e) {
            LOG(ERROR) << "Exception caught: " << e.what();
            continue;
        }
        if (!component) {
            LOG(ERROR) << "Component is null!";
            continue;
        }
        TObject* component_clone;
        try {
            component_clone = component->Clone();
        } catch (const std::exception& e) {
            LOG(ERROR) << "Exception caught: " << e.what();
            continue;
        }
        // TObject* component = primitives->At(i);
        // TObject* component_clone = component->Clone();
        if (!component_clone) {
            LOG(ERROR) << "Component clone is null!";
            continue;
        }

        if (component_clone->InheritsFrom("TH1")) {
            TH1* hist = (TH1*) component_clone;
            hist->SetStats(0);
            if (is_first) {
                hist->Draw();
                is_first = false;
            } else {
                hist->Draw("same");
            }
        } else if (component_clone->InheritsFrom("TH2")) {
            TH2* hist = (TH2*) component_clone;
            hist->SetStats(0);
            // auto _draw_option = hist->GetDrawOption();
            if (is_first) {
                hist->Draw("colz");
                is_first = false;
            } else {
                hist->Draw("colz same");
            }
        } else if (component_clone->InheritsFrom("TGraphErrors")) {
            TGraphErrors* graph = (TGraphErrors*) component_clone;
            if (is_first) {
                graph->Draw("AP");
                is_first = false;
            } else {
                graph->Draw("P same");
            }
        } else if (component_clone->InheritsFrom("TGraph")) {
            TGraph* graph = (TGraph*) component_clone;
            if (is_first) {
                graph->Draw("AP");
                is_first = false;
            } else {
                graph->Draw("P same");
            }
        } else if (component_clone->InheritsFrom("TLine")) {
            TLine* line = (TLine*) component_clone;
            line->Draw();
        } else if (component_clone->InheritsFrom("TFrame")) {
            TFrame* frame = (TFrame*) component_clone;
            frame->Draw();
        } else if (component_clone->InheritsFrom("TBox")) {
            TBox* box = (TBox*) component_clone;
            box->Draw();
        } else if (component_clone->InheritsFrom("TPaveText")) {
            TPaveText* pave_text = (TPaveText*) component_clone;
            pave_text->Draw();
        } else if (component_clone->InheritsFrom("TF1")) {
            TF1* func = (TF1*) component_clone;
            if (is_first) {
                func->Draw();
                is_first = false;
            } else {
                func->Draw("same");
            }
        } else if (component_clone->InheritsFrom("TLatex")) {
            TLatex* latex = (TLatex*) component_clone;
            latex->Draw();
        } else if (component_clone->InheritsFrom("TLegend")) {
            TLegend* legend = (TLegend*) component_clone;
            legend->Draw();
        } else {
            LOG(ERROR) << "Unknown object type: " << component_clone->ClassName();
        }
    }
    target_pad->Update();
}

TCanvas* GlobalChannelPainter::draw_module_channel_canvas(std::vector <TCanvas*> canvas_list,  std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, int module_index) {
    if (is_EEEMCal_mapping) {
        return nullptr;
    } else {
        if (is_FoCal_mapping) {
            TCanvas* _module_canvas = new TCanvas(hist_name.c_str(), hist_title.c_str(), 800, 600);
            int _pad_divide_x_y;
            if (module_index == 5) {
                _pad_divide_x_y = 7;
            } else {
                _pad_divide_x_y = 5;
            }
            _module_canvas->Divide(_pad_divide_x_y, _pad_divide_x_y, 0.0, 0.0);
            auto _pad_module_board = this->focal_module_board_list[module_index - 1];
            auto _pad_module_channel = this->focal_module_channel_list[module_index - 1];
            for (int _sub_sub_canvas_index = 1; _sub_sub_canvas_index <= _pad_divide_x_y*_pad_divide_x_y; _sub_sub_canvas_index++) {
                auto _sub_pad = (TPad*) _module_canvas->GetPad(_sub_sub_canvas_index);
                _sub_pad->SetLeftMargin(0);
                _sub_pad->SetRightMargin(0);
                _sub_pad->SetTopMargin(0);
                _sub_pad->SetBottomMargin(0);

                auto _board_number = _pad_module_board[_sub_sub_canvas_index - 1];
                auto _channel_number = _pad_module_channel[_sub_sub_canvas_index - 1];
                if (this->focal_channel_map.find(_channel_number) == this->focal_channel_map.end()) {
                    LOG(ERROR) << "Channel number " << _channel_number << " not found in the channel map!";
                    exit(1);
                }
                auto _h2gcroc_channel_number = this->focal_channel_map[_channel_number];
                if (this->focal_fpga_map.find(_board_number) == this->focal_fpga_map.end()) {
                    LOG(ERROR) << "Board number " << _board_number << " not found in the fpga map!";
                    exit(1);
                }
                auto _h2gcroc_board_number = this->focal_fpga_map[_board_number];

                auto _unified_channel = _h2gcroc_board_number * (FPGA_CHANNEL_NUMBER_VALID / 2)  + _h2gcroc_channel_number;
                auto _channel_index = hists_channel_map[_unified_channel];
                auto _canvas = canvas_list[_channel_index];

                draw_canvas_components(_canvas, _sub_pad);
            }
            return _module_canvas;
        }
    }
    return nullptr;
}

void GlobalChannelPainter::draw_global_channel_canvas(std::vector <TCanvas*> canvas_list,  std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title) {
    this->clear_canvas();
    this->painter_canvas = new TCanvas(hist_name.c_str(), hist_title.c_str(), 3200, 3200);
    if (is_EEEMCal_mapping) {
        this->painter_canvas->Divide(5, 5, 0.0, 0.0);
        for (int _sub_canvas_index = 1; _sub_canvas_index <= 25; _sub_canvas_index++) {
            auto _pad = (TPad*) this->painter_canvas->GetPad(_sub_canvas_index);
            _pad->SetLeftMargin(0);
            _pad->SetRightMargin(0);
            _pad->SetTopMargin(0);
            _pad->SetBottomMargin(0);

            auto _pad_module_fpga = this->module_fpga_list[_sub_canvas_index - 1];
            auto _pad_module_asic = this->module_asic_list[_sub_canvas_index - 1];
            auto _pad_module_connector = this->module_connector_list[_sub_canvas_index - 1];
            auto _pad_connector = this->connector_list_list[_pad_module_connector-1];

            auto _pad_canvas = new TCanvas(("sub_canvas_" + std::to_string(_sub_canvas_index)).c_str(), ("Sub Canvas " + std::to_string(_sub_canvas_index)).c_str(), 800, 600);
            this->sub_canvas_list.push_back(_pad_canvas);
            _pad_canvas->Divide(4, 4, 0.0, 0.0);
            for (int _sub_sub_canvas_index = 1; _sub_sub_canvas_index <= 16; _sub_sub_canvas_index++) {
                auto _sub_pad = (TPad*) _pad_canvas->GetPad(_sub_sub_canvas_index);
                _sub_pad->SetLeftMargin(0);
                _sub_pad->SetRightMargin(0);
                _sub_pad->SetTopMargin(0);
                _sub_pad->SetBottomMargin(0);

                auto _unified_channel = _pad_connector[_sub_sub_canvas_index - 1] + _pad_module_fpga * FPGA_CHANNEL_NUMBER_VALID + _pad_module_asic * int(FPGA_CHANNEL_NUMBER_VALID / 2);
                auto _channel_index = hists_channel_map[_unified_channel];
                auto _canvas = canvas_list[_channel_index];
                draw_canvas_components(_canvas, _sub_pad);
            }
            _pad->cd();
            _pad_canvas->DrawClonePad();
            _pad_canvas->Update();
        }
    } else {
        if (is_FoCal_mapping) {
            this->painter_canvas->Divide(3, 3, 0.0, 0.0);
            for (int _sub_canvas_index = 1; _sub_canvas_index <= 9; _sub_canvas_index++) {
                auto _pad = (TPad*) this->painter_canvas->GetPad(_sub_canvas_index);
                _pad->SetLeftMargin(0);
                _pad->SetRightMargin(0);
                _pad->SetTopMargin(0);
                _pad->SetBottomMargin(0);

                auto _pad_module_board = this->focal_module_board_list[_sub_canvas_index - 1];
                auto _pad_module_channel = this->focal_module_channel_list[_sub_canvas_index - 1];

                auto _pad_canvas = new TCanvas(("sub_canvas_" + std::to_string(_sub_canvas_index)).c_str(), ("Sub Canvas " + std::to_string(_sub_canvas_index)).c_str(), 800, 600);
                this->sub_canvas_list.push_back(_pad_canvas);
                int _pad_divide_x_y;
                if (_sub_canvas_index == 5) {
                    _pad_divide_x_y = 7;
                } else {
                    _pad_divide_x_y = 5;
                }
                _pad_canvas->Divide(_pad_divide_x_y, _pad_divide_x_y, 0.0, 0.0);
                for (int _sub_sub_canvas_index = 1; _sub_sub_canvas_index <= _pad_divide_x_y*_pad_divide_x_y; _sub_sub_canvas_index++) {
                    auto _sub_pad = (TPad*) _pad_canvas->GetPad(_sub_sub_canvas_index);
                    _sub_pad->SetLeftMargin(0);
                    _sub_pad->SetRightMargin(0);
                    _sub_pad->SetTopMargin(0);
                    _sub_pad->SetBottomMargin(0);

                    auto _board_number = _pad_module_board[_sub_sub_canvas_index - 1];
                    auto _channel_number = _pad_module_channel[_sub_sub_canvas_index - 1];
                    if (this->focal_channel_map.find(_channel_number) == this->focal_channel_map.end()) {
                        LOG(ERROR) << "Channel number " << _channel_number << " not found in the channel map!";
                        exit(1);
                    }
                    auto _h2gcroc_channel_number = this->focal_channel_map[_channel_number];

                    if (this->focal_fpga_map.find(_board_number) == this->focal_fpga_map.end()) {
                        LOG(ERROR) << "Board number " << _board_number << " not found in the fpga map!";
                        exit(1);
                    }

                    auto _h2gcroc_board_number = this->focal_fpga_map[_board_number];

                    auto _unified_channel = _h2gcroc_board_number * (FPGA_CHANNEL_NUMBER_VALID / 2)  + _h2gcroc_channel_number;
                    auto _channel_index = hists_channel_map[_unified_channel];
                    auto _canvas = canvas_list[_channel_index];
                    draw_canvas_components(_canvas, _sub_pad);
                }
                _pad->cd();
                _pad_canvas->DrawClonePad();
                _pad_canvas->Update();
            }
        }
    }
    this->painter_canvas->Update();
}
