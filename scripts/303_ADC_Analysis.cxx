#include "H2GCROC_Common.hxx"
#include "H2GCROC_Lib.hxx"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv) {
    ScriptOptions opts = parse_arguments_single_root(argc, argv, "1.0");

    std::vector<Color_t> fpga_colors = {
        kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kPink, kViolet
    };

    gROOT->SetBatch(kTRUE);
    const double sample_time = 25.0; // unit: ns
    const double phase_shift_time = 1.5625; // unit: ns
    // get run number from the input file name
    size_t run_number_pos = opts.input_file.find("Run");
    std::string run_info_str = "Run ";
    if (run_number_pos != std::string::npos) {
        run_info_str += opts.input_file.substr(run_number_pos + 3, 3);
    } else {
        run_info_str += "Unknown";
    }

    // * --- Read the root file ---------------------------------------------------------
    // * --------------------------------------------------------------------------------
    int fpga_count = -1;
    int machine_gun_samples = -1;
    int entry_max = -1;

    TFile *input_root;
    TTree *input_tree;
    std::vector<UShort_t> legal_fpga_id_list;

    if (!readRootMetaData(opts.input_file.c_str(), input_root, input_tree, fpga_count, machine_gun_samples, entry_max, legal_fpga_id_list)) {
        LOG(ERROR) << "Failed to read metadata from input file " << opts.input_file;
        return 1;
    }

    if (entry_max == 0) {
        LOG(ERROR) << "No events in the input file!";
        return 1;
    }
    if (opts.n_events > 0 && opts.n_events < entry_max) {
        entry_max = opts.n_events;
    } else {
        if (opts.n_events > entry_max) {
            LOG(WARNING) << "Requested number of events " << opts.n_events << " is larger than the number of events in the input file " << entry_max;
        }
    }

    // std::vector<int> interested_channels = {68, 72, 62, 58, 54, 50, 46, 42, 74, 70, 64, 60, 52, 48, 44, 40, 71, 67 ,63, 59, 55, 49, 43, 39, 73, 69, 65, 61, 53, 51, 45, 41};
    std::vector<int> interested_but_covered_channels = {54, 58};

    std::vector<int> interested_channels = {50};
    // * --- Create the output file -----------------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *output_root = new TFile(opts.output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << opts.output_file;
        return 1;
    }

    // std::unordered_map <int, int> unifiedToHistIndex;
    // std::unordered_map <int, int> histIndexToUnified;

    // int channel_array_index = 0;
    // int sample_time_bins = int((machine_gun_samples * sample_time) / phase_shift_time);
    // for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
    //     auto _fpga_id = legal_fpga_id_list[_fpga_index];
    //     for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
    //         auto _channel_valid = get_valid_fpga_channel(_channel_index);
    //         if (_channel_valid == -1){
    //             continue;
    //         }
    //         auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
    //         unifiedToHistIndex[_unified_valid_channel_number] = channel_array_index;
    //         histIndexToUnified[channel_array_index] = _unified_valid_channel_number;

    //         channel_array_index++;
    //     } // end of channel loop
    // } // end of fpga loop
   
    LOG(INFO) << "FPGA count: " << fpga_count << " Machine gun samples: " << machine_gun_samples << " Entry max: " << entry_max;

    std::vector <ULong64_t*> branch_timestamps_list;
    std::vector <UInt_t*> branch_daqh_list_list;
    std::vector <Bool_t*> branch_tc_list_list;
    std::vector <Bool_t*> branch_tp_list_list;
    std::vector <UInt_t*> branch_val0_list_list;
    std::vector <UInt_t*> branch_val1_list_list;
    std::vector <UInt_t*> branch_val2_list_list;
    std::vector <UInt_t*> branch_crc32_list_list;
    std::vector <UInt_t*> branch_last_heartbeat_list;

    for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
        auto _fpga_id = legal_fpga_id_list[_fpga_index];
        auto *branch_timestamps = new ULong64_t[machine_gun_samples];  // 64 bits
        auto *branch_daqh_list  = new UInt_t[4 * machine_gun_samples];    // 32 bits
        auto *branch_tc_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_tp_list    = new Bool_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_val0_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_val1_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_val2_list  = new UInt_t[FPGA_CHANNEL_NUMBER * machine_gun_samples];
        auto *branch_crc32_list = new UInt_t[4 * machine_gun_samples];
        auto *branch_last_heartbeat = new UInt_t[machine_gun_samples];

        input_tree->SetBranchAddress(("timestamps_"    + std::to_string(_fpga_id)).c_str(), branch_timestamps);
        input_tree->SetBranchAddress(("daqh_list_"     + std::to_string(_fpga_id)).c_str(), branch_daqh_list);
        input_tree->SetBranchAddress(("tc_list_"       + std::to_string(_fpga_id)).c_str(), branch_tc_list);
        input_tree->SetBranchAddress(("tp_list_"       + std::to_string(_fpga_id)).c_str(), branch_tp_list);
        input_tree->SetBranchAddress(("val0_list_"     + std::to_string(_fpga_id)).c_str(), branch_val0_list);
        input_tree->SetBranchAddress(("val1_list_"     + std::to_string(_fpga_id)).c_str(), branch_val1_list);
        input_tree->SetBranchAddress(("val2_list_"     + std::to_string(_fpga_id)).c_str(), branch_val2_list);
        input_tree->SetBranchAddress(("crc32_list_"    + std::to_string(_fpga_id)).c_str(), branch_crc32_list);
        input_tree->SetBranchAddress(("last_heartbeat_"+ std::to_string(_fpga_id)).c_str(), branch_last_heartbeat);

        branch_timestamps_list.push_back(branch_timestamps);
        branch_daqh_list_list.push_back(branch_daqh_list);
        branch_tc_list_list.push_back(branch_tc_list);
        branch_tp_list_list.push_back(branch_tp_list);
        branch_val0_list_list.push_back(branch_val0_list);
        branch_val1_list_list.push_back(branch_val1_list);
        branch_val2_list_list.push_back(branch_val2_list);
        branch_crc32_list_list.push_back(branch_crc32_list);
        branch_last_heartbeat_list.push_back(branch_last_heartbeat);
    }

    long counter_total_events = 0;
    long counter_hamming_code_error = 0;
    long counter_bad_daqh_start_end = 0;

    long counter_double_tot = 0;
    long counter_double_toa = 0;
    long counter_valid_tot  = 0;
    long counter_valid_toa  = 0;

    std::vector<TH2D*> h2d_adc_channel_sample_list;
    for (int i=0; i<FPGA_CHANNEL_NUMBER * fpga_count; i++) {
        std::string hist_name = "h2d_adc_samples_channel_" + std::to_string(i);
        TH2D *h2d_adc_channel_sample = new TH2D(
            hist_name.c_str(), (hist_name + ";Sample Index;ADC Value").c_str(),
            machine_gun_samples, 0, machine_gun_samples,
            512, 0, 1024);
        h2d_adc_channel_sample->SetDirectory(nullptr);
        h2d_adc_channel_sample_list.push_back(h2d_adc_channel_sample);
    }

    std::vector<TH2D*> h2d_tot_channel_sample_list;
    for (int i=0; i<FPGA_CHANNEL_NUMBER * fpga_count; i++) {
        std::string hist_name = "h2d_tot_samples_channel_" + std::to_string(i);
        TH2D *h2d_tot_channel_sample = new TH2D(
            hist_name.c_str(), (hist_name + ";Sample Index;TOT Value").c_str(),
            machine_gun_samples, 0, machine_gun_samples,
            512, 0, 1024);
        h2d_tot_channel_sample->SetDirectory(nullptr);
        h2d_tot_channel_sample_list.push_back(h2d_tot_channel_sample);
    }

    std::vector<TH2D*> h2d_toa_channel_sample_list;
    for (int i=0; i<FPGA_CHANNEL_NUMBER * fpga_count; i++) {
        std::string hist_name = "h2d_toa_samples_channel_" + std::to_string(i);
        TH2D *h2d_toa_channel_sample = new TH2D(
            hist_name.c_str(), (hist_name + ";Sample Index;TOA Value").c_str(),
            machine_gun_samples, 0, machine_gun_samples,
            512, 0, 1024);
        h2d_toa_channel_sample->SetDirectory(nullptr);
        h2d_toa_channel_sample_list.push_back(h2d_toa_channel_sample);
    }

    std::vector<TH2D*> h2d_adc_toa_correlation_channel_list;
    for (int i=0; i<FPGA_CHANNEL_NUMBER * fpga_count; i++) {
        std::string hist_name = "h2d_adc_toa_correlation_channel_" + std::to_string(i);
        TH2D *h2d_adc_toa_correlation_channel = new TH2D(
            hist_name.c_str(), (hist_name + ";TOA Value;ADC Value").c_str(),
            (machine_gun_samples+4)*16, -4.0 * sample_time, (machine_gun_samples) * sample_time,
            512, 0, 1024);
        h2d_adc_toa_correlation_channel->SetDirectory(nullptr);
        h2d_adc_toa_correlation_channel_list.push_back(h2d_adc_toa_correlation_channel);
    }

    // * Th2Ds with ToA
    std::vector<TH2D*> h2d_waveform_adc_list;
    for (int i=0; i<FPGA_CHANNEL_NUMBER * fpga_count; i++) {
        std::string hist_name = "h2d_waveform_adc_channel_" + std::to_string(i);
        TH2D *h2d_waveform_adc = new TH2D(
            hist_name.c_str(), (hist_name + ";TOA (ns);ADC Value").c_str(),
            (machine_gun_samples+4)*16, -4.0 * sample_time, (machine_gun_samples) * sample_time,
            512, 0, 1024);
        h2d_waveform_adc->SetDirectory(nullptr);
        h2d_waveform_adc_list.push_back(h2d_waveform_adc);
    }

    // * TH1Ds for ADC Peak distribution
    std::vector<TH1D*> h1d_adc_peak_list;
    for (int i=0; i<FPGA_CHANNEL_NUMBER * fpga_count; i++) {
        std::string hist_name = "h1d_adc_peak_channel_" + std::to_string(i);
        TH1D *h1d_adc_peak = new TH1D(
            hist_name.c_str(), (hist_name + ";ADC Peak Value;Counts").c_str(),
            256, 0, 1024);
        h1d_adc_peak->SetDirectory(nullptr);
        h1d_adc_peak_list.push_back(h1d_adc_peak);
    }

    TH1D *h1d_adc_average = new TH1D(
        "h1d_adc_average", "Average ADC Value;ADC Value;Counts",
        256, 0, 1024 * 1.2);
    h1d_adc_average->SetDirectory(nullptr);

    TH1D *h1d_toa_average = new TH1D(
        "h1d_toa_average", "Average TOA Value;TOA Value (ns);Counts",
        256, 0, 1024 * 1.2);
    h1d_toa_average->SetDirectory(nullptr);

    TH1D *h1d_tot_average = new TH1D(
        "h1d_tot_average", "Average TOT Value;TOT Value;Counts",
        256, 0, 4096 * 1.2);
    h1d_toa_average->SetDirectory(nullptr);

    // * TH1Ds for TOA raw value distribution
    std::vector<TH1D*> h1d_toa_list;
    for (int i=0; i<FPGA_CHANNEL_NUMBER * fpga_count; i++) {
        std::string hist_name = "h1d_toa_channel_" + std::to_string(i);
        TH1D *h1d_toa = new TH1D(
            hist_name.c_str(), (hist_name + ";TOA Value;Counts").c_str(),
            256, 0, 1024);
        h1d_toa->SetDirectory(nullptr);
        h1d_toa_list.push_back(h1d_toa);
    }

    // * TH1Ds for TOT raw value distribution
    std::vector<TH1D*> h1d_tot_list;
    for (int i=0; i<FPGA_CHANNEL_NUMBER * fpga_count; i++) {
        std::string hist_name = "h1d_tot_channel_" + std::to_string(i);
        TH1D *h1d_tot = new TH1D(
            hist_name.c_str(), (hist_name + ";TOT Value;Counts").c_str(),
            256, 0, 4096);
        h1d_tot->SetDirectory(nullptr);
        h1d_tot_list.push_back(h1d_tot);
    }

    std::vector<TH1D*> h1d_toa_sample_index_channel_list;
    for (int i=0; i<FPGA_CHANNEL_NUMBER * fpga_count; i++) {
        std::string hist_name = "h1d_toa_sample_index_channel_" + std::to_string(i);
        TH1D *h1d_toa_sample_index_channel = new TH1D(
            hist_name.c_str(), (hist_name + ";Sample Index of First ToA;Counts").c_str(),
            machine_gun_samples, 0, machine_gun_samples);
        h1d_toa_sample_index_channel->SetDirectory(nullptr);
        h1d_toa_sample_index_channel_list.push_back(h1d_toa_sample_index_channel);
    }

    std::vector<TH1D*> h1d_adc_peak_index_channel_list;
    for (int i=0; i<FPGA_CHANNEL_NUMBER * fpga_count; i++) {
        std::string hist_name = "h1d_adc_peak_index_channel_" + std::to_string(i);
        TH1D *h1d_adc_peak_index_channel = new TH1D(
            hist_name.c_str(), (hist_name + ";Sample Index of ADC Peak;Counts").c_str(),
            machine_gun_samples, 0, machine_gun_samples);
        h1d_adc_peak_index_channel->SetDirectory(nullptr);
        h1d_adc_peak_index_channel_list.push_back(h1d_adc_peak_index_channel);
    }

    for (int _entry = 0; _entry < entry_max; _entry++) {
        input_tree->GetEntry(_entry);
        counter_total_events++;

        if (_entry % 5000 == 0) {
            LOG(INFO) << "Processing entry " << _entry << " / " << entry_max;
        }

        double event_adc_sum = 0.0;
        double event_adc_count = 0.0;
        double event_tot_sum = 0.0;
        double event_tot_count = 0.0;
        double event_toa_count = 0.0;
        double event_toa_sum = 0.0;

        for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
            auto _fpga_id    = legal_fpga_id_list[_fpga_index];
            auto _timestamp  = branch_timestamps_list[_fpga_index][0];
            auto _daqh_list  = branch_daqh_list_list[_fpga_index];
            auto _tc_list    = branch_tc_list_list[_fpga_index];
            auto _tp_list    = branch_tp_list_list[_fpga_index];
            auto _val0_list  = branch_val0_list_list[_fpga_index];
            auto _val1_list  = branch_val1_list_list[_fpga_index];
            auto _val2_list  = branch_val2_list_list[_fpga_index];

            double _fpga_event_adc_sum_a0 = -1.0;
            double _fpga_event_adc_sum_a1 = -1.0;
            double _fpga_event_tot_sum_a0 = -1.0;
            double _fpga_event_tot_sum_a1 = -1.0;
            double _fpga_event_tot_count_a0 = -1.0;
            double _fpga_event_tot_count_a1 = -1.0;
            double _fpga_event_toa_count_a0 = -1.0;
            double _fpga_event_toa_count_a1 = -1.0;

            double _hamming_code_error_count_a0 = 0.0;
            double _hamming_code_error_count_a1 = 0.0;

            double _daqh_error_count_a0 = 0.0;
            double _daqh_error_count_a1 = 0.0;

            // -- Check Hamming code and header bytes --
            // -----------------------------------------
            std::vector <bool> _hamming_code_pass_list;
            bool skip_event_hamming_code = false;
            bool good_header_bytes = true;
            for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                bool _hamming_code_pass = true;
                for (int _daqh_index = 0; _daqh_index < 4; _daqh_index++){
                    auto _daqh = _daqh_list[_sample_index * 4 + _daqh_index];
                    auto _h1h2h3 = (_daqh >> 4) & 0x7;
                    if (_h1h2h3 != 0x00){
                        _hamming_code_pass = false;
                        if (_daqh_index < 2) {
                            _hamming_code_error_count_a0++;
                        } else {
                            _hamming_code_error_count_a1++;
                        }
                    }
                    auto _daqh_first_half_byte = (_daqh >> 28) & 0xf;
                    auto _daqh_last_half_byte  = _daqh & 0xf;
                    if (_daqh_first_half_byte != 5 ||  _daqh_last_half_byte != 5) {
                        good_header_bytes = false;
                        if (_daqh_index < 2) {
                            _daqh_error_count_a0++;
                        } else {
                            _daqh_error_count_a1++;
                        }
                    }
                }
                _hamming_code_pass_list.push_back(_hamming_code_pass);
                if (!_hamming_code_pass){
                    skip_event_hamming_code = true;
                    break;
                }
            }
            if (skip_event_hamming_code) {
                counter_hamming_code_error++;
            }
            if (!good_header_bytes) {
                counter_bad_daqh_start_end++;
            }

            std::vector <int> _channel_val0_max_list;
            std::vector <int> _channel_val1_max_list;
            std::vector <int> _channel_val2_max_list;
            std::vector <int> _channel_val0_max_index_list;
            std::vector <int> _channel_val1_max_index_list;
            std::vector <int> _channel_val2_max_index_list;

            for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                auto _channel_valid = get_valid_fpga_channel(_channel_index);
                if (_channel_valid == -1){
                    
                    _channel_val0_max_list.push_back(-1);
                    _channel_val1_max_list.push_back(-1);
                    _channel_val2_max_list.push_back(-1);
                    _channel_val0_max_index_list.push_back(-1);
                    _channel_val1_max_index_list.push_back(-1);
                    _channel_val2_max_index_list.push_back(-1);
                    continue;
                }
                auto _unified_valid_channel_number = get_unified_valid_fpga_channel(_fpga_id, _channel_valid);
                (void)_unified_valid_channel_number;
                const int hist_index = _fpga_index * FPGA_CHANNEL_NUMBER + _channel_index;

                int _val0_max = -1;
                int _val0_max_index = -1;
                int _val1_max = -1;
                int _val1_max_index = -1;
                int _val2_max = -1;
                int _val2_max_index = -1;

                int _adc_pedestal = 0; // TODO: get the pedestal value from the metadata

                std::vector<int> _adc_samples;
                int _toa_first = 0;
                int _toa_first_index = -1;
                int _tot_first = 0;
                int _tot_first_index = -1;

                for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                    const int _val0 = static_cast<int>(_val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER] & 0x3FFu);
                    const int _val1 = static_cast<int>(_val1_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER] & 0x3FFu);
                    const int _val2 = static_cast<int>(_val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER] & 0x3FFu);
                    if (_sample_index == 0) {
                        _adc_pedestal = _val0; // assuming the first sample is the pedestal
                    }
                    _adc_samples.push_back(_val0);
                    

                    if (_val0 > 0 && _val0 > _val0_max){
                        _val0_max = _val0;
                        _val0_max_index = _sample_index;
                    }
                    if (_val1 > 0){
                        if (_tot_first == 0) {
                            _tot_first = _val1;
                            _tot_first_index = _sample_index;
                        }
                        if (_val1_max == -1) {
                            _val1_max = _val1;
                            _val1_max_index = _sample_index;
                            counter_valid_tot++;
                            if (_channel_index < 76) { // for the first ASIC
                                if (_fpga_event_tot_sum_a0 < 0) {
                                    _fpga_event_tot_sum_a0 = 0;
                                }
                                if (_fpga_event_tot_count_a0 < 0) {
                                    _fpga_event_tot_count_a0 = 0;
                                }
                            } else {
                                if (_fpga_event_tot_sum_a1 < 0) {
                                    _fpga_event_tot_sum_a1 = 0;
                                }
                                if (_fpga_event_tot_count_a1 < 0) {
                                    _fpga_event_tot_count_a1 = 0;
                                }
                            }
                            double _val1_decoded = static_cast<double>(_val1);
                            if (_val1_decoded > 512) {
                                _val1_decoded -= 512;
                                _val1_decoded *= 8.0;
                            }
                            if (_channel_index < 76) { // for the first ASIC
                                _fpga_event_tot_sum_a0 += _val1_decoded;
                                _fpga_event_tot_count_a0 += 1.0;
                            } else {
                                _fpga_event_tot_sum_a1 += _val1_decoded;
                                _fpga_event_tot_count_a1 += 1.0;
                            }

                        } else
                            counter_double_tot++;
                    }
                    if (_val2 > 0) {
                        if (_toa_first == 0) {
                            _toa_first = _val2;
                            _toa_first_index = _sample_index;
                        }
                        if (_val2_max == -1) {
                            _val2_max = _val2;
                            _val2_max_index = _sample_index;
                            counter_valid_toa++;
                            if (_channel_index < 76) { // for the first ASIC
                                if (_fpga_event_toa_count_a0 < 0) {
                                    _fpga_event_toa_count_a0 = 0;
                                }
                                _fpga_event_toa_count_a0 += 1.0;
                            } else {
                                if (_fpga_event_toa_count_a1 < 0) {
                                    _fpga_event_toa_count_a1 = 0;
                                }
                                _fpga_event_toa_count_a1 += 1.0;
                            }
                        } else
                            counter_double_toa++;
                    }
                } // end of sample loop
                _channel_val0_max_list.push_back(_val0_max);
                _channel_val1_max_list.push_back(_val1_max);
                _channel_val2_max_list.push_back(_val2_max);
                // LOG(DEBUG) << "toa code raw: " << _val2_max;
                _channel_val0_max_index_list.push_back(_val0_max_index);
                _channel_val1_max_index_list.push_back(_val1_max_index);
                _channel_val2_max_index_list.push_back(_val2_max_index);

                if (_val0_max >= 0) {
                    double _adc_peak = static_cast<double>(_val0_max - _adc_pedestal);
                    h1d_adc_peak_index_channel_list[hist_index]->Fill(static_cast<double>(_val0_max_index));
                    if (_adc_peak > 0 && _toa_first_index <= 7 && (_val0_max_index >= 4 && _val0_max_index <= 7)) {
                        h1d_adc_peak_list[hist_index]->Fill(_adc_peak);
                        // if the channel is interested and not covered, add to the event sum
                        if (std::find(interested_channels.begin(), interested_channels.end(), _channel_index) != interested_channels.end() &&
                            std::find(interested_but_covered_channels.begin(), interested_but_covered_channels.end(), _channel_index) == interested_but_covered_channels.end()) {
                            event_adc_sum += static_cast<double>(_adc_peak);
                            event_adc_count += 1.0;
                        }
                    }
                }

                if (_toa_first > 0) {
                    h1d_toa_list[hist_index]->Fill(static_cast<double>(_toa_first));
                    h1d_toa_sample_index_channel_list[hist_index]->Fill(static_cast<double>(_toa_first_index));
                    // if the channel is interested and not covered, add to the event sum
                    if (std::find(interested_channels.begin(), interested_channels.end(), _channel_index) != interested_channels.end() &&
                        std::find(interested_but_covered_channels.begin(), interested_but_covered_channels.end(), _channel_index) == interested_but_covered_channels.end()) {
                        event_toa_sum += static_cast<double>(_toa_first);
                        event_toa_count += 1.0;
                    }
                }

                if (_tot_first > 0 && _toa_first_index <= 7) {
                    double _tot_first_decoded = static_cast<double>(_tot_first);
                    if (_tot_first_decoded > 512) {
                        _tot_first_decoded -= 512;
                        _tot_first_decoded *= 8.0;
                    }
                    h1d_tot_list[hist_index]->Fill(_tot_first_decoded);
                    // if the channel is interested and not covered, add to the event sum
                    if (std::find(interested_channels.begin(), interested_channels.end(), _channel_index) != interested_channels.end() &&
                        std::find(interested_but_covered_channels.begin(), interested_but_covered_channels.end(), _channel_index) == interested_but_covered_channels.end()) {
                        event_tot_sum += _tot_first_decoded;
                        event_tot_count += 1.0;
                    }
                }


                for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                    h2d_adc_channel_sample_list[hist_index]->Fill(_sample_index, _adc_samples[_sample_index]);
                    if (_toa_first > 0) {
                        double _toa_first_ns = static_cast<double>(_toa_first) * 0.025; // assuming the ToA unit is 25 ps
                        _toa_first_ns += static_cast<double>(_toa_first_index) * sample_time - 50.0;
                        h2d_waveform_adc_list[hist_index]->Fill(_sample_index * sample_time - _toa_first_ns, _adc_samples[_sample_index]);
                    }
                }

                if (_val0_max >= 0 && _toa_first_index <= 7) {
                    if (_channel_index < 76) { // for the first ASIC
                        if (_fpga_event_adc_sum_a0 < 0) {
                            _fpga_event_adc_sum_a0 = 0;
                        }
                        _fpga_event_adc_sum_a0 += _val0_max - _adc_pedestal;
                    } else {
                        if (_fpga_event_adc_sum_a1 < 0) {
                            _fpga_event_adc_sum_a1 = 0;
                        }
                        _fpga_event_adc_sum_a1 += _val0_max - _adc_pedestal;
                    }
                    if (_toa_first > 0){
                        // fill the ADC-TOA correlation histogram
                        double _toa_first_ns = static_cast<double>(_toa_first) * 0.025; // assuming the ToA unit is 25 ps
                        _toa_first_ns += static_cast<double>(_toa_first_index) * sample_time - 50.0;
                        h2d_adc_toa_correlation_channel_list[hist_index]->Fill(_toa_first_ns, _val0_max - _adc_pedestal);
                    }
                }

            } // end of channel loop
        } // end of fpga loop
        if (event_adc_count > 0) {
            h1d_adc_average->Fill(event_adc_sum / event_adc_count);
        }
        if (event_toa_count > 0) {
            h1d_toa_average->Fill(event_toa_sum / event_toa_count);
        }
        if (event_tot_count > 0) {
            h1d_tot_average->Fill(event_tot_sum / event_tot_count);
        }
    } // end of entry loop

    LOG(INFO) << "=== Data Reading Summary ===============================";
    LOG(INFO) << "Total events: " << counter_total_events;
    // LOG(INFO) << "Hamming code errors: " << counter_hamming_code_error << " (" << (double(counter_hamming_code_error) / counter_total_events * 100.0) << "%)";
    // LOG(INFO) << "Bad DAQH start/end: " << counter_bad_daqh_start_end << " (" << (double(counter_bad_daqh_start_end) / counter_total_events * 100.0) << "%)";
    LOG(INFO) << "Valid TOT: " << counter_valid_tot << " (" << (double(counter_valid_tot) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "Valid TOA: " << counter_valid_toa << " (" << (double(counter_valid_toa) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "Double TOT: " << counter_double_tot << " (" << (double(counter_double_tot) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "Double TOA: " << counter_double_toa << " (" << (double(counter_double_toa) / counter_total_events / machine_gun_samples / FPGA_CHANNEL_NUMBER / fpga_count * 100.0) << "%)";
    LOG(INFO) << "=== End of Data Reading Summary ========================";

    input_root->Close();

    std::string mapping_json_file = "config/mapping_Feb2026_re.json";
    std::ifstream mapping_json_ifs(mapping_json_file);
    if (!mapping_json_ifs.is_open()) {
        LOG(ERROR) << "Failed to open mapping JSON file: " << mapping_json_file;
        return 1;
    }

    json mapping_json;
    mapping_json_ifs >> mapping_json;
    const auto& sipm_board      = mapping_json.at("SiPM_Board");
    const auto& board_loc       = mapping_json.at("Board_Loc");
    const auto& board_rotation  = mapping_json.at("Board_Rotation");
    const auto& board_flip      = mapping_json.at("Board_Flip");

    const int NX = 16, NY = 4;
    const int board_cols = 8, board_rows = 4;
    const int NPIX = NX * NY;
    const int vldb_number = fpga_count;
    const int TOTAL_CH = FPGA_CHANNEL_NUMBER * vldb_number;

    std::vector<int> ch2pid(TOTAL_CH, -1);
    std::vector<int> pid2ch(NPIX, -1);

    auto chan2pad = build_chan2pad_LUT(
        vldb_number, FPGA_CHANNEL_NUMBER,
        NX, NY, board_cols, board_rows,
        sipm_board, board_loc, board_rotation, board_flip
    );

    int used = 0;
    std::vector<int> row_count(NY, 0);
    for (int i = 0; i < (int)chan2pad.size(); ++i) {
        int pad = chan2pad[i];
        if (pad < 0) continue;
        used++;
        int row = pad / NX;
        if (row >= 0 && row < NY) row_count[row]++;
    }
    LOG(INFO) << "Mapped pads: " << used << " / " << chan2pad.size();
    for (int r = 0; r < NY; ++r) {
        LOG(INFO) << "Row " << r << ": " << row_count[r];
    }

    MosaicTopology topo_wave;
    topo_wave.NX = NX;
    topo_wave.NY = NY;
    topo_wave.vldb_number = vldb_number;
    topo_wave.channels_per_vldb = FPGA_CHANNEL_NUMBER;
    topo_wave.reverse_row = true;
    topo_wave.minimalist_axis = true;
    topo_wave.th2_logz = true;
    topo_wave.chan2pad = chan2pad;

    MosaicTopology topo_ped_median = topo_wave;
    topo_ped_median = topo_wave;
    topo_ped_median.th2_logz = true;

    // save each of the channel sample histograms to the output root file's subdirectory
    output_root->mkdir("Channel_Samples");
    output_root->cd("Channel_Samples");
    for (size_t i = 0; i < h2d_adc_channel_sample_list.size(); ++i) {
        h2d_adc_channel_sample_list[i]->Write();
    }

    // save each of the waveform histograms to the output root file's subdirectory
    output_root->mkdir("Waveforms");
    output_root->cd("Waveforms");
    for (size_t i = 0; i < h2d_waveform_adc_list.size(); ++i) {
        h2d_waveform_adc_list[i]->Write();
    }

    // save each of the ADC peak histograms to the output root file's subdirectory
    output_root->mkdir("ADC_Peaks");
    output_root->cd("ADC_Peaks");
    for (size_t i = 0; i < h1d_adc_peak_list.size(); ++i) {
        h1d_adc_peak_list[i]->Write();
    }

    output_root->cd();

    TCanvas *canvas_channel_samples = new TCanvas("canvas_channel_samples", "ADC Samples for Each Channel", 1200, 800);
    draw_mosaic_fixed(*canvas_channel_samples, h2d_adc_channel_sample_list, topo_wave);
    canvas_channel_samples->Modified();
    canvas_channel_samples->Update();
    canvas_channel_samples->Write();
    // save as png as well for quick checking
    // std::string canvas_sample_png = opts.output_folder + "/adc_samples_mosaic_" + run_info_str + ".png";
    // canvas_channel_samples->SaveAs(canvas_sample_png.c_str());
    canvas_channel_samples->Close();

    TCanvas *canvas_waveform_adc = new TCanvas("canvas_waveform_adc", "ADC Waveforms for Each Channel", 1200, 800);
    draw_mosaic_fixed(*canvas_waveform_adc, h2d_waveform_adc_list, topo_wave);
    canvas_waveform_adc->Modified();
    canvas_waveform_adc->Update();
    canvas_waveform_adc->Write();
    // std::string canvas_waveform_png = opts.output_folder + "/adc_waveforms_mosaic_" + run_info_str + ".png";
    // canvas_waveform_adc->SaveAs(canvas_waveform_png.c_str());
    canvas_waveform_adc->Close();

    // draw the ADC Th1Ds in a mosaic layout as well
    TCanvas *canvas_adc_peak = new TCanvas("canvas_adc_peak", "ADC Peak Values for Interested Channels", 1200, 800);
    draw_mosaic_fixed(*canvas_adc_peak, h1d_adc_peak_list, topo_ped_median);
    canvas_adc_peak->Modified();
    canvas_adc_peak->Update();
    canvas_adc_peak->Write();
    // std::string canvas_peak_png = opts.output_folder + "/adc_peaks_mosaic_" + run_info_str + ".png";
    // canvas_adc_peak->SaveAs(canvas_peak_png.c_str());
    canvas_adc_peak->Close();

    // draw the ADC-ToA correlation Th2Ds in a mosaic layout as well
    TCanvas *canvas_adc_toa_correlation = new TCanvas("canvas_adc_toa_correlation", "ADC-ToA Correlation for Each Channel", 1200, 800);
    draw_mosaic_fixed(*canvas_adc_toa_correlation, h2d_adc_toa_correlation_channel_list, topo_wave);
    canvas_adc_toa_correlation->Modified();
    canvas_adc_toa_correlation->Update();
    canvas_adc_toa_correlation->Write();
    // std::string canvas_adc_toa_correlation_png = opts.output_folder + "/adc_toa_correlation_mosaic_" + run_info_str + ".png";
    // canvas_adc_toa_correlation->SaveAs(canvas_adc_toa_correlation_png.c_str());
    canvas_adc_toa_correlation->Close();

    // draw the ToA sample index distribution Th1Ds in a mosaic layout as well
    TCanvas *canvas_toa_sample_index = new TCanvas("canvas_toa_sample_index", "ToA Sample Index Distribution for Each Channel", 1200, 800);
    draw_mosaic_fixed(*canvas_toa_sample_index, h1d_toa_sample_index_channel_list, topo_wave);
    canvas_toa_sample_index->Modified();
    canvas_toa_sample_index->Update();
    canvas_toa_sample_index->Write();
    // std::string canvas_toa_sample_index_png = opts.output_folder + "/toa_sample_index_mosaic_" + run_info_str + ".png";
    // canvas_toa_sample_index->SaveAs(canvas_toa_sample_index_png.c_str());
    canvas_toa_sample_index->Close();

    // draw the ADC peak sample index distribution Th1Ds in a mosaic layout as well
    TCanvas *canvas_adc_peak_sample_index = new TCanvas("canvas_adc_peak_sample_index", "ADC Peak Sample Index Distribution for Each Channel", 1200, 800);
    draw_mosaic_fixed(*canvas_adc_peak_sample_index, h1d_adc_peak_index_channel_list, topo_wave);
    canvas_adc_peak_sample_index->Modified();
    canvas_adc_peak_sample_index->Update();
    canvas_adc_peak_sample_index->Write();
    // std::string canvas_adc_peak_sample_index_png = opts.output_folder + "/adc_peak_sample_index_mosaic_" + run_info_str + ".png";
    // canvas_adc_peak_sample_index->SaveAs(canvas_adc_peak_sample_index_png.c_str());
    canvas_adc_peak_sample_index->Close();

    // draw the interested channels's ADC peak on top of each other for comparison
    TCanvas *canvas_adc_peak_overlay = new TCanvas("canvas_adc_peak_overlay", "ADC Peak Overlay for Interested Channels", 1200, 800);
    canvas_adc_peak_overlay->cd();
    TLegend *legend = new TLegend(0.7, 0.5, 0.89, 0.89);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);


    for (size_t i = 0; i < interested_channels.size(); ++i) {
        int ch = interested_channels[i];
        if (ch < 0 || ch >= (int)h1d_adc_peak_list.size()) {
            LOG(WARNING) << "Interested channel " << ch << " is out of range, skipping.";
            continue;
        }
        h1d_adc_peak_list[ch]->SetLineColor(fpga_colors[i % fpga_colors.size()]);
        h1d_adc_peak_list[ch]->SetLineWidth(2);
        // remove the statistics box
        h1d_adc_peak_list[ch]->SetStats(0);
        if (i == 0) {
            gPad->SetLogy();
            // set axis label
            h1d_adc_peak_list[ch]->GetXaxis()->SetTitle("ADC Peak Value");
            h1d_adc_peak_list[ch]->GetYaxis()->SetTitle("Counts");
            h1d_adc_peak_list[ch]->Draw("HIST");
        } else {
            h1d_adc_peak_list[ch]->Draw("HIST SAME");
        }
        legend->AddEntry(h1d_adc_peak_list[ch], ("Channel " + std::to_string(ch)).c_str(), "l");
    }
    legend->Draw();
    canvas_adc_peak_overlay->Modified();
    canvas_adc_peak_overlay->Update();
    canvas_adc_peak_overlay->Write();
    std::string canvas_peak_overlay_png = opts.output_folder + "/adc_peaks_overlay_" + run_info_str + ".png";
    canvas_adc_peak_overlay->SaveAs(canvas_peak_overlay_png.c_str());
    canvas_adc_peak_overlay->Close();

    // similarly, draw the interested channels's ToA distribution on top of each other for comparison
    TCanvas *canvas_toa_overlay = new TCanvas("canvas_toa_overlay", "ToA Overlay for Interested Channels", 1200, 800);
    canvas_toa_overlay->cd();
    TLegend *legend_toa = new TLegend(0.7, 0.5, 0.89, 0.89);
    legend_toa->SetFillStyle(0);
    legend_toa->SetBorderSize(0);

    for (size_t i = 0; i < interested_channels.size(); ++i) {
        int ch = interested_channels[i];
        if (ch < 0 || ch >= (int)h1d_toa_list.size()) {
            LOG(WARNING) << "Interested channel " << ch << " is out of range, skipping.";
            continue;
        }
        h1d_toa_list[ch]->SetLineColor(fpga_colors[i % fpga_colors.size()]);
        h1d_toa_list[ch]->SetLineWidth(2);
        // remove the statistics box
        h1d_toa_list[ch]->SetStats(0);
        if (i == 0) {
            gPad->SetLogy();
            // set axis label
            h1d_toa_list[ch]->GetXaxis()->SetTitle("ToA Value");
            h1d_toa_list[ch]->GetYaxis()->SetTitle("Counts");
            h1d_toa_list[ch]->Draw("HIST");
        } else {
            h1d_toa_list[ch]->Draw("HIST SAME");
        }
        legend_toa->AddEntry(h1d_toa_list[ch], ("Channel " + std::to_string(ch)).c_str(), "l");
    }
    legend_toa->Draw();
    canvas_toa_overlay->Modified();
    canvas_toa_overlay->Update();
    canvas_toa_overlay->Write();
    std::string canvas_toa_overlay_png = opts.output_folder + "/toa_overlay_" + run_info_str + ".png";
    canvas_toa_overlay->SaveAs(canvas_toa_overlay_png.c_str());
    canvas_toa_overlay->Close();


    // similarly, draw the interested channels's TOT distribution on top of each other for comparison
    TCanvas *canvas_tot_overlay = new TCanvas("canvas_tot_overlay", "TOT Overlay for Interested Channels", 1200, 800);
    canvas_tot_overlay->cd();
    TLegend *legend_tot = new TLegend(0.7, 0.5, 0.89, 0.89);
    legend_tot->SetFillStyle(0);
    legend_tot->SetBorderSize(0);

    for (size_t i = 0; i < interested_channels.size(); ++i) {
        int ch = interested_channels[i];
        if (ch < 0 || ch >= (int)h1d_tot_list.size()) {
            LOG(WARNING) << "Interested channel " << ch << " is out of range, skipping.";
            continue;
        }
        h1d_tot_list[ch]->SetLineColor(fpga_colors[i % fpga_colors.size()]);
        h1d_tot_list[ch]->SetLineWidth(2);
        // remove the statistics box
        h1d_tot_list[ch]->SetStats(0);
        if (i == 0) {
            gPad->SetLogy();
            // set axis label
            h1d_tot_list[ch]->GetXaxis()->SetTitle("TOT Value");
            h1d_tot_list[ch]->GetYaxis()->SetTitle("Counts");
            h1d_tot_list[ch]->Draw("HIST");
        } else {
            h1d_tot_list[ch]->Draw("HIST SAME");
        }
        legend_tot->AddEntry(h1d_tot_list[ch], ("Channel " + std::to_string(ch)).c_str(), "l");
    }
    legend_tot->Draw();
    canvas_tot_overlay->Modified();
    canvas_tot_overlay->Update();
    canvas_tot_overlay->Write();
    std::string canvas_tot_overlay_png = opts.output_folder + "/tot_overlay_" + run_info_str + ".png";
    canvas_tot_overlay->SaveAs(canvas_tot_overlay_png.c_str());
    canvas_tot_overlay->Close();

    // draw the average ADC, ToA, and TOT distribution as well
    TCanvas *canvas_average = new TCanvas("canvas_average", "Average ADC, ToA, and TOT Distributions", 1200, 800);
    canvas_average->Divide(1, 3);

    h1d_adc_average->SetLineColor(kRed);
    h1d_adc_average->SetLineWidth(2);
    h1d_adc_average->SetStats(0);

    h1d_toa_average->SetLineColor(kBlue);
    h1d_toa_average->SetLineWidth(2);
    h1d_toa_average->SetStats(0);

    h1d_tot_average->SetLineColor(kGreen+2);
    h1d_tot_average->SetLineWidth(2);
    h1d_tot_average->SetStats(0);

    canvas_average->cd(1);
    gPad->SetLogy();
    h1d_adc_average->GetXaxis()->SetTitle("ADC Value");
    h1d_adc_average->GetYaxis()->SetTitle("Counts");
    h1d_adc_average->Draw("HIST");

    canvas_average->cd(2);
    gPad->SetLogy();
    h1d_toa_average->GetXaxis()->SetTitle("TOA Value (ns)");
    h1d_toa_average->GetYaxis()->SetTitle("Counts");
    h1d_toa_average->Draw("HIST");

    canvas_average->cd(3);
    gPad->SetLogy();
    h1d_tot_average->GetXaxis()->SetTitle("TOT Value");
    h1d_tot_average->GetYaxis()->SetTitle("Counts");
    h1d_tot_average->Draw("HIST");

    canvas_average->Modified();
    canvas_average->Update();
    canvas_average->Write();
    std::string canvas_average_png = opts.output_folder + "/average_distributions_" + run_info_str + ".png";
    canvas_average->SaveAs(canvas_average_png.c_str());
    canvas_average->Close();

    output_root->cd();

    output_root->Close();

    return 0;
}