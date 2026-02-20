#include "H2GCROC_Common.hxx"
#include "H2GCROC_Lib.hxx"

#include "TH2D.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TString.h"
#include "TObject.h"
#include <limits>

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
    

    std::vector<int> interested_channels = {50, 52};

    std::vector<int> interested_but_covered_channels = {54, 58};

    const int pedestal_index_max = 1; // think the pedestal is stable, and the first two samples are enough to calculate the pedestal
    const int peak_index_min = 5;
    const int peak_index_max = 9; // the sample index of the signal peak should be within this range
    const int first_toa_index_min = 2;
    const int first_toa_index_max = 5; // the sample index of the first ToA should be within this range

    TFile *output_root = new TFile(opts.output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << opts.output_file;
        return 1;
    }

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

    Long64_t processed_entries = 0;
    Long64_t hamming_error_entries = 0;
    bool skip_hamming_error_entries = true;

    const double adc_hist_min = 0.0;
    const double adc_hist_max = 1024.0;
    const int adc_hist_bins = 1024;

    const double toa_peak_window_min = 115.0; // unit: ns
    const double toa_peak_window_max = 125.0; // unit: ns

    const double toa_ns_hist_min = 0.0;
    const double toa_ns_hist_max = machine_gun_samples * sample_time;
    const int toa_ns_hist_bins = machine_gun_samples * 16; // 16 bins per sample

    const double toa_shifted_ns_hist_min = -sample_time;
    const double toa_shifted_ns_hist_max = machine_gun_samples * sample_time;
    const int toa_shifted_ns_hist_bins = (machine_gun_samples + 1) * 16; // 16 bins per sample, and add one sample before the first sample to cover the negative shifted ToA

    const double sample_hist_min = 0.0;
    const double sample_hist_max = machine_gun_samples;
    const int sample_hist_bins = machine_gun_samples;

    // th1d for pedestal distribution for each channel
    std::vector<TH1D*> h1d_adc_pedestal_channel_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::string hist_name = "h1d_adc_pedestal_channel_" + std::to_string(_chn);
        auto *h1d_adc_pedestal_channel = new TH1D(
            hist_name.c_str(),
            ("ADC Pedestal for Channel " + std::to_string(_chn) + ";ADC;Count").c_str(),
            adc_hist_bins/4,
            adc_hist_min,
            adc_hist_max/4
        );
        h1d_adc_pedestal_channel->SetDirectory(nullptr);
        h1d_adc_pedestal_channel_list.push_back(h1d_adc_pedestal_channel);
    }

    // th2d for adc samples for each channel
    std::vector<TH2D*> h2d_adc_samples_channel_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::string hist_name = "h2d_adc_samples_channel_" + std::to_string(_chn);
        auto *h2d_adc_samples_channel = new TH2D(
            hist_name.c_str(),
            ("ADC Samples for Channel " + std::to_string(_chn) + ";Sample Index;ADC").c_str(),
            sample_hist_bins,
            sample_hist_min,
            sample_hist_max,
            adc_hist_bins,
            adc_hist_min,
            adc_hist_max
        );
        h2d_adc_samples_channel->SetDirectory(nullptr);
        h2d_adc_samples_channel_list.push_back(h2d_adc_samples_channel);
    }

    // th1d for adc peak index distribution for each channel
    std::vector<TH1D*> h1d_adc_peak_index_channel_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::string hist_name = "h1d_adc_peak_index_channel_" + std::to_string(_chn);
        auto *h1d_adc_peak_index_channel = new TH1D(
            hist_name.c_str(),
            ("ADC Peak Sample Index for Channel " + std::to_string(_chn) + ";Sample Index;Count").c_str(),
            sample_hist_bins,
            sample_hist_min,
            sample_hist_max
        );
        h1d_adc_peak_index_channel->SetDirectory(nullptr);
        h1d_adc_peak_index_channel_list.push_back(h1d_adc_peak_index_channel);
    }

    // th1d for limited adc peak (pedestal subtracted) distribution for each channel
    std::vector<TH1D*> h1d_adc_peak_limited_minus_pedestal_channel_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::string hist_name = "h1d_adc_peak_limited_minus_pedestal_channel_" + std::to_string(_chn);
        auto *h1d_adc_peak_limited_minus_pedestal_channel = new TH1D(
            hist_name.c_str(),
            ("ADC Peak (Limited to Index " + std::to_string(peak_index_max) + ") Minus Pedestal for Channel " + std::to_string(_chn) + ";ADC Peak - Pedestal;Count").c_str(),
            adc_hist_bins,
            adc_hist_min,
            adc_hist_max
        );
        h1d_adc_peak_limited_minus_pedestal_channel->SetDirectory(nullptr);
        h1d_adc_peak_limited_minus_pedestal_channel_list.push_back(h1d_adc_peak_limited_minus_pedestal_channel);
    }

    // th1d for limited, ToA filtered adc peak (pedestal subtracted) distribution for each channel
    std::vector<TH1D*> h1d_adc_peak_limited_toa_filtered_minus_pedestal_channel_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::string hist_name = "h1d_adc_peak_limited_toa_filtered_minus_pedestal_channel_" + std::to_string(_chn);
        auto *h1d_adc_peak_limited_toa_filtered_minus_pedestal_channel = new TH1D(
            hist_name.c_str(),
            ("ADC Peak (Limited to Index " + std::to_string(peak_index_max) + ") with ToA in [" + std::to_string(toa_peak_window_min) + ", " + std::to_string(toa_peak_window_max) + "] ns Minus Pedestal for Channel " + std::to_string(_chn) + ";ADC Peak - Pedestal;Count").c_str(),
            adc_hist_bins,
            adc_hist_min,
            adc_hist_max
        );
        h1d_adc_peak_limited_toa_filtered_minus_pedestal_channel->SetDirectory(nullptr);
        h1d_adc_peak_limited_toa_filtered_minus_pedestal_channel_list.push_back(h1d_adc_peak_limited_toa_filtered_minus_pedestal_channel);
    }

    // th2d for limited adc peak v.s. toa correlation for each channel
    std::vector<TH2D*> h2d_adc_peak_limited_v_toa_channel_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::string hist_name = "h2d_adc_peak_limited_v_toa_channel_" + std::to_string(_chn);
        auto *h2d_adc_peak_limited_v_toa_channel = new TH2D(
            hist_name.c_str(),
            ("Correlation between ADC Peak (Limited to Index " + std::to_string(peak_index_max) + ") and ToA for Channel " + std::to_string(_chn) + ";ToA (ns);ADC Peak").c_str(),
            toa_ns_hist_bins,
            toa_ns_hist_min,
            toa_ns_hist_max,
            adc_hist_bins,
            adc_hist_min,
            adc_hist_max
        );
        h2d_adc_peak_limited_v_toa_channel->SetDirectory(nullptr);
        h2d_adc_peak_limited_v_toa_channel_list.push_back(h2d_adc_peak_limited_v_toa_channel);
    } 

    // th2d for toa shifted adc waveform for each channel
    std::vector<TH2D*> h2d_toa_shifted_adc_waveform_channel_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::string hist_name = "h2d_toa_shifted_adc_waveform_channel_" + std::to_string(_chn);
        auto *h2d_toa_shifted_adc_waveform_channel = new TH2D(
            hist_name.c_str(),
            ("ToA Shifted ADC Waveform for Channel " + std::to_string(_chn) + ";Sample Index Shifted by ToA;ADC").c_str(),
            toa_shifted_ns_hist_bins,
            toa_shifted_ns_hist_min,
            toa_shifted_ns_hist_max,
            adc_hist_bins,
            adc_hist_min,
            adc_hist_max
        );
        h2d_toa_shifted_adc_waveform_channel->SetDirectory(nullptr);
        h2d_toa_shifted_adc_waveform_channel_list.push_back(h2d_toa_shifted_adc_waveform_channel);
    }

    std::vector<double> channel_valid_toa_count_list(FPGA_CHANNEL_NUMBER * fpga_count, 0.0); // to count the number of events with valid ToA for each channel, which will decide whether to use toa for max searching
    std::vector<std::vector<double>> channel_adc_peak_sample_values_list(FPGA_CHANNEL_NUMBER * fpga_count);
    for (auto& vec : channel_adc_peak_sample_values_list) {
        vec.reserve(entry_max);
    }
    std::vector<std::vector<double>> channel_pedestal_list(FPGA_CHANNEL_NUMBER * fpga_count);
    for (auto& vec : channel_pedestal_list) {
        vec.reserve(entry_max);
    }

    for (int _entry = 0; _entry < entry_max; _entry++) {
        input_tree->GetEntry(_entry);
        processed_entries++;
        if (_entry % 5000 == 0) {
            LOG(INFO) << "Processing entry " << _entry << " / " << entry_max;
        }

        bool _hamming_code_passed = true;

        for (int _fpga_index = 0; _fpga_index < fpga_count; _fpga_index++) {
            auto _fpga_id    = legal_fpga_id_list[_fpga_index];
            auto _timestamp  = branch_timestamps_list[_fpga_index][0];
            auto _daqh_list  = branch_daqh_list_list[_fpga_index];
            auto _tc_list    = branch_tc_list_list[_fpga_index];
            auto _tp_list    = branch_tp_list_list[_fpga_index];
            auto _val0_list  = branch_val0_list_list[_fpga_index];
            auto _val1_list  = branch_val1_list_list[_fpga_index];
            auto _val2_list  = branch_val2_list_list[_fpga_index];

            for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                for (int _daqh_index = 0; _daqh_index < 4; _daqh_index++){
                    auto _daqh = _daqh_list[_sample_index * 4 + _daqh_index];
                    auto _h1h2h3 = (_daqh >> 4) & 0x7;
                    if (_h1h2h3 != 0x00){
                        _hamming_code_passed = false;
                    }
                }
                if (skip_hamming_error_entries && !_hamming_code_passed) {
                    break;
                }
            } // sample loop

            if (skip_hamming_error_entries && !_hamming_code_passed) {
                break;
            }

            for (int _channel_index = 0; _channel_index < FPGA_CHANNEL_NUMBER; _channel_index++) {
                auto _channel_valid = get_valid_fpga_channel(_channel_index);
                if (_channel_valid == -1){ // if it is CM or Calib channel, skip it
                    continue;
                }
                // ! --- values for one event in one channel --------------------------------------
                double _adc_pedestal = 0.0;
                std::vector<int> _adc_samples;
                _adc_samples.reserve(machine_gun_samples);
                std::vector<double> _adc_pedestal_samples;
                _adc_pedestal_samples.reserve(pedestal_index_max + 1);

                double _adc_peak = 0.0;
                std::vector<int> _adc_peak_sample_index_list; // in case there are multiple samples with the same peak value
                double _adc_peak_limited = 0.0; // the peak value limited to the samples within the expected peak index range [peak_index_min, peak_index_max]

                int _toa_first = 0;
                int _toa_first_sample_index = 0;

                for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                    auto _adc_raw = static_cast<int>(_val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER] & 0x3FFu);
                    auto _tot_raw = static_cast<int>(_val1_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER] & 0x3FFu);
                    auto _toa_raw = static_cast<int>(_val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER] & 0x3FFu);

                    _adc_samples.push_back(_adc_raw);
                    if (_sample_index <= pedestal_index_max) {
                        _adc_pedestal_samples.push_back(static_cast<double>(_adc_raw));
                    }

                    if (_toa_raw != 0 && _toa_first == 0 && _sample_index >= first_toa_index_min && _sample_index <= first_toa_index_max) {
                        _toa_first = _toa_raw;
                        _toa_first_sample_index = _sample_index;
                    }
                } // sample loop
                _adc_pedestal = std::accumulate(_adc_pedestal_samples.begin(), _adc_pedestal_samples.end(), 0.0) / _adc_pedestal_samples.size();
                h1d_adc_pedestal_channel_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index]->Fill(_adc_pedestal);
                for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                    h2d_adc_samples_channel_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index]->Fill(_sample_index, _adc_samples[_sample_index]);
                }
                channel_pedestal_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index].push_back(_adc_pedestal);

                _adc_peak = *std::max_element(_adc_samples.begin(), _adc_samples.end());
                for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                    if (_adc_samples[_sample_index] == _adc_peak) {
                        _adc_peak_sample_index_list.push_back(_sample_index);
                        h1d_adc_peak_index_channel_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index]->Fill(_sample_index);
                    }
                }

                _adc_peak_limited = 0.0;
                for (int _sample_index = peak_index_min; _sample_index <= peak_index_max; _sample_index++) {
                    if (_adc_samples[_sample_index] > _adc_peak_limited) {
                        _adc_peak_limited = _adc_samples[_sample_index];
                    }
                }
                // double _adc_peak_limited_minus_pedestal = _adc_peak_limited - _adc_pedestal;
                double _adc_peak_limited_minus_pedestal = _adc_peak_limited;
                h1d_adc_peak_limited_minus_pedestal_channel_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index]->Fill(_adc_peak_limited_minus_pedestal);
                channel_adc_peak_sample_values_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index].push_back(_adc_peak_limited);

                if (_toa_first != 0) {
                    double _toa_first_ns = static_cast<double>(_toa_first) * 0.025 + static_cast<double>(_toa_first_sample_index) * sample_time; // convert ToA to ns, and add the time of the first sample index
                    if (_toa_first >= 268){
                        _toa_first_ns -= 25.0;
                    }
                    h2d_adc_peak_limited_v_toa_channel_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index]->Fill(_toa_first_ns, _adc_peak_limited);
                    for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                        double _sample_index_shifted_by_toa = static_cast<double>(_sample_index*sample_time) - _toa_first_ns+75.0; // shift the sample index by the ToA, and add 75 ns to make the shifted sample index positive
                        h2d_toa_shifted_adc_waveform_channel_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index]->Fill(_sample_index_shifted_by_toa, _adc_samples[_sample_index]);
                    }
                    if (_toa_first_ns >= toa_peak_window_min && _toa_first_ns <= toa_peak_window_max) {
                        double _adc_peak_limited_toa_filtered_minus_pedestal = _adc_peak_limited_minus_pedestal;
                        h1d_adc_peak_limited_toa_filtered_minus_pedestal_channel_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index]->Fill(_adc_peak_limited_toa_filtered_minus_pedestal);
                    }
                    channel_valid_toa_count_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index] += 1.0;
                }
            } // channel loop
        } // fpga loop
        if (skip_hamming_error_entries && !_hamming_code_passed) {
            hamming_error_entries++;
            continue;
        }
    } // event loop

    LOG(INFO) << "Processed entries: " << processed_entries << " (" << (float)processed_entries / entry_max * 100 << "%)";
    LOG(INFO) << "Hamming code error: " << hamming_error_entries << " (" << (float)hamming_error_entries / processed_entries * 100 << "%)";
    // print the valid toa ratio for example channels
    for (int _channel : interested_channels) {
        double valid_toa_ratio = channel_valid_toa_count_list[_channel] / processed_entries;
        LOG(INFO) << "Channel " << _channel << " valid ToA ratio: " << valid_toa_ratio * 100 << "%";
    }

    input_root->Close();

    // Initialize mosaic topology for visualization
    std::string mapping_json_file = "config/mapping_Feb2026_re.json";
    MosaicTopoSetup mosaic_setup = initialize_mosaic_topology(fpga_count, mapping_json_file, FPGA_CHANNEL_NUMBER);

    // calculate the mean pedestal for each channel
    std::vector<double> channel_pedestal_mean_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        double pedestal_mean = std::accumulate(channel_pedestal_list[_chn].begin(), channel_pedestal_list[_chn].end(), 0.0) / channel_pedestal_list[_chn].size();
        channel_pedestal_mean_list.push_back(pedestal_mean);
    }

    output_root->cd();

    auto canvas_adc_pedestal = new TCanvas("canvas_adc_pedestal", "ADC Pedestal Distribution;Channel;ADC", 1200, 800);
    draw_mosaic_fixed(*canvas_adc_pedestal, h1d_adc_pedestal_channel_list, mosaic_setup.topo_ped_median);
    canvas_adc_pedestal->Modified();
    canvas_adc_pedestal->Update();
    canvas_adc_pedestal->Write();

    auto canvas_adc_samples = new TCanvas("canvas_adc_samples", "ADC Samples;Channel;ADC", 1200, 800);
    draw_mosaic_fixed(*canvas_adc_samples, h2d_adc_samples_channel_list, mosaic_setup.topo_ped_median);
    canvas_adc_samples->Modified();
    canvas_adc_samples->Update();
    canvas_adc_samples->Write();

    auto canvas_adc_peak_index = new TCanvas("canvas_adc_peak_index", "ADC Peak Sample Index Distribution;Channel;Sample Index", 1200, 800);
    draw_mosaic_fixed(*canvas_adc_peak_index, h1d_adc_peak_index_channel_list, mosaic_setup.topo_ped_median);
    canvas_adc_peak_index->Modified();
    canvas_adc_peak_index->Update();
    canvas_adc_peak_index->Write();

    auto canvas_adc_peak_limited_minus_pedestal = new TCanvas("canvas_adc_peak_limited_minus_pedestal", ("ADC Peak (Limited to Sample Index " + std::to_string(peak_index_max) + ") Minus Pedestal Distribution;Channel;ADC Peak - Pedestal").c_str(), 1200, 800);
    draw_mosaic_fixed(*canvas_adc_peak_limited_minus_pedestal, h1d_adc_peak_limited_minus_pedestal_channel_list, mosaic_setup.topo_ped_median);
    canvas_adc_peak_limited_minus_pedestal->Modified();
    canvas_adc_peak_limited_minus_pedestal->Update();
    canvas_adc_peak_limited_minus_pedestal->Write();

    auto canvas_adc_peak_limited_v_toa = new TCanvas("canvas_adc_peak_limited_v_toa", ("Correlation between ADC Peak (Limited to Sample Index " + std::to_string(peak_index_max) + ") and ToA;Channel;Correlation").c_str(), 1200, 800);
    draw_mosaic_fixed(*canvas_adc_peak_limited_v_toa, h2d_adc_peak_limited_v_toa_channel_list, mosaic_setup.topo_ped_median);
    canvas_adc_peak_limited_v_toa->Modified();
    canvas_adc_peak_limited_v_toa->Update();
    canvas_adc_peak_limited_v_toa->Write();

    auto canvas_toa_shifted_adc_waveform = new TCanvas("canvas_toa_shifted_adc_waveform", "ToA Shifted ADC Waveform;Channel;ToA Shifted ADC Waveform", 1200, 800);
    draw_mosaic_fixed(*canvas_toa_shifted_adc_waveform, h2d_toa_shifted_adc_waveform_channel_list, mosaic_setup.topo_ped_median);
    canvas_toa_shifted_adc_waveform->Modified();
    canvas_toa_shifted_adc_waveform->Update();
    canvas_toa_shifted_adc_waveform->Write();

    auto canvas_adc_peak_limited_toa_filtered_minus_pedestal = new TCanvas("canvas_adc_peak_limited_toa_filtered_minus_pedestal", ("ADC Peak (Limited to Sample Index " + std::to_string(peak_index_max) + ") with ToA in [" + std::to_string(toa_peak_window_min) + ", " + std::to_string(toa_peak_window_max) + "] ns Minus Pedestal;Channel;ADC Peak - Pedestal").c_str(), 1200, 800);
    draw_mosaic_fixed(*canvas_adc_peak_limited_toa_filtered_minus_pedestal, h1d_adc_peak_limited_toa_filtered_minus_pedestal_channel_list, mosaic_setup.topo_ped_median);
    canvas_adc_peak_limited_toa_filtered_minus_pedestal->Modified();
    canvas_adc_peak_limited_toa_filtered_minus_pedestal->Update();
    canvas_adc_peak_limited_toa_filtered_minus_pedestal->Write();

    const double saturation_threshold = 0.2; // if the last bin has more than 10% of counts, consider it as saturation
    const double saturation_ignore_threshold = 0.01; // if the last bin has more than 1% of counts, ignore it in fitting

    const double threshold_toa_ratio_valid = 0.2;
    const int sliding_x_window_bins = 6; // unit: 1.5625 ns
    const int sliding_x_window_step_bins = 1;
    TDirectory* dir_interested = output_root->mkdir("Interested_Channels");
    dir_interested->cd();
    for (int _channel : interested_channels) {
        if (_channel < 0 || _channel >= FPGA_CHANNEL_NUMBER * fpga_count) {
            LOG(WARNING) << "Interested but covered channel " << _channel << " is out of range, skipping";
            continue;
        }
        double valid_toa_ratio = channel_valid_toa_count_list[_channel] / processed_entries;
        if (valid_toa_ratio >= threshold_toa_ratio_valid) {
            LOG(INFO) << "Channel " << _channel << " has valid ToA ratio " << valid_toa_ratio * 100 << "%, which is above the threshold " << threshold_toa_ratio_valid * 100 << "%, so it will be included in the ToA filtered ADC peak analysis";
            auto& th2d_shifted_waveform = h2d_toa_shifted_adc_waveform_channel_list[_channel];
            auto sliding_result = FindBestXWindowByYMean(th2d_shifted_waveform, sliding_x_window_bins, sliding_x_window_step_bins);
            double toa_window_min = sliding_result.xmin;
            double toa_window_max = sliding_result.xmax;

            AutoZoomByQuantile(sliding_result.hY, 0.01, 0.99, 0.2);
            
            if (sliding_result.xbin1 != -1 && sliding_result.hY) {
                // save the sliding result histogram
                auto canvas_sliding = new TCanvas(("canvas_sliding_channel_" + std::to_string(_channel)).c_str(), ("Finding the Best ToA Window by Sliding for Channel " + std::to_string(_channel)).c_str(), 800, 600);
                canvas_sliding->SetLeftMargin(0.12);
                canvas_sliding->SetBottomMargin(0.12);
                canvas_sliding->SetRightMargin(0.05);
                canvas_sliding->SetTopMargin(0.08);
                
                TLegend* legend = new TLegend(0.6, 0.7, 0.89, 0.89);
                legend->SetFillStyle(0);
                legend->SetBorderSize(0);

                // Set histogram properties
                sliding_result.hY->SetStats(kFALSE);
                sliding_result.hY->SetTitle(("ADC Distribution in Best ToA Window for Channel " + std::to_string(_channel) + ";ADC;Count").c_str());
                sliding_result.hY->Draw("hist");
                
                // Set axis properties after drawing
                sliding_result.hY->GetXaxis()->SetTitleSize(0.045);
                sliding_result.hY->GetXaxis()->SetLabelSize(0.04);
                sliding_result.hY->GetYaxis()->SetTitleSize(0.045);
                sliding_result.hY->GetYaxis()->SetLabelSize(0.04);
                sliding_result.hY->GetXaxis()->SetNdivisions(505);
                sliding_result.hY->GetYaxis()->SetNdivisions(505);
                legend->AddEntry(sliding_result.hY, "Y Projection in Sliding X Window", "l");
                
                // Check for saturation: if last bin has > 10% of counts, skip fit
                int last_bin = sliding_result.hY->GetNbinsX();
                double last_bin_content = sliding_result.hY->GetBinContent(last_bin);
                double total_entries = sliding_result.hY->GetEntries();
                double last_bin_fraction = last_bin_content / total_entries;
                
                double fit_mean, fit_sigma;
                if (last_bin_fraction > saturation_threshold) {
                    // Saturation detected
                    fit_mean = 1023.0;
                    fit_sigma = 1.0 / sqrt(12.0); // assuming uniform distribution between 1023 and 1024 for the saturated bin
                    LOG(WARNING) << "Channel " << _channel << " shows saturation (last bin has " << last_bin_fraction * 100 << "% of counts). Assigning mean=1023, sigma=" << fit_sigma << " and skipping Gaussian fit.";
                    // Create dummy TF1 with fixed parameters and use Fit with N option (no actual fit, just store)
                    TF1* dummy_fit = new TF1("fit_func", "gaus", 0, 1023);
                    dummy_fit->AddToGlobalList(false);
                    dummy_fit->SetParameter(0, sliding_result.hY->GetMaximum()); // amplitude
                    dummy_fit->SetParameter(1, fit_mean); // mean = 1023
                    dummy_fit->SetParameter(2, fit_sigma); // sigma
                    dummy_fit->FixParameter(0, sliding_result.hY->GetMaximum());
                    dummy_fit->FixParameter(1, fit_mean);
                    dummy_fit->FixParameter(2, fit_sigma);
                    // Add to histogram's function list without fitting
                    sliding_result.hY->GetListOfFunctions()->Add(dummy_fit);
                    legend->AddEntry(dummy_fit, ("Saturation detected: mean = " + std::to_string(fit_mean) + " ADC, sigma = " + std::to_string(fit_sigma) + " ADC").c_str(), "");
                    // set the last bin to 0 count
                    sliding_result.hY->SetBinContent(last_bin, 0);
                } else {
                    // Two-stage Gaussian fit
                    // double hist_mean = sliding_result.hY->GetMean();
                    // rebin
                    sliding_result.hY->Rebin(4); // rebin by a factor of 4 to reduce statistical fluctuation for fitting
                    // force all the bins under 80 ADC to be 0
                    for (int bin = 1; bin <= sliding_result.hY->GetNbinsX(); ++bin) {
                        if (sliding_result.hY->GetXaxis()->GetBinCenter(bin) < 80) {
                            sliding_result.hY->SetBinContent(bin, 0);
                        }
                    }
                    int max_value_bin = sliding_result.hY->GetMaximumBin();
                    double hist_mean = sliding_result.hY->GetXaxis()->GetBinCenter(max_value_bin);
                    double hist_rms = sliding_result.hY->GetRMS();
                    double prefit_range_min = hist_mean - 4.0 * hist_rms;
                    double prefit_range_max = hist_mean + 4.0 * hist_rms;
                    if (prefit_range_min < 0) {
                        prefit_range_min = 0;
                    }
                    if (prefit_range_max > 1023) {
                        prefit_range_max = 1023;
                    }
                    
                    // If last bin has 1-10% of counts, exclude it from fit range
                    if (last_bin_fraction > saturation_ignore_threshold) {
                        double last_bin_low_edge = sliding_result.hY->GetXaxis()->GetBinLowEdge(last_bin);
                        if (prefit_range_max > last_bin_low_edge) {
                            prefit_range_max = last_bin_low_edge;
                            LOG(INFO) << "Channel " << _channel << " has " << last_bin_fraction * 100 << "% in last bin. Excluding last bin from fit (max range: " << prefit_range_max << ")";
                        }
                    }
                    
                    // Stage 1: Pre-fit to get rough parameters
                    TF1* prefit_func = new TF1("prefit_func", "gaus", prefit_range_min, prefit_range_max);
                    prefit_func->SetParameters(sliding_result.hY->GetMaximum(), hist_mean, hist_rms);
                    sliding_result.hY->Fit(prefit_func, "RQN");  // N option: don't store in histogram's function list
                    double prefit_mean = prefit_func->GetParameter(1);
                    double prefit_sigma = prefit_func->GetParameter(2);
                    delete prefit_func;
                    
                    // Stage 2: Final fit using pre-fit results
                    double fit_range_min = prefit_mean - 3.0 * prefit_sigma;
                    double fit_range_max = prefit_mean + 3.0 * prefit_sigma;
                    if (fit_range_min < 0) {
                        fit_range_min = 0;
                    }
                    if (fit_range_max > 1023) {
                        fit_range_max = 1023;
                    }
                    
                    TF1* fit_func = new TF1("fit_func", "gaus", fit_range_min, fit_range_max);
                    fit_func->SetParameters(sliding_result.hY->GetMaximum(), prefit_mean, prefit_sigma);
                    // limit the sigma parameter to be within [0.5*prefit_sigma, 2*prefit_sigma] to avoid unphysical fit results due to remaining statistical fluctuation
                    fit_func->SetParLimits(2, 0.5 * prefit_sigma, 2.0 * prefit_sigma);
                    sliding_result.hY->Fit(fit_func, "RQ");
                    fit_func->SetLineColorAlpha(kRed, 0.7);
                    fit_func->Draw("same");
                    fit_mean = fit_func->GetParameter(1);
                    fit_sigma = fit_func->GetParameter(2);
                    legend->AddEntry(fit_func, ("Gaussian Fit: mean = " + std::to_string(fit_mean) + " ADC, sigma = " + std::to_string(fit_sigma) + " ADC").c_str(), "l");
                }
                legend->Draw();
                
                canvas_sliding->Modified();
                canvas_sliding->Update();
                canvas_sliding->Write();
                // Don't close canvas explicitly - let ROOT handle cleanup when file closes
                // canvas_sliding->Close();
                
                LOG(INFO) << "For channel " << _channel << ", the best ToA window is [" << sliding_result.xmin << ", " << sliding_result.xmax << "] ns with mean ADC " << sliding_result.meanY << ", Gaussian fit mean = " << fit_mean << " ADC, sigma = " << fit_sigma << " ADC";

            } else {
                LOG(WARNING) << "Failed to find a valid ToA window for channel " << _channel;
            }


            LOG(INFO) << "For channel " << _channel << ", the best ToA window is [" << toa_window_min << ", " << toa_window_max << "] ns with mean ADC " << sliding_result.meanY;
            // draw the 2d histogram of the waveform and highlight the best ToA window
            auto canvas_waveform = new TCanvas(("canvas_waveform_channel_" + std::to_string(_channel)).c_str(), ("ToA Shifted ADC Waveform for Channel " + std::to_string(_channel)).c_str(), 1000, 600);
            canvas_waveform->SetLeftMargin(0.12);
            canvas_waveform->SetBottomMargin(0.12);
            canvas_waveform->SetRightMargin(0.15);
            canvas_waveform->SetTopMargin(0.08);
            
            th2d_shifted_waveform->SetStats(kFALSE);
            th2d_shifted_waveform->SetTitle(";Time [ns];ADC");
            // set logz
            canvas_waveform->SetLogz();
            // set axis ticks
            th2d_shifted_waveform->GetXaxis()->SetNdivisions(505);
            th2d_shifted_waveform->GetYaxis()->SetNdivisions(505);
            th2d_shifted_waveform->Draw("colz");

            // Enable and make ticks visible
            gPad->SetTicks(1, 1);
            
            // Explicitly set tick length so they're actually visible
            th2d_shifted_waveform->GetXaxis()->SetTickLength(0.03);
            th2d_shifted_waveform->GetYaxis()->SetTickLength(0.03);
            th2d_shifted_waveform->GetZaxis()->SetTickLength(0);  // Hide Z-axis ticks on right edge

            // Set divisions and sizes
            th2d_shifted_waveform->GetXaxis()->SetNdivisions(510);
            th2d_shifted_waveform->GetYaxis()->SetNdivisions(510);
            th2d_shifted_waveform->GetXaxis()->SetTitleSize(0.045);
            th2d_shifted_waveform->GetXaxis()->SetLabelSize(0.04);
            th2d_shifted_waveform->GetYaxis()->SetTitleSize(0.045);
            th2d_shifted_waveform->GetYaxis()->SetLabelSize(0.04);
            th2d_shifted_waveform->GetZaxis()->SetTitleSize(0.045);
            th2d_shifted_waveform->GetZaxis()->SetLabelSize(0.04);
            
            TLine* line_min = new TLine(toa_window_min, adc_hist_min, toa_window_min, adc_hist_max);
            line_min->SetLineColor(kRed+2);
            line_min->SetLineWidth(2);
            // dashed line style for the window boundary
            line_min->SetLineStyle(2);
            line_min->Draw("same");
            TLine* line_max = new TLine(toa_window_max, adc_hist_min, toa_window_max, adc_hist_max);
            line_max->SetLineColor(kRed+2);
            line_max->SetLineWidth(2);
            line_max->SetLineStyle(2);
            line_max->Draw("same");

            // write latex info
            const double x_text = 0.82;
            double y_start = 0.89;
            double y_step = 0.04;
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.04);
            latex.SetTextFont(62);
            latex.SetTextAlign(33);
            latex.DrawLatex(x_text, y_start, "Laser Test with H2GCROC");
            latex.SetTextSize(0.03);
            latex.SetTextFont(42);
            latex.DrawLatex(x_text, y_start - y_step, ("ADC Waveform Shifted by ToA for Channel " + std::to_string(_channel)).c_str());
            latex.DrawLatex(x_text, y_start - 2 * y_step, ("Best ToA Window: [" + std::to_string(std::round(toa_window_min * 1000) / 1000).substr(0, 5) + ", " + std::to_string(std::round(toa_window_max * 1000) / 1000).substr(0, 5) + "] ns").c_str());
            latex.DrawLatex(x_text, y_start - 3 * y_step, "CERN, February 2026");
            
            canvas_waveform->Modified();
            canvas_waveform->Update();
            canvas_waveform->Write();
            // save as a separte pdf file
            canvas_waveform->SaveAs((opts.output_file + "_channel_" + std::to_string(_channel) + "_toa_shifted_waveform.pdf").c_str());
            // Don't close canvas explicitly - let ROOT handle cleanup when file closes
            // canvas_waveform->Close();

        } else {
            LOG(INFO) << "Channel " << _channel << " has valid ToA ratio " << valid_toa_ratio * 100 << "%, which is below the threshold " << threshold_toa_ratio_valid * 100 << "%, so it will NOT be included in the ToA filtered ADC peak analysis"; 
            auto& channel_adc_peak_list = channel_adc_peak_sample_values_list[_channel];
            std::sort(channel_adc_peak_list.begin(), channel_adc_peak_list.end());
            auto canvas_peak = new TCanvas(("canvas_peak_channel_" + std::to_string(_channel)).c_str(), ("ADC Peak (Limited to Sample Index " + std::to_string(peak_index_max) + ") Distribution for Channel " + std::to_string(_channel)).c_str(), 800, 600);
            TLegend* legend_sample = new TLegend(0.6, 0.7, 0.89, 0.89);
            legend_sample->SetFillStyle(0);
            legend_sample->SetBorderSize(0);
            auto* h1d_peak = new TH1D(("h1d_peak_channel_" + std::to_string(_channel)).c_str(), ("ADC Peak (Limited to Sample Index " + std::to_string(peak_index_max) + ") Distribution for Channel " + std::to_string(_channel) + ";ADC Peak;Count").c_str(), adc_hist_bins, adc_hist_min, adc_hist_max);
            int start_index = static_cast<int>(channel_adc_peak_list.size() * 0.9);
            for (size_t i = start_index; i < channel_adc_peak_list.size(); i++) {
                h1d_peak->Fill(channel_adc_peak_list[i]);
            }
            h1d_peak->SetStats(kFALSE);
            h1d_peak->Draw("hist");
            legend_sample->AddEntry(h1d_peak, ("ADC Peak (Limited to Sample Index " + std::to_string(peak_index_max) + ")").c_str(), "l");
            
            // Check for saturation: if last bin has > 10% of counts, skip fit
            int last_bin = h1d_peak->GetNbinsX();
            double last_bin_content = h1d_peak->GetBinContent(last_bin);
            double total_entries = h1d_peak->GetEntries();
            double last_bin_fraction = last_bin_content / total_entries;
            LOG(INFO) << "Channel " << _channel << ": last bin content = " << last_bin_content << ", total entries = " << total_entries << ", last bin fraction = " << last_bin_fraction * 100 << "%";
            
            double fit_mean, fit_sigma;
            if (last_bin_fraction > saturation_threshold) {
                // Saturation detected
                fit_mean = 1023.0;
                fit_sigma = 1.0 / sqrt(12.0); // assuming uniform distribution between 1023 and 1024 for the saturated bin
                LOG(WARNING) << "Channel " << _channel << " shows saturation (last bin has " << last_bin_fraction * 100 << "% of counts). Assigning mean=1023, sigma=" << fit_sigma << " and skipping Gaussian fit. ";
                // Create dummy TF1 with fixed parameters and use Fit with N option (no actual fit, just store)
                TF1* dummy_fit = new TF1("fit_func", "gaus", 0, 1023);
                dummy_fit->AddToGlobalList(false);
                dummy_fit->SetParameter(0, h1d_peak->GetMaximum()); // amplitude
                dummy_fit->SetParameter(1, fit_mean); // mean = 1023
                dummy_fit->SetParameter(2, fit_sigma); // sigma
                dummy_fit->FixParameter(0, h1d_peak->GetMaximum());
                dummy_fit->FixParameter(1, fit_mean);
                dummy_fit->FixParameter(2, fit_sigma);
                // Add to histogram's function list without fitting
                h1d_peak->GetListOfFunctions()->Add(dummy_fit);
                legend_sample->AddEntry(dummy_fit, ("Saturation detected: mean = " + std::to_string(fit_mean) + " ADC, sigma = " + std::to_string(fit_sigma) + " ADC").c_str(), "l");
                // set the last bin to 0 count
                h1d_peak->SetBinContent(last_bin, 0);
            } else {
                // Two-stage Gaussian fit
                // double hist_mean = h1d_peak->GetMean();
                // rebin
                h1d_peak->Rebin(4);
                // force all the bins under 80 ADC to be 0
                for (int bin = 1; bin <= h1d_peak->GetNbinsX(); ++bin) {
                    if (h1d_peak->GetXaxis()->GetBinCenter(bin) < 80) {
                        h1d_peak->SetBinContent(bin, 0);
                    }
                }
                int max_value_bin = h1d_peak->GetMaximumBin();
                double hist_mean = h1d_peak->GetXaxis()->GetBinCenter(max_value_bin);
                double hist_rms = h1d_peak->GetRMS();
                double prefit_range_min = hist_mean - 4.0 * hist_rms;
                double prefit_range_max = hist_mean + 4.0 * hist_rms;
                if (prefit_range_min < 0) {
                    prefit_range_min = 0;
                }
                if (prefit_range_max > 1023) {
                    prefit_range_max = 1023;
                }
                
                // If last bin has 1-10% of counts, exclude it from fit range
                if (last_bin_fraction > saturation_ignore_threshold) {
                    double last_bin_low_edge = h1d_peak->GetXaxis()->GetBinLowEdge(last_bin);
                    if (prefit_range_max > last_bin_low_edge) {
                        prefit_range_max = last_bin_low_edge;
                        LOG(INFO) << "Channel " << _channel << " has " << last_bin_fraction * 100 << "% in last bin. Excluding last bin from fit (max range: " << prefit_range_max << ")";
                    }
                }
                
                // Stage 1: Pre-fit to get rough parameters
                TF1* prefit_func = new TF1("prefit_func", "gaus", prefit_range_min, prefit_range_max);
                prefit_func->SetParameters(h1d_peak->GetMaximum(), hist_mean, hist_rms);
                // limit the sigma parameter to be within [0.5*hist_rms, 2*hist_rms] to avoid unphysical pre-fit results due to statistical fluctuation
                prefit_func->SetParLimits(2, 0.5 * hist_rms, 2.0 * hist_rms);
                h1d_peak->Fit(prefit_func, "RQN");  // N option: don't store in histogram's function list
                double prefit_mean = prefit_func->GetParameter(1);
                double prefit_sigma = prefit_func->GetParameter(2);
                delete prefit_func;
                
                // Stage 2: Final fit using pre-fit results
                double fit_range_min = prefit_mean - 3.0 * prefit_sigma;
                double fit_range_max = prefit_mean + 3.0 * prefit_sigma;
                if (fit_range_min < 0) {
                    fit_range_min = 0;
                }
                if (fit_range_max > 1023) {
                    fit_range_max = 1023;
                }
                
                TF1* fit_func = new TF1("fit_func", "gaus", fit_range_min, fit_range_max);
                fit_func->SetParameters(h1d_peak->GetMaximum(), prefit_mean, prefit_sigma);
                h1d_peak->Fit(fit_func, "RQ");
                fit_func->SetLineColor(kRed);
                fit_func->Draw("same");
                fit_mean = fit_func->GetParameter(1);
                fit_sigma = fit_func->GetParameter(2);
                legend_sample->AddEntry(fit_func, ("Gaussian Fit: mean = " + std::to_string(fit_mean) + " ADC, sigma = " + std::to_string(fit_sigma) + " ADC").c_str(), "l");
            }
            legend_sample->Draw();
            canvas_peak->Modified();
            canvas_peak->Update();
            canvas_peak->Write();
        }

    }

    output_root->Close();

    return 0;
}