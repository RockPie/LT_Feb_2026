#include "H2GCROC_Common.hxx"
#include "H2GCROC_Lib.hxx"
#include "H2GCROC_ToT.hxx"

#include "TH2D.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TString.h"
#include "TObject.h"
#include <limits>

INITIALIZE_EASYLOGGINGPP

double adc_tot_combine(double adc_raw, double adc_pedestal, double tot, const TotToAdcLUT& lut) {
    // If ToT is zero, return the ADC value directly
    if (tot == 0) {
        return adc_raw - adc_pedestal;
    }
    if (adc_raw < 1023){
        return adc_raw - adc_pedestal;
    }
    double adc_interp = lut.Eval(tot);
    if (adc_interp <= adc_raw - adc_pedestal) {
        // If the LUT returns zero or negative, fallback to raw ADC
        return adc_raw - adc_pedestal;
    }

    return adc_interp;
}

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

    std::string example_tot_lut_file = "dump/405_ToT_Scan/ToTScan5.root_LUT_Channel_50.txt";
    TotToAdcLUT lut = LoadTotToAdcLUT(example_tot_lut_file);

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

    const double tot_hist_min = 0.0;
    const double tot_hist_max = 4096.0;
    const int tot_hist_bins = 512;

    const double sample_hist_min = 0.0;
    const double sample_hist_max = machine_gun_samples;
    const int sample_hist_bins = machine_gun_samples;

    std::vector<TH1D*> h1d_tot_channel_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::string hist_name = "h1d_tot_channel_" + std::to_string(_chn);
        auto* h1d_tot_channel = new TH1D(
            hist_name.c_str(),
            ("ToT distribution for channel " + std::to_string(_chn)).c_str(),
            tot_hist_bins,
            tot_hist_min,
            tot_hist_max
        );
        h1d_tot_channel->SetDirectory(nullptr);
        h1d_tot_channel_list.push_back(h1d_tot_channel);
    }

    std::vector<TH1D*> h1d_tot_sample_index_channel_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::string hist_name = "h1d_tot_sample_index_channel_" + std::to_string(_chn);
        auto* h1d_tot_sample_index_channel = new TH1D(
            hist_name.c_str(),
            ("ToT sample index distribution for channel " + std::to_string(_chn)).c_str(),
            sample_hist_bins,
            sample_hist_min,
            sample_hist_max
        );
        h1d_tot_sample_index_channel->SetDirectory(nullptr);
        h1d_tot_sample_index_channel_list.push_back(h1d_tot_sample_index_channel);
    }

    std::vector<TH2D*> h2d_tot_vs_toa_channel_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::string hist_name = "h2d_tot_vs_toa_channel_" + std::to_string(_chn);
        auto* h2d_tot_vs_toa_channel = new TH2D(
            hist_name.c_str(),
            ("ToT vs ToA for channel " + std::to_string(_chn)).c_str(),
            tot_hist_bins,
            tot_hist_min,
            tot_hist_max,
            toa_ns_hist_bins,
            toa_ns_hist_min,
            toa_ns_hist_max
        );
        h2d_tot_vs_toa_channel->SetDirectory(nullptr);
        h2d_tot_vs_toa_channel_list.push_back(h2d_tot_vs_toa_channel);
    }
    
    std::vector<double> channel_valid_tot_count_list(FPGA_CHANNEL_NUMBER * fpga_count, 0.0); // to count the number of events with valid ToT for each channel, which will decide whether to use tot for max searching
    std::vector<std::vector<double>> channel_adc_peak_sample_values_list(FPGA_CHANNEL_NUMBER * fpga_count);
    for (auto& vec : channel_adc_peak_sample_values_list) {
        vec.reserve(entry_max);
    }
    std::vector<std::vector<double>> channel_pedestal_list(FPGA_CHANNEL_NUMBER * fpga_count);
    for (auto& vec : channel_pedestal_list) {
        vec.reserve(entry_max);
    }

    std::vector<std::vector<double>> channel_combined_peak_sample_values_list(FPGA_CHANNEL_NUMBER * fpga_count);
    for (auto& vec : channel_combined_peak_sample_values_list) {
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

                double _tot_first = 0.0;
                int _tot_first_sample_index = 0;

                double _toa_first = 0.0;
                int _toa_first_sample_index = 0;

                for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                    auto _adc_raw = static_cast<int>(_val0_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER] & 0x3FFu);
                    auto _tot_raw = static_cast<int>(_val1_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER] & 0x3FFu);
                    auto _toa_raw = static_cast<int>(_val2_list[_channel_index + _sample_index * FPGA_CHANNEL_NUMBER] & 0x3FFu);

                    _adc_samples.push_back(_adc_raw);
                    if (_sample_index <= pedestal_index_max) {
                        _adc_pedestal_samples.push_back(static_cast<double>(_adc_raw));
                    }

                    if (_tot_raw > 0 && _tot_first == 0.0) {
                        _tot_first = static_cast<double>(_tot_raw);
                        _tot_first_sample_index = _sample_index;
                    }

                    if (_toa_raw > 0 && _toa_first == 0.0) {
                        _toa_first = static_cast<double>(_toa_raw);
                        _toa_first_sample_index = _sample_index;
                    }

                } // sample loop
                _adc_pedestal = std::accumulate(_adc_pedestal_samples.begin(), _adc_pedestal_samples.end(), 0.0) / _adc_pedestal_samples.size();

                channel_pedestal_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index].push_back(_adc_pedestal);

                _adc_peak = *std::max_element(_adc_samples.begin(), _adc_samples.end());
                for (int _sample_index = 0; _sample_index < machine_gun_samples; _sample_index++) {
                    if (_adc_samples[_sample_index] == _adc_peak) {
                        _adc_peak_sample_index_list.push_back(_sample_index);
                    }
                }

                _adc_peak_limited = 0.0;
                for (int _sample_index = peak_index_min; _sample_index <= peak_index_max; _sample_index++) {
                    if (_adc_samples[_sample_index] > _adc_peak_limited) {
                        _adc_peak_limited = _adc_samples[_sample_index];
                    }
                }

                channel_adc_peak_sample_values_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index].push_back(_adc_peak_limited);

                if (_tot_first > 0.0) {
                    double _tot_decoded = _tot_first;
                    if (_tot_decoded > 512.0) {
                        _tot_decoded -= 512.0;
                        _tot_decoded = _tot_decoded * 8.0; // ToT extension for large ToT values
                    }
                    h1d_tot_channel_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index]->Fill(_tot_decoded);
                    h1d_tot_sample_index_channel_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index]->Fill(_tot_first_sample_index);
                    channel_valid_tot_count_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index] += 1.0;
                    if (_toa_first > 0.0) {
                        double _toa_ns = _toa_first * 0.025 + static_cast<double>(_toa_first_sample_index) * sample_time;
                        if (_toa_first >= 268.0) {
                            _toa_ns -= 25.0; // ToA extension for large ToA values
                        }
                        h2d_tot_vs_toa_channel_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index]->Fill(_tot_decoded, _toa_ns);
                    }
                }

                double _tot_decoded_for_combination = _tot_first;
                if (_tot_decoded_for_combination > 512.0) {
                    _tot_decoded_for_combination -= 512.0;
                    _tot_decoded_for_combination = _tot_decoded_for_combination * 8.0; // ToT extension for large ToT values
                }
                double _combined_value = adc_tot_combine(_adc_peak, _adc_pedestal, _tot_decoded_for_combination, lut);
                channel_combined_peak_sample_values_list[_fpga_index * FPGA_CHANNEL_NUMBER + _channel_index].push_back(_combined_value);
                // if is the first 100 events, print the details for channel 50
                if (processed_entries <= 100 && (_channel_index == 50 || _channel_index == 52)) {
                    LOG(INFO) << "Entry " << _entry << " FPGA "<< _fpga_id << " Channel " << _channel_index << " - ADC samples: " << _adc_samples[0] << ", " << _adc_samples[1] << ", " << _adc_samples[2] << ", " << _adc_samples[3] << ", " << _adc_samples[4] << ", ... ToT first: " << _tot_decoded_for_combination << " ToA first: " << _toa_first * 0.025 + static_cast<double>(_toa_first_sample_index) * sample_time << " ns, Combined value: " << _combined_value;
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
    // print the tot valid ratio for each example channel
    for (int _chn : interested_channels) {
        double tot_valid_ratio = channel_valid_tot_count_list[_chn] / processed_entries;
        LOG(INFO) << "Channel " << _chn << " ToT valid ratio: " << tot_valid_ratio * 100 << "%";
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

    // create combined value histograms from the bin range defined by the 5% min and 95% max of the combined values, to avoid the influence of outliers
    double combined_hist_min = std::numeric_limits<double>::max();
    double combined_hist_max = std::numeric_limits<double>::min();
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::vector<double>& combined_values = channel_combined_peak_sample_values_list[_chn];
        std::sort(combined_values.begin(), combined_values.end());
        double hist_min = combined_values[static_cast<int>(combined_values.size() * 0.05)];
        double hist_max = combined_values[static_cast<int>(combined_values.size() * 0.95)];
        if (hist_min < combined_hist_min) {
            combined_hist_min = hist_min;
        }
        if (hist_max > combined_hist_max) {
            combined_hist_max = hist_max;
        }
    }
    combined_hist_min -= 0.1 * (combined_hist_max - combined_hist_min);
    combined_hist_max += 0.1 * (combined_hist_max - combined_hist_min);
    int combined_hist_bins = 256;

    std::vector<TH1D*> h1d_combined_channel_list;
    for (int _chn=0; _chn<FPGA_CHANNEL_NUMBER * fpga_count; _chn++) {
        std::string hist_name = "h1d_combined_channel_" + std::to_string(_chn);
        auto* h1d_combined_channel = new TH1D(
            hist_name.c_str(),
            ("Combined ADC-ToT distribution for channel " + std::to_string(_chn)).c_str(),
            combined_hist_bins,
            combined_hist_min,
            combined_hist_max
        );
        h1d_combined_channel->SetDirectory(nullptr);
        for (double value : channel_combined_peak_sample_values_list[_chn]) {
            h1d_combined_channel->Fill(value);
        }
        h1d_combined_channel_list.push_back(h1d_combined_channel);
    }

    auto canvas_combined_channel_distribution = new TCanvas("canvas_combined_channel_distribution", "Combined ADC-ToT distribution for all channels", 1200, 800);
    draw_mosaic_fixed(*canvas_combined_channel_distribution, h1d_combined_channel_list, mosaic_setup.topo_ped_median);
    canvas_combined_channel_distribution->Modified();
    canvas_combined_channel_distribution->Update();
    canvas_combined_channel_distribution->Write();
    canvas_combined_channel_distribution->Close();

    auto canvas_channel_tot_distribution = new TCanvas("canvas_channel_tot_distribution", "ToT distribution for all channels", 1200, 800);
    draw_mosaic_fixed(*canvas_channel_tot_distribution, h1d_tot_channel_list, mosaic_setup.topo_ped_median);
    canvas_channel_tot_distribution->Modified();
    canvas_channel_tot_distribution->Update();
    canvas_channel_tot_distribution->Write();
    canvas_channel_tot_distribution->Close();

    auto canvas_channel_tot_sample_index_distribution = new TCanvas("canvas_channel_tot_sample_index_distribution", "ToT sample index distribution for all channels", 1200, 800);
    draw_mosaic_fixed(*canvas_channel_tot_sample_index_distribution, h1d_tot_sample_index_channel_list, mosaic_setup.topo_ped_median);
    canvas_channel_tot_sample_index_distribution->Modified();
    canvas_channel_tot_sample_index_distribution->Update();
    canvas_channel_tot_sample_index_distribution->Write();
    canvas_channel_tot_sample_index_distribution->Close();

    auto canvas_channel_tot_vs_toa = new TCanvas("canvas_channel_tot_vs_toa", "ToT vs ToA for all channels", 1200, 800);
    draw_mosaic_fixed(*canvas_channel_tot_vs_toa, h2d_tot_vs_toa_channel_list, mosaic_setup.topo_ped_median);
    canvas_channel_tot_vs_toa->Modified();
    canvas_channel_tot_vs_toa->Update();
    canvas_channel_tot_vs_toa->Write();
    canvas_channel_tot_vs_toa->Close();

    TDirectory *dir_channel = output_root->mkdir("Interested_Channels");
    dir_channel->cd();
    for (int _chn : interested_channels) {
        // save the tot distribution histogram for this channel
        auto canvas_tot_distribution = new TCanvas(("canvas_tot_distribution_channel_" + std::to_string(_chn)).c_str(), ("ToT distribution for channel " + std::to_string(_chn)).c_str(), 800, 600);
        h1d_tot_channel_list[_chn]->Draw();
        canvas_tot_distribution->Modified();
        canvas_tot_distribution->Update();
        canvas_tot_distribution->Write();
        canvas_tot_distribution->Close();
    }

    output_root->Close();

    return 0;
}