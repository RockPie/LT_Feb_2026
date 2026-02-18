#include "H2GCROC_Common.hxx"
#include "H2GCROC_Lib.hxx"
#include "TKey.h"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv) {
    ScriptOptions opts = parse_arguments_single_json(argc, argv, "1.0");

    gROOT->SetBatch(kTRUE);
    const double sample_time = 25.0; // unit: ns
    const double phase_shift_time = 1.5625; // unit: ns

    std::string input_scan_json = opts.input_file;

    LOG(INFO) << "Input scan JSON file: " << input_scan_json;

    json scan_json;
    try {
        std::ifstream ifs(input_scan_json);
        if (!ifs.is_open()) {
            LOG(ERROR) << "Failed to open input scan JSON file " << input_scan_json;
            return 1;
        }
        ifs >> scan_json;
    } catch (const std::exception& e) {
        LOG(ERROR) << "Failed to parse input scan JSON file " << input_scan_json << ": " << e.what();
        return 1;
    }
    const auto& scan_brief = scan_json["scan_brief"].get<std::string>();
    const auto& scan_bias = scan_json["scan_bias"].get<double>();
    const auto& scan_CC = scan_json["scan_CC"].get<double>();
    const auto& scan_Cf = scan_json["scan_Cf"].get<double>();
    const auto& scan_Cfcomp = scan_json["scan_Cfcomp"].get<double>();
    const auto& run_numbers = scan_json["run_numbers"].get<std::vector<int>>();
    const auto& laser_intensities = scan_json["laser_intensities"].get<std::vector<double>>();

    LOG(INFO) << "Scan brief: " << scan_brief;
    LOG(INFO) << "Scan bias: " << scan_bias;
    LOG(INFO) << "Scan CC: " << scan_CC;
    LOG(INFO) << "Scan Cf: " << scan_Cf;
    LOG(INFO) << "Scan Cfcomp: " << scan_Cfcomp;

    size_t scan_number_pos = input_scan_json.find("scan_number_");
    std::string scan_info_str = "Scan";
    if (scan_number_pos != std::string::npos) {
        scan_info_str += " " + input_scan_json.substr(scan_number_pos + 12, 1);
    } else {
        scan_info_str += " Unknown";
    }

    TFile *output_root = new TFile(opts.output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << opts.output_file;
        return 1;
    }

    std::string data_file_prefix = "dump/401_ADC_Analysis/Run";
    std::vector<int> interested_channels;
    std::vector<std::vector<TH1D*>> channel_histograms; // indexed by channel, then by run

    for (int run_number_index = 0; run_number_index < run_numbers.size(); run_number_index++) {
        auto& run_number = run_numbers[run_number_index];
        auto& laser_intensity = laser_intensities[run_number_index];
        LOG(INFO) << "Run number: " << run_number << ", Laser intensity: " << laser_intensity;
        std::string input_data_file = data_file_prefix + std::to_string(run_number) + ".root";
        TFile *input_root = TFile::Open(input_data_file.c_str(), "READ");
        if (!input_root || input_root->IsZombie()) {
            LOG(ERROR) << "Failed to open input data file " << input_data_file;
            continue;
        }

        // open the directory "Interested_Channels"
        TDirectory *interested_channels_dir = input_root->GetDirectory("Interested_Channels");
        if (!interested_channels_dir) {
            LOG(ERROR) << "Failed to get directory Interested_Channels from input data file " << input_data_file;
            input_root->Close();
            continue;
        }

        // print all the Canvas in the directory "Interested_Channels"
        TIter next(interested_channels_dir->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*) next())) {
            if (strcmp(key->GetClassName(), "TCanvas") == 0) {
                TCanvas *canvas = (TCanvas*) key->ReadObj();
                if (canvas) {
                    std::string canvas_name = canvas->GetName();
                    LOG(INFO) << "Saved canvas " << canvas_name << " to output file " << opts.output_file;
                    // if it is canvas_peak_channel_50, extract the channel number and save the histograms in this canvas to the channel_histograms vector for later drawing
                    if (canvas_name.find("canvas_peak_channel_") != std::string::npos) {
                        size_t channel_pos = canvas_name.rfind('_');
                        int channel = (channel_pos != std::string::npos) ? std::stoi(canvas_name.substr(channel_pos + 1)) : -1;
                        LOG(INFO) << "Extracted channel " << channel << " from canvas name " << canvas_name;
                        if (channel < 0) {
                            LOG(WARNING) << "Failed to parse channel from canvas name " << canvas_name << ". Skipping.";
                            continue;
                        }
                        if (run_number_index == 0) {
                            // add if not in the interested channels
                            if (std::find(interested_channels.begin(), interested_channels.end(), channel) == interested_channels.end()) {
                                interested_channels.push_back(channel);
                                channel_histograms.push_back(std::vector<TH1D*>());
                            }
                        }
                        auto channel_it = std::find(interested_channels.begin(), interested_channels.end(), channel);
                        if (channel_it == interested_channels.end()) {
                            // Skip channels not in the interested list
                            continue;
                        }
                        int channel_index = std::distance(interested_channels.begin(), channel_it);
                        TIter next_primitive(canvas->GetListOfPrimitives());
                        while (auto primitive = next_primitive()) {
                            if (strcmp(primitive->ClassName(), "TH1D") == 0) {
                                TH1D *hist = (TH1D*) primitive;
                                TH1D *hist_clone = (TH1D*) hist->Clone((canvas_name + "_Run" + std::to_string(run_number)).c_str());
                                hist_clone->SetDirectory(output_root);
                                channel_histograms[channel_index].push_back(hist_clone);
                            }
                        }
                    }
                    // canvas_sliding_channel_50
                    else if (canvas_name.find("canvas_sliding_channel_") != std::string::npos) {
                        size_t channel_pos = canvas_name.rfind('_');
                        int channel = (channel_pos != std::string::npos) ? std::stoi(canvas_name.substr(channel_pos + 1)) : -1;
                        LOG(INFO) << "Extracted channel " << channel << " from canvas name " << canvas_name;
                        if (channel < 0) {
                            LOG(WARNING) << "Failed to parse channel from canvas name " << canvas_name << ". Skipping.";
                            continue;
                        }
                        if (run_number_index == 0) {
                            // add if not in the interested channels
                            if (std::find(interested_channels.begin(), interested_channels.end(), channel) == interested_channels.end()) {
                                interested_channels.push_back(channel);
                                channel_histograms.push_back(std::vector<TH1D*>());
                            }
                        }
                        auto channel_it = std::find(interested_channels.begin(), interested_channels.end(), channel);
                        if (channel_it == interested_channels.end()) {
                            // Skip channels not in the interested list
                            continue;
                        }
                        int channel_index = std::distance(interested_channels.begin(), channel_it);
                        TIter next_primitive(canvas->GetListOfPrimitives());
                        while (auto primitive = next_primitive()) {
                            if (strcmp(primitive->ClassName(), "TH1D") == 0) {
                                LOG(INFO) << "Extracted TH1D " << primitive->GetName() << " from canvas " << canvas_name;
                                TH1D *hist = (TH1D*) primitive;
                                TH1D *hist_clone = (TH1D*) hist->Clone((canvas_name + "_Run" + std::to_string(run_number)).c_str());
                                hist_clone->SetDirectory(output_root);
                                channel_histograms[channel_index].push_back(hist_clone);
                            }
                        }
                    }
                } else {
                    LOG(ERROR) << "Failed to read canvas " << key->GetName() << " from input data file " << input_data_file;
                }
                canvas->Close();
            }
        }


        input_root->Close();

    }

    LOG(INFO) << "Interested channels: ";
    for (size_t i = 0; i < interested_channels.size(); i++) {
        LOG(INFO) << "Channel " << interested_channels[i] << ": " << channel_histograms[i].size() << " histograms";
    }   

    output_root->cd();

    std::vector<std::vector<double>> channel_adc_peak_mean_list(interested_channels.size());
    std::vector<std::vector<double>> channel_adc_peak_sigma_list(interested_channels.size());
    std::vector<double> channel_laser_intensity_list;
    std::vector<double> channel_laser_intensity_error_list;

    double max_y_global = 0.0;
    for (size_t channel_index = 0; channel_index < interested_channels.size(); channel_index++) {
        for (size_t run_number_index = 0; run_number_index < channel_histograms[channel_index].size(); run_number_index++) {
            TH1D *hist = channel_histograms[channel_index][run_number_index];
            if (hist) {
                double hist_max = hist->GetMaximum();
                if (hist_max > max_y_global) {
                    max_y_global = hist_max;
                }
            }
        }
    }

    // draw the th1d in the same canvas for each interested channel
    for (size_t channel_index = 0; channel_index < interested_channels.size(); channel_index++) {
        int channel = interested_channels[channel_index];
        TCanvas *canvas = new TCanvas(("canvas_peak_channel_" + std::to_string(channel)).c_str(), ("Channel " + std::to_string(channel)).c_str(), 800, 600);
        canvas->cd();
        for (size_t run_number_index = 0; run_number_index < run_numbers.size(); run_number_index++) {
            if (channel_index < channel_histograms.size() && run_number_index < channel_histograms[channel_index].size()) {
                TH1D *hist = channel_histograms[channel_index][run_number_index];
                // get the gaussian fit with it
                TF1 *gaus_fit = (TF1*) hist->GetFunction("fit_func");
                if (gaus_fit) {
                    double mean = gaus_fit->GetParameter(1);
                    double sigma = gaus_fit->GetParameter(2);
                    LOG(INFO) << "Run " << run_numbers[run_number_index] << ", Channel " << channel << ": mean = " << mean << ", sigma = " << sigma;
                    channel_adc_peak_mean_list[channel_index].push_back(mean);
                    channel_adc_peak_sigma_list[channel_index].push_back(sigma);
                } else {
                    LOG(WARNING) << "No gaussian fit found for Run " << run_numbers[run_number_index] << ", Channel " << channel;
                    channel_adc_peak_mean_list[channel_index].push_back(0);
                    channel_adc_peak_sigma_list[channel_index].push_back(0);
                }
                if (channel_index == 0) {
                    channel_laser_intensity_list.push_back(laser_intensities[run_number_index]);
                    channel_laser_intensity_error_list.push_back(0.01); // assume 1% error for the laser intensity
                }
                if (hist) {
                    // set x range
                    hist->GetXaxis()->SetRangeUser(0, 1023);
                    hist->GetYaxis()->SetRangeUser(0, max_y_global * 1.2);
                    hist->SetLineColor(run_number_index + 1);
                    // Remove all functions from the histogram before drawing to hide fit curves
                    // hist->GetListOfFunctions()->Clear();
                    hist->Draw(run_number_index == 0 ? "" : "HIST SAME");
                }
            }
        }
        canvas->Write();
        canvas->Close();
    }

    // print the mean ADC peak and sigma for each channel and each run number
    for (size_t channel_index = 0; channel_index < interested_channels.size(); channel_index++) {
        int channel = interested_channels[channel_index];
        LOG(INFO) << "Channel " << channel << ":";
        for (size_t run_number_index = 0; run_number_index < run_numbers.size(); run_number_index++) {
            LOG(INFO) << "  Run " << run_numbers[run_number_index] << ": mean ADC peak = " << channel_adc_peak_mean_list[channel_index][run_number_index] << ", sigma = " << channel_adc_peak_sigma_list[channel_index][run_number_index];
        }
    }

    // draw the error graph of mean ADC peak vs laser intensity for each channel
    for (size_t channel_index = 0; channel_index < interested_channels.size(); channel_index++) {
        int channel = interested_channels[channel_index];
        TCanvas *canvas = new TCanvas(("canvas_mean_peak_vs_laser_channel_" + std::to_string(channel)).c_str(), ("Mean ADC Peak vs Laser Intensity - Channel " + std::to_string(channel)).c_str(), 800, 600);
        canvas->cd();
        TGraphErrors *graph = new TGraphErrors(channel_laser_intensity_list.size());
        for (size_t i = 0; i < channel_laser_intensity_list.size(); i++) {
            // filter out nan value
            if (std::isnan(channel_adc_peak_mean_list[channel_index][i]) || std::isnan(channel_adc_peak_sigma_list[channel_index][i])) {
                LOG(WARNING) << "NaN value found for Channel " << channel << ", Laser Intensity " << channel_laser_intensity_list[i] << ": mean = " << channel_adc_peak_mean_list[channel_index][i] << ", sigma = " << channel_adc_peak_sigma_list[channel_index][i] << ". Skipping this point.";
                continue;
            }
            graph->SetPoint(i, channel_laser_intensity_list[i], channel_adc_peak_mean_list[channel_index][i]);
            graph->SetPointError(i, channel_laser_intensity_error_list[i], channel_adc_peak_sigma_list[channel_index][i]);
        }
        graph->SetTitle(("Mean ADC Peak vs Laser Intensity - Channel " + std::to_string(channel)).c_str());
        graph->GetXaxis()->SetTitle("Laser Intensity");
        graph->GetYaxis()->SetTitle("Mean ADC Peak");
        // set axis range
        graph->GetXaxis()->SetRangeUser(0, *std::max_element(channel_laser_intensity_list.begin(), channel_laser_intensity_list.end()) * 1.2);
        graph->GetYaxis()->SetRangeUser(0, *std::max_element(channel_adc_peak_mean_list[channel_index].begin(), channel_adc_peak_mean_list[channel_index].end()) * 1.2);
        graph->SetMarkerStyle(20);
        graph->Draw("AP");
        canvas->Write();
        canvas->Close();
    }

    output_root->Close();

    return 0;
}
