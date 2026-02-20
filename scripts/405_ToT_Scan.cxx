#include "H2GCROC_Common.hxx"
#include "H2GCROC_Lib.hxx"
#include "H2GCROC_ToT.hxx"
#include "TKey.h"
#include <regex>

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

    std::string scan_info_str = "Scan";
    std::smatch scan_match;
    if (std::regex_search(input_scan_json, scan_match, std::regex("scan_number_(\\d+)"))) {
        scan_info_str += " " + scan_match[1].str();
    } else {
        scan_info_str += " Unknown";
    }

    TFile *output_root = new TFile(opts.output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << opts.output_file;
        return 1;
    }

    std::string data_adc_file_prefix = "dump/401_ADC_Analysis/Run";
    std::string data_tot_file_prefix = "dump/404_ToT_Analysis/Run";
    std::vector<int> interested_channels;
    std::vector<std::vector<TH1D*>> channel_adc_th1ds; // indexed by channel, then by run
    std::vector<std::vector<TH1D*>> channel_tot_th1ds; // indexed by channel, then by run
    std::vector<TGraphErrors*> adc_laser_intensity_graphs; // indexed by channel
    std::vector<TGraphErrors*> tot_laser_intensity_graphs; // indexed by

    for (int run_number_index = 0; run_number_index < run_numbers.size(); run_number_index++) {
        auto& run_number = run_numbers[run_number_index];
        auto& laser_intensity = laser_intensities[run_number_index];
        LOG(INFO) << "Run number: " << run_number << ", Laser intensity: " << laser_intensity;
        std::string input_adc_data_file = data_adc_file_prefix + std::to_string(run_number) + ".root";
        std::string input_tot_data_file = data_tot_file_prefix + std::to_string(run_number) + ".root";

        // * --- Read the ADC data file ---
        TFile *input_adc_root = TFile::Open(input_adc_data_file.c_str(), "READ");
        if (!input_adc_root || input_adc_root->IsZombie()) {
            LOG(ERROR) << "Failed to open input data file " << input_adc_data_file;
            continue;   
        }

        // open the directory "Interested_Channels"
        TDirectory *adc_interested_channels_dir = input_adc_root->GetDirectory("Interested_Channels");
        if (!adc_interested_channels_dir) {
            LOG(ERROR) << "Failed to get directory Interested_Channels from input data file " << input_adc_data_file;
            input_adc_root->Close();
            continue;
        }

        // print all the Canvas in the directory "Interested_Channels"
        TIter next(adc_interested_channels_dir->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*) next())) {
            if (strcmp(key->GetClassName(), "TCanvas") == 0) {
                TCanvas *canvas = (TCanvas*) key->ReadObj();
                if (canvas) {
                    std::string canvas_name = canvas->GetName();
                    LOG(INFO) << "Saved canvas " << canvas_name << " to output file " << opts.output_file;
                    // if it is canvas_peak_channel_50, extract the channel number and save the histograms in this canvas to the channel_adc_th1ds vector for later drawing
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
                                channel_adc_th1ds.push_back(std::vector<TH1D*>());
                                channel_tot_th1ds.push_back(std::vector<TH1D*>());
                            }
                        }
                        auto channel_it = std::find(interested_channels.begin(), interested_channels.end(), channel);
                        if (channel_it == interested_channels.end()) {
                            // Skip channels not in the interested list
                            continue;
                        }
                        int channel_index = std::distance(interested_channels.begin(), channel_it);
                        if (channel_index < 0 || channel_index >= (int)channel_adc_th1ds.size()) {
                            LOG(ERROR) << "Invalid channel_index " << channel_index << " for channel " << channel << ". Skipping.";
                            continue;
                        }
                        TIter next_primitive(canvas->GetListOfPrimitives());
                        while (auto primitive = next_primitive()) {
                            if (strcmp(primitive->ClassName(), "TH1D") == 0) {
                                TH1D *hist = (TH1D*) primitive;
                                TH1D *hist_clone = (TH1D*) hist->Clone((canvas_name + "_Run" + std::to_string(run_number)).c_str());
                                if (!hist_clone) {
                                    LOG(WARNING) << "Failed to clone histogram " << primitive->GetName() << ". Skipping.";
                                    continue;
                                }
                                hist_clone->SetDirectory(output_root);
                                channel_adc_th1ds[channel_index].push_back(hist_clone);
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
                                channel_adc_th1ds.push_back(std::vector<TH1D*>());
                                channel_tot_th1ds.push_back(std::vector<TH1D*>());
                            }
                        }
                        auto channel_it = std::find(interested_channels.begin(), interested_channels.end(), channel);
                        if (channel_it == interested_channels.end()) {
                            // Skip channels not in the interested list
                            continue;
                        }
                        int channel_index = std::distance(interested_channels.begin(), channel_it);
                        if (channel_index < 0 || channel_index >= (int)channel_adc_th1ds.size()) {
                            LOG(ERROR) << "Invalid channel_index " << channel_index << " for channel " << channel << ". Skipping.";
                            continue;
                        }
                        TIter next_primitive(canvas->GetListOfPrimitives());
                        while (auto primitive = next_primitive()) {
                            if (strcmp(primitive->ClassName(), "TH1D") == 0) {
                                LOG(INFO) << "Extracted TH1D " << primitive->GetName() << " from canvas " << canvas_name;
                                TH1D *hist = (TH1D*) primitive;
                                TH1D *hist_clone = (TH1D*) hist->Clone((canvas_name + "_Run" + std::to_string(run_number)).c_str());
                                if (!hist_clone) {
                                    LOG(WARNING) << "Failed to clone histogram " << primitive->GetName() << ". Skipping.";
                                    continue;
                                }
                                hist_clone->SetDirectory(output_root);
                                channel_adc_th1ds[channel_index].push_back(hist_clone);
                            }
                        }
                    }
                } else {
                    LOG(ERROR) << "Failed to read canvas " << key->GetName() << " from input data file " << input_adc_data_file;
                }
            }
        }
        input_adc_root->Close();

        // * --- Read the ToT data file ---
        TFile *input_tot_root = TFile::Open(input_tot_data_file.c_str(), "READ");
        if (!input_tot_root || input_tot_root->IsZombie()) {
            LOG(ERROR) << "Failed to open input data file " << input_tot_data_file;
            continue;
        }
        // open the directory "Interested_Channels"
        TDirectory *tot_interested_channels_dir = input_tot_root->GetDirectory("Interested_Channels");
        if (!tot_interested_channels_dir) {
            LOG(ERROR) << "Failed to get directory Interested_Channels from input data file " << input_tot_data_file;
            input_tot_root->Close();
            continue;
        }

        // print all the Canvas in the directory "Interested_Channels"
        TIter next_tot(tot_interested_channels_dir->GetListOfKeys());
        while ((key = (TKey*) next_tot())) {
            if (strcmp(key->GetClassName(), "TCanvas") == 0) {
                TCanvas *canvas = (TCanvas*) key->ReadObj();
                if (canvas) {
                    std::string canvas_name = canvas->GetName();
                    LOG(INFO) << "Saved canvas " << canvas_name << " to output file " << opts.output_file;
                    // if it is canvas_tot_distribution_channel_50, extract the channel number and save the histogram in this canvas to the channel_tot_th1ds vector for later drawing
                    if (canvas_name.find("canvas_tot_distribution_channel_") != std::string::npos) {
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
                                channel_adc_th1ds.push_back(std::vector<TH1D*>());
                                channel_tot_th1ds.push_back(std::vector<TH1D*>());
                            }
                        }
                        auto channel_it = std::find(interested_channels.begin(), interested_channels.end(), channel);
                        if (channel_it == interested_channels.end()) {
                            // Skip channels not in the interested list
                            continue;
                        }
                        int channel_index = std::distance(interested_channels.begin(), channel_it);
                        if (channel_index < 0 || channel_index >= (int)channel_tot_th1ds.size()) {
                            LOG(ERROR) << "Invalid channel_index " << channel_index << " for channel " << channel << ". Skipping.";
                            continue;
                        }
                        TIter next_primitive(canvas->GetListOfPrimitives());
                        while (auto primitive = next_primitive()) {
                            if (strcmp(primitive->ClassName(), "TH1D") == 0) {
                                TH1D *hist = (TH1D*) primitive;
                                TH1D *hist_clone = (TH1D*) hist->Clone((canvas_name + "_Run" + std::to_string(run_number)).c_str());
                                if (!hist_clone) {
                                    LOG(WARNING) << "Failed to clone histogram " << primitive->GetName() << ". Skipping.";
                                    continue;
                                }
                                hist_clone->SetDirectory(output_root);
                                channel_tot_th1ds[channel_index].push_back(hist_clone);
                            }
                        }
                    }
                } else {
                    LOG(ERROR) << "Failed to read canvas " << key->GetName() << " from input data file " << input_tot_data_file;
                }
            }
        }

        input_tot_root->Close();
    }

    LOG(INFO) << "Interested channels: ";
    for (size_t i = 0; i < interested_channels.size(); i++) {
        LOG(INFO) << "Channel " << interested_channels[i] << ": " << channel_adc_th1ds[i].size() << " histograms";
    }   

    output_root->cd();

    std::vector<std::vector<double>> channel_adc_peak_mean_list(interested_channels.size());
    std::vector<std::vector<double>> channel_adc_peak_sigma_list(interested_channels.size());
    std::vector<double> channel_laser_intensity_list;
    std::vector<double> channel_laser_intensity_error_list;

    double max_y_global = 0.0;
    for (size_t channel_index = 0; channel_index < interested_channels.size(); channel_index++) {
        for (size_t run_number_index = 0; run_number_index < channel_adc_th1ds[channel_index].size(); run_number_index++) {
            TH1D *hist = channel_adc_th1ds[channel_index][run_number_index];
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
            if (channel_index < channel_adc_th1ds.size() && run_number_index < channel_adc_th1ds[channel_index].size()) {
                TH1D *hist = channel_adc_th1ds[channel_index][run_number_index];
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
        graph->SetTitle("");
        graph->GetXaxis()->SetTitle("Laser Intensity");
        graph->GetYaxis()->SetTitle("Mean ADC Peak");
        // set axis range
        if (!channel_laser_intensity_list.empty()) {
            graph->GetXaxis()->SetRangeUser(0, *std::max_element(channel_laser_intensity_list.begin(), channel_laser_intensity_list.end()) * 1.2);
        }
        if (!channel_adc_peak_mean_list[channel_index].empty()) {
            graph->GetYaxis()->SetRangeUser(0, *std::max_element(channel_adc_peak_mean_list[channel_index].begin(), channel_adc_peak_mean_list[channel_index].end()) * 1.4);
        }
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.0);
        graph->Draw("AEP");

        TGraphErrors* graph_clone = (TGraphErrors*) graph->Clone(("graph_mean_peak_vs_laser_channel_" + std::to_string(channel)).c_str());
        adc_laser_intensity_graphs.push_back(graph_clone);

        
        // write latex info
        TLatex latex;
        latex.SetNDC();
        latex.SetTextAlign(13);
        latex.SetTextSize(0.04);
        latex.SetTextFont(62);
        latex.DrawLatex(0.13, 0.88, "Laser Test for H2GCROC");
        latex.SetTextSize(0.03);
        latex.SetTextFont(42);
        latex.DrawLatex(0.13, 0.84, ("ADC laser intensity scan, Channel " + std::to_string(channel)).c_str());
        latex.DrawLatex(0.13, 0.80, "Hamamatsu S14160-6010PS");
        latex.DrawLatex(0.13, 0.76, "CERN, February 2026");
        
        canvas->Update();
        canvas->Write();
        // save as a seprate pdf file
        std::string pdf_output_file = opts.output_file;
        pdf_output_file.replace(pdf_output_file.find(".root"), 5, "_mean_peak_vs_laser_channel_" + std::to_string(channel) + ".pdf");
        canvas->SaveAs(pdf_output_file.c_str());
        canvas->Close();
    }

    for (size_t channel_index = 0; channel_index < interested_channels.size(); channel_index++) {
        int channel = interested_channels[channel_index];
        TCanvas *canvas = new TCanvas(("canvas_tot_distribution_channel_" + std::to_string(channel)).c_str(), ("ToT Distribution - Channel " + std::to_string(channel)).c_str(), 1000, 600);
        canvas->cd();
        std::vector<double> tot_means;
        std::vector<double> tot_mean_errors;
        std::vector<double> tot_laser_intensities;
        std::vector<double> tot_laser_intensity_errors;
        std::vector<TH1D*> hists_to_draw;
        std::vector<TGraphErrors*> graphs_to_draw;

        // get the 1% min and 99% max for all hist to set the global x range for the histograms
        double global_x_min = 4096;
        double global_x_max = 0;
        for (size_t run_number_index = 0; run_number_index < run_numbers.size(); run_number_index++) {
            if (channel_index < channel_tot_th1ds.size() && run_number_index < channel_tot_th1ds[channel_index].size()) {
                TH1D *hist = channel_tot_th1ds[channel_index][run_number_index];
                if (hist) {
                    // get the 1% min and 99% max quantiles
                    double quantile_result_min = 0;
                    double quantile_prob_min = 0.05;
                    hist->GetQuantiles(1, &quantile_result_min, &quantile_prob_min);
                    double x_min = quantile_result_min;
                    
                    double quantile_result_max = 0;
                    double quantile_prob_max = 0.95;
                    hist->GetQuantiles(1, &quantile_result_max, &quantile_prob_max);
                    double x_max = quantile_result_max;
                    
                    if (x_min < global_x_min) global_x_min = x_min;
                    if (x_max > global_x_max) global_x_max = x_max;
                }
            }
        }
        
        
        // First pass: collect all histograms and compute statistics
        double max_hist_value = 0;
        double max_tot_mean = 0;
        for (size_t run_number_index = 0; run_number_index < run_numbers.size(); run_number_index++) {
            double laser_intensity = laser_intensities[run_number_index];
            if (channel_index < channel_tot_th1ds.size() && run_number_index < channel_tot_th1ds[channel_index].size()) {
                TH1D *hist = channel_tot_th1ds[channel_index][run_number_index];
                if (hist) {
                    hist->SetLineColor(run_number_index + 1);
                    hist->GetXaxis()->SetRangeUser(global_x_min, global_x_max);
                    double hist_max = hist->GetMaximum();
                    if (hist_max > max_hist_value) max_hist_value = hist_max;
                    hists_to_draw.push_back(hist);
                }
                // use static value to get mean and error
                // double mean = hist ? hist->GetMean() : 0.0;
                // double error = hist ? hist->GetMeanError() : 0.0;
                double probs[3] = {0.16, 0.50, 0.84};
                double quantiles[3];
                if (hist) {
                    hist->GetQuantiles(3, quantiles, probs);
                } else {
                    quantiles[0] = quantiles[1] = quantiles[2] = 0.0;
                }
                double x16 = quantiles[0];
                double x50 = quantiles[1];
                double x84 = quantiles[2];

                double err_low = x50 - x16;
                double err_high = x84 - x50;
                double mean = x50;
                double error = (err_low + err_high) / 2.0; // use

                if (mean > max_tot_mean) max_tot_mean = mean;
                LOG(INFO) << "Run " << run_numbers[run_number_index] << ", Channel " << channel << ": mean ToT = " << mean << ", error = " << error;
                // create the mean value with error bar
                TGraphErrors *graph = new TGraphErrors(1);
                graph->SetPoint(0, mean, 100); // x is mean ToT, y is fixed at 100 for visualization
                graph->SetPointError(0, error, 0); // x error is the error of mean ToT, y error is 0
                graph->SetMarkerStyle(20);
                graph->SetMarkerSize(1.0);
                graph->SetMarkerColor(run_number_index + 1);
                graph->SetLineColor(run_number_index + 1);
                graphs_to_draw.push_back(graph);

                tot_means.push_back(mean);
                tot_mean_errors.push_back(error);
                tot_laser_intensities.push_back(laser_intensity);
                tot_laser_intensity_errors.push_back(0.01); // assume 1% error for the laser intensity
            }
        }
        
        // Determine appropriate Y range to fit both histograms and data points
        double y_max = std::max(max_hist_value * 1.2, 120.0); // add some margin above the max histogram value, and ensure it's above the data points at y=100
        
        // Second pass: draw histograms with proper range
        for (size_t i = 0; i < hists_to_draw.size(); i++) {
            TH1D *hist = hists_to_draw[i];
            hist->GetYaxis()->SetRangeUser(0, y_max);
            hist->GetXaxis()->SetRangeUser(global_x_min, global_x_max);
            // add axis label and tick
            hist->GetXaxis()->SetTitle("ToT Value");
            hist->GetYaxis()->SetTitle("Counts");
            hist->GetXaxis()->SetTitleSize(0.05);
            hist->GetYaxis()->SetTitleSize(0.05);
            hist->GetXaxis()->SetLabelSize(0.04);
            hist->GetYaxis()->SetLabelSize(0.04);
            hist->GetXaxis()->SetNdivisions(505);
            hist->GetYaxis()->SetNdivisions(505);
            hist->Draw(i == 0 ? "HIST" : "HIST SAME");
        }
        
        // Third pass: draw all graphs on top
        for (auto graph : graphs_to_draw) {
            graph->Draw("|> SAME");
        }
        
        canvas->Update();
        canvas->Write();
        canvas->Close();

        TCanvas *canvas_mean_tot_vs_laser = new TCanvas(("canvas_mean_tot_vs_laser_channel_" + std::to_string(channel)).c_str(), ("Mean ToT vs Laser Intensity - Channel " + std::to_string(channel)).c_str(), 800, 600);
        canvas_mean_tot_vs_laser->cd();
        TGraphErrors *graph_mean_tot = new TGraphErrors(tot_means.size());
        for (size_t i = 0; i < tot_means.size(); i++) {
            graph_mean_tot->SetPoint(i, tot_laser_intensities[i], tot_means[i]);
            graph_mean_tot->SetPointError(i, tot_laser_intensity_errors[i], tot_mean_errors[i]);
        }
        graph_mean_tot->SetTitle("");
        graph_mean_tot->GetXaxis()->SetTitle("Laser Intensity");
        graph_mean_tot->GetYaxis()->SetTitle("Mean ToT");
        // set axis range
        if (!tot_laser_intensities.empty()) {
            graph_mean_tot->GetXaxis()->SetRangeUser(0, *std::max_element(tot_laser_intensities.begin(), tot_laser_intensities.end()) * 1.2);
        }
        if (!tot_means.empty()) {
            graph_mean_tot->GetYaxis()->SetRangeUser(0, *std::max_element(tot_means.begin(), tot_means.end()) * 1.2);
        }
        graph_mean_tot->SetMarkerStyle(20);
        graph_mean_tot->SetMarkerSize(1.0);
        graph_mean_tot->Draw("AEP");

        TGraphErrors* graph_mean_tot_clone = (TGraphErrors*)graph_mean_tot->Clone(("graph_tot_vs_laser_channel_" + std::to_string(channel)).c_str());
        tot_laser_intensity_graphs.push_back(graph_mean_tot_clone);

        // write latex info
        TLatex latex;
        latex.SetNDC();
        latex.SetTextAlign(13);
        latex.SetTextSize(0.04);
        latex.SetTextFont(62);
        latex.DrawLatex(0.13, 0.88, "Laser Test for H2GCROC");
        latex.SetTextSize(0.03);
        latex.SetTextFont(42);
        latex.DrawLatex(0.13, 0.84, ("ToT laser intensity scan, Channel " + std::to_string(channel)).c_str());
        latex.DrawLatex(0.13, 0.80, "Hamamatsu S14160-6010PS");
        latex.DrawLatex(0.13, 0.76, "CERN, February 2026");

        canvas_mean_tot_vs_laser->Update();
        canvas_mean_tot_vs_laser->Write();
        // save as a seprate pdf file
        std::string pdf_output_file = opts.output_file;
        pdf_output_file.replace(pdf_output_file.find(".root"), 5, "_mean_tot_vs_laser_channel_" + std::to_string(channel) + ".pdf");
        canvas_mean_tot_vs_laser->SaveAs(pdf_output_file.c_str());
        canvas_mean_tot_vs_laser->Close();
    }


    // ! start building LUT
    for (size_t channel_index = 0; channel_index < interested_channels.size(); channel_index++) {
        int channel = interested_channels[channel_index];
        auto& graph_adc_laser = adc_laser_intensity_graphs[channel_index];
        auto& graph_tot_laser = tot_laser_intensity_graphs[channel_index];
        if (graph_adc_laser && graph_tot_laser) {
            auto r = BuildTotToAdcLUT_FromGraphs(
                graph_adc_laser,
                graph_tot_laser,
                4096,          // tot bins
                150.0,         // linear fit min RAW ADC
                950.0,         // linear fit max RAW ADC
                4000,          // laser samples (int!!)
                0.5,           // pit epsilon
                "samples_channel" + std::to_string(channel), // string
                100.0,         // baseline
                65535.0,       // adc_out_max (不要再用1023!)
                1023.0         // pit raw override
            );
            LOG(INFO) << "Channel " << channel << " - ADC(L) fit: ADC = " << r.alpha << " * L + " << r.beta;
            if (r.has_pit) {
                LOG(WARNING) << "Detected ToT pit in L interval: [" << r.pit_L_start << ", " << r.pit_L_end
                            << "], mapped to ADC=1023.";
            }
            // std::ofstream lut_out("LUT_Channel_" + std::to_string(channel) + ".txt");
            std::ofstream lut_out((opts.output_file + "_LUT_Channel_" + std::to_string(channel) + ".txt").c_str());
            for (int t = 0; t < (int)r.lut.size(); t++) {
                lut_out << t << " " << r.lut[t] << "\n";
            }
            lut_out.close();

        } else {
            LOG(WARNING) << "Missing graph for Channel " << channel << ": " 
                         << (graph_adc_laser ? "" : "ADC graph ") 
                         << (graph_tot_laser ? "" : "ToT graph ") 
                         << ". Skipping LUT building for this channel.";
        }
        
    }
    output_root->Close();


    return 0;
}
