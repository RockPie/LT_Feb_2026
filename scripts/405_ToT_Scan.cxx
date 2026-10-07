#include "H2GCROC_Common.hxx"
#include "H2GCROC_Lib.hxx"
#include "H2GCROC_ToT.hxx"
#include "TKey.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include <regex>

INITIALIZE_EASYLOGGINGPP

struct TrimmedMeanStats {
    double mean = 0.0;
    double standard_error = 0.0;
    double q_low = 0.0;
    double q_high = 0.0;
    double retained_entries = 0.0;
    bool valid = false;
};

static TrimmedMeanStats CalculateTrimmedMean(const TH1D* hist, double trim_fraction) {
    TrimmedMeanStats result;
    if (!hist || trim_fraction < 0.0 || trim_fraction >= 0.5) {
        return result;
    }

    const double total_weight = hist->Integral(1, hist->GetNbinsX());
    if (total_weight <= 1.0) {
        return result;
    }

    const double probabilities[2] = {trim_fraction, 1.0 - trim_fraction};
    double quantiles[2] = {0.0, 0.0};
    const_cast<TH1D*>(hist)->GetQuantiles(2, quantiles, probabilities);
    result.q_low = quantiles[0];
    result.q_high = quantiles[1];

    const double retained_begin = trim_fraction * total_weight;
    const double retained_end = (1.0 - trim_fraction) * total_weight;
    double cumulative_weight = 0.0;
    double trimmed_sum = 0.0;
    for (int bin = 1; bin <= hist->GetNbinsX(); bin++) {
        const double bin_weight = hist->GetBinContent(bin);
        const double bin_begin = cumulative_weight;
        const double bin_end = cumulative_weight + bin_weight;
        const double retained_weight = std::max(
            0.0,
            std::min(bin_end, retained_end) - std::max(bin_begin, retained_begin));
        trimmed_sum += retained_weight * hist->GetBinCenter(bin);
        result.retained_entries += retained_weight;
        cumulative_weight = bin_end;
    }
    if (result.retained_entries <= 0.0) {
        return result;
    }
    result.mean = trimmed_sum / result.retained_entries;

    double winsorized_sum = 0.0;
    double winsorized_sum_squares = 0.0;
    for (int bin = 1; bin <= hist->GetNbinsX(); bin++) {
        const double bin_weight = hist->GetBinContent(bin);
        const double winsorized_value = std::clamp(
            hist->GetBinCenter(bin), result.q_low, result.q_high);
        winsorized_sum += bin_weight * winsorized_value;
        winsorized_sum_squares += bin_weight * winsorized_value * winsorized_value;
    }
    const double winsorized_mean = winsorized_sum / total_weight;
    const double winsorized_variance = std::max(
        0.0,
        (winsorized_sum_squares - total_weight * winsorized_mean * winsorized_mean)
            / (total_weight - 1.0));
    result.standard_error = std::sqrt(winsorized_variance / total_weight)
        / (1.0 - 2.0 * trim_fraction);
    result.valid = true;
    return result;
}

int main(int argc, char **argv) {
    ScriptOptions opts = parse_arguments_single_json(argc, argv, "1.0");

    gROOT->SetBatch(kTRUE);
    const double sample_time = 25.0; // unit: ns
    const double phase_shift_time = 1.5625; // unit: ns

    const auto execution_time = std::time(nullptr);
    const auto local_execution_time = *std::localtime(&execution_time);
    char execution_time_buffer[64];
    std::strftime(execution_time_buffer, sizeof(execution_time_buffer), "%d-%m-%Y %H:%M:%S", &local_execution_time);
    const std::string execution_time_label = "CERN, " + std::string(execution_time_buffer);

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
    const bool has_configured_example_channels = scan_json.contains("example_channels");
    const std::vector<int> example_channels = has_configured_example_channels
        ? scan_json["example_channels"].get<std::vector<int>>()
        : std::vector<int>();

    const int min_entrys_threshold = 1;

    LOG(INFO) << "Scan brief: " << scan_brief;
    LOG(INFO) << "Scan bias: " << scan_bias;
    LOG(INFO) << "Scan CC: " << scan_CC;
    LOG(INFO) << "Scan Cf: " << scan_Cf;
    LOG(INFO) << "Scan Cfcomp: " << scan_Cfcomp;
    LOG(INFO) << "Example channels source: "
              << (has_configured_example_channels ? "scan config" : "ADC/ToT analysis output");

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
    if (input_scan_json.find("LT_Sep_2026") != std::string::npos) {
        data_adc_file_prefix = "dump/401_ADC_Analysis/LT_Sep_2026/Run";
        data_tot_file_prefix = "dump/404_ToT_Analysis/LT_Sep_2026/Run";
    }
    std::vector<int> interested_channels = example_channels;
    std::vector<std::vector<TH1D*>> channel_adc_th1ds(
        interested_channels.size(), std::vector<TH1D*>(run_numbers.size(), nullptr));
    std::vector<std::vector<TH1D*>> channel_tot_th1ds(
        interested_channels.size(), std::vector<TH1D*>(run_numbers.size(), nullptr));
    std::vector<TGraphErrors*> adc_laser_intensity_graphs; // indexed by channel
    std::vector<TGraphAsymmErrors*> tot_laser_intensity_graphs; // indexed by channel

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
                        if (!has_configured_example_channels && run_number_index == 0
                            && std::find(interested_channels.begin(), interested_channels.end(), channel) == interested_channels.end()) {
                            interested_channels.push_back(channel);
                            channel_adc_th1ds.push_back(std::vector<TH1D*>(run_numbers.size(), nullptr));
                            channel_tot_th1ds.push_back(std::vector<TH1D*>(run_numbers.size(), nullptr));
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
                                channel_adc_th1ds[channel_index][run_number_index] = hist_clone;
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
                        if (!has_configured_example_channels && run_number_index == 0
                            && std::find(interested_channels.begin(), interested_channels.end(), channel) == interested_channels.end()) {
                            interested_channels.push_back(channel);
                            channel_adc_th1ds.push_back(std::vector<TH1D*>(run_numbers.size(), nullptr));
                            channel_tot_th1ds.push_back(std::vector<TH1D*>(run_numbers.size(), nullptr));
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
                                if (!channel_adc_th1ds[channel_index][run_number_index]) {
                                    channel_adc_th1ds[channel_index][run_number_index] = hist_clone;
                                }
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
                        if (!has_configured_example_channels && run_number_index == 0
                            && std::find(interested_channels.begin(), interested_channels.end(), channel) == interested_channels.end()) {
                            interested_channels.push_back(channel);
                            channel_adc_th1ds.push_back(std::vector<TH1D*>(run_numbers.size(), nullptr));
                            channel_tot_th1ds.push_back(std::vector<TH1D*>(run_numbers.size(), nullptr));
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
                                if (!channel_tot_th1ds[channel_index][run_number_index]) {
                                    channel_tot_th1ds[channel_index][run_number_index] = hist_clone;
                                } else {
                                    LOG(WARNING) << "Multiple ToT histograms found for Run " << run_number
                                                 << ", Channel " << channel << ". Keeping the first one.";
                                }
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
                if (!hist) {
                    LOG(WARNING) << "No ADC histogram found for Run " << run_numbers[run_number_index]
                                 << ", Channel " << channel << ".";
                    channel_adc_peak_mean_list[channel_index].push_back(0);
                    channel_adc_peak_sigma_list[channel_index].push_back(0);
                    if (channel_index == 0) {
                        channel_laser_intensity_list.push_back(laser_intensities[run_number_index]);
                        channel_laser_intensity_error_list.push_back(0.01);
                    }
                    continue;
                }
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
        if (channel_adc_peak_mean_list[channel_index].size() != run_numbers.size()
            || channel_adc_peak_sigma_list[channel_index].size() != run_numbers.size()) {
            LOG(WARNING) << "Skipping ADC summary for Channel " << channel
                         << ": expected " << run_numbers.size() << " runs, found "
                         << channel_adc_peak_mean_list[channel_index].size() << " ADC values.";
            continue;
        }
        LOG(INFO) << "Channel " << channel << ":";
        for (size_t run_number_index = 0; run_number_index < run_numbers.size(); run_number_index++) {
            LOG(INFO) << "  Run " << run_numbers[run_number_index] << ": mean ADC peak = " << channel_adc_peak_mean_list[channel_index][run_number_index] << ", sigma = " << channel_adc_peak_sigma_list[channel_index][run_number_index];
        }
    }

    // draw the error graph of mean ADC peak vs laser intensity for each channel
    for (size_t channel_index = 0; channel_index < interested_channels.size(); channel_index++) {
        int channel = interested_channels[channel_index];
        if (channel_adc_peak_mean_list[channel_index].size() != channel_laser_intensity_list.size()
            || channel_adc_peak_sigma_list[channel_index].size() != channel_laser_intensity_list.size()) {
            LOG(WARNING) << "Skipping ADC graph for Channel " << channel
                         << ": ADC and laser-intensity value counts differ.";
            adc_laser_intensity_graphs.push_back(nullptr);
            continue;
        }
        TCanvas *canvas = new TCanvas(("canvas_mean_peak_vs_laser_channel_" + std::to_string(channel)).c_str(), ("Mean ADC Peak vs Laser Intensity - Channel " + std::to_string(channel)).c_str(), 800, 600);
        canvas->cd();
        TGraphErrors *graph = new TGraphErrors();
        int point_index = 0;
        for (size_t i = 0; i < channel_laser_intensity_list.size(); i++) {
            double mean = channel_adc_peak_mean_list[channel_index][i];
            double sigma = channel_adc_peak_sigma_list[channel_index][i];
            // filter out nan and zero value
            if (std::isnan(mean) || std::isnan(sigma)) {
                LOG(WARNING) << "NaN value found for Channel " << channel << ", Laser Intensity " << channel_laser_intensity_list[i] << ": mean = " << channel_adc_peak_mean_list[channel_index][i] << ", sigma = " << channel_adc_peak_sigma_list[channel_index][i] << ". Skipping this point.";
                continue;
            }
            graph->SetPoint(point_index, channel_laser_intensity_list[i], mean);
            graph->SetPointError(point_index, channel_laser_intensity_error_list[i], sigma);
            point_index++;
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
        latex.DrawLatex(0.13, 0.76, execution_time_label.c_str());
        
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
        std::vector<double> tot_trimmed_means;
        std::vector<double> tot_trimmed_mean_errors;
        std::vector<double> tot_laser_intensities;
        std::vector<double> tot_laser_intensity_errors;
        std::vector<TH1D*> hists_to_draw;
        std::vector<TGraphAsymmErrors*> graphs_to_draw;
        
        
        // First pass: compute statistics from counts, then normalize for drawing.
        double max_hist_value = 0;
        for (size_t run_number_index = 0; run_number_index < run_numbers.size(); run_number_index++) {
            double laser_intensity = laser_intensities[run_number_index];
            if (channel_index < channel_tot_th1ds.size() && run_number_index < channel_tot_th1ds[channel_index].size()) {
                TH1D *hist = channel_tot_th1ds[channel_index][run_number_index];
                if (hist && hist->GetEntries() < min_entrys_threshold) {
                    LOG(WARNING) << "Histogram for Run " << run_numbers[run_number_index] << ", Channel " << channel << " has less than " << min_entrys_threshold << " entries (" << hist->GetEntries() << "). Skipping this histogram.";
                    continue;
                }
                if (hist) {
                    hist->SetLineColor(run_number_index + 1);
                    hist->GetXaxis()->SetRangeUser(0, 4096);
                }
                constexpr double trim_fraction = 0.10;
                const TrimmedMeanStats stats = CalculateTrimmedMean(hist, trim_fraction);
                const double histogram_integral = hist
                    ? hist->Integral(1, hist->GetNbinsX())
                    : 0.0;
                if (histogram_integral <= 0.0) {
                    LOG(WARNING) << "Histogram for Run " << run_numbers[run_number_index]
                                 << ", Channel " << channel
                                 << " has zero integral. Skipping this histogram.";
                    continue;
                }
                hist->Scale(1.0 / histogram_integral);
                max_hist_value = std::max(max_hist_value, hist->GetMaximum());
                hists_to_draw.push_back(hist);
                if (!stats.valid) {
                    LOG(WARNING) << "Failed to calculate trimmed mean for Run "
                                 << run_numbers[run_number_index] << ", Channel " << channel
                                 << ". Skipping this histogram.";
                    continue;
                }
                LOG(INFO) << "Run " << run_numbers[run_number_index] << ", Channel " << channel
                          << ": 10% trimmed mean ToT = " << stats.mean
                          << ", Winsorized standard error = " << stats.standard_error
                          << ", q10 = " << stats.q_low << ", q90 = " << stats.q_high;
                TGraphAsymmErrors *graph = new TGraphAsymmErrors(1);
                graph->SetPoint(0, stats.mean, 0);
                graph->SetPointError(
                    0, stats.standard_error, stats.standard_error, 0, 0);
                graph->SetMarkerStyle(20);
                graph->SetMarkerSize(1.0);
                graph->SetMarkerColor(run_number_index + 1);
                graph->SetLineColor(run_number_index + 1);
                graphs_to_draw.push_back(graph);

                tot_trimmed_means.push_back(stats.mean);
                tot_trimmed_mean_errors.push_back(stats.standard_error);
                tot_laser_intensities.push_back(laser_intensity);
                tot_laser_intensity_errors.push_back(0.01); // assume 1% error for the laser intensity
            }
        }
        
        const double y_max = max_hist_value * 1.5;
        
        // Second pass: draw histograms with proper range
        for (size_t i = 0; i < hists_to_draw.size(); i++) {
            TH1D *hist = hists_to_draw[i];
            hist->GetYaxis()->SetRangeUser(0, y_max);
            hist->GetXaxis()->SetRangeUser(0, 4096);
            // add axis label and tick
            hist->GetXaxis()->SetTitle("ToT Value");
            hist->GetYaxis()->SetTitle("Normalized entries");
            hist->GetXaxis()->SetTitleSize(0.05);
            hist->GetYaxis()->SetTitleSize(0.05);
            hist->GetXaxis()->SetLabelSize(0.04);
            hist->GetYaxis()->SetLabelSize(0.04);
            hist->GetXaxis()->SetNdivisions(505);
            hist->GetYaxis()->SetNdivisions(505);
            hist->Draw(i == 0 ? "HIST" : "HIST SAME");
        }
        
        // Third pass: draw trimmed-mean markers near the top of the normalized distributions.
        for (auto graph : graphs_to_draw) {
            double graph_x = 0.0;
            double graph_y = 0.0;
            graph->GetPoint(0, graph_x, graph_y);
            graph->SetPoint(0, graph_x, y_max * 0.9);
            graph->Draw("|> SAME");
        }
        
        canvas->Update();
        canvas->Write();
        canvas->Close();

        TCanvas *canvas_mean_tot_vs_laser = new TCanvas(("canvas_mean_tot_vs_laser_channel_" + std::to_string(channel)).c_str(), ("10% Trimmed Mean ToT vs Laser Intensity - Channel " + std::to_string(channel)).c_str(), 800, 600);
        canvas_mean_tot_vs_laser->cd();
        TGraphAsymmErrors *graph_mean_tot = new TGraphAsymmErrors();
        int point_index = 0;
        for (size_t i = 0; i < tot_trimmed_means.size(); i++) {
            // if (!std::isfinite(tot_trimmed_means[i])) {
            //     continue;
            // }
            graph_mean_tot->SetPoint(point_index, tot_laser_intensities[i], tot_trimmed_means[i]);
            graph_mean_tot->SetPointError(
                point_index,
                tot_laser_intensity_errors[i],
                tot_laser_intensity_errors[i],
                tot_trimmed_mean_errors[i],
                tot_trimmed_mean_errors[i]);
            point_index++;
        }
        graph_mean_tot->SetTitle("");
        graph_mean_tot->GetXaxis()->SetTitle("Laser Intensity");
        graph_mean_tot->GetYaxis()->SetTitle("10% Trimmed Mean ToT");
        // set axis range
        if (!tot_laser_intensities.empty()) {
            graph_mean_tot->GetXaxis()->SetRangeUser(0, *std::max_element(tot_laser_intensities.begin(), tot_laser_intensities.end()) * 1.2);
        }
        if (!tot_trimmed_means.empty()) {
            graph_mean_tot->GetYaxis()->SetRangeUser(0, *std::max_element(tot_trimmed_means.begin(), tot_trimmed_means.end()) * 1.2);
        }
        graph_mean_tot->SetMarkerStyle(20);
        graph_mean_tot->SetMarkerSize(1.0);
        graph_mean_tot->Draw("AEP");

        TGraphAsymmErrors* graph_mean_tot_clone = (TGraphAsymmErrors*)graph_mean_tot->Clone(("graph_tot_vs_laser_channel_" + std::to_string(channel)).c_str());
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
        latex.DrawLatex(0.13, 0.76, execution_time_label.c_str());

        canvas_mean_tot_vs_laser->Update();
        canvas_mean_tot_vs_laser->Write();
        // save as a seprate pdf file
        std::string pdf_output_file = opts.output_file;
        pdf_output_file.replace(pdf_output_file.find(".root"), 5, "_mean_tot_vs_laser_channel_" + std::to_string(channel) + ".pdf");
        canvas_mean_tot_vs_laser->SaveAs(pdf_output_file.c_str());
        canvas_mean_tot_vs_laser->Close();
    }

    TCanvas *canvas_mean_tot_vs_laser_all_channels = new TCanvas(
        "canvas_mean_tot_vs_laser_all_channels",
        "10% Trimmed Mean ToT vs Laser Intensity - All Example Channels",
        1000,
        700);
    canvas_mean_tot_vs_laser_all_channels->cd();
    gPad->SetTicks(1, 1);

    TMultiGraph *mean_tot_vs_laser_all_channels = new TMultiGraph(
        "mean_tot_vs_laser_all_channels",
        ";Laser Intensity;10% Trimmed Mean ToT");
    TLegend *mean_tot_legend = new TLegend(0.14, 0.14, 0.31, 0.35);
    mean_tot_legend->SetFillStyle(0);
    mean_tot_legend->SetBorderSize(0);

    const std::vector<int> channel_colors = {
        kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1,
        kOrange + 7, kCyan + 2, kViolet + 2, kBlack
    };
    int plotted_channel_count = 0;
    for (size_t channel_index = 0; channel_index < interested_channels.size(); channel_index++) {
        TGraphAsymmErrors *graph = tot_laser_intensity_graphs[channel_index];
        if (!graph || graph->GetN() == 0) {
            LOG(WARNING) << "No ToT vs laser-intensity points for Channel "
                         << interested_channels[channel_index] << ". Skipping it in the combined plot.";
            continue;
        }

        const int color = channel_colors[channel_index % channel_colors.size()];
        graph->SetMarkerColor(color);
        graph->SetLineColor(color);
        graph->SetMarkerStyle(20 + static_cast<int>(channel_index % 8));
        graph->SetMarkerSize(1.2);
        mean_tot_vs_laser_all_channels->Add(graph, "P");
        mean_tot_legend->AddEntry(
            graph,
            ("Channel " + std::to_string(interested_channels[channel_index])).c_str(),
            "p");
        plotted_channel_count++;
    }

    if (plotted_channel_count > 0) {
        mean_tot_vs_laser_all_channels->SetMinimum(0);
        mean_tot_vs_laser_all_channels->SetMaximum(4096);
        mean_tot_vs_laser_all_channels->Draw("A P");
        mean_tot_vs_laser_all_channels->GetYaxis()->SetRangeUser(0, 4096);
        mean_tot_vs_laser_all_channels->GetXaxis()->SetLimits(200, 1100);
        mean_tot_legend->Draw();
    } else {
        TLatex empty_label;
        empty_label.SetNDC();
        empty_label.SetTextAlign(22);
        empty_label.SetTextSize(0.05);
        empty_label.DrawLatex(0.5, 0.5, "No valid ToT vs laser-intensity points");
    }

    std::vector<std::string> scan_brief_lines;
    std::istringstream scan_brief_stream(scan_brief);
    std::string scan_brief_part;
    std::string scan_brief_line;
    constexpr size_t scan_brief_line_length = 70;
    while (std::getline(scan_brief_stream, scan_brief_part, ',')) {
        if (!scan_brief_part.empty() && scan_brief_part.front() == ' ') {
            scan_brief_part.erase(0, 1);
        }
        const std::string candidate = scan_brief_line.empty()
            ? scan_brief_part
            : scan_brief_line + ", " + scan_brief_part;
        if (!scan_brief_line.empty() && candidate.size() > scan_brief_line_length) {
            scan_brief_lines.push_back(scan_brief_line);
            scan_brief_line = scan_brief_part;
        } else {
            scan_brief_line = candidate;
        }
    }
    if (!scan_brief_line.empty()) {
        scan_brief_lines.push_back(scan_brief_line);
    }

    TLatex combined_latex;
    combined_latex.SetNDC();
    combined_latex.SetTextAlign(13);
    combined_latex.SetTextFont(62);
    combined_latex.SetTextSize(0.035);
    combined_latex.DrawLatex(0.13, 0.86, "Laser Test with H2GCROC");

    combined_latex.SetTextFont(42);
    combined_latex.SetTextSize(0.024);
    double combined_text_y = 0.82;
    for (const auto& brief_line : scan_brief_lines) {
        combined_latex.DrawLatex(0.13, combined_text_y, brief_line.c_str());
        combined_text_y -= 0.032;
    }

    if (!run_numbers.empty() && !laser_intensities.empty()) {
        const auto intensity_range = std::minmax_element(
            laser_intensities.begin(), laser_intensities.end());
        std::ostringstream scan_range_label;
        scan_range_label << "Runs " << run_numbers.front() << "-" << run_numbers.back()
                         << " (" << run_numbers.size() << "), laser intensity "
                         << *intensity_range.first << "-" << *intensity_range.second;
        combined_latex.DrawLatex(0.13, combined_text_y, scan_range_label.str().c_str());
        combined_text_y -= 0.032;
    }
    combined_latex.DrawLatex(0.13, combined_text_y, execution_time_label.c_str());

    canvas_mean_tot_vs_laser_all_channels->Modified();
    canvas_mean_tot_vs_laser_all_channels->Update();
    canvas_mean_tot_vs_laser_all_channels->Write();
    std::string combined_tot_pdf_file = opts.output_file;
    combined_tot_pdf_file.replace(
        combined_tot_pdf_file.find(".root"),
        5,
        "_mean_tot_vs_laser_all_channels.pdf");
    canvas_mean_tot_vs_laser_all_channels->SaveAs(combined_tot_pdf_file.c_str());
    canvas_mean_tot_vs_laser_all_channels->Close();

    TCanvas *canvas_adc_tot_vs_laser_example_channels = new TCanvas(
        "canvas_adc_tot_vs_laser_example_channels",
        "ADC and ToT vs Laser Intensity - Example Channels",
        1000,
        700);
    canvas_adc_tot_vs_laser_example_channels->cd();
    gPad->SetTicks(1, 0);
    gPad->SetRightMargin(0.14);

    constexpr double adc_axis_min = 0.0;
    constexpr double adc_axis_max = 1100.0;
    constexpr double tot_axis_min = 0.0;
    constexpr double tot_axis_max = 4096.0;
    const double tot_to_adc_scale =
        (adc_axis_max - adc_axis_min) / (tot_axis_max - tot_axis_min);

    TMultiGraph *adc_tot_comparison = new TMultiGraph(
        "adc_tot_vs_laser_example_channels",
        ";Laser Intensity;Mean ADC Peak");
    TLegend *channel_legend = new TLegend(0.14, 0.66, 0.31, 0.86);
    channel_legend->SetFillStyle(0);
    channel_legend->SetBorderSize(0);
    TLegend *quantity_legend = new TLegend(0.70, 0.78, 0.84, 0.89);
    quantity_legend->SetFillStyle(0);
    quantity_legend->SetBorderSize(0);

    TGraph *adc_marker_example = new TGraph();
    adc_marker_example->SetMarkerStyle(20);
    adc_marker_example->SetMarkerColor(kBlack);
    TGraph *tot_marker_example = new TGraph();
    tot_marker_example->SetMarkerStyle(25);
    tot_marker_example->SetMarkerColor(kBlack);
    quantity_legend->AddEntry(adc_marker_example, "ADC (left axis)", "p");
    quantity_legend->AddEntry(tot_marker_example, "ToT (right axis)", "p");

    int comparison_channel_count = 0;
    for (size_t channel_index = 0; channel_index < interested_channels.size(); channel_index++) {
        TGraphErrors *adc_graph = channel_index < adc_laser_intensity_graphs.size()
            ? adc_laser_intensity_graphs[channel_index]
            : nullptr;
        TGraphAsymmErrors *tot_graph = channel_index < tot_laser_intensity_graphs.size()
            ? tot_laser_intensity_graphs[channel_index]
            : nullptr;
        if ((!adc_graph || adc_graph->GetN() == 0) && (!tot_graph || tot_graph->GetN() == 0)) {
            LOG(WARNING) << "No ADC or ToT points for Channel "
                         << interested_channels[channel_index]
                         << ". Skipping it in the ADC/ToT comparison plot.";
            continue;
        }

        const int color = channel_colors[channel_index % channel_colors.size()];
        if (adc_graph && adc_graph->GetN() > 0) {
            adc_graph->SetMarkerColor(color);
            adc_graph->SetLineColor(color);
            adc_graph->SetMarkerStyle(20);
            adc_graph->SetMarkerSize(1.2);
            adc_tot_comparison->Add(adc_graph, "P");
            channel_legend->AddEntry(
                adc_graph,
                ("Channel " + std::to_string(interested_channels[channel_index])).c_str(),
                "p");
        }

        if (tot_graph && tot_graph->GetN() > 0) {
            TGraphAsymmErrors *scaled_tot_graph = new TGraphAsymmErrors();
            scaled_tot_graph->SetName(
                ("graph_scaled_tot_vs_laser_channel_"
                    + std::to_string(interested_channels[channel_index])).c_str());
            for (int point = 0; point < tot_graph->GetN(); point++) {
                double laser_intensity = 0.0;
                double mean_tot = 0.0;
                tot_graph->GetPoint(point, laser_intensity, mean_tot);
                scaled_tot_graph->SetPoint(
                    point,
                    laser_intensity,
                    adc_axis_min + (mean_tot - tot_axis_min) * tot_to_adc_scale);
                scaled_tot_graph->SetPointError(
                    point,
                    tot_graph->GetErrorXlow(point),
                    tot_graph->GetErrorXhigh(point),
                    tot_graph->GetErrorYlow(point) * tot_to_adc_scale,
                    tot_graph->GetErrorYhigh(point) * tot_to_adc_scale);
            }
            scaled_tot_graph->SetMarkerColor(color);
            scaled_tot_graph->SetLineColor(color);
            scaled_tot_graph->SetMarkerStyle(25);
            scaled_tot_graph->SetMarkerSize(1.2);
            adc_tot_comparison->Add(scaled_tot_graph, "P");
            if (!adc_graph || adc_graph->GetN() == 0) {
                channel_legend->AddEntry(
                    scaled_tot_graph,
                    ("Channel " + std::to_string(interested_channels[channel_index])).c_str(),
                    "p");
            }
        }
        comparison_channel_count++;
    }

    if (comparison_channel_count > 0 && !laser_intensities.empty()) {
        adc_tot_comparison->SetMinimum(adc_axis_min);
        adc_tot_comparison->SetMaximum(adc_axis_max);
        adc_tot_comparison->Draw("A P");
        adc_tot_comparison->GetXaxis()->SetLimits(200, 1100);
        adc_tot_comparison->GetYaxis()->SetRangeUser(adc_axis_min, adc_axis_max);
        gPad->Modified();
        gPad->Update();

        TGaxis *tot_axis = new TGaxis(
            gPad->GetUxmax(),
            adc_axis_min,
            gPad->GetUxmax(),
            adc_axis_max,
            tot_axis_min,
            tot_axis_max,
            510,
            "+L");
        tot_axis->SetTitle("10% Trimmed Mean ToT");
        tot_axis->SetTitleOffset(1.25);
        tot_axis->SetLabelFont(adc_tot_comparison->GetYaxis()->GetLabelFont());
        tot_axis->SetLabelSize(adc_tot_comparison->GetYaxis()->GetLabelSize());
        tot_axis->SetTitleFont(adc_tot_comparison->GetYaxis()->GetTitleFont());
        tot_axis->SetTitleSize(adc_tot_comparison->GetYaxis()->GetTitleSize());
        tot_axis->Draw();
        channel_legend->Draw();
        quantity_legend->Draw();
    } else {
        TLatex empty_label;
        empty_label.SetNDC();
        empty_label.SetTextAlign(22);
        empty_label.SetTextSize(0.05);
        empty_label.DrawLatex(0.5, 0.5, "No valid ADC or ToT comparison points");
    }

    canvas_adc_tot_vs_laser_example_channels->Modified();
    canvas_adc_tot_vs_laser_example_channels->Update();
    canvas_adc_tot_vs_laser_example_channels->Write();
    std::string adc_tot_comparison_pdf_file = opts.output_file;
    adc_tot_comparison_pdf_file.replace(
        adc_tot_comparison_pdf_file.find(".root"),
        5,
        "_adc_tot_vs_laser_example_channels.pdf");
    canvas_adc_tot_vs_laser_example_channels->SaveAs(adc_tot_comparison_pdf_file.c_str());
    canvas_adc_tot_vs_laser_example_channels->Close();


    const int overview_columns = std::max(1, static_cast<int>(std::ceil(std::sqrt(interested_channels.size()))));
    const int overview_rows = std::max(1, static_cast<int>((interested_channels.size() + overview_columns - 1) / overview_columns));
    TCanvas *canvas_tot_distributions_all_channels = new TCanvas(
        "canvas_tot_distributions_all_channels",
        "ToT Distributions - All Channels",
        800 * overview_columns,
        600 * overview_rows);
    canvas_tot_distributions_all_channels->Divide(overview_columns, overview_rows);

    for (size_t channel_index = 0; channel_index < interested_channels.size(); channel_index++) {
        const int channel = interested_channels[channel_index];
        canvas_tot_distributions_all_channels->cd(static_cast<int>(channel_index) + 1);
        gPad->SetTicks(1, 1);

        double global_y_max = 0;
        std::vector<std::pair<size_t, TH1D*>> valid_histograms;
        for (size_t run_index = 0; run_index < channel_tot_th1ds[channel_index].size(); run_index++) {
            TH1D *hist = channel_tot_th1ds[channel_index][run_index];
            if (!hist || hist->GetEntries() < min_entrys_threshold) {
                continue;
            }

            global_y_max = std::max(global_y_max, hist->GetMaximum());
            valid_histograms.emplace_back(run_index, hist);
        }

        if (valid_histograms.empty()) {
            TLatex empty_label;
            empty_label.SetNDC();
            empty_label.SetTextAlign(22);
            empty_label.SetTextSize(0.06);
            empty_label.DrawLatex(0.5, 0.5, ("Channel " + std::to_string(channel) + ": no ToT distributions").c_str());
            continue;
        }

        const int legend_columns = std::max(
            1,
            static_cast<int>((valid_histograms.size() + 5) / 6));
        TLegend *legend = new TLegend(0.16, 0.58, 0.88, 0.90);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetNColumns(legend_columns);
        legend->SetTextSize(valid_histograms.size() > 24 ? 0.017 : 0.021);

        for (size_t histogram_index = 0; histogram_index < valid_histograms.size(); histogram_index++) {
            const size_t run_index = valid_histograms[histogram_index].first;
            TH1D *hist = valid_histograms[histogram_index].second;
            hist->SetStats(false);
            hist->SetTitle(("Channel " + std::to_string(channel) + ";ToT;Normalized entries").c_str());
            hist->SetLineColor(static_cast<int>(run_index) + 1);
            hist->SetLineWidth(2);
            hist->GetXaxis()->SetRangeUser(0, 4096);
            hist->GetYaxis()->SetRangeUser(0, global_y_max * 1.5);
            hist->Draw(histogram_index == 0 ? "HIST" : "HIST SAME");

            if (run_index < run_numbers.size() && run_index < laser_intensities.size()) {
                legend->AddEntry(
                    hist,
                    ("R" + std::to_string(run_numbers[run_index]) + ", L="
                        + std::to_string(static_cast<int>(laser_intensities[run_index]))).c_str(),
                    "l");
            }
        }
        legend->Draw();
    }

    output_root->cd();
    canvas_tot_distributions_all_channels->Modified();
    canvas_tot_distributions_all_channels->Update();
    canvas_tot_distributions_all_channels->Write();
    std::string overview_pdf_file = opts.output_file;
    overview_pdf_file.replace(overview_pdf_file.find(".root"), 5, "_tot_distributions_all_channels.pdf");
    canvas_tot_distributions_all_channels->SaveAs(overview_pdf_file.c_str());
    canvas_tot_distributions_all_channels->Close();


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
