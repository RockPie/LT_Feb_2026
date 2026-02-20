#include "H2GCROC_Common.hxx"
#include "H2GCROC_Lib.hxx"
#include "TKey.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include <iomanip>
#include <regex>
#include <sstream>

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
    const auto& scan_data = scan_json["scan_data"].get<std::string>();
    const auto& scan_configs = scan_json["scan_configs"].get<std::vector<std::string>>();
    const auto& scan_config_labels = scan_json["scan_config_labels"].get<std::vector<std::string>>();

    // laser_scan_0.json
    std::string scan_info_str = "Laser Scan";
    std::smatch scan_match;
    if (std::regex_search(input_scan_json, scan_match, std::regex("laser_scan_(\\d+)"))) {
        scan_info_str += " " + scan_match[1].str();
    } else {
        scan_info_str += " Unknown";
    }

    TFile *output_root = new TFile(opts.output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << opts.output_file;
        return 1;
    }

    std::string data_file_prefix = "";
    if (scan_data == "ADC"){
        data_file_prefix = "dump/402_ADC_Scan/Scan";
    }

    std::vector<TGraphErrors*> graph_list;
    std::vector<int> graph_sub_config_index_list;
    std::vector<TF1*> fit_list;

    for (int sub_config_index = 0; sub_config_index < scan_configs.size(); sub_config_index++) {
        const auto& sub_config_str = scan_configs[sub_config_index];
        const auto& sub_config_label = scan_config_labels[sub_config_index];
        LOG(INFO) << "Processing sub-config " << sub_config_index << ": " << sub_config_str << " (" << sub_config_label << ")";

        // config/scan_number_0.json

        int sub_config_number = -1;
        size_t sub_config_number_pos = sub_config_str.find("scan_number_");
        if (sub_config_number_pos != std::string::npos) {
            sub_config_number = std::stoi(sub_config_str.substr(sub_config_number_pos + std::string("scan_number_").size()));
        } else {
            LOG(WARNING) << "Failed to parse sub-config number from " << sub_config_str << ". Expected format: scan_number_X.json where X is the sub-config number. Skipping this sub-config.";
            continue;
        }

        std::string sub_config_result_root_file = data_file_prefix + std::to_string(sub_config_number) + ".root";

        LOG(INFO) << "Opening sub-config result file: " << sub_config_result_root_file;
        TFile *sub_config_result_root = TFile::Open(sub_config_result_root_file.c_str(), "READ");
        if (!sub_config_result_root || sub_config_result_root->IsZombie()) {
            LOG(ERROR) << "Failed to open sub-config result file " << sub_config_result_root_file;
            continue;
        }

        // get the keys and find "canvas_mean_peak_vs_laser_channel_"
        TIter next_key(sub_config_result_root->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)next_key())) {
            std::string key_name = key->GetName();
            if (key_name.find("canvas_mean_peak_vs_laser_channel_") != std::string::npos) {
                LOG(INFO) << "Found canvas key: " << key_name;
                TCanvas *canvas = (TCanvas*)sub_config_result_root->Get(key_name.c_str());
                if (canvas) {
                    std::string canvas_name = canvas->GetName();
                    if (canvas_name.find("canvas_mean_peak_vs_laser_channel_") != std::string::npos) {
                        // get the TErrorGraph from the canvas
                        TIter next_primitive(canvas->GetListOfPrimitives());
                        TGraphErrors *graph = nullptr;
                        while (auto primitive = next_primitive()) {
                            if (std::string(primitive->ClassName()) == "TGraphErrors") {
                                graph = (TGraphErrors*)primitive;
                                TGraphErrors *graph_clone = (TGraphErrors*)graph->Clone(canvas_name.c_str());
                                graph_list.push_back(graph_clone);
                                graph_sub_config_index_list.push_back(sub_config_index);
                                break;
                            }
                        }
                    } else {
                        LOG(WARNING) << "Canvas name " << canvas_name << " does not match expected pattern. Skipping.";
                    }
                } else {
                    LOG(WARNING) << "Failed to retrieve canvas from key " << key_name;
                }
                canvas->Close();
            }
        }

        sub_config_result_root->Close();
    }

    output_root->cd();
    auto format_value = [](double value) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3) << value;
        return oss.str();
    };
    TCanvas *combined_canvas = new TCanvas("combined_canvas", "Combined Laser Scan Results", 1000, 600);
    combined_canvas->cd();
    TMultiGraph *multi_graph = new TMultiGraph("multi_graph", ";Laser Intensity;Fitted ADC Mean");
    const std::vector<int> color_wheel = {
        kRed, kBlue, kGreen + 2, kMagenta, kCyan + 2, kOrange + 7, kViolet + 2, kTeal + 2
    };
    for (size_t i = 0; i < graph_list.size(); i++) {
        TGraphErrors *graph = graph_list[i];
        graph->SetMarkerStyle(20 + i);
        int color = color_wheel[i % color_wheel.size()];
        graph->SetMarkerColor(color);
        graph->SetLineColor(color);
        multi_graph->Add(graph, "P");
        // do the linear fit over the range of y > 200 and y < 950 for each graph
        const int min_fit_points = 3;
        const double min_x_span = 1e-6;
        double fit_min_x = std::numeric_limits<double>::max();
        double fit_max_x = std::numeric_limits<double>::lowest();
        int fit_point_count = 0;
        for (int point_index = 0; point_index < graph->GetN(); point_index++) {
            double x, y;
            graph->GetPoint(point_index, x, y);
            if (std::isfinite(x) && std::isfinite(y) && y > 200 && y < 950) {
                fit_min_x = std::min(fit_min_x, x);
                fit_max_x = std::max(fit_max_x, x);
                fit_point_count++;
            }
        }
        if (fit_point_count >= min_fit_points && (fit_max_x - fit_min_x) > min_x_span) {
            TF1 *fit = new TF1((std::string(graph->GetName()) + "_fit").c_str(), "pol1", fit_min_x, fit_max_x);
            fit->SetLineColor(color);
            graph->Fit(fit, "R");
            fit_list.push_back(fit);
        } else {
            LOG(WARNING) << "Graph " << graph->GetName() << " does not have enough points or x-span for fit (y > 200 and y < 950). Skipping fit.";
            fit_list.push_back(nullptr);
        }
    }

    std::vector<double> linear_fit_slopes;
    std::vector<double> linear_fit_x_intercepts_at_100;
    std::vector<double> linear_fit_slope_errors;
    std::vector<double> linear_fit_x_intercept_at_100_errors;
    std::vector<int> linear_fit_sub_config_indices;
    std::vector<int> linear_fit_channel_numbers;
    multi_graph->Draw("A");
    TLegend *legend = new TLegend(0.55, 0.12, 0.90, 0.35);
    for (size_t i = 0; i < graph_list.size(); i++) {
        TGraphErrors *graph = graph_list[i];
        std::string graph_name = graph->GetName();
        // expected graph name format:
        // canvas_mean_peak_vs_laser_channel_<chn>
        const std::string channel_prefix = "canvas_mean_peak_vs_laser_channel_";
        size_t channel_pos = graph_name.find(channel_prefix);
        // if channel is found, extract it; otherwise, use graph name
        if (channel_pos != std::string::npos) {
            size_t channel_start = channel_pos + channel_prefix.size();
            int channel = std::stoi(graph_name.substr(channel_start));
            int sub_config_index = (i < graph_sub_config_index_list.size()) ? graph_sub_config_index_list[i] : -1;
            std::string legend_entry = "Channel " + std::to_string(channel);
            if (sub_config_index >= 0 && sub_config_index < static_cast<int>(scan_config_labels.size())) {
                // legend_entry += " Sub-config " + std::to_string(sub_config_number);
                legend_entry += " (" + scan_config_labels[sub_config_index] + ")";
            } else if (sub_config_index >= 0) {
                LOG(WARNING) << "Sub-config index " << sub_config_index << " exceeds scan_config_labels size. Using default legend entry.";
            }
            if (i < fit_list.size() && fit_list[i]) {
                double slope = fit_list[i]->GetParameter(1);
                double intercept = fit_list[i]->GetParameter(0);
                double slope_error = fit_list[i]->GetParError(1);
                double intercept_error = fit_list[i]->GetParError(0);
                // x_intercept_at_100
                double x_intercept_at_100 = (slope != 0.0) ? ((100.0 - intercept) / slope) : 0.0;
                double x_intercept_at_100_error = (slope != 0.0) ? std::sqrt(std::pow(intercept_error / slope, 2) + std::pow((intercept * slope_error) / (slope * slope), 2)) : 0.0;
                // double x_intercept = (slope != 0.0) ? (-intercept / slope) : 0.0;
                // legend_entry += ", slope=" + format_value(slope) + ", x-intercept=" + format_value(x_intercept_at_100);
            }
            legend->AddEntry(graph, legend_entry.c_str(), "P");
            linear_fit_slopes.push_back((fit_list[i]) ? fit_list[i]->GetParameter(1) : std::numeric_limits<double>::quiet_NaN());
            linear_fit_x_intercepts_at_100.push_back((fit_list[i]) ? ((fit_list[i]->GetParameter(1) != 0.0) ? ((100.0 - fit_list[i]->GetParameter(0)) / fit_list[i]->GetParameter(1)) : 0.0) : std::numeric_limits<double>::quiet_NaN());
            linear_fit_slope_errors.push_back((fit_list[i]) ? fit_list[i]->GetParError(1) : std::numeric_limits<double>::quiet_NaN());
            linear_fit_x_intercept_at_100_errors.push_back((fit_list[i]) ? ((fit_list[i]->GetParameter(1) != 0.0) ? std::sqrt(std::pow(fit_list[i]->GetParError(0) / fit_list[i]->GetParameter(1), 2) + std::pow((fit_list[i]->GetParameter(0) * fit_list[i]->GetParError(1)) / (fit_list[i]->GetParameter(1) * fit_list[i]->GetParameter(1)), 2)) : 0.0) : std::numeric_limits<double>::quiet_NaN());
            linear_fit_sub_config_indices.push_back(sub_config_index);
            linear_fit_channel_numbers.push_back(channel);
        } else {
            LOG(WARNING) << "Graph name " << graph_name << " does not match expected pattern. Using graph name as legend entry.";
            std::string legend_entry = graph_name;
            if (i < fit_list.size() && fit_list[i]) {
                double slope = fit_list[i]->GetParameter(1);
                double intercept = fit_list[i]->GetParameter(0);
                double x_intercept_at_100 = (slope != 0.0) ? ((100.0 - intercept) / slope) : 0.0;
                // double x_intercept = (slope != 0.0) ? (-intercept / slope) : 0.0;
                legend_entry += ", slope=" + format_value(slope) + ", x-intercept=" + format_value(x_intercept_at_100);
            }
            legend->AddEntry(graph, legend_entry.c_str(), "P");
        }
    }
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetNColumns(2);
    legend->Draw();

    combined_canvas->Modified();
    combined_canvas->Update();
    combined_canvas->Write();
    // save as a seprate pdf file
    combined_canvas->SaveAs((opts.output_file + "_combined.pdf").c_str());
    combined_canvas->Close();

    auto canvas_slopes = new TCanvas("canvas_slopes", "Linear Fit Slopes", 1000, 600);
    canvas_slopes->cd();
    // adjust margin
    canvas_slopes->SetLeftMargin(0.15);
    canvas_slopes->SetRightMargin(0.15);
    TGraphErrors *slope_graph = new TGraphErrors(linear_fit_slopes.size());
    bool has_cc_label = false;
    bool has_bias_label = false;
    for (const auto& label : scan_config_labels) {
        has_cc_label = has_cc_label || (label.find("CC") != std::string::npos);
        has_bias_label = has_bias_label || (label.find("V") != std::string::npos);
    }
    std::string slope_x_title = "Bias Voltage (V)";
    std::string slope_scan_label = "Bias Voltage Scan";
    std::string right_axis_title = "Relative slope to lowest bias";
    if (has_cc_label && !has_bias_label) {
        slope_x_title = "Current Conveyor (CC)";
        slope_scan_label = "Current Conveyor Scan";
        right_axis_title = "Relative slope to lowest CC";
    }
    auto parse_label_value = [](const std::string& label, const std::string& key) -> double {
        size_t key_pos = label.find(key);
        if (key_pos == std::string::npos) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        std::string value_str = label.substr(key_pos + key.size());
        std::istringstream iss(value_str);
        double value = std::numeric_limits<double>::quiet_NaN();
        iss >> value;
        return value;
    };

    for (size_t i = 0; i < linear_fit_slopes.size(); i++) {
        int sub_config_index = linear_fit_sub_config_indices[i];
        int channel_number = linear_fit_channel_numbers[i];
        //  "config/scan_number_0.json"
        if (sub_config_index >= 0 && sub_config_index < static_cast<int>(scan_config_labels.size())) {
            std::string sub_label = scan_config_labels[sub_config_index];
            // if it is bias scan, then it is "XX V"
            if (sub_label.find("V") != std::string::npos) {
                double bias_voltage = parse_label_value(sub_label, "");
                if (!std::isfinite(bias_voltage)) {
                    LOG(WARNING) << "Failed to parse bias voltage from label: " << sub_label;
                    continue;
                }
                slope_graph->SetPoint(i, bias_voltage, linear_fit_slopes[i]);
                slope_graph->SetPointError(i, 0.01, linear_fit_slope_errors[i]);
            } 
    //             "scan_config_labels":[
    //     "CC 15",
    //     "CC 8",
    //     "CC 2"
    // ]
            if (sub_label.find("CC") != std::string::npos) {
                double cc_value = parse_label_value(sub_label, "CC");
                if (!std::isfinite(cc_value)) {
                    LOG(WARNING) << "Failed to parse CC value from label: " << sub_label;
                    continue;
                }
                slope_graph->SetPoint(i, cc_value, linear_fit_slopes[i]);
                slope_graph->SetPointError(i, 0.01, linear_fit_slope_errors[i]);
            }
        } else {
            LOG(WARNING) << "Sub-config index " << sub_config_index << " not found in scan_config_labels. Skipping.";
            continue;
        }
    }
    slope_graph->SetTitle((";" + slope_x_title + ";Slope (ADC / Laser Intensity)").c_str());

    // zoom by setting y-axis range to better visualize the slope distribution
    slope_graph->GetYaxis()->SetRangeUser(0, *std::max_element(linear_fit_slopes.begin(), linear_fit_slopes.end()) * 1.5);
    slope_graph->SetMarkerStyle(20);
    // add plot information
    slope_graph->Draw("AP");
    canvas_slopes->Update();

    // add right axis with relative change to the lowest slope
    double min_slope = std::numeric_limits<double>::max();
    for (double slope : linear_fit_slopes) {
        if (std::isfinite(slope) && slope > 0.0) {
            min_slope = std::min(min_slope, slope);
        }
    }
    if (min_slope != std::numeric_limits<double>::max()) {
        double y_min = slope_graph->GetYaxis()->GetXmin();
        double y_max = slope_graph->GetYaxis()->GetXmax();
        double rel_min = (y_min / min_slope);
        double rel_max = (y_max / min_slope);
        TGaxis *right_axis = new TGaxis(gPad->GetUxmax(), y_min, gPad->GetUxmax(), y_max,
                                        rel_min, rel_max, 510, "+L");
        right_axis->SetTitle(right_axis_title.c_str());
        right_axis->SetTitleOffset(1.2);
        right_axis->SetLabelSize(slope_graph->GetYaxis()->GetLabelSize());
        right_axis->SetTitleSize(slope_graph->GetYaxis()->GetTitleSize());
        // set font
        right_axis->SetLabelFont(slope_graph->GetYaxis()->GetLabelFont());
        right_axis->SetTitleFont(slope_graph->GetYaxis()->GetTitleFont());
        right_axis->Draw();
    }

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextAlign(13);
    latex->SetTextSize(0.04);
    latex->SetTextFont(62);
    latex->DrawLatex(0.18, 0.88, "Laser Test Setup");
    latex->SetTextSize(0.03);
    latex->SetTextFont(42);
    latex->DrawLatex(0.18, 0.84, slope_scan_label.c_str());
    latex->DrawLatex(0.18, 0.80, "Hamamatsu S14160-6010PS");
    latex->DrawLatex(0.18, 0.76, "CERN, February 2026");
    canvas_slopes->Modified();
    canvas_slopes->Update();
    canvas_slopes->Write();
    // save as a seprate pdf file
    canvas_slopes->SaveAs((opts.output_file + "_slopes.pdf").c_str());
    canvas_slopes->Close();

    output_root->Close();

    return 0;
}