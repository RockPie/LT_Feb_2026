#include "H2GCROC_Common.hxx"
#include "H2GCROC_Lib.hxx"
#include "TKey.h"
#include "TMultiGraph.h"
#include <iomanip>
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
    size_t scan_number_pos = scan_brief.find("laser_scan_");
    int scan_number = -1;
    std::string scan_info_str = "Laser Scan";
    if (scan_number_pos != std::string::npos) {
        scan_number = std::stoi(input_scan_json.substr(scan_number_pos + 1 + std::string("laser_scan_").size(), 1));
        scan_info_str += " " + std::to_string(scan_number);
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
                                TGraphErrors *graph_clone = (TGraphErrors*)graph->Clone((std::string(graph->GetName()) + "_clone").c_str());
                                graph_list.push_back(graph_clone);
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
    TCanvas *combined_canvas = new TCanvas("combined_canvas", "Combined Laser Scan Results", 800, 600);
    combined_canvas->cd();
    TMultiGraph *multi_graph = new TMultiGraph("multi_graph", "Mean ADC Peak vs Laser Channel;Laser Channel;Mean ADC Peak");
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
        // do the linear fit over the range of y > 150 and y < 950 for each graph
        double fit_min_x = std::numeric_limits<double>::max();
        double fit_max_x = std::numeric_limits<double>::lowest();
        for (int point_index = 0; point_index < graph->GetN(); point_index++) {
            double x, y;
            graph->GetPoint(point_index, x, y);
            if (y > 200 && y < 950) {
                fit_min_x = std::min(fit_min_x, x);
                fit_max_x = std::max(fit_max_x, x);
            }
        }
        if (fit_min_x != std::numeric_limits<double>::max() && fit_max_x != std::numeric_limits<double>::lowest()) {
            TF1 *fit = new TF1((std::string(graph->GetName()) + "_fit").c_str(), "pol1", fit_min_x, fit_max_x);
            fit->SetLineColor(color);
            graph->Fit(fit, "R");
            fit_list.push_back(fit);
        } else {
            LOG(WARNING) << "Graph " << graph->GetName() << " does not have points in the fit range (y > 150 and y < 950). Skipping fit.";
            fit_list.push_back(nullptr);
        }
    }
    multi_graph->Draw("A");
    TLegend *legend = new TLegend(0.7, 0.5, 0.89, 0.89);
    for (size_t i = 0; i < graph_list.size(); i++) {
        TGraphErrors *graph = graph_list[i];
        std::string graph_name = graph->GetName();
        // expected graph name format: graph_mean_peak_vs_laser_channel_scan_number_X_sub_config_Y_clone
        size_t scan_number_pos = graph_name.find("scan_number_");
        size_t sub_config_pos = graph_name.find("sub_config_");
        if (scan_number_pos != std::string::npos && sub_config_pos != std::string::npos) {
            int scan_number = std::stoi(graph_name.substr(scan_number_pos + std::string("scan_number_").size(), sub_config_pos - (scan_number_pos + std::string("scan_number_").size())));
            int sub_config_number = std::stoi(graph_name.substr(sub_config_pos + std::string("sub_config_").size(), graph_name.find("_clone") - (sub_config_pos + std::string("sub_config_").size())));
            std::string legend_entry = "Scan " + std::to_string(scan_number) + " Sub-config " + std::to_string(sub_config_number);
            if (sub_config_number < scan_config_labels.size()) {
                legend_entry += " (" + scan_config_labels[sub_config_number] + ")";
            } else {
                LOG(WARNING) << "Sub-config number " << sub_config_number << " exceeds scan_config_labels size. Using default legend entry.";
            }
            if (i < fit_list.size() && fit_list[i]) {
                double slope = fit_list[i]->GetParameter(1);
                double intercept = fit_list[i]->GetParameter(0);
                double x_intercept = (slope != 0.0) ? (-intercept / slope) : 0.0;
                legend_entry += ", slope=" + format_value(slope) + ", x-intercept=" + format_value(x_intercept);
            }
            legend->AddEntry(graph, legend_entry.c_str(), "P");
        } else {
            LOG(WARNING) << "Graph name " << graph_name << " does not match expected pattern. Using graph name as legend entry.";
            std::string legend_entry = graph_name;
            if (i < fit_list.size() && fit_list[i]) {
                double slope = fit_list[i]->GetParameter(1);
                double intercept = fit_list[i]->GetParameter(0);
                double x_intercept = (slope != 0.0) ? (-intercept / slope) : 0.0;
                legend_entry += ", slope=" + format_value(slope) + ", x-intercept=" + format_value(x_intercept);
            }
            legend->AddEntry(graph, legend_entry.c_str(), "P");
        }
    }
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    combined_canvas->Modified();
    combined_canvas->Update();
    combined_canvas->Write();
    combined_canvas->Close();

    output_root->Close();

    return 0;
}