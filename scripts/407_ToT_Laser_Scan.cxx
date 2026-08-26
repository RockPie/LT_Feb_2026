#include "H2GCROC_Common.hxx"
#include "H2GCROC_Lib.hxx"
#include "H2GCROC_ToT.hxx"
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
    if (scan_data == "ToT"){
        data_file_prefix = "dump/405_ToT_Scan/ToTScan";
    }

    std::vector<TGraphErrors*> graph_list;
    std::vector<int> graph_sub_config_index_list;
    std::vector<int> graph_channel_list;

    for (int sub_config_index = 0; sub_config_index < scan_configs.size(); sub_config_index++) {
        const auto& sub_config_str = scan_configs[sub_config_index];
        const auto& sub_config_label = scan_config_labels[sub_config_index];
        LOG(INFO) << "Processing sub-config " << sub_config_index << ": " << sub_config_str << " (" << sub_config_label << ")";

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

        TIter next_key(sub_config_result_root->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)next_key())) {
            std::string key_name = key->GetName();
            // canvas_mean_tot_vs_laser_channel_
            if (key_name.find("canvas_mean_tot_vs_laser_channel_") != std::string::npos) {
                TCanvas *canvas = (TCanvas*)sub_config_result_root->Get(key_name.c_str());
                if (canvas) {
                    std::string canvas_name = canvas->GetName();
                    if (canvas_name.find("canvas_mean_tot_vs_laser_channel_") != std::string::npos) {
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
                        std::string channel_str = canvas_name.substr(std::string("canvas_mean_tot_vs_laser_channel_").size());
                        int channel = std::stoi(channel_str);
                        graph_channel_list.push_back(channel);
                        // Do something with the canvas, e.g. save it to the output root file
                        output_root->cd();
                        canvas->Write((sub_config_label + "_" + key_name).c_str());
                    }
                } // end of processing tot plot
            }
        } // end of loop over all objects in the sub-config result root file

        sub_config_result_root->Close();
    } // end of loop over sub-configs

    // ! === Data analysis ===
    const std::vector<int> color_wheel = {
        kRed, kBlue, kGreen + 2, kMagenta, kCyan + 2, kOrange + 7, kViolet + 2, kTeal + 2
    };
    output_root->cd();

    if (graph_list.empty()) {
        LOG(WARNING) << "No graphs were loaded for ToT scan. Check input scan configs and files.";
        output_root->Close();
        return 0;
    }

    auto canvas_multiple_graph = new TCanvas("canvas_multiple_graph", "Multiple Graphs", 1000, 600);
    auto multigraph = new TMultiGraph("", ";Laser Intensity [a.u.];Mean ToT");
    auto multigraph_legend = new TLegend(0.5, 0.77, 0.89, 0.89);
    multigraph_legend->SetFillStyle(0);
    multigraph_legend->SetBorderSize(0);
    for (size_t i = 0; i < graph_list.size(); i++) {
        TGraphErrors *graph = graph_list[i];
        int sub_config_index = graph_sub_config_index_list[i];
        int color = color_wheel[sub_config_index % color_wheel.size()];
        graph->SetMarkerColor(color);
        graph->SetLineColor(color);
        multigraph->Add(graph, "P");
        std::string legend_entry = "Chn " + std::to_string(graph_channel_list[i]) + " (" + scan_config_labels[sub_config_index] + ")";
        multigraph_legend->AddEntry(graph, legend_entry.c_str(), "P");
    }
    // set y axis range to [0, 4096]
    multigraph->SetMinimum(0);
    multigraph->SetMaximum(4096);
    multigraph->Draw("A");
    // draw the legend in three columns
    multigraph_legend->SetNColumns(4);
    multigraph_legend->Draw();
    // Write info
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextFont(62);
    double text_x = 0.13;
    double text_y = 0.85;
    double text_y_step = 0.045;
    latex.DrawLatex(text_x, text_y, "Laser Test with H2GCROC");
    latex.SetTextSize(0.03);
    latex.SetTextFont(42);
    latex.DrawLatex(text_x, text_y - text_y_step, scan_brief.c_str());
    latex.DrawLatex(text_x, text_y - 2 * text_y_step, "February 2026");
    // save as a separate pdf file
    canvas_multiple_graph->SaveAs((opts.output_file + "_multiple_graph.pdf").c_str());
    canvas_multiple_graph->Write();
    canvas_multiple_graph->Close();

    // Do the scale similarity analysis
    const double ssa_x0 = 5.65;
    const double ssa_sMin = 0.005;
    const double ssa_sMax = 1;
    const int ssa_nPoints = 10000;
    auto channel_number_set = std::set<int>(graph_channel_list.begin(), graph_channel_list.end());

    std::vector<double> root_save_channel_numbers;
    std::vector<double> root_save_bias_voltages;
    std::vector<double> root_save_similarity_values;

    for (int channel : channel_number_set) {
        auto canvas_ssa = new TCanvas(("canvas_ssa_channel_" + std::to_string(channel)).c_str(), ("Scale Similarity Analysis - Channel " + std::to_string(channel)).c_str(), 1000, 600);
        std::vector<TGraph*> ssa_graphs;
        std::vector<TGraphErrors*> ssa_toa_graphs;
        std::vector<std::string> ssa_graph_scan_variable_list;
        // first round to find the graph with the maximum average y value as the reference graph
        int reference_graph_index = -1;
        double max_average_y = -1;
        for (size_t i = 0; i < graph_list.size(); i++){
            if (graph_channel_list[i] == channel) {
                TGraphErrors *graph = graph_list[i];
                double average_y = 0;
                int n_points = graph->GetN();
                for (int j = 0; j < n_points; j++) {
                    average_y += graph->GetY()[j];
                }
                average_y /= n_points;
                if (average_y > max_average_y) {
                    max_average_y = average_y;
                    reference_graph_index = i;
                }
            }
        } // end of loop to find reference graph
        if (reference_graph_index == -1) {
            LOG(WARNING) << "Failed to find reference graph for channel " << channel << ". Skipping scale similarity analysis for this channel.";
            continue;
        }
        TGraphErrors *reference_graph = (TGraphErrors*)graph_list[reference_graph_index]->Clone();
        // second round to calculate the scale factor for each graph by comparing with the reference graph
        double ssa_graphs_y_max = -1;
        for (size_t i = 0; i < graph_list.size(); i++){
            if (graph_channel_list[i] == channel) {
                if (i == reference_graph_index) {
                    LOG(INFO) << "Graph " << graph_list[i]->GetName() << " is the reference graph for channel " << channel << ". Scale factor = 1.";
                    continue;
                }
                // Skip 43V early to avoid index mismatch
                if (scan_config_labels[graph_sub_config_index_list[i]].find("43") != std::string::npos) {
                    LOG(INFO) << "Skipping 43 V graph early to avoid alignment issues.";
                    continue;
                }
                TGraphErrors *graph = graph_list[i];
                TGraphErrors *graph_clone = (TGraphErrors*)graph->Clone();
                ssa_toa_graphs.push_back(graph_clone);
                auto graphSim = ScanScaleSimilarity(reference_graph, graph, ssa_x0, ssa_sMin, ssa_sMax, ssa_nPoints, true);
                if (!graphSim) {
                    LOG(WARNING) << "Failed to compute scale similarity for graph " << graph->GetName() << ". Skipping this graph.";
                    continue;
                }
                for (int j = 0; j < graphSim->GetN(); j++) {
                    if (graphSim->GetY()[j] > ssa_graphs_y_max) {
                        ssa_graphs_y_max = graphSim->GetY()[j];
                    }
                }
                TGraph *graphSim_clone = (TGraph*)graphSim->Clone();
                graphSim_clone->SetName((std::string(graph->GetName()) + "_sim").c_str());
                graphSim_clone->SetTitle("");
                graphSim_clone->GetXaxis()->SetTitle("Scale Factor");
                graphSim_clone->GetYaxis()->SetTitle("Similarity");
                ssa_graphs.push_back(graphSim_clone);
                ssa_graph_scan_variable_list.push_back(scan_config_labels[graph_sub_config_index_list[i]]);
            }
        } // end of loop to calculate scale factor

        // set log x and log y
        canvas_ssa->SetLogx();
        canvas_ssa->SetLogy();
        auto ssa_legend = TLegend(0.45, 0.75, 0.89, 0.89);
        ssa_legend.SetFillStyle(0);
        ssa_legend.SetBorderSize(0);
        ssa_legend.SetTextSize(0.02);
        std::vector<double> ssa_optimal_scale_factors;
        for (size_t i = 0; i < ssa_graphs.size(); i++) {
            ssa_graphs[i]->SetMarkerColor(color_wheel[i % color_wheel.size()]);
            ssa_graphs[i]->SetLineColor(color_wheel[i % color_wheel.size()]);
            ssa_graphs[i]->SetMinimum(1);
            ssa_graphs[i]->SetMaximum(ssa_graphs_y_max * 10);
            ssa_graphs[i]->Draw(i == 0 ? "AP" : "P SAME");
            // find the minmum point of the graph and print the corresponding scale factor and similarity value
            double min_y = 1e9;
            double min_x = 1;
            for (int j = 0; j < ssa_graphs[i]->GetN(); j++) {
                if (ssa_graphs[i]->GetY()[j] < min_y) {
                    min_y = ssa_graphs[i]->GetY()[j];
                    min_x = ssa_graphs[i]->GetX()[j];
                }
            }
            ssa_optimal_scale_factors.push_back(min_x);
            LOG(INFO) << "Graph " << ssa_graphs[i]->GetName() << ": minimum similarity = " << min_y << " at scale factor = " << min_x;
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(3) << (1.0 / min_x);
            std::string legend_entry = ssa_graph_scan_variable_list[i] + " (X_{Smin}=1/" + oss.str() + ")";
            ssa_legend.AddEntry(ssa_graphs[i], legend_entry.c_str(), "L");

            int bias_voltage = -1;
            // label might have V and "43 V"
            std::smatch voltage_match;
            if (std::regex_search(ssa_graph_scan_variable_list[i], voltage_match, std::regex("(\\d+)\\s*V"))) {
                bias_voltage = std::stoi(voltage_match[1].str());
                root_save_channel_numbers.push_back(channel);
                root_save_bias_voltages.push_back(bias_voltage);
                root_save_similarity_values.push_back(1.0/min_x);
            } else {
                LOG(WARNING) << "Failed to parse bias voltage from graph scan variable label: " << ssa_graph_scan_variable_list[i] << ". Expected format: something like '43 V'. Skipping saving bias voltage for this graph.";
            }
        }
        ssa_legend.SetNColumns(3);
        ssa_legend.Draw();
        TLatex latex_ssa;
        latex_ssa.SetNDC();
        latex_ssa.SetTextSize(0.04);
        latex_ssa.SetTextFont(62);
        double text_x_ssa = 0.13;
        double text_y_ssa = 0.85;
        double text_y_step_ssa = 0.045;
        latex_ssa.DrawLatex(text_x_ssa, text_y_ssa, "Laser Test with H2GCROC");
        latex_ssa.SetTextSize(0.03);
        latex_ssa.SetTextFont(42);
        latex_ssa.DrawLatex(text_x_ssa, text_y_ssa - text_y_step_ssa, scan_brief.c_str());
        latex_ssa.DrawLatex(text_x_ssa, text_y_ssa - 2 * text_y_step_ssa, "February 2026"); 
        canvas_ssa->Write();
        // save as a separate pdf file
        canvas_ssa->SaveAs((opts.output_file + "_ssa_channel_" + std::to_string(channel) + ".pdf").c_str());
        canvas_ssa->Close();

        auto canvas_ssa_scaled = new TCanvas(("canvas_ssa_scaled_channel_" + std::to_string(channel)).c_str(), ("Scale Similarity Analysis with Optimal Scaling - Channel " + std::to_string(channel)).c_str(), 1000, 600);
        auto ssa_scaled_legend = TLegend(0.5, 0.11, 0.89, 0.4);
        ssa_scaled_legend.SetFillStyle(0);
        ssa_scaled_legend.SetBorderSize(0);
        // draw the reference graph
        reference_graph->SetMarkerColor(kBlack);
        reference_graph->SetMarkerStyle(20);
        reference_graph->SetMarkerSize(1.0);
        reference_graph->SetLineColor(kBlack);
        reference_graph->SetLineWidth(2);
        reference_graph->Draw("APE");
        ssa_scaled_legend.AddEntry(reference_graph, "Reference (54 V)", "PE");
        for (size_t i = 0; i < ssa_toa_graphs.size(); i++) {
            auto tot_graph = ssa_toa_graphs[i];
            double optimal_scale_factor = ssa_optimal_scale_factors[i];
            auto graph_scaled = ScaleX_TGraphErrors(tot_graph, optimal_scale_factor, ssa_x0, (std::string(tot_graph->GetName()) + "_scaled").c_str());
            graph_scaled->SetMarkerColor(color_wheel[i % color_wheel.size()]);
            graph_scaled->SetMarkerStyle(20);
            graph_scaled->SetMarkerSize(1.0);
            graph_scaled->SetLineColor(color_wheel[i % color_wheel.size()]);
            graph_scaled->SetLineWidth(2);
            graph_scaled->Draw("P SAME");
            std::ostringstream oss_scaled;
            oss_scaled << std::fixed << std::setprecision(3) << (1.0 / optimal_scale_factor);
            std::string legend_entry = ssa_graph_scan_variable_list[i] + " scaled (X_{Smin}=1/" + oss_scaled.str() + ")";
            ssa_scaled_legend.AddEntry(graph_scaled, legend_entry.c_str(), "PE");
        }
        ssa_scaled_legend.SetNColumns(2);
        ssa_scaled_legend.Draw();
        TLatex latex_ssa_scaled;
        latex_ssa_scaled.SetNDC();
        latex_ssa_scaled.SetTextSize(0.04);
        latex_ssa_scaled.SetTextFont(62);
        double text_x_ssa_scaled = 0.13;
        double text_y_ssa_scaled = 0.85;
        double text_y_step_ssa_scaled = 0.045;
        latex_ssa_scaled.DrawLatex(text_x_ssa_scaled, text_y_ssa_scaled, "Laser Test with H2GCROC");
        latex_ssa_scaled.SetTextSize(0.03);
        latex_ssa_scaled.SetTextFont(42);
        std::string scaled_scan_brief = "Scaled ToT - Laser Intensity, Channel " + std::to_string(channel);
        latex_ssa_scaled.DrawLatex(text_x_ssa_scaled, text_y_ssa_scaled - text_y_step_ssa_scaled, scaled_scan_brief.c_str());
        latex_ssa_scaled.DrawLatex(text_x_ssa_scaled, text_y_ssa_scaled - 2 * text_y_step_ssa_scaled, "February 2026"); 
        canvas_ssa_scaled->Write();
        // save as a separate pdf file
        canvas_ssa_scaled->SaveAs((opts.output_file + "_ssa_scaled_channel_" + std::to_string(channel) + ".pdf").c_str());
        canvas_ssa_scaled->Close();

    }

    // save the SSA results to TVectorD in the output root file
    TVectorD channel_numbers(root_save_channel_numbers.size());
    TVectorD bias_voltages(root_save_bias_voltages.size());
    TVectorD similarity_values(root_save_similarity_values.size());
    for (size_t i = 0; i < root_save_channel_numbers.size(); i++) {
        channel_numbers[i] = root_save_channel_numbers[i];
        bias_voltages[i] = root_save_bias_voltages[i];
        similarity_values[i] = root_save_similarity_values[i];
    }
    channel_numbers.Write("ssa_channel_numbers");
    bias_voltages.Write("ssa_bias_voltages");
    similarity_values.Write("ssa_similarity_values");

    output_root->Close();
    return 0;
}