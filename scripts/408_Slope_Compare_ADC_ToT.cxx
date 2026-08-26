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
    gROOT->SetBatch(kTRUE);

    std::string input_laser_adc_scan_file = "dump/403_Laser_Scan/LaserScan0.root";
    std::string input_laser_tot_scan_file = "dump/407_ToT_Laser_Scan/ToTLaserScan3.root";

    std::string output_file = "dump/408_Slope_Compare_ADC_ToT/Compare0.root";

    TFile *laser_adc_scan_root = TFile::Open(input_laser_adc_scan_file.c_str(), "READ");
    if (!laser_adc_scan_root || laser_adc_scan_root->IsZombie()) {
        LOG(ERROR) << "Failed to open laser ADC scan file " << input_laser_adc_scan_file;
        return 1;
    }

    // Read ADC slope TVectors
    TVectorD *adc_slope_values = (TVectorD*)laser_adc_scan_root->Get("slope_values");
    TVectorD *adc_slope_errors = (TVectorD*)laser_adc_scan_root->Get("slope_errors");
    TVectorD *adc_channel_numbers = (TVectorD*)laser_adc_scan_root->Get("channel_numbers");
    TVectorD *adc_bias_voltages = (TVectorD*)laser_adc_scan_root->Get("bias_voltages");

    if (!adc_slope_values || !adc_slope_errors || !adc_channel_numbers) {
        LOG(ERROR) << "Failed to read ADC slope vectors from " << input_laser_adc_scan_file;
        laser_adc_scan_root->Close();
        return 1;
    }

    const double bias_error = 0.01; // unit: volt

    std::vector<double> adc_slopes, adc_slope_errs, adc_channels, adc_biases;
    for (int i = 0; i < adc_slope_values->GetNrows(); i++) {
        adc_slopes.push_back((*adc_slope_values)[i]);
        adc_slope_errs.push_back((*adc_slope_errors)[i]);
        adc_channels.push_back((*adc_channel_numbers)[i]);
        if (adc_bias_voltages && i < adc_bias_voltages->GetNrows()) {
            adc_biases.push_back((*adc_bias_voltages)[i]);
        }
    }

    LOG(INFO) << "Read " << adc_slopes.size() << " ADC slope values";
    laser_adc_scan_root->Close();

    TFile *laser_tot_scan_root = TFile::Open(input_laser_tot_scan_file.c_str(), "READ");
    if (!laser_tot_scan_root || laser_tot_scan_root->IsZombie()) {
        LOG(ERROR) << "Failed to open laser ToT scan file " << input_laser_tot_scan_file;
        return 1;
    }

    // Read ToT SSA TVectors
    TVectorD *tot_ssa_channel_numbers = (TVectorD*)laser_tot_scan_root->Get("ssa_channel_numbers");
    TVectorD *tot_ssa_bias_voltages = (TVectorD*)laser_tot_scan_root->Get("ssa_bias_voltages");
    TVectorD *tot_ssa_similarity_values = (TVectorD*)laser_tot_scan_root->Get("ssa_similarity_values");

    if (!tot_ssa_channel_numbers || !tot_ssa_bias_voltages || !tot_ssa_similarity_values) {
        LOG(ERROR) << "Failed to read ToT SSA vectors from " << input_laser_tot_scan_file;
        laser_tot_scan_root->Close();
        return 1;
    }

    std::vector<double> tot_ssa_channels, tot_ssa_biases, tot_ssa_similarities;
    for (int i = 0; i < tot_ssa_channel_numbers->GetNrows(); i++) {
        tot_ssa_channels.push_back((*tot_ssa_channel_numbers)[i]);
        tot_ssa_biases.push_back((*tot_ssa_bias_voltages)[i]);
        tot_ssa_similarities.push_back((*tot_ssa_similarity_values)[i]);
    }

    LOG(INFO) << "Read " << tot_ssa_channels.size() << " ToT SSA values";
    laser_tot_scan_root->Close();

    TFile *output_root = new TFile(output_file.c_str(), "RECREATE");
    if (!output_root || output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << output_file;
        return 1;
    }

    // arrage the slopes by each channel
    // set of the channel
    std::vector<TGraphErrors*> adc_ratio_graphs; // indexed by channel
    std::vector<TGraphErrors*> tot_ratio_graphs; // indexed by channel
    std::vector<std::string> adc_legend_entries;
    std::vector<std::string> tot_legend_entries;
    double global_y_max = 0.0;
    auto unique_channels = std::set<double>(adc_channels.begin(), adc_channels.end());
    for (size_t i = 0; i < unique_channels.size(); i++) {
        double channel_number = adc_channels[i];
        std::vector<double> channel_adc_slopes, channel_adc_slope_errs, channel_adc_biases;
        std::vector<double> channel_tot_similarities, channel_tot_biases;
        for (size_t j = 0; j < adc_channels.size(); j++) {
            if (adc_channels[j] == channel_number) {
                channel_adc_slopes.push_back(adc_slopes[j]);
                channel_adc_slope_errs.push_back(adc_slope_errs[j]);
                if (j < adc_biases.size()) {
                    channel_adc_biases.push_back(adc_biases[j]);
                }
            }
        }
        // Get ToT SSA data for this channel
        for (size_t j = 0; j < tot_ssa_channels.size(); j++) {
            if (tot_ssa_channels[j] == channel_number) {
                channel_tot_similarities.push_back(tot_ssa_similarities[j]);
                channel_tot_biases.push_back(tot_ssa_biases[j]);
            }
        }
        
        LOG(INFO) << "Channel " << channel_number << ": " << channel_adc_slopes.size() 
                  << " ADC slopes, " << channel_tot_similarities.size() << " ToT SSA values";

        // calculate the ratio of slope by referencing 54 V as the baseline
        std::vector<double> adc_slope_ratios;
        std::vector<double> adc_slope_ratio_errs;
        std::vector<double> adc_slope_ratio_biases;
        double baseline_slope = 0.0;
        double baseline_slope_err = 0.0;
        for (size_t j = 0; j < channel_adc_slopes.size(); j++) {
            if (channel_adc_biases[j] == 54.0) {
                baseline_slope = channel_adc_slopes[j];
                baseline_slope_err = channel_adc_slope_errs[j];
                break;
            }
        }
        TGraphErrors *adc_slope_graph = new TGraphErrors();
        int point_index = 0;
        for (size_t j = 0; j < channel_adc_slopes.size(); j++) {
            // skip the 54 V and 43 V
            if (channel_adc_biases[j] == 54.0 || channel_adc_biases[j] == 43.0) {
                continue;
            }
            if (baseline_slope != 0.0) {
                double slope_ratio = baseline_slope / channel_adc_slopes[j];
                double relative_err_squared = std::pow(baseline_slope_err / baseline_slope, 2) + std::pow(channel_adc_slope_errs[j] / channel_adc_slopes[j], 2);
                double slope_ratio_err = slope_ratio * std::sqrt(relative_err_squared);
                adc_slope_ratios.push_back(slope_ratio);
                adc_slope_ratio_errs.push_back(slope_ratio_err);
                adc_slope_ratio_biases.push_back(channel_adc_biases[j]);
                adc_slope_graph->SetPoint(point_index, channel_adc_biases[j], slope_ratio);
                adc_slope_graph->SetPointError(point_index, bias_error, slope_ratio_err);
                point_index++;
            } else {
                adc_slope_ratios.push_back(0.0);
                adc_slope_ratio_errs.push_back(0.0);
                adc_slope_ratio_biases.push_back(channel_adc_biases[j]);
            }
        }
        adc_legend_entries.push_back("ADC Ratio, Channel " + std::to_string((int)channel_number));
        adc_ratio_graphs.push_back(adc_slope_graph);
        // go through the graph to find the maximum y value for setting the same y axis range for all graphs
        for (size_t j = 0; j < channel_adc_slopes.size(); j++) {
            if (channel_adc_biases[j] != 54.0) {
                double slope_ratio = baseline_slope / channel_adc_slopes[j];
                if (slope_ratio > global_y_max) {
                    global_y_max = slope_ratio;
                }
            }
        }

        TGraphErrors *tot_similarity_graph = new TGraphErrors();
        int tot_point_index = 0;
        for (size_t j = 0; j < channel_tot_similarities.size(); j++) {
            tot_similarity_graph->SetPoint(tot_point_index, channel_tot_biases[j], channel_tot_similarities[j]);
            tot_similarity_graph->SetPointError(tot_point_index, bias_error, 0.0); // assuming no error for similarity for now
            if (channel_tot_similarities[j] > global_y_max) {
                global_y_max = channel_tot_similarities[j];
            }
            tot_point_index++;
        }
        tot_legend_entries.push_back("ToT Ratio, Channel " + std::to_string((int)channel_number));
        tot_ratio_graphs.push_back(tot_similarity_graph);
    }

    // draw all the graphs in the same canvas
    TCanvas *canvas = new TCanvas("canvas", "Slope Ratio and ToT SSA Comparison", 1000, 600);
    TLegend *legend = new TLegend(0.7, 0.6, 0.89, 0.89);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    for (size_t i = 0; i < adc_ratio_graphs.size(); i++) {
        adc_ratio_graphs[i]->SetMarkerStyle(20+i);
        adc_ratio_graphs[i]->SetMarkerColor(kCyan+2+i); // different color for each channel
        adc_ratio_graphs[i]->SetLineColor(kCyan+2+i);
        adc_ratio_graphs[i]->SetMaximum(global_y_max * 1.0);
        tot_ratio_graphs[i]->SetMarkerStyle(22+i);
        tot_ratio_graphs[i]->SetMarkerColor(kPink+2+i); // different color for each channel
        tot_ratio_graphs[i]->SetLineColor(kPink+2+i);

        adc_ratio_graphs[i]->GetYaxis()->SetTitle("Gain @ 54 V / Gain @ Bias Voltage");
        adc_ratio_graphs[i]->GetXaxis()->SetTitle("Bias Voltage (V)");

        legend->AddEntry(adc_ratio_graphs[i], adc_legend_entries[i].c_str(), "PE");
        legend->AddEntry(tot_ratio_graphs[i], tot_legend_entries[i].c_str(), "PE");

        if (i == 0) {
            adc_ratio_graphs[i]->Draw("AP");
        } else {
            adc_ratio_graphs[i]->Draw("P");
        }
        tot_ratio_graphs[i]->Draw("P SAME");
    }
    legend->Draw();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->SetTextFont(62);
    double text_x = 0.13;
    double text_y = 0.85;
    double text_y_step = 0.05;
    latex->DrawLatex(text_x, text_y, "Laser Test with H2GCROC");
    text_y -= text_y_step;
    latex->SetTextSize(0.03);
    latex->SetTextFont(42);
    latex->DrawLatex(text_x, text_y, ("ADC and ToT Ratio Relative to 54 V, Channel " + std::to_string((int)adc_channels[0])).c_str());
    text_y -= text_y_step;
    latex->DrawLatex(text_x, text_y, ("Channel count: " + std::to_string(unique_channels.size())).c_str());
    text_y -= text_y_step;
    latex->DrawLatex(text_x, text_y, "February 2026");
    // save as a separate pdf file
    canvas->SaveAs("dump/408_Slope_Compare_ADC_ToT/Compare0.pdf");
    canvas->Write();
    canvas->Close();

    output_root->Close();

    return 0;
}