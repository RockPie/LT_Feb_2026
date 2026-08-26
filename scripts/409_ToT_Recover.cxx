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

// Michaelis-Menten fit function
double michaelis_menten(double *x, double *par) {
    double y_0 = par[0]; // baseline offset
    double A = par[1];   // maximum amplitude
    double x_0 = par[2]; // half-maximum point (Km)
    double x_1_2 = par[3]; // slope factor (Hill coefficient)
    double x_val = x[0];
    double u = x_val - x_0;
    if (u < 0) {
        return y_0; // For x < x_0, return the baseline offset
    } else {
        return y_0 + A * u / (u + x_1_2);
    }
}

// Hill function
double hill_function(double *x, double *par) {
    double y_0 = par[0]; // baseline offset
    double A = par[1];   // maximum amplitude
    double x_0 = par[2]; // half-maximum point (Km)
    double n = par[3];   // Hill coefficient
    double x_1_2 = par[4]; // slope factor for x < x_0
    double x_val = x[0];
    double u = x_val - x_0;
    if (u < 0) {
        return y_0;
    } else {
        return y_0 + A * std::pow(u, n) / (std::pow(u, n) + std::pow(x_1_2, n));
    }
}

// Scaled Hill function for fitting ToT data
double scaled_hill_function(double *x, double *par) {
    double y_0 = par[0]; // baseline offset
    double A = par[1];   // maximum amplitude
    double x_0 = par[2]; // half-maximum point (Km)
    double n = par[3];   // Hill coefficient
    double x_1_2 = par[4]; // slope factor for x < x_0
    double x_scale = par[5]; // scale factor for x
    double x_val = x[0]; // scale the input x value
    double u = (x_val - x_0) * x_scale; // adjust x_0 accordingly
    if (u < 0) {
        return y_0;
    } else {
        return y_0 + A * std::pow(u, n) / (std::pow(u, n) + std::pow(x_1_2 * x_scale, n));
    }
}

double saturating_linear_function(double *x, double *par) {
    double slope = par[0]; // slope of the linear part
    double intercept = par[1]; // y-intercept
    double x_value = x[0];
    double linear_part = slope * x_value + intercept;
    return std::min(linear_part, 1023.0); // Cap the output at 1023
}

int main(int argc, char **argv) {
    gROOT->SetBatch(kTRUE);

    if (argc != 4) {
        LOG(ERROR) << "Usage: " << argv[0] << " <ToT_scan_file> <ADC_scan_file> <output_file>";
        LOG(ERROR) << "Example: " << argv[0] << " dump/405_ToT_Scan/ToTScan5.root dump/402_ADC_Scan/Scan5.root dump/409_ToT_Recover/RecoveredToTScan5.root";
        return 1;
    }

    std::string input_tot_scan_file = argv[1];
    std::string input_adc_scan_file = argv[2];
    std::string output_file = argv[3];

    // Extract scan numbers from input files to verify they match
    auto extract_scan_number = [](const std::string& filename) -> int {
        std::regex scan_regex(R"(Scan(\d+)|ToTScan(\d+))");
        std::smatch match;
        if (std::regex_search(filename, match, scan_regex)) {
            // match[1] for Scan, match[2] for ToTScan
            std::string num_str = match[1].matched ? match[1].str() : match[2].str();
            return std::stoi(num_str);
        }
        return -1;
    };

    int tot_scan_num = extract_scan_number(input_tot_scan_file);
    int adc_scan_num = extract_scan_number(input_adc_scan_file);

    if (tot_scan_num == -1 || adc_scan_num == -1) {
        LOG(ERROR) << "Failed to extract scan numbers from input files";
        LOG(ERROR) << "ToT file: " << input_tot_scan_file << " (scan number: " << tot_scan_num << ")";
        LOG(ERROR) << "ADC file: " << input_adc_scan_file << " (scan number: " << adc_scan_num << ")";
        return 1;
    }

    if (tot_scan_num != adc_scan_num) {
        LOG(ERROR) << "Scan number mismatch!";
        LOG(ERROR) << "ToT file: " << input_tot_scan_file << " (scan " << tot_scan_num << ")";
        LOG(ERROR) << "ADC file: " << input_adc_scan_file << " (scan " << adc_scan_num << ")";
        return 1;
    }

    LOG(INFO) << "Processing scan " << tot_scan_num;
    LOG(INFO) << "ToT input: " << input_tot_scan_file;
    LOG(INFO) << "ADC input: " << input_adc_scan_file;
    LOG(INFO) << "Output: " << output_file;

    // Extract output folder from output file path
    std::string output_folder = output_file.substr(0, output_file.find_last_of("/\\") + 1);

    auto color_adc = TColor::GetColor(255, 62, 155);
    auto color_tot = TColor::GetColor(58, 139, 149);
    auto color_tot_variation = TColor::GetColor(102,208,188);
    auto color_adc_variation = TColor::GetColor(255, 136, 186);

    TFile *adc_scan_root = TFile::Open(input_adc_scan_file.c_str(), "READ");
    if (!adc_scan_root || adc_scan_root->IsZombie()) {
        LOG(ERROR) << "Failed to open ADC scan file " << input_adc_scan_file;
        return 1;
    }

    auto adc_laser_scan_ch50_canvas = (TCanvas*)adc_scan_root->Get("canvas_mean_peak_vs_laser_channel_50");
    TGraphErrors *adc_laser_scan_ch50_graph = nullptr;
    {
        TIter next_primitive(adc_laser_scan_ch50_canvas->GetListOfPrimitives());
        while (auto primitive = next_primitive()) {
            if (std::string(primitive->ClassName()) == "TGraphErrors") {
                adc_laser_scan_ch50_graph = (TGraphErrors*)primitive;
                break;
            }
        }
    }
    auto adc_laser_scan_ch50_graph_clone = (TGraphErrors*)adc_laser_scan_ch50_graph->Clone("canvas_mean_peak_vs_laser_channel_50");

    auto adc_laser_scan_ch52_canvas = (TCanvas*)adc_scan_root->Get("canvas_mean_peak_vs_laser_channel_52");
    TGraphErrors *adc_laser_scan_ch52_graph = nullptr;
    {
        TIter next_primitive(adc_laser_scan_ch52_canvas->GetListOfPrimitives());
        while (auto primitive = next_primitive()) {
            if (std::string(primitive->ClassName()) == "TGraphErrors") {
                adc_laser_scan_ch52_graph = (TGraphErrors*)primitive;
                break;
            }
        }
    }
    auto adc_laser_scan_ch52_graph_clone = (TGraphErrors*)adc_laser_scan_ch52_graph->Clone("canvas_mean_peak_vs_laser_channel_52");

    adc_scan_root->Close();

    TFile *tot_scan_root = TFile::Open(input_tot_scan_file.c_str(), "READ");
    if (!tot_scan_root || tot_scan_root->IsZombie()) {
        LOG(ERROR) << "Failed to open ToT scan file " << input_tot_scan_file;
        return 1;
    }

    auto tot_laser_scan_ch50_canvas = (TCanvas*)tot_scan_root->Get("canvas_mean_tot_vs_laser_channel_50");
    TGraphErrors *tot_laser_scan_ch50_graph = nullptr;
    {
        TIter next_primitive(tot_laser_scan_ch50_canvas->GetListOfPrimitives());
        while (auto primitive = next_primitive()) {
            if (std::string(primitive->ClassName()) == "TGraphErrors") {
                tot_laser_scan_ch50_graph = (TGraphErrors*)primitive;
                break;
            }
        }
    }
    auto tot_laser_scan_ch50_graph_clone = (TGraphErrors*)tot_laser_scan_ch50_graph->Clone("canvas_mean_tot_vs_laser_channel_50");


    auto tot_laser_scan_ch52_canvas = (TCanvas*)tot_scan_root->Get("canvas_mean_tot_vs_laser_channel_52");
    TGraphErrors *tot_laser_scan_ch52_graph = nullptr;
    {
        TIter next_primitive(tot_laser_scan_ch52_canvas->GetListOfPrimitives());
        while (auto primitive = next_primitive()) {
            if (std::string(primitive->ClassName()) == "TGraphErrors") {
                tot_laser_scan_ch52_graph = (TGraphErrors*)primitive;
                break;
            }
        }
    }
    auto tot_laser_scan_ch52_graph_clone = (TGraphErrors*)tot_laser_scan_ch52_graph->Clone("canvas_mean_tot_vs_laser_channel_52");

    tot_scan_root->Close();

    TFile *output_root = new TFile(output_file.c_str(), "RECREATE");
    if (!output_root || output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << output_file;
        return 1;
    }

    output_root->cd();

    auto canvas_adc_tot_ch50 = new TCanvas("canvas_adc_tot_ch50", "ADC vs ToT for Channel 50", 1000, 600);
    double y_scale_factor = 1.4;
    const double tot_to_adc_scale = 1024.0 / 4096.0;
    const double axis_text_size = 0.04;
    const double x_axis_min = 5.4;
    const double x_axis_max = 10.2;

    auto legend_adc_tot_ch50 = new TLegend(0.5, 0.7, 0.89, 0.89);
    legend_adc_tot_ch50->SetFillStyle(0);
    legend_adc_tot_ch50->SetBorderSize(0);
    legend_adc_tot_ch50->SetTextFont(42);
    legend_adc_tot_ch50->SetTextSize(0.025);

    adc_laser_scan_ch50_graph_clone->SetMarkerColor(color_adc);
    adc_laser_scan_ch50_graph_clone->SetLineColor(color_adc);
    adc_laser_scan_ch50_graph_clone->SetTitle(";Laser Intensity [a.u.];Mean ADC Peak");
    adc_laser_scan_ch50_graph_clone->GetYaxis()->SetRangeUser(0, 1024*y_scale_factor);
    adc_laser_scan_ch50_graph_clone->GetXaxis()->SetLimits(x_axis_min, x_axis_max);
    adc_laser_scan_ch50_graph_clone->GetXaxis()->SetRangeUser(x_axis_min, x_axis_max);
    adc_laser_scan_ch50_graph_clone->Draw("AP");

    double fit_x_min = 5.73;
    double fit_x_max = x_axis_max;
    // TF1 *adc_fit = new TF1("adc_fit", "pol1", fit_x_min, fit_x_max);
    // adc_laser_scan_ch50_graph_clone->Fit(adc_fit, "QR", "", fit_x_min, fit_x_max);
    // double adc_fit_slope = adc_fit->GetParameter(1);
    // double adc_fit_intercept = adc_fit->GetParameter(0);
    TF1 *adc_fit = new TF1("adc_fit", saturating_linear_function, fit_x_min, fit_x_max, 2);
    adc_fit->SetParameters(1000.0, -5000.0); // Initial guess
    adc_fit->SetParLimits(0, 0, 1e5);
    adc_fit->SetParLimits(1, -1e6, -1000.0);
    adc_laser_scan_ch50_graph_clone->Fit(adc_fit, "QR", "", fit_x_min, fit_x_max);
    double adc_fit_slope = adc_fit->GetParameter(0);
    double adc_fit_intercept = adc_fit->GetParameter(1);
    adc_fit->SetLineColorAlpha(color_adc_variation, 0.7);
    adc_fit->Draw("SAME");
    LOG(INFO) << "ADC fit parameters: slope = " << adc_fit_slope << ", intercept = " << adc_fit_intercept;

    legend_adc_tot_ch50->AddEntry(adc_laser_scan_ch50_graph_clone, "Mean ADC Peak", "PE");
    if (adc_fit_intercept > 0)
        legend_adc_tot_ch50->AddEntry(adc_fit, Form("y = %.2f x + %.2f", adc_fit_slope, adc_fit_intercept), "L");
    else
        legend_adc_tot_ch50->AddEntry(adc_fit, Form("y = %.2f x - %.2f", adc_fit_slope, -adc_fit_intercept), "L");
    canvas_adc_tot_ch50->Update();


    const int axis_label_font = 42; // Helvetica
    auto tot_laser_scan_ch50_graph_scaled = (TGraphErrors*)tot_laser_scan_ch50_graph_clone->Clone("canvas_mean_tot_vs_laser_channel_50_scaled");
    auto tot_interpolated_ch50 = InterpolateWithUncertainty_AkimaMC(tot_laser_scan_ch50_graph, 100, 2000, false);
    auto tot_interpolated_ch50_scaled = (TGraphErrors*)tot_interpolated_ch50->Clone("canvas_mean_tot_vs_laser_channel_50_interpolated_scaled");
    for (int i = 0; i < tot_laser_scan_ch50_graph_scaled->GetN(); i++) {
        double x = 0.0;
        double y = 0.0;
        tot_laser_scan_ch50_graph_scaled->GetPoint(i, x, y);
        tot_laser_scan_ch50_graph_scaled->SetPoint(i, x, y * tot_to_adc_scale);
        tot_laser_scan_ch50_graph_scaled->SetPointError(i,
            tot_laser_scan_ch50_graph_scaled->GetErrorX(i),
            tot_laser_scan_ch50_graph_scaled->GetErrorY(i) * tot_to_adc_scale);
    }
    for (int i = 0; i < tot_interpolated_ch50_scaled->GetN(); i++) {
        double x = 0.0;
        double y = 0.0;
        tot_interpolated_ch50_scaled->GetPoint(i, x, y);
        tot_interpolated_ch50_scaled->SetPoint(i, x, y * tot_to_adc_scale);
        tot_interpolated_ch50_scaled->SetPointError(i,
            tot_interpolated_ch50_scaled->GetErrorX(i),
            tot_interpolated_ch50_scaled->GetErrorY(i) * tot_to_adc_scale);
    }
    
    tot_laser_scan_ch50_graph_scaled->SetMarkerColor(color_tot);
    tot_laser_scan_ch50_graph_scaled->SetLineColor(color_tot);
    tot_laser_scan_ch50_graph_scaled->SetTitle(";Laser Intensity [a.u.];Mean ToT");
    tot_laser_scan_ch50_graph_scaled->GetYaxis()->SetRangeUser(0, 4096*y_scale_factor*tot_to_adc_scale);
    tot_laser_scan_ch50_graph_scaled->GetXaxis()->SetLimits(x_axis_min, x_axis_max);
    tot_laser_scan_ch50_graph_scaled->GetXaxis()->SetRangeUser(x_axis_min, x_axis_max);
    tot_laser_scan_ch50_graph_scaled->Draw("P SAME");
    legend_adc_tot_ch50->AddEntry(tot_laser_scan_ch50_graph_scaled, "Mean ToT", "PE");

    // do the Michaelis-Menten fit
    // TF1 *tot_fit = new TF1("tot_fit", michaelis_menten, 5.2, 10.5, 4);
    // tot_fit->SetParameters(0, 4096*tot_to_adc_scale, 6.0, 0.5);
    // tot_fit->SetParLimits(0, 0, 4096*tot_to_adc_scale);
    // tot_fit->SetParLimits(1, 0, 4096*tot_to_adc_scale*2);
    // tot_fit->SetParLimits(2, 5.0, 10.0);
    // tot_fit->SetParLimits(3, 0.1, 5.0);
    // tot_laser_scan_ch50_graph_scaled->Fit(tot_fit, "QR", "", 5.2, 10.5);
    // tot_fit->SetLineColorAlpha(color_tot_variation, 0.7);
    // tot_fit->Draw("SAME");
    // double tot_fit_y_max = tot_fit->Eval(x_axis_max);
    // LOG(INFO) << "ToT fit parameters: y_0 = " << tot_fit->GetParameter(0) << ", A = " << tot_fit->GetParameter(1)
    //           << ", x_0 = " << tot_fit->GetParameter(2) << ", x_1_2 = " << tot_fit->GetParameter(3);
    // legend_adc_tot_ch50->AddEntry(tot_fit, Form("ToT Fit: y = %.1f + %.1f * (x - %.1f) / ((x - %.1f) + %.1f)", tot_fit->GetParameter(0), tot_fit->GetParameter(1), tot_fit->GetParameter(2), tot_fit->GetParameter(2), tot_fit->GetParameter(3)), "L");

    // do the Hill fit
    // TF1 *tot_fit = new TF1("tot_fit", hill_function, 5.2, 10.5, 5);
    // tot_fit->SetParameters(0, 4096*tot_to_adc_scale, 6.0, 2.0, 0.5);
    // tot_fit->FixParameter(0, 0); // fix baseline offset to 0
    // tot_fit->SetParLimits(1, 0, 4096*tot_to_adc_scale*2);
    // tot_fit->SetParLimits(2, 5.9, 6.0);
    // tot_fit->SetParLimits(3, 0.01, 5.0);
    // tot_fit->SetParLimits(4, 0.01, 5.0);

    // // set initial parameter values
    // tot_fit->SetParameter(0, 0); // baseline offset
    // tot_fit->SetParameter(1, 2500*tot_to_adc_scale); // maximum amplitude
    // tot_fit->SetParameter(2, 5.70);
    // tot_fit->SetParameter(3, 0.7); // Hill coefficient
    // tot_fit->SetParameter(4, 0.5); // slope factor for x < x_0

    
    // tot_fit->SetLineColorAlpha(color_tot_variation, 0.7);

    // tot_laser_scan_ch50_graph_scaled->Fit(tot_fit, "QR", "", 5.2, 10.5);
    // double fitted_y_0 = tot_fit->GetParameter(0) / tot_to_adc_scale;
    // double fitte_d_A = tot_fit->GetParameter(1) / tot_to_adc_scale;
    // double fitted_x_0 = tot_fit->GetParameter(2);
    // double fitted_n = tot_fit->GetParameter(3);
    // double fitted_x_1_2 = tot_fit->GetParameter(4);
    // double fitted_chi2 = tot_fit->GetChisquare();
    // double fitted_ndf = tot_fit->GetNDF();
    // tot_fit->SetLineColorAlpha(color_tot_variation, 0.7);
    // tot_fit->Draw("SAME");
    // legend_adc_tot_ch50->AddEntry(tot_fit, Form("y = %.0f * (x - %.1f)^{%.2f} / ((x - %.1f)^{%.2f} + %.1f^{%.2f})", fitte_d_A, fitted_x_0, fitted_n, fitted_x_0, fitted_n, fitted_x_1_2, fitted_n), "L");
    // legend_adc_tot_ch50->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF = %.1f / %d", fitted_chi2, (int)fitted_ndf), "");

    // do the scaled Hill fit
    TF1 *tot_fit = new TF1("tot_fit", scaled_hill_function, 5.2, 10.5, 6);
    tot_fit->SetParameters(0, 4096*tot_to_adc_scale, 6.0, 2.0, 0.5, 1.0);
    tot_fit->FixParameter(0, 0); // fix baseline offset to 0
    tot_fit->SetParLimits(1, 0, 4096*tot_to_adc_scale*2);
    tot_fit->SetParLimits(2, 5.0, 7.7);
    tot_fit->SetParLimits(3, 0.01, 5.0);
    tot_fit->SetParLimits(4, 0.01, 5.0);
    tot_fit->SetParLimits(5, 0.01, 10.0);

    tot_fit->SetParameter(0, 0); // baseline offset
    tot_fit->SetParameter(1, 2500*tot_to_adc_scale); // maximum amplitude
    tot_fit->SetParameter(2, 5.70);
    tot_fit->SetParameter(3, 0.7); // Hill coefficient
    tot_fit->SetParameter(4, 0.5); // slope factor for x < x_0
    tot_fit->SetParameter(5, 1.0); // x scale factor

    // tot_fit->FixParameter(4, 0.2); // fix baseline offset to 0
    // tot_fit->FixParameter(3, 0.76); // fix baseline offset to 0

    tot_laser_scan_ch50_graph_scaled->Fit(tot_fit, "QR", "", 5.2, 10.5);
    double fitted_y_0 = tot_fit->GetParameter(0) / tot_to_adc_scale;
    double fitted_A = tot_fit->GetParameter(1) / tot_to_adc_scale;
    double fitted_x_0 = tot_fit->GetParameter(2);
    double fitted_n = tot_fit->GetParameter(3);
    double fitted_x_1_2 = tot_fit->GetParameter(4);
    double fitted_x_scale = tot_fit->GetParameter(5);
    double fitted_chi2 = tot_fit->GetChisquare();
    double fitted_ndf = tot_fit->GetNDF();

    tot_fit->SetLineColor(color_tot_variation);
    tot_fit->SetLineWidth(2);
    tot_fit->Draw("SAME");

    legend_adc_tot_ch50->AddEntry(tot_fit, Form("y = %.0f #times u^{%.2f} / (u^{%.2f} + %.2f^{%.2f})", fitted_A, fitted_n, fitted_n, fitted_x_1_2, fitted_n), "L");
    legend_adc_tot_ch50->AddEntry((TObject*)nullptr, Form("u = (x - %.1f) #times %.2f,   #chi^{2}/NDF = %.1f / %d", fitted_x_0, fitted_x_scale, fitted_chi2, (int)fitted_ndf), "");

    // tot_interpolated_ch50_scaled->SetMarkerColor(color_tot_variation);
    // tot_interpolated_ch50_scaled->SetMarkerStyle(0);
    // tot_interpolated_ch50_scaled->SetLineColor(color_tot_variation);
    // // tot_interpolated_ch50_scaled->SetLineWidth(0.2);
    // tot_interpolated_ch50_scaled->Draw("LE SAME");
    // legend_adc_tot_ch50->AddEntry(tot_interpolated_ch50_scaled, "Interpolated ToT", "LE");

    auto right_axis = new TGaxis(x_axis_max, 0, x_axis_max, 1024*y_scale_factor, 0, 4096, 510, "+L");
    right_axis->SetTitle("Mean ToT");
    right_axis->SetTitleSize(axis_text_size);
    right_axis->SetTitleFont(axis_label_font);
    right_axis->SetLabelSize(axis_text_size);
    right_axis->SetLabelFont(axis_label_font);
    right_axis->Draw();
    
    legend_adc_tot_ch50->Draw();
    TLatex latex_adc_tot_ch50;
    latex_adc_tot_ch50.SetNDC();
    latex_adc_tot_ch50.SetTextSize(0.04);
    latex_adc_tot_ch50.SetTextFont(62);
    double text_x = 0.13;
    double text_y = 0.85;
    double text_y_step = 0.045;
    latex_adc_tot_ch50.DrawLatex(text_x, text_y, "Laser Test with H2GCROC");
    latex_adc_tot_ch50.SetTextSize(0.03);
    latex_adc_tot_ch50.SetTextFont(42);
    text_y -= text_y_step;
    latex_adc_tot_ch50.DrawLatex(text_x, text_y, "Channel 50");
    text_y -= text_y_step;
    latex_adc_tot_ch50.DrawLatex(text_x, text_y, "ADC and ToT response comparison");
    text_y -= text_y_step;
    latex_adc_tot_ch50.DrawLatex(text_x, text_y, "February 2026, CERN");
    canvas_adc_tot_ch50->Write();
    canvas_adc_tot_ch50->SaveAs((output_folder + "adc_vs_tot_channel_50_Scan" + std::to_string(tot_scan_num) + ".pdf").c_str());
    canvas_adc_tot_ch50->Close();

    output_root->Close();
    return 0;
}