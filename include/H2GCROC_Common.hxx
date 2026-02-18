#ifndef COMMON_HPP
#define COMMON_HPP

// ============================================================================
// INCLUDES
// ============================================================================
#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
#include <string>
#include <regex>

#include "TCanvas.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TF1.h"
#include "TF1Convolution.h"
#include "TH2.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TVirtualFitter.h"
#include "TPad.h"
#include "TStyle.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadExecutor.hxx>

#include "easylogging++.h"
#include "nlohmann/json.hpp"
#include "argparse/argparse.hpp"

// ============================================================================
// DEFINES
// ============================================================================
#ifndef FPGA_CHANNEL_NUMBER
#define FPGA_CHANNEL_NUMBER 152
#endif

#ifndef FPGA_CHANNEL_NUMBER_VALID
#define FPGA_CHANNEL_NUMBER_VALID 144
#endif

using json = nlohmann::json;

// ============================================================================
// LOGGING
// ============================================================================
void set_easylogger();

// ============================================================================
// DATA STRUCTURES
// ============================================================================

struct ScriptOptions {
    std::string input_file;
    std::string output_file;
    std::string output_folder;
    int n_events;
    bool verbose;
    bool focal;
    bool timewalk;
    std::string script_name;
    std::string script_version;
    std::string pedestal_file;
    std::string csv_file;
    std::string timewalk_file;
    std::string fitting_file;
};

// ============================================================================
// ARGUMENT PARSING
// ============================================================================
ScriptOptions parse_arguments_single_root(int argc, char **argv, const std::string& version = "0.1");
ScriptOptions parse_arguments_single_json(int argc, char **argv, const std::string& version = "0.1");
ScriptOptions parse_arguments_single_root_single_csv(int argc, char **argv, const std::string& version = "0.1");

// ============================================================================
// CHANNEL UTILITIES
// ============================================================================
int get_valid_fpga_channel(int fpga_channel);
int get_total_fpga_channel(int fpga_channel);

inline int get_unified_fpga_channel(int fpga_id, int fpga_channel){
    return fpga_id * FPGA_CHANNEL_NUMBER + fpga_channel;
}

inline int get_unified_valid_fpga_channel(int fpga_id, int fpga_channel){
    return fpga_id * FPGA_CHANNEL_NUMBER_VALID + fpga_channel;
}

inline UInt_t decode_tot_value(UInt_t val1) {
    UInt_t mask = -(val1 >= 512);  // mask is 0xFFFFFFFF if val1 >= 512, else 0
    return (val1 & ~mask) | ((val1 - 512) * 8 & mask);
}

inline double decode_toa_value_ns(UInt_t val2) {
    constexpr double scale0 = 0.025;
    constexpr double scale1 = 0.2;
    constexpr double scale2 = 6.25;

    UInt_t part0 = val2 & 0x07;
    UInt_t part1 = (val2 >> 3) & 0x1F;
    UInt_t part2 = (val2 >> 8) & 0x1F;

    return part0 * scale0 + part1 * scale1 + part2 * scale2;
}

// ============================================================================
// GLOBAL CHANNEL PAINTER
// ============================================================================
// ============================================================================
// GLOBAL CHANNEL PAINTER
// ============================================================================
class GlobalChannelPainter {
public:
    GlobalChannelPainter(const std::string& mapping_file);
    GlobalChannelPainter(const std::string& mapping_file, const std::string& channel_mapping_file);
    ~GlobalChannelPainter();

    TCanvas* get_canvas() { return painter_canvas; }
    
    void draw_global_channel_hists2D(std::vector <TH2D*> hists, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title);
    void draw_global_channel_hists1D(std::vector <TH1D*> hists, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title);
    void draw_global_channel_hists1D_group(std::vector <std::vector <TH1D*>> hists_list, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, std::vector <EColor> colors, std::vector <std::string> legend_labels);
    void draw_global_channel_hists1D_run_group(std::vector <std::vector <TH1D*>> hists_list, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, std::vector <EColor> colors, std::vector <std::string> legend_labels);
    void draw_global_channel_canvas(std::vector <TCanvas*> canvas_list,  std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title);
    void draw_canvas_components(TCanvas* source_canvas, TPad* target_pad);

    TCanvas* draw_module_channel_canvas(std::vector <TCanvas*> canvas_list, std::unordered_map <int, int> hists_channel_map, const std::string& hist_name, const std::string& hist_title, int module_index);

private:
    void clear_canvas();

private:
    TCanvas *painter_canvas;
    json mapping_json;

    json focal_mapping_json;
    json focal_channel_mapping_json;

    bool is_EEEMCal_mapping;
    bool is_FoCal_mapping;

    std::vector <int> module_fpga_list;
    std::vector <int> module_asic_list;
    std::vector <int> module_connector_list;
    std::vector <std::vector <int>> connector_list_list;

    std::vector <std::vector <int>> focal_module_board_list;
    std::vector <std::vector <int>> focal_module_channel_list;
    std::unordered_map <int, int> focal_channel_map;
    std::unordered_map <int, int> focal_fpga_map;

    std::vector <TCanvas*> sub_canvas_list;
};

// ============================================================================
// DRAWING UTILITIES
// ============================================================================

inline void draw_on_pad(TPad* pad, TObject* obj, bool minimalist_axis, bool th2_logz, TF1* fit_func = nullptr)
{
    if (!pad || !obj) return;
    pad->cd();

    pad->SetMargin(0, 0, 0, 0);
    pad->SetBorderMode(0);
    pad->SetFrameBorderMode(0);
    pad->SetFillStyle(0);

    // TH2* histograms
    if (obj->InheritsFrom(TH2::Class())) {
        auto* h2 = static_cast<TH2*>(obj);
        h2->SetStats(0);
        pad->SetLogz(th2_logz ? 1 : 0);
        h2->Draw("COLZ");

        if (minimalist_axis) {
            auto *xa = h2->GetXaxis(), *ya = h2->GetYaxis(), *za = h2->GetZaxis();
            if (xa) { xa->SetLabelSize(0); xa->SetTitleSize(0); xa->SetTickLength(0); }
            if (ya) { ya->SetLabelSize(0); ya->SetTitleSize(0); ya->SetTickLength(0); }
            if (za) { za->SetLabelSize(0); za->SetTitleSize(0); }
        }
        pad->Modified();
        return;
    }

    // TH1* histograms
    if (obj->InheritsFrom(TH1::Class())) {
        auto* h1 = static_cast<TH1*>(obj);
        h1->SetStats(0);
        pad->SetLogz(0);
        h1->Draw("HIST");
        if (minimalist_axis) {
            auto *xa = h1->GetXaxis(), *ya = h1->GetYaxis();
            if (xa) { xa->SetLabelSize(0); xa->SetTitleSize(0); xa->SetTickLength(0); }
            if (ya) { ya->SetLabelSize(0); ya->SetTitleSize(0); ya->SetTickLength(0); }
        }
        if (fit_func) {
            fit_func->SetLineColor(kRed);
            fit_func->Draw("SAME");
        }
        pad->Modified();
        return;
    }

    // TGraphErrors* and TGraph*
    if (obj->InheritsFrom(TGraphErrors::Class()) || obj->InheritsFrom(TGraph::Class())) {
        auto* gr = static_cast<TGraph*>(obj);
        pad->SetLogz(0);
        gr->Draw("AP");
        if (minimalist_axis) {
            if (auto* hframe = gr->GetHistogram()) {
                auto *xa = hframe->GetXaxis(), *ya = hframe->GetYaxis();
                if (xa) { xa->SetLabelSize(0); xa->SetTitleSize(0); xa->SetTickLength(0); }
                if (ya) { ya->SetLabelSize(0); ya->SetTitleSize(0); ya->SetTickLength(0); }
            }
        }
        pad->Modified();
        return;
    }

    // Default fallback
    pad->SetLogz(0);
    obj->Draw();
    pad->Modified();
}

inline void draw_on_pad(TPad* pad, TObject* obj, TObject* fit_obj, bool minimalist_axis, bool th2_logz, TF1* fit_func = nullptr){
    draw_on_pad(pad, obj, minimalist_axis, th2_logz, fit_func);
    if (fit_obj) {
        pad->cd();
        fit_obj->Draw("L SAME");
        pad->Modified();
    }
}

// ============================================================================
// GRID BUILDER
// ============================================================================

inline void build_gapless_grid(TCanvas& canvas, int NX, int NY,
                               std::vector<std::vector<TPad*>>& pads)
{
    assert(NX > 0 && NY > 0);
    pads.assign(NY, std::vector<TPad*>(NX, nullptr));

    // Setup global style
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);

    const double dx  = 1.0 / NX;
    const double dy  = 1.0 / NY;
    const double eps = 1e-6;

    for (int col = 0; col < NX; ++col) {
        for (int row = 0; row < NY; ++row) {
            const double x1 = col * dx;
            const double x2 = (col + 1) * dx;
            const double y1 = 1.0 - (row + 1) * dy;
            const double y2 = 1.0 - row * dy;

            const double xl = (col == 0)    ? 0.0 : x1 - eps;
            const double xr = (col == NX-1) ? 1.0 : x2 + eps;
            const double yb = (row == NY-1) ? 0.0 : y1 - eps;
            const double yt = (row == 0)    ? 1.0 : y2 + eps;

            TString nm = Form("pad_r%d_c%d", row, col);
            auto* p = new TPad(nm, nm, xl, yb, xr, yt);
            p->SetMargin(0, 0, 0, 0);
            p->SetBorderMode(0);
            p->SetFrameBorderMode(0);
            p->SetFillStyle(0);
            p->Draw();
            pads[row][col] = p;
        }
    }
}

    // Mapping structure for placing histograms on canvas pads.
    // Uses chan2pad LUT to map channel indices to pad positions.
    // Linear index definitions:
    //   chan_linear = vldb * channels_per_vldb + channel
    //   pad_linear  = row * NX + col (row-major)
struct MosaicTopology {
    int NX{16};
    int NY{12};
    int vldb_number{0};
    int channels_per_vldb{0};
    bool reverse_row{true};     // 是否最后做行反转（与你当前图片一致）
    bool minimalist_axis{true};
    bool th2_logz{true};

    // 固定顺序 LUT：长度应为 vldb_number * channels_per_vldb。
    // 若某通道不想画，置为 -1。
    std::vector<int> chan2pad;  // pad 的线性索引（0..NX*NY-1）或 -1

    inline int pad_linear_of(int col, int row) const {
        if (reverse_row) row = (NY - 1) - row;
        return row * NX + col;
    }

    // 安全性检查
    bool valid() const {
        const int need = vldb_number * channels_per_vldb;
        if (NX <= 0 || NY <= 0 || vldb_number <= 0 || channels_per_vldb <= 0) return false;
        if ((int)chan2pad.size() != need) return false;
        return true;
    }
};

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TObject*>& items,
                              const std::vector<TF1*>& fits,
                              const MosaicTopology& topo)
{
    if (!topo.valid()) return;
    if (items.size() != fits.size()) return;

    std::vector<std::vector<TPad*>> pads;
    build_gapless_grid(canvas, topo.NX, topo.NY, pads);

    const int pad_count = topo.NX * topo.NY;

    // 预构建 pad* 的线性数组，O(1) 访问
    std::vector<TPad*> pad_linear(pad_count, nullptr);
    for (int r = 0; r < topo.NY; ++r) {
        for (int c = 0; c < topo.NX; ++c) {
            int pr = r;
            if (topo.reverse_row) pr = (topo.NY - 1) - r; // 因为 pad_linear 用的就是最终的行序
            pad_linear[pr * topo.NX + c] = pads[r][c];
        }
    }

    const int n_items_expected = topo.vldb_number * topo.channels_per_vldb;
    const int n_items = std::min<int>( (int)items.size(), n_items_expected );

    for (int i = 0; i < n_items; ++i) {
        TObject* obj = items[i];
        TF1* fit = fits[i];
        if (!obj || !fit) continue;

        // Skip empty histograms
        if (auto* h1 = dynamic_cast<TH1*>(obj)) {
            if (h1->GetEntries() == 0) continue;
        } else if (auto* gr = dynamic_cast<TGraph*>(obj)) {
            if (gr->GetN() == 0) continue;
        }

        int pad_lnr = topo.chan2pad[i];
        if (pad_lnr < 0 || pad_lnr >= pad_count) continue;

        TPad* pad = pad_linear[pad_lnr];
        if (!pad) continue;

        draw_on_pad(pad, obj, topo.minimalist_axis, topo.th2_logz, fit);
    }

    canvas.Update();
}

// ============================================================================
// MOSAIC DRAWING FUNCTIONS (OVERLOADS)
// ============================================================================

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TObject*>& items,
                              const MosaicTopology& topo)
{
    if (!topo.valid()) return;

    std::vector<std::vector<TPad*>> pads;
    build_gapless_grid(canvas, topo.NX, topo.NY, pads);

    const int pad_count = topo.NX * topo.NY;

    // 预构建 pad* 的线性数组，O(1) 访问
    std::vector<TPad*> pad_linear(pad_count, nullptr);
    for (int r = 0; r < topo.NY; ++r) {
        for (int c = 0; c < topo.NX; ++c) {
            int pr = r;
            if (topo.reverse_row) pr = (topo.NY - 1) - r; // 因为 pad_linear 用的就是最终的行序
            pad_linear[pr * topo.NX + c] = pads[r][c];
        }
    }

    const int n_items_expected = topo.vldb_number * topo.channels_per_vldb;
    const int n_items = std::min<int>( (int)items.size(), n_items_expected );

    for (int i = 0; i < n_items; ++i) {
        TObject* obj = items[i];
        if (!obj) continue;

        // Skip empty histograms
        if (auto* h1 = dynamic_cast<TH1*>(obj)) {
            if (h1->GetEntries() == 0) continue;
        } else if (auto* gr = dynamic_cast<TGraph*>(obj)) {
            if (gr->GetN() == 0) continue;
        }

        int pad_lnr = topo.chan2pad[i];
        if (pad_lnr < 0 || pad_lnr >= pad_count) continue;

        TPad* pad = pad_linear[pad_lnr];
        if (!pad) continue;

        draw_on_pad(pad, obj, topo.minimalist_axis, topo.th2_logz);
    }

    canvas.Update();
}

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TObject*>& items,
                              const std::vector<TObject*>& fits,
                              const MosaicTopology& topo)
{
    if (!topo.valid()) return;

    std::vector<std::vector<TPad*>> pads;
    build_gapless_grid(canvas, topo.NX, topo.NY, pads);

    const int pad_count = topo.NX * topo.NY;

    // 预构建 pad* 的线性数组，O(1) 访问
    std::vector<TPad*> pad_linear(pad_count, nullptr);
    for (int r = 0; r < topo.NY; ++r) {
        for (int c = 0; c < topo.NX; ++c) {
            int pr = r;
            if (topo.reverse_row) pr = (topo.NY - 1) - r; // 因为 pad_linear 用的就是最终的行序
            pad_linear[pr * topo.NX + c] = pads[r][c];
        }
    }

    const int n_items_expected = topo.vldb_number * topo.channels_per_vldb;
    const int n_items = std::min<int>( (int)items.size(), n_items_expected );

    for (int i = 0; i < n_items; ++i) {
        TObject* obj = items[i];
        TObject* fit_obj = fits[i];
        if (!obj) continue;
        if (!fit_obj) continue;

        // Skip empty histograms
        if (auto* h1 = dynamic_cast<TH1*>(obj)) {
            if (h1->GetEntries() == 0) continue;
        } else if (auto* gr = dynamic_cast<TGraph*>(obj)) {
            if (gr->GetN() == 0) continue;
        }

        int pad_lnr = topo.chan2pad[i];
        if (pad_lnr < 0 || pad_lnr >= pad_count) continue;

        TPad* pad = pad_linear[pad_lnr];
        if (!pad) continue;

        draw_on_pad(pad, obj, fit_obj, topo.minimalist_axis, topo.th2_logz);
    }

    canvas.Update();
}

inline void AutoZoomByQuantile(TH1D* h,
                        double qlow = 0.01,
                        double qhigh = 0.99,
                        double expand_frac = 0.10)
{
    if (!h || h->GetEntries() == 0) return;

    double probs[2] = {qlow, qhigh};
    double quantiles[2];

    h->GetQuantiles(2, quantiles, probs);

    double qmin = quantiles[0];
    double qmax = quantiles[1];

    if (qmax <= qmin) return;

    double width = qmax - qmin;
    double xmin = qmin - expand_frac * width;
    double xmax = qmax + expand_frac * width;

    h->GetXaxis()->SetRangeUser(xmin, xmax);
}

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TH2D*>& h2,
                              const std::vector<TGraph*>& fits,
                              const MosaicTopology& topo)
{
    std::vector<TObject*> objs; objs.reserve(h2.size());
    for (auto* p : h2) objs.push_back(static_cast<TObject*>(p));
    std::vector<TObject*> fit_objs; fit_objs.reserve(fits.size());
    for (auto* p : fits) fit_objs.push_back(static_cast<TObject*>(p));
    draw_mosaic_fixed(canvas, objs, fit_objs, topo);
}

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TH2D*>& h2,
                              const MosaicTopology& topo)
{
    std::vector<TObject*> objs; objs.reserve(h2.size());
    for (auto* p : h2) objs.push_back(static_cast<TObject*>(p));
    draw_mosaic_fixed(canvas, objs, topo);
}

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TH1D*>& h1,
                              const MosaicTopology& topo)
{
    std::vector<TObject*> objs; objs.reserve(h1.size());
    for (auto* p : h1) objs.push_back(static_cast<TObject*>(p));
    // 对 TH1*，如果你不想 LogZ，在构建 topo 时把 th2_logz=false 即可（只影响 TH2）
    draw_mosaic_fixed(canvas, objs, topo);
}

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TH1I*>& h1,
                              const MosaicTopology& topo)
{
    std::vector<TObject*> objs; objs.reserve(h1.size());
    for (auto* p : h1) objs.push_back(static_cast<TObject*>(p));
    draw_mosaic_fixed(canvas, objs, topo);
}

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TH1D*>& h1,
                              const std::vector<TF1*>& fits,
                              const MosaicTopology& topo)
{
    std::vector<TObject*> objs; objs.reserve(h1.size());
    for (auto* p : h1) objs.push_back(static_cast<TObject*>(p));
    draw_mosaic_fixed(canvas, objs, fits, topo);
}

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TGraph*>& gr,
                              const MosaicTopology& topo)
{
    std::vector<TObject*> objs; objs.reserve(gr.size());
    for (auto* p : gr) objs.push_back(static_cast<TObject*>(p));
    draw_mosaic_fixed(canvas, objs, topo);
}

inline void draw_mosaic_fixed(TCanvas& canvas,
                              const std::vector<TGraphErrors*>& gr,
                              const MosaicTopology& topo)
{
    std::vector<TObject*> objs; objs.reserve(gr.size());
    for (auto* p : gr) objs.push_back(static_cast<TObject*>(p));
    draw_mosaic_fixed(canvas, objs, topo);
}

// ============================================================================
// LUT AND TOPOLOGY INITIALIZATION
// ============================================================================

inline std::vector<int> build_chan2pad_LUT(
    int vldb_number,
    int channels_per_vldb,
    int NX, int NY,
    int board_cols, int board_rows,
    const json& sipm_board,          // 通道到(row,col)
    const json& board_loc,           // 板到(row,col)
    const json& board_rotation,      // 板是否旋转
    const json& board_flip           // 板是否翻转（可不使用）
) {
    std::vector<int> lut(vldb_number * channels_per_vldb, -1);

    for (int vldb_id = 0; vldb_id < vldb_number; ++vldb_id) {
        for (int ch = 0; ch < channels_per_vldb; ++ch) {

            // ---------------- 通道拓扑逻辑 ----------------
            int channel_index_in_h2g       = ch % 76;
            int channel_index_in_h2g_half  = channel_index_in_h2g % 38;
            int asic_id                    = ch / 76;
            int half_id                    = channel_index_in_h2g / 38;

            if (channel_index_in_h2g_half == 19) continue;
            if (channel_index_in_h2g_half == 0)  continue;

            int channel_index_no_CM_Calib = channel_index_in_h2g_half;
            if (channel_index_in_h2g_half > 19)      channel_index_no_CM_Calib -= 2;
            else if (channel_index_in_h2g_half > 0)  channel_index_no_CM_Calib -= 1;

            std::string chn_key = std::to_string(channel_index_no_CM_Calib);
            int col = -1, row = -1;
            if (sipm_board.contains(chn_key)) {
                col = sipm_board.at(chn_key).at("col").get<int>();
                row = sipm_board.at(chn_key).at("row").get<int>();
            } else {
                continue;
            }

            std::string board_key = std::to_string(half_id + asic_id * 2 + vldb_id * 4);
            int board_col = -1, board_row = -1;
            int board_rotated = 0;
            int board_flipped = 0;

            if (board_loc.contains(board_key)) {
                board_col = board_loc.at(board_key).at("col").get<int>();
                board_row = board_loc.at(board_key).at("row").get<int>();
            }
            if (board_rotation.contains(board_key)) {
                board_rotated = board_rotation.at(board_key).get<int>();
            }
            if (board_flip.contains(board_key)) {
                board_flipped = board_flip.at(board_key).get<int>();
            }

            int uni_col = -1, uni_row = -1;
            if (board_rotated == 0) {
                uni_col = board_col * board_cols + col;
                uni_row = board_row * board_rows + row;
            } else {
                // 旋转180°
                uni_col = board_col * board_cols + (board_cols - 1 - col);
                uni_row = board_row * board_rows + (board_rows - 1 - row);
            }

            if (uni_col < 0 || uni_row < 0) continue;
            if (uni_col >= NX || uni_row >= NY) continue;

            int pad_linear = uni_row * NX + uni_col;
            int chan_linear = vldb_id * channels_per_vldb + ch;

            lut[chan_linear] = pad_linear;
        }
    }
    return lut;
}

// ============================================================================
// MOSAIC TOPOLOGY SETUP
// ============================================================================

struct MosaicTopoSetup {
    MosaicTopology topo_wave;
    MosaicTopology topo_ped_median;
};

inline MosaicTopoSetup initialize_mosaic_topology(int fpga_count, const std::string& mapping_json_file, int fpga_channel_number) {
    MosaicTopoSetup setup;
    
    std::ifstream mapping_json_ifs(mapping_json_file);
    if (!mapping_json_ifs.is_open()) {
        LOG(ERROR) << "Failed to open mapping JSON file: " << mapping_json_file;
        throw std::runtime_error("Cannot open mapping JSON file");
    }

    json mapping_json;
    mapping_json_ifs >> mapping_json;
    const auto& sipm_board      = mapping_json.at("SiPM_Board");
    const auto& board_loc       = mapping_json.at("Board_Loc");
    const auto& board_rotation  = mapping_json.at("Board_Rotation");
    const auto& board_flip      = mapping_json.at("Board_Flip");

    const int NX = 16, NY = 4;
    const int board_cols = 8, board_rows = 4;

    auto chan2pad = build_chan2pad_LUT(
        fpga_count, fpga_channel_number,
        NX, NY, board_cols, board_rows,
        sipm_board, board_loc, board_rotation, board_flip
    );

    setup.topo_wave.NX = NX;
    setup.topo_wave.NY = NY;
    setup.topo_wave.vldb_number = fpga_count;
    setup.topo_wave.channels_per_vldb = fpga_channel_number;
    setup.topo_wave.reverse_row = true;
    setup.topo_wave.minimalist_axis = true;
    setup.topo_wave.th2_logz = true;
    setup.topo_wave.chan2pad = chan2pad;

    setup.topo_ped_median = setup.topo_wave;
    setup.topo_ped_median.th2_logz = true;

    return setup;
}

#endif // COMMON_HPP