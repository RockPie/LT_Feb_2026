#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <numeric>

#include "TGraphErrors.h"
#include "TGraph.h"
#include "TF1.h"
#include "TSpline.h"
#include "TRandom3.h"
#include "Math/Interpolator.h"

// =============================================================================
// Data containers
// =============================================================================

struct TotToAdcBuildResult {
    // lut[t] = baseline-subtracted IDEAL ADC (can be >1023), stored as uint16_t
    std::vector<uint16_t> lut;   // size = tot_bins, index = decoded ToT code (0..tot_bins-1)
    double alpha = 0.0;          // ADC_raw_ideal = alpha * L + beta  (linear fit in RAW ADC domain)
    double beta  = 0.0;
    bool has_pit = false;
    double pit_L_start = 0.0;
    double pit_L_end   = 0.0;
};

struct TotToAdcLUT {
    // lut[t] = ADC_out (baseline-subtracted ideal ADC), stored as uint16_t
    std::vector<uint16_t> lut;
    int tot_bins = 0;

    // Optional clamp range for Eval output (kept wide by default)
    uint16_t adc_min = 0;
    uint16_t adc_max = 65535;

    // Clamp + lookup (integer ToT code)
    inline uint16_t Eval(int tot_code) const {
        if (lut.empty()) return 0;
        if (tot_code < 0) return lut.front();
        if (tot_code >= (int)lut.size()) return lut.back();
        uint16_t v = lut[(size_t)tot_code];
        if (v < adc_min) v = adc_min;
        if (v > adc_max) v = adc_max;
        return v;
    }

    // Optional: linear interpolation for non-integer ToT (double)
    inline double Eval(double tot) const {
        if (lut.empty()) return 0.0;
        if (!std::isfinite(tot)) return 0.0;

        if (tot <= 0.0) {
            double v = (double)lut.front();
            if (v < adc_min) v = adc_min;
            if (v > adc_max) v = adc_max;
            return v;
        }

        const double max_idx = (double)lut.size() - 1.0;
        if (tot >= max_idx) {
            double v = (double)lut.back();
            if (v < adc_min) v = adc_min;
            if (v > adc_max) v = adc_max;
            return v;
        }

        int i0 = (int)std::floor(tot);
        int i1 = i0 + 1;
        double frac = tot - (double)i0;

        double y0 = (double)lut[(size_t)i0];
        double y1 = (double)lut[(size_t)i1];
        double y  = y0 + frac * (y1 - y0);

        if (y < (double)adc_min) y = (double)adc_min;
        if (y > (double)adc_max) y = (double)adc_max;
        return y;
    }
};

// =============================================================================
// Helpers
// =============================================================================

static inline std::string TrimCopy(std::string s) {
    auto is_space = [](unsigned char c) { return std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [&](char c){ return !is_space((unsigned char)c); }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [&](char c){ return !is_space((unsigned char)c); }).base(), s.end());
    return s;
}

// Helper: sort points by x, return vectors
static void ExtractSortedXY(const TGraphErrors& g,
                            std::vector<double>& xs,
                            std::vector<double>& ys,
                            std::vector<double>& exs,
                            std::vector<double>& eys)
{
    const int N = g.GetN();
    xs.resize(N); ys.resize(N); exs.resize(N); eys.resize(N);
    for (int i = 0; i < N; i++) {
        double x, y;
        g.GetPoint(i, x, y);
        xs[i] = x; ys[i] = y;
        exs[i] = g.GetErrorX(i);
        eys[i] = g.GetErrorY(i);
    }

    std::vector<int> idx(N);
    for (int i = 0; i < N; i++) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b){ return xs[a] < xs[b]; });

    auto xs0 = xs, ys0 = ys, ex0 = exs, ey0 = eys;
    for (int k = 0; k < N; k++) {
        xs[k]  = xs0[idx[k]];
        ys[k]  = ys0[idx[k]];
        exs[k] = ex0[idx[k]];
        eys[k] = ey0[idx[k]];
    }
}

static void ExtractSortedXY(const TGraph& g,
                            std::vector<double>& xs,
                            std::vector<double>& ys,
                            std::vector<double>& exs,
                            std::vector<double>& eys) {
    const int n = g.GetN();
    xs.resize(n);
    ys.resize(n);
    exs.assign(n, 0.0);
    eys.assign(n, 0.0);

    for (int i = 0; i < n; i++) {
        g.GetPoint(i, xs[i], ys[i]);
    }

    std::vector<int> indices(n);
    for (int i = 0; i < n; i++) indices[i] = i;
    std::sort(indices.begin(), indices.end(), [&](int a, int b){ return xs[a] < xs[b]; });

    auto original_xs = xs;
    auto original_ys = ys;
    for (int i = 0; i < n; i++) {
        xs[i] = original_xs[indices[i]];
        ys[i] = original_ys[indices[i]];
    }
}

// Detect a “pit” interval in ToT(L): find first significant decrease, then end when it recovers above previous max.
// epsTot is in ToT units; tune if needed.
static bool DetectTotPitInterval(const std::vector<double>& L,
                                 const std::vector<double>& T,
                                 double epsTot,
                                 double& L_start,
                                 double& L_end)
{
    const int N = (int)L.size();
    if (N < 4) return false;

    double running_max = T[0];
    int pit_start_idx = -1;
    int pit_end_idx = -1;

    for (int i = 1; i < N; i++) {
        running_max = std::max(running_max, T[i-1]);
        // pit start: ToT drops below previous point AND below running max by eps
        if (pit_start_idx < 0) {
            if ((T[i] < T[i-1] - epsTot) && (T[i] < running_max - epsTot)) {
                pit_start_idx = i - 1; // start around the turning point
            }
        } else {
            // pit end: ToT recovers to (previous running max - eps)
            if (T[i] >= running_max - epsTot) {
                pit_end_idx = i;
                break;
            }
        }
    }

    if (pit_start_idx < 0) return false;
    if (pit_end_idx < 0) pit_end_idx = N - 1;

    L_start = L[pit_start_idx];
    L_end   = L[pit_end_idx];
    if (L_end < L_start) std::swap(L_start, L_end);
    return true;
}

// =============================================================================
// LUT builder
// =============================================================================
//
// IMPORTANT: 输出 LUT 存的是 “baseline-subtracted 的 ideal ADC”（可 >1023）
// - 线性模型是在 RAW ADC 域拟合：ADC_raw_ideal = alpha*L + beta
// - 最终输出：ADC_out = max(0, ADC_raw_ideal - baseline)
// - 不再对 1023 做上限截断（除非你主动给 adc_out_max）
//
// Pit：如果你仍希望 pit 映射到某个固定 RAW 值（比如 1023），用 pit_raw_override 控制。
// 如果你不想特殊处理 pit，把 pit_raw_override 设成 NaN 即可（默认 NaN）。
//
TotToAdcBuildResult BuildTotToAdcLUT_FromGraphs(
    const TGraphErrors* gAdcVsLaser,
    const TGraph* gTotVsLaser,
    int tot_bins = 4096,
    double adc_lin_min = 150.0,            // RAW ADC range used to select linear region
    double adc_lin_max = 950.0,            // RAW ADC range used to select linear region
    int laser_samples = 4000,
    double pit_eps_tot = 0.5,              // ToT decrease threshold to declare pit (units: ToT code)
    const std::string& dump_map_txt = "",  // if non-empty: write "ToT ADCraw_ideal ADCout L" (debug)
    double adc_baseline = 0.0,             // baseline in RAW ADC units, e.g. 100
    double adc_out_max = 65535.0,          // output clamp max (baseline-subtracted ideal ADC)
    double pit_raw_override = std::numeric_limits<double>::quiet_NaN() // e.g. 1023.0 if you want
) {
    TotToAdcBuildResult res;
    res.lut.assign((size_t)tot_bins, 0);

    if (!gAdcVsLaser || !gTotVsLaser || gAdcVsLaser->GetN() < 2 || gTotVsLaser->GetN() < 2) {
        return res;
    }

    if (adc_out_max < 0) adc_out_max = 0;
    if (adc_out_max > 65535.0) adc_out_max = 65535.0;

    // -------- 1) Fit ADCraw(L) linearly in [adc_lin_min, adc_lin_max] using points filtered by RAW ADC range --------
    std::vector<double> L_adc, A_adc, ex_adc, ey_adc;
    ExtractSortedXY(*gAdcVsLaser, L_adc, A_adc, ex_adc, ey_adc);

    TGraphErrors g_adc_fit;
    for (int i = 0; i < (int)L_adc.size(); i++) {
        const double Araw = A_adc[i];
        const double L    = L_adc[i];
        if (std::isfinite(Araw) && std::isfinite(L) && Araw >= adc_lin_min && Araw <= adc_lin_max) {
            int p = g_adc_fit.GetN();
            g_adc_fit.SetPoint(p, L, Araw);
            g_adc_fit.SetPointError(p, ex_adc[i], ey_adc[i]);
        }
    }

    TF1 f_lin("f_lin", "pol1");
    if (g_adc_fit.GetN() < 2) {
        TGraphErrors g_adc_all = *gAdcVsLaser;
        g_adc_all.Fit(&f_lin, "Q");
    } else {
        g_adc_fit.Fit(&f_lin, "Q");
    }

    // ADCraw_ideal = alpha * L + beta
    res.beta  = f_lin.GetParameter(0);
    res.alpha = f_lin.GetParameter(1);

    // -------- 2) Build spline for ToT(L): L -> ToT --------
    std::vector<double> L_tot, T_tot, ex_tot, ey_tot;
    ExtractSortedXY(*gTotVsLaser, L_tot, T_tot, ex_tot, ey_tot);

    TGraph g_LT((int)L_tot.size(), L_tot.data(), T_tot.data());
    TSpline3 spline_LT("spline_LT", &g_LT);

    // -------- 3) Detect pit interval in ToT(L) (non-monotonic) --------
    double pit_L_start = 0.0, pit_L_end = 0.0;
    res.has_pit = DetectTotPitInterval(L_tot, T_tot, pit_eps_tot, pit_L_start, pit_L_end);
    res.pit_L_start = pit_L_start;
    res.pit_L_end   = pit_L_end;

    // -------- 4) Sample laser densely, generate (ToT, ADCout) samples --------
    const double L_min = *std::min_element(L_tot.begin(), L_tot.end());
    const double L_max = *std::max_element(L_tot.begin(), L_tot.end());
    if (!(L_max > L_min)) return res;

    std::vector<std::pair<double, double>> samples; // (ToT, ADCout)
    samples.reserve((size_t)laser_samples + 1);

    std::ofstream fout;
    if (!dump_map_txt.empty()) fout.open(dump_map_txt);

    for (int i = 0; i <= laser_samples; i++) {
        const double u = (double)i / (double)laser_samples;
        const double L = L_min + (L_max - L_min) * u;

        const double T = spline_LT.Eval(L);
        if (!std::isfinite(T)) continue;

        // Ideal unsaturated RAW ADC from linear model
        double Araw_ideal = res.alpha * L + res.beta;
        if (!std::isfinite(Araw_ideal)) continue;
        if (Araw_ideal < 0) Araw_ideal = 0;

        bool in_pit = false;
        if (res.has_pit && (L >= res.pit_L_start && L <= res.pit_L_end)) {
            in_pit = true;
        }

        // Optional pit override in RAW domain (ONLY if you pass a finite value)
        if (in_pit && std::isfinite(pit_raw_override)) {
            Araw_ideal = pit_raw_override;
            if (Araw_ideal < 0) Araw_ideal = 0;
        }

        // Baseline subtraction on output (output is baseline-subtracted ideal ADC)
        double Aout = Araw_ideal - adc_baseline;
        if (Aout < 0) Aout = 0;
        if (Aout > adc_out_max) Aout = adc_out_max;

        samples.emplace_back(T, Aout);

        if (fout.is_open()) {
            // ToT  Araw_ideal  Aout  L  [PIT]
            fout << T << " " << Araw_ideal << " " << Aout << " " << L << (in_pit ? " PIT" : "") << "\n";
        }
    }

    if (fout.is_open()) fout.close();
    if (samples.empty()) return res;

    std::sort(samples.begin(), samples.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    // -------- 5) Build LUT: for each integer ToT code, pick ADCout of nearest sample in ToT --------
    for (int t = 0; t < tot_bins; t++) {
        const double Tq = (double)t;

        auto it = std::lower_bound(samples.begin(), samples.end(), std::make_pair(Tq, -1e300),
                                   [](const auto& a, const auto& b){ return a.first < b.first; });

        double best_adc = 0.0;
        double best_dT  = std::numeric_limits<double>::infinity();

        // check neighbors around lower_bound
        for (int k = -2; k <= 2; k++) {
            auto it2 = it;

            if (k < 0) {
                for (int j = 0; j < -k && it2 != samples.begin(); j++) --it2;
            } else if (k > 0) {
                for (int j = 0; j < k && it2 != samples.end(); j++) ++it2;
            }
            if (it2 == samples.end()) continue;

            const double Ts = it2->first;
            const double As = it2->second; // already baseline-subtracted ideal ADC
            const double dT = std::fabs(Ts - Tq);

            if (dT < best_dT - 1e-12) {
                best_dT  = dT;
                best_adc = As;
            } else if (std::fabs(dT - best_dT) < 1e-12) {
                // tie -> take larger ADC
                best_adc = std::max(best_adc, As);
            }
        }

        if (best_adc < 0) best_adc = 0;
        if (best_adc > adc_out_max) best_adc = adc_out_max;

        res.lut[(size_t)t] = (uint16_t)std::lround(best_adc);
    }

    return res;
}

// =============================================================================
// LUT loader (reads "ToT ADC" pairs OR ADC-only list)
// =============================================================================
//
// CRITICAL FIX vs your old version:
// - 不再把 ADC clamp 到 1023（默认 adc_max=65535）
// - 不再 “把 0 forward-fill” （0 可能是真实值）
// - 如果文件是 sparse pairs，你可以选择 fill_missing_with_last=true 才做 forward fill
//
inline TotToAdcLUT LoadTotToAdcLUT(const std::string& filepath,
                                  int expected_tot_bins = 4096,
                                  uint16_t adc_min = 0,
                                  uint16_t adc_max = 65535,
                                  bool fill_missing_with_last = false)
{
    std::ifstream fin(filepath);
    if (!fin.is_open()) {
        throw std::runtime_error("LoadTotToAdcLUT: failed to open file: " + filepath);
    }

    TotToAdcLUT out;
    out.adc_min = adc_min;
    out.adc_max = adc_max;

    std::string line;
    std::vector<std::pair<int, long long>> pairs; // (tot, adc)
    std::vector<long long> adc_only;

    while (std::getline(fin, line)) {
        // strip comments (# or //)
        auto pos_hash = line.find('#');
        if (pos_hash != std::string::npos) line = line.substr(0, pos_hash);
        auto pos_slash = line.find("//");
        if (pos_slash != std::string::npos) line = line.substr(0, pos_slash);

        line = TrimCopy(line);
        if (line.empty()) continue;

        std::istringstream iss(line);

        long long a = 0, b = 0;
        if (!(iss >> a)) continue;

        if (iss >> b) {
            // tot adc
            pairs.emplace_back((int)a, b);
        } else {
            // adc only
            adc_only.push_back(a);
        }
    }

    if (!pairs.empty()) {
        int max_tot = -1;
        for (const auto& kv : pairs) max_tot = std::max(max_tot, kv.first);

        int bins = expected_tot_bins > 0 ? expected_tot_bins : (max_tot + 1);
        bins = std::max(bins, max_tot + 1);

        out.lut.assign((size_t)bins, 0);

        // fill provided points
        for (const auto& kv : pairs) {
            const int tot = kv.first;
            long long adc = kv.second;
            if (tot < 0 || tot >= bins) continue;

            if (adc < (long long)adc_min) adc = (long long)adc_min;
            if (adc > (long long)adc_max) adc = (long long)adc_max;

            out.lut[(size_t)tot] = (uint16_t)adc;
        }

        // optional forward fill for sparse LUT
        if (fill_missing_with_last && !out.lut.empty()) {
            uint16_t last = out.lut.front();
            for (size_t i = 0; i < out.lut.size(); i++) {
                // treat "unset" as exactly 0 AND i!=0 AND last!=0 (conservative)
                if (out.lut[i] == 0 && i != 0 && last != 0) out.lut[i] = last;
                last = out.lut[i];
            }
        }

        out.tot_bins = (int)out.lut.size();
        return out;
    }

    if (!adc_only.empty()) {
        out.lut.resize(adc_only.size());
        for (size_t i = 0; i < adc_only.size(); i++) {
            long long adc = adc_only[i];
            if (adc < (long long)adc_min) adc = (long long)adc_min;
            if (adc > (long long)adc_max) adc = (long long)adc_max;
            out.lut[i] = (uint16_t)adc;
        }
        out.tot_bins = (int)out.lut.size();
        return out;
    }

    throw std::runtime_error("LoadTotToAdcLUT: no valid data found in: " + filepath);
}

static double LinInterp(const std::vector<double>& xs,
                        const std::vector<double>& ys,
                        double xq)
{
    if (xs.size() < 2) return NAN;
    if (xq < xs.front() || xq > xs.back()) return NAN;

    auto it = std::upper_bound(xs.begin(), xs.end(), xq);
    size_t j = std::max<size_t>(1, std::distance(xs.begin(), it)) - 1;

    double x0 = xs[j],   x1 = xs[j+1];
    double y0 = ys[j],   y1 = ys[j+1];
    double t  = (xq - x0) / (x1 - x0);
    return y0 + t * (y1 - y0);
}

static void GraphToSortedArrays(const TGraphErrors* g,
                                std::vector<double>& xs,
                                std::vector<double>& ys,
                                std::vector<double>& exs,
                                std::vector<double>& eys)
{
    xs.clear(); ys.clear(); exs.clear(); eys.clear();
    const int n = g ? g->GetN() : 0;
    xs.reserve(n); ys.reserve(n); exs.reserve(n); eys.reserve(n);

    std::vector<std::tuple<double,double,double,double>> tmp;
    tmp.reserve(n);

    for (int i = 0; i < n; ++i) {
        double x, y;
        g->GetPoint(i, x, y);
        double ex = g->GetErrorX(i);
        double ey = g->GetErrorY(i);
        tmp.emplace_back(x, y, ex, ey);
    }

    std::sort(tmp.begin(), tmp.end(),
              [](auto& a, auto& b){ return std::get<0>(a) < std::get<0>(b); });

    for (auto& t : tmp) {
        xs.push_back(std::get<0>(t));
        ys.push_back(std::get<1>(t));
        exs.push_back(std::get<2>(t));
        eys.push_back(std::get<3>(t));
    }
}

static std::pair<double,int> Similarity_Chi2NDF_ScaledX(
    const TGraphErrors* gRef,
    const TGraphErrors* gMov,
    double X0,
    double s,
    bool useSplineForRefY = true
){
    if (!gRef || !gMov) return {NAN, 0};
    if (gRef->GetN() < 2 || gMov->GetN() < 1) return {NAN, 0};

    std::vector<double> xR, yR, exR, eyR;
    GraphToSortedArrays(gRef, xR, yR, exR, eyR);

    double xMin = xR.front();
    double xMax = xR.back();

    std::unique_ptr<TSpline3> sp;
    if (useSplineForRefY) {
        TGraph gTmp((int)xR.size(), xR.data(), yR.data());
        sp.reset(new TSpline3("sp_ref_y", &gTmp));
    }

    auto evalRefY = [&](double xq)->double {
        if (xq < xMin || xq > xMax) return NAN;
        if (useSplineForRefY) return sp->Eval(xq);
        return LinInterp(xR, yR, xq);
    };

    auto evalRefEy = [&](double xq)->double {
        return LinInterp(xR, eyR, xq);
    };
    double chi2 = 0.0;
    int nUsed = 0;

    const int nM = gMov->GetN();
    for (int i = 0; i < nM; ++i) {
        double x, y;
        gMov->GetPoint(i, x, y);
        double eyM = gMov->GetErrorY(i);

        double xScaled = X0 + s * (x - X0);
        if (xScaled < xMin || xScaled > xMax) continue;

        double yRef = evalRefY(xScaled);
        if (!std::isfinite(yRef)) continue;

        double eyR_i = evalRefEy(xScaled);
        if (!std::isfinite(eyR_i)) eyR_i = 0.0;

        double sigma2 = eyM*eyM + eyR_i*eyR_i;
        if (sigma2 <= 0) continue;

        double r = (y - yRef) / std::sqrt(sigma2);
        chi2 += r*r;
        nUsed++;
    }

    if (nUsed < 2) return {NAN, nUsed};
    double chi2ndf = chi2 / (nUsed - 1);
    return {chi2ndf, nUsed};
}

TGraph* ScanScaleSimilarity(
    const TGraphErrors* gRef,
    const TGraphErrors* gMov,
    double X0,
    double sMin,
    double sMax,
    int nSteps,
    bool useSplineForRefY = true
){
    if (nSteps < 2) return nullptr;

    auto out = new TGraph();
    out->SetName("g_scale_similarity");
    out->SetTitle("Similarity vs x-scale;scale s;#chi^{2}/NDF");

    for (int k = 0; k < nSteps; ++k) {
        // Use logarithmic spacing instead of linear
        double s = sMin * std::pow(sMax / sMin, (double)k / (double)(nSteps - 1));
        auto [chi2ndf, nUsed] = Similarity_Chi2NDF_ScaledX(gRef, gMov, X0, s, useSplineForRefY);

        if (!std::isfinite(chi2ndf)) continue;

        out->SetPoint(out->GetN(), s, chi2ndf);
    }
    return out;
}

TGraphErrors* ScaleX_TGraphErrors(
    const TGraphErrors* g,
    double s,
    double X0,
    const char* newName = "g_scaled"
){
    if (!g) return nullptr;

    int n = g->GetN();
    if (n <= 0) return nullptr;

    auto* gOut = new TGraphErrors(n);
    gOut->SetName(newName);

    for (int i = 0; i < n; ++i) {
        double x, y;
        g->GetPoint(i, x, y);

        double ex = g->GetErrorX(i);
        double ey = g->GetErrorY(i);

        // 进行放缩
        double xScaled  = X0 + s * (x - X0);
        double exScaled = std::abs(s) * ex;

        gOut->SetPoint(i, xScaled, y);
        gOut->SetPointError(i, exScaled, ey);
    }

    gOut->SetTitle(g->GetTitle());
    return gOut;
}

static inline void DedupByX(std::vector<double>& xs,
                     std::vector<double>& ys,
                     std::vector<double>& exs,
                     std::vector<double>& eys,
                     double eps = 1e-12)
{
    if (xs.empty()) return;
    std::vector<double> nx, ny, nex, ney;
    nx.reserve(xs.size()); ny.reserve(xs.size());
    nex.reserve(xs.size()); ney.reserve(xs.size());

    nx.push_back(xs[0]);
    ny.push_back(ys[0]);
    nex.push_back(exs[0]);
    ney.push_back(eys[0]);

    for (size_t i=1;i<xs.size();++i){
        if (std::abs(xs[i]-nx.back()) < eps){
            // skip duplicate-x point
            continue;
        }
        nx.push_back(xs[i]);
        ny.push_back(ys[i]);
        nex.push_back(exs[i]);
        ney.push_back(eys[i]);
    }
    xs.swap(nx); ys.swap(ny); exs.swap(nex); eys.swap(ney);
}

inline TGraphErrors* InterpolateWithUncertainty_AkimaMC(const TGraphErrors* gin,
                                                 int nFine = 400,
                                                 int nMC   = 800,
                                                 bool includeXUnc = true,
                                                 UInt_t seed = 0)
{
    if (!gin || gin->GetN() < 3) {
        throw std::runtime_error("Need at least 3 points for Akima interpolation.");
    }
    if (nFine < 2) nFine = 2;
    if (nMC < 50)  nMC = 50;

    // 1) 取点并按 x 排序
    std::vector<std::tuple<double,double,double,double>> pts;
    pts.reserve(gin->GetN());
    for (int i=0;i<gin->GetN();++i){
        double x,y; gin->GetPoint(i,x,y);
        double ex = gin->GetErrorX(i);
        double ey = gin->GetErrorY(i);
        pts.emplace_back(x,y,ex,ey);
    }
    std::sort(pts.begin(), pts.end(),
              [](auto& a, auto& b){ return std::get<0>(a) < std::get<0>(b); });

    std::vector<double> xs, ys, exs, eys;
    xs.reserve(pts.size()); ys.reserve(pts.size());
    exs.reserve(pts.size()); eys.reserve(pts.size());
    for (auto &t: pts){
        xs.push_back(std::get<0>(t));
        ys.push_back(std::get<1>(t));
        exs.push_back(std::get<2>(t));
        eys.push_back(std::get<3>(t));
    }

    DedupByX(xs, ys, exs, eys);
    if (xs.size() < 3) throw std::runtime_error("Too few unique-x points after dedup.");

    // 2) 定义输出 x 网格（均匀）
    const double xmin = xs.front();
    const double xmax = xs.back();
    std::vector<double> xq(nFine);
    for (int i=0;i<nFine;++i){
        xq[i] = xmin + (xmax-xmin)*i/(nFine-1.0);
    }

    // 3) MC 累积：对每个 xq 累加 y 和 y^2
    std::vector<double> sumY(nFine, 0.0), sumY2(nFine, 0.0);

    TRandom3 rng(seed);

    // MC loop
    for (int k=0;k<nMC;++k){
        // 3.1) 采样一组“伪数据”
        std::vector<std::pair<double,double>> spts;
        spts.reserve(xs.size());
        for (size_t i=0;i<xs.size();++i){
            double xk = xs[i];
            double yk = ys[i];

            const double ex = exs[i];
            const double ey = eys[i];

            if (includeXUnc && ex > 0) xk = rng.Gaus(xs[i], ex);
            if (ey > 0) yk = rng.Gaus(ys[i], ey);

            spts.push_back({xk, yk});
        }

        // 3.2) 排序 + 去重（MC 后 x 可能乱序/重复）
        std::sort(spts.begin(), spts.end(),
                  [](auto&a, auto&b){ return a.first < b.first; });

        std::vector<double> sx, sy;
        sx.reserve(spts.size()); sy.reserve(spts.size());
        for (auto &p: spts){
            if (!sx.empty() && std::abs(p.first - sx.back()) < 1e-12) continue;
            sx.push_back(p.first);
            sy.push_back(p.second);
        }
        if (sx.size() < 3) continue; // 极端情况下跳过

        // 3.3) 构造 Akima 插值器并评估
        ROOT::Math::Interpolator itp(sx, sy, ROOT::Math::Interpolation::kAKIMA);

        for (int i=0;i<nFine;++i){
            double yv = itp.Eval(xq[i]);
            sumY[i]  += yv;
            sumY2[i] += yv*yv;
        }
    }

    // 4) 汇总均值与方差 -> 输出 TGraphErrors
    auto gout = new TGraphErrors(nFine);
    for (int i=0;i<nFine;++i){
        const double mean = sumY[i] / nMC;
        double var = sumY2[i] / nMC - mean*mean;
        if (var < 0) var = 0; // 数值保护
        const double ey = std::sqrt(var);

        gout->SetPoint(i, xq[i], mean);

        // 输出 ex 怎么设：你可以设 0，或者设为网格半步长（代表“bin宽”概念）
        double ex_out = 0.0;
        // double ex_out = (xmax-xmin)/(nFine-1.0)/2.0;
        gout->SetPointError(i, ex_out, ey);
    }

    return gout;
}