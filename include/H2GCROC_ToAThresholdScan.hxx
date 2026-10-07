#ifndef H2GCROC_TOA_THRESHOLD_SCAN_HXX
#define H2GCROC_TOA_THRESHOLD_SCAN_HXX

// Header-only, C++11 or newer. No changes to the original ToA conversion units.
// Define H2GCROC_TOA_SCAN_NO_ROOT before including this file to test the core
// algorithm without ROOT. Do NOT define it in the H2GCROC analysis program.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace h2g_toa {

struct Config {
    int firstSample = 2;
    int lastSample = 5;
    double sampleNs = 25.0;
    double lsbNs = 0.025;
    int binWidthTicks = 10;       // 10 * 0.025 = 0.25 ns; identical for all thresholds
    double smoothSigmaNs = 1.0;  // fixed Gaussian sigma, NOT chosen separately per threshold
    std::uint64_t lowStatsEntries = 200; // warning only; still scan every nonempty channel
    double tieRelativeTolerance = 1e-10;
    double tieAbsoluteTolerance = 1e-20;
};

// Sufficient statistics: one first valid ToA (sample, code) per channel/event.
// The sample index must be retained; a raw-code histogram alone is insufficient.
class Counts {
public:
    Counts(int firstSample, int lastSample)
        : first_(firstSample), last_(lastSample), entries_(0) {
        if (first_ < 0 || last_ < first_ || last_ > 100000)
            throw std::invalid_argument("Invalid ToA sample range");
        data_.resize(static_cast<std::size_t>(last_ - first_ + 1));
        for (auto& row : data_) row.fill(0);
    }
    void add(int sample, int raw, std::uint64_t weight = 1) {
        // Preserve the input program's convention: raw code zero means no ToA.
        if (raw == 0) return;
        if (sample < first_ || sample > last_ || raw < 1 || raw > 1023)
            throw std::out_of_range("Invalid first-ToA sample or raw code");
        if (weight > std::numeric_limits<std::uint64_t>::max() - entries_)
            throw std::overflow_error("ToA count overflow");
        data_[static_cast<std::size_t>(sample - first_)][raw] += weight;
        entries_ += weight;
    }
    int firstSample() const { return first_; }
    int lastSample() const { return last_; }
    std::uint64_t entries() const { return entries_; }
    std::uint64_t at(int sample, int raw) const {
        return data_.at(static_cast<std::size_t>(sample - first_)).at(raw);
    }
private:
    int first_, last_;
    std::uint64_t entries_;
    std::vector<std::array<std::uint64_t, 1024>> data_;
};

struct Result {
    int threshold = -1; // -1 = no usable ToA; NEVER a fixed-code fallback
    std::uint64_t entries = 0;
    bool lowStatistics = false;
    bool allThresholdsEquivalent = false;
    int equivalentMinima = 0;
    int plateauFirst = -1;
    int plateauLast = -1;
    double minimumScore = std::numeric_limits<double>::quiet_NaN();
    double chosenScore = std::numeric_limits<double>::quiet_NaN();
    double scoreContrast = 0.0;
    double xminNs = 0.0, xmaxNs = 0.0, binWidthNs = 0.0;
    std::vector<double> scores;
    std::vector<double> before, after, smoothBefore, smoothAfter;
};

namespace detail {
inline long long floorDiv(long long a, long long b) {
    long long q = a / b;
    if (a % b < 0) --q;
    return q;
}

struct Geometry {
    long long ticksPerSample = 0;
    long long firstBinNumber = 0;
    int nbins = 0, radius = 0, binWidthTicks = 0;
    double binWidthNs = 0.0, xminNs = 0.0, xmaxNs = 0.0;
    std::vector<double> kernel;

    int bin(int sample, int raw, bool subtractCycle) const {
        const long long t = static_cast<long long>(sample) * ticksPerSample
                          + raw - (subtractCycle ? ticksPerSample : 0);
        const long long i = floorDiv(t, binWidthTicks) - firstBinNumber;
        if (i < radius + 2 || i >= nbins - radius - 2)
            throw std::logic_error("ToA lies outside the padded scan axis");
        return static_cast<int>(i);
    }
};

inline Geometry geometry(const Config& c) {
    if (c.firstSample < 0 || c.lastSample < c.firstSample || c.lastSample > 100000 ||
        !(c.sampleNs > 0.0) || !std::isfinite(c.sampleNs) ||
        !(c.lsbNs > 0.0) || !std::isfinite(c.lsbNs) || c.binWidthTicks < 1 ||
        !(c.smoothSigmaNs > 0.0) || !std::isfinite(c.smoothSigmaNs) ||
        !(c.tieRelativeTolerance >= 0.0) || !std::isfinite(c.tieRelativeTolerance) ||
        !(c.tieAbsoluteTolerance >= 0.0) || !std::isfinite(c.tieAbsoluteTolerance))
        throw std::invalid_argument("Invalid ToA scan configuration");
    Geometry g;
    const double ratio = c.sampleNs / c.lsbNs;
    if (!std::isfinite(ratio) || ratio < 1.0 || ratio > 1e9)
        throw std::invalid_argument("Invalid sampleNs / lsbNs");
    g.ticksPerSample = std::llround(ratio);
    if (std::abs(ratio - static_cast<double>(g.ticksPerSample)) > 1e-8)
        throw std::invalid_argument("sampleNs must be an integer multiple of lsbNs");
    // 25 ns / 0.025 ns = 1000 ticks, NOT 1024 ticks.
    if (g.ticksPerSample % c.binWidthTicks != 0)
        throw std::invalid_argument("binWidthTicks must divide the ticks in one sample");
    g.binWidthTicks = c.binWidthTicks;
    g.binWidthNs = c.binWidthTicks * c.lsbNs;
    const double radius = std::ceil(4.0 * c.smoothSigmaNs / g.binWidthNs);
    if (!std::isfinite(radius) || radius > 100000)
        throw std::invalid_argument("Smoothing kernel is too large");
    g.radius = std::max(1, static_cast<int>(radius));
    const int padding = g.radius + 3;
    // Include ALL possible converted and uncorrected times and the entire kernel.
    const long long lo = (static_cast<long long>(c.firstSample) - 1) * g.ticksPerSample;
    const long long hi = static_cast<long long>(c.lastSample) * g.ticksPerSample + 1023;
    g.firstBinNumber = floorDiv(lo, c.binWidthTicks) - padding;
    const long long last = floorDiv(hi, c.binWidthTicks) + padding;
    const long long bins = last - g.firstBinNumber + 1;
    if (bins < 5 || bins > 1000000)
        throw std::invalid_argument("Invalid or excessively large ToA histogram");
    g.nbins = static_cast<int>(bins);
    g.xminNs = static_cast<double>(g.firstBinNumber) * g.binWidthNs;
    g.xmaxNs = g.xminNs + g.nbins * g.binWidthNs;
    g.kernel.resize(2 * g.radius + 1);
    for (int d = -g.radius; d <= g.radius; ++d) {
        const double z = d * g.binWidthNs / c.smoothSigmaNs;
        g.kernel[d + g.radius] = std::exp(-0.5 * z * z);
    }
    const double norm = std::accumulate(g.kernel.begin(), g.kernel.end(), 0.0);
    for (double& v : g.kernel) v /= norm;
    return g;
}

inline void addSmoothed(std::vector<double>& smoothed, int bin, double weight,
                        const Geometry& g) {
    for (int d = -g.radius; d <= g.radius; ++d) {
        const int i = bin + d;
        if (i >= 0 && i < g.nbins)
            smoothed[i] += weight * g.kernel[d + g.radius];
    }
}

inline std::vector<double> smooth(const std::vector<double>& hist, const Geometry& g) {
    std::vector<double> out(g.nbins, 0.0);
    for (int b = 0; b < g.nbins; ++b)
        if (hist[b] != 0.0) addSmoothed(out, b, hist[b], g);
    return out;
}

// S(T) = sum_i [p_s(i-1) - 2 p_s(i) + p_s(i+1)]^2 / sum_i p_s(i)^2.
// p_s is a Gaussian-smoothed probability MASS per fixed-width bin.
// The L2 denominator avoids rewarding mere peak-height dilution/splitting.
// No Gaussian fit, variable range, peak-window cut, or width minimization is used.
inline double roughness(const std::vector<double>& smoothed, std::uint64_t n) {
    if (n == 0) return std::numeric_limits<double>::quiet_NaN();
    const double invN = 1.0 / static_cast<double>(n);
    double value = 0.0, norm2 = 0.0;
    for (double bin : smoothed) {
        const double p = bin * invN;
        norm2 += p * p;
    }
    for (std::size_t i = 1; i + 1 < smoothed.size(); ++i) {
        const double d2 = (smoothed[i-1] - 2.0 * smoothed[i] + smoothed[i+1]) * invN;
        value += d2 * d2;
    }
    if (!(norm2 > 0.0)) throw std::runtime_error("Empty smoothed ToA distribution");
    return value / norm2;
}

inline std::vector<double> histogram(const Counts& counts, const Geometry& g, int t) {
    // t=1024 is the uncorrected diagnostic baseline, NOT an optimization candidate.
    std::vector<double> h(g.nbins, 0.0);
    for (int s = counts.firstSample(); s <= counts.lastSample(); ++s)
        for (int raw = 1; raw < 1024; ++raw)
            h[g.bin(s, raw, raw >= t)] += static_cast<double>(counts.at(s, raw));
    return h;
}
} // namespace detail

using ScanObserver = std::function<void(int, const std::vector<double>&)>;

inline Result scan(const Counts& counts, const Config& config,
                   const ScanObserver& observer = ScanObserver()) {
    if (counts.firstSample() != config.firstSample || counts.lastSample() != config.lastSample)
        throw std::invalid_argument("Counts/config sample windows do not match");
    const detail::Geometry g = detail::geometry(config);
    Result r;
    r.entries = counts.entries();
    r.lowStatistics = r.entries < config.lowStatsEntries;
    r.xminNs = g.xminNs; r.xmaxNs = g.xmaxNs; r.binWidthNs = g.binWidthNs;
    r.scores.assign(1024, std::numeric_limits<double>::quiet_NaN());
    r.before = detail::histogram(counts, g, 1024);
    r.smoothBefore = detail::smooth(r.before, g);
    if (r.entries == 0) return r;

    std::vector<double> current = detail::histogram(counts, g, 0);
    std::vector<double> smoothed = detail::smooth(current, g);
    // Scan exactly 0,1,...,1023. Changing T to T+1 only moves raw code T
    // from the shifted location back to the unshifted location, for each sample.
    for (int t = 0; t < 1024; ++t) {
        if (t > 0) {
            const int raw = t - 1;
            for (int s = config.firstSample; s <= config.lastSample; ++s) {
                const double n = static_cast<double>(counts.at(s, raw));
                if (n == 0.0) continue;
                const int from = g.bin(s, raw, true), to = g.bin(s, raw, false);
                current[from] -= n;
                current[to] += n;
                detail::addSmoothed(smoothed, from, -n, g);
                detail::addSmoothed(smoothed, to, n, g);
            }
        }
        r.scores[t] = detail::roughness(smoothed, r.entries);
        if (observer) observer(t, current);
    }
    r.minimumScore = *std::min_element(r.scores.begin(), r.scores.end());
    const double maxScore = *std::max_element(r.scores.begin(), r.scores.end());
    const double tol = config.tieAbsoluteTolerance
                     + config.tieRelativeTolerance * std::abs(r.minimumScore);
    // Resolve numerical/empty-code plateaus without using any prior threshold.
    // Choose the middle of the longest contiguous minimum plateau. If two
    // plateaus have the same length, choose the lower-code one deterministically.
    int runStart = -1, bestLength = 0;
    for (int t = 0; t <= 1024; ++t) {
        const bool atMin = t < 1024 && r.scores[t] <= r.minimumScore + tol;
        if (atMin) {
            ++r.equivalentMinima;
            if (runStart < 0) runStart = t;
        } else if (runStart >= 0) {
            const int length = t - runStart;
            if (length > bestLength) {
                bestLength = length;
                r.plateauFirst = runStart;
                r.plateauLast = t - 1;
            }
            runStart = -1;
        }
    }
    if (bestLength == 0) throw std::runtime_error("No finite threshold score");
    r.threshold = r.plateauFirst + (r.plateauLast - r.plateauFirst) / 2;
    r.chosenScore = r.scores[r.threshold];
    r.allThresholdsEquivalent = r.equivalentMinima == 1024;
    r.scoreContrast = maxScore > 0.0 ? (maxScore - r.minimumScore) / maxScore : 0.0;
    r.after = detail::histogram(counts, g, r.threshold);
    r.smoothAfter = detail::smooth(r.after, g);
    return r;
}

// raw==0 is not a valid measurement. Invalid calibration (-1) explicitly leaves
// the time uncorrected; it never substitutes a fixed threshold.
inline double correctedNs(int sample, int raw, int threshold, const Config& c) {
    if (sample < 0 || raw < 1 || raw > 1023 || threshold < -1 || threshold > 1023)
        throw std::invalid_argument("Invalid ToA conversion argument");
    return sample * c.sampleNs + raw * c.lsbNs
         - ((threshold >= 0 && raw >= threshold) ? c.sampleNs : 0.0);
}

} // namespace h2g_toa

#ifndef H2GCROC_TOA_SCAN_NO_ROOT
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TPad.h"
#include "TROOT.h"
#include <memory>
#include <sstream>

namespace h2g_toa {

inline bool passesHamming(const std::vector<UInt_t*>& daqh, int samples) {
    if (samples <= 0) throw std::invalid_argument("Invalid number of samples");
    for (const UInt_t* words : daqh) {
        if (!words) throw std::invalid_argument("Null DAQ-header buffer");
        for (int s = 0; s < samples; ++s)
            for (int i = 0; i < 4; ++i)
                if (((words[4*s + i] >> 4) & 0x7u) != 0u) return false;
    }
    return true;
}

namespace detail {
inline std::unique_ptr<TH1D> rootHist(const char* name, const char* title,
                                    const Result& r, const std::vector<double>& values,
                                    bool isSmoothed = false) {
    std::unique_ptr<TH1D> h(new TH1D(name, title, static_cast<int>(values.size()),
                                    r.xminNs, r.xmaxNs));
    h->SetDirectory(nullptr);
    h->SetStats(kFALSE);
    for (std::size_t i = 0; i < values.size(); ++i) {
        h->SetBinContent(static_cast<int>(i) + 1, values[i]);
        // Smoothed bins are correlated; do not attach Poisson errors to them.
        h->SetBinError(static_cast<int>(i) + 1,
                       isSmoothed ? 0.0 : std::sqrt(std::max(0.0, values[i])));
    }
    h->ResetStats();
    h->SetEntries(static_cast<double>(r.entries));
    return h;
}
inline void writeOrThrow(TObject& o) {
    if (o.Write(o.GetName(), TObject::kOverwrite) <= 0)
        throw std::runtime_error(std::string("Cannot write ROOT object: ") + o.GetName());
}
} // namespace detail

// All channels get ROOT diagnostics. A nonempty pdfPrefix additionally exports
// three separate PDF figures for this channel. Directory ownership is restored.
inline Result scanAndWrite(const Counts& counts, const Config& config,
                           TDirectory& parent, int globalChannel,
                           const std::string& pdfPrefix = std::string()) {
    const std::string channelName = "Channel_" + std::to_string(globalChannel);
    TDirectory* directory = parent.GetDirectory(channelName.c_str());
    if (!directory) directory = parent.mkdir(channelName.c_str());
    if (!directory) throw std::runtime_error("Cannot create channel scan directory");
    TDirectory::TContext context(directory);
    const detail::Geometry geom = detail::geometry(config);
    std::unique_ptr<TH2D> map;
    if (counts.entries() > 0) {
        map.reset(new TH2D("toa_vs_threshold",
            "Converted ToA for every threshold;Threshold code;Converted ToA [ns];Events / bin",
            1024, -0.5, 1023.5, geom.nbins, geom.xminNs, geom.xmaxNs));
        map->SetDirectory(nullptr);
        map->SetStats(kFALSE);
    }
    const ScanObserver observer = [&](int threshold, const std::vector<double>& hist) {
        for (int b = 0; b < geom.nbins; ++b)
            map->SetBinContent(threshold + 1, b + 1, hist[b]);
    };
    Result r = scan(counts, config, map ? observer : ScanObserver());

    TParameter<int> threshold("best_threshold", r.threshold);
    TParameter<Long64_t> entries("n_valid_toa", static_cast<Long64_t>(r.entries));
    TParameter<int> lowStats("low_statistics", r.lowStatistics ? 1 : 0);
    TParameter<int> equivalent("equivalent_minima", r.equivalentMinima);
    TParameter<int> plateauFirst("minimum_plateau_first", r.plateauFirst);
    TParameter<int> plateauLast("minimum_plateau_last", r.plateauLast);
    TParameter<int> flat("all_thresholds_equivalent", r.allThresholdsEquivalent ? 1 : 0);
    TParameter<double> minimum("minimum_score", r.minimumScore);
    TParameter<double> selected("selected_score", r.chosenScore);
    TParameter<double> contrast("score_contrast", r.scoreContrast);
    std::ostringstream desc;
    desc.precision(17);
    desc << "t=sample*" << config.sampleNs << "+raw*" << config.lsbNs
         << "-(raw>=T ? " << config.sampleNs << " : 0); raw=0 excluded; sample=["
         << config.firstSample << "," << config.lastSample << "]; bin_ns=" << r.binWidthNs
         << "; Gaussian_sigma_ns=" << config.smoothSigmaNs
         << "; kernel=+/-4sigma; score=sum(second_difference(p)^2)/sum(p^2),p=smoothed_counts/N"
         << "; tie_relative_tolerance=" << config.tieRelativeTolerance
         << "; tie_absolute_tolerance=" << config.tieAbsoluteTolerance
         << "; selected=middle_of_longest_minimum_plateau; raw histograms NOT smoothed"
         << "; low_statistics_warning_below=" << config.lowStatsEntries;
    TNamed metadata("scan_definition", desc.str().c_str());
    detail::writeOrThrow(threshold); detail::writeOrThrow(entries);
    detail::writeOrThrow(lowStats); detail::writeOrThrow(equivalent);
    detail::writeOrThrow(plateauFirst); detail::writeOrThrow(plateauLast);
    detail::writeOrThrow(flat); detail::writeOrThrow(minimum);
    detail::writeOrThrow(selected); detail::writeOrThrow(contrast);
    detail::writeOrThrow(metadata);
    if (r.entries == 0) {
        TNamed status("status", "No valid ToA: threshold=-1; no scan/fit/conversion inferred.");
        detail::writeOrThrow(status);
        return r;
    }

    TH2D joint("first_toa_raw_vs_sample", "First valid ToA;Raw ToA code;Sample index;Events",
               1024, -0.5, 1023.5, config.lastSample - config.firstSample + 1,
               config.firstSample - 0.5, config.lastSample + 0.5);
    joint.SetDirectory(nullptr); joint.SetStats(kFALSE);
    for (int sample = config.firstSample; sample <= config.lastSample; ++sample)
        for (int raw = 1; raw <= 1023; ++raw)
            joint.SetBinContent(raw + 1, sample - config.firstSample + 1,
                                static_cast<double>(counts.at(sample, raw)));
    joint.ResetStats(); joint.SetEntries(static_cast<double>(r.entries));
    detail::writeOrThrow(joint);

    // Every x-bin is one threshold distribution, not an independent population.
    map->ResetStats();
    map->SetEntries(1024.0 * static_cast<double>(r.entries));
    detail::writeOrThrow(*map);
    std::array<double, 1024> x;
    for (int t = 0; t < 1024; ++t) x[t] = t;
    TGraph score(1024, x.data(), r.scores.data());
    score.SetName("smoothness_score");
    score.SetTitle("ToA threshold scan;Threshold code;Roughness score (lower is smoother)");
    score.SetLineWidth(2);
    detail::writeOrThrow(score);

    auto before = detail::rootHist("toa_before", "Before wrap correction;ToA [ns];Events / bin", r, r.before);
    auto after = detail::rootHist("toa_after", "After best-threshold correction;ToA [ns];Events / bin", r, r.after);
    auto smoothBefore = detail::rootHist("toa_before_smoothed", "Smoothed before;ToA [ns];Smoothed events / bin", r, r.smoothBefore, true);
    auto smoothAfter = detail::rootHist("toa_after_smoothed", "Smoothed after;ToA [ns];Smoothed events / bin", r, r.smoothAfter, true);
    detail::writeOrThrow(*before); detail::writeOrThrow(*after);
    detail::writeOrThrow(*smoothBefore); detail::writeOrThrow(*smoothAfter);

    const std::string suffix = "_channel_" + std::to_string(globalChannel);
    {
        TCanvas canvas(("canvas_score" + suffix).c_str(), "Threshold smoothness", 1000, 650);
        canvas.SetLeftMargin(0.14); canvas.SetBottomMargin(0.12);
        score.Draw("AL");
        score.GetXaxis()->SetLimits(-0.5, 1023.5);
        const double ymax = *std::max_element(r.scores.begin(), r.scores.end());
        score.SetMinimum(0.0); score.SetMaximum(ymax > 0.0 ? 1.15 * ymax : 1.0);
        TLine line(r.threshold, 0.0, r.threshold, ymax > 0.0 ? 1.15 * ymax : 1.0);
        line.SetLineColor(kRed + 1); line.SetLineStyle(2); line.SetLineWidth(2); line.Draw();
        TLegend legend(0.48, 0.74, 0.89, 0.89);
        legend.SetFillStyle(0); legend.SetBorderSize(0);
        const std::string label = "Best T=" + std::to_string(r.threshold) + ", plateau ["
            + std::to_string(r.plateauFirst) + "," + std::to_string(r.plateauLast) + "]";
        legend.AddEntry(&line, label.c_str(), "l");
        const std::string nlabel = "Channel " + std::to_string(globalChannel)
            + ", N=" + std::to_string(r.entries) + (r.lowStatistics ? " (LOW STATISTICS)" : "");
        legend.AddEntry(&score, nlabel.c_str(), "l"); legend.Draw();
        canvas.Modified(); canvas.Update(); detail::writeOrThrow(canvas);
        if (!pdfPrefix.empty()) canvas.SaveAs((pdfPrefix + "_score.pdf").c_str());
        canvas.Clear(); // remove pointers to stack graphics objects before destruction
    }
    {
        TCanvas canvas(("canvas_distribution" + suffix).c_str(), "ToA before and after", 1000, 650);
        canvas.SetLeftMargin(0.12); canvas.SetBottomMargin(0.12);
        before->SetTitle(("ToA correction, channel " + std::to_string(globalChannel)
                        + ";ToA [ns];Events / bin").c_str());
        before->SetLineColor(kBlue + 1); after->SetLineColor(kRed + 1);
        smoothBefore->SetLineColor(kBlue + 1); smoothAfter->SetLineColor(kRed + 1);
        smoothBefore->SetLineStyle(2); smoothAfter->SetLineStyle(2);
        smoothBefore->SetLineWidth(2); smoothAfter->SetLineWidth(2);
        before->SetMinimum(0.0);
        before->SetMaximum(1.30 * std::max(before->GetMaximum(), after->GetMaximum()));
        before->Draw("HIST"); after->Draw("HIST SAME");
        smoothBefore->Draw("HIST SAME"); smoothAfter->Draw("HIST SAME");
        TLegend legend(0.53, 0.69, 0.89, 0.89);
        legend.SetFillStyle(0); legend.SetBorderSize(0);
        legend.AddEntry(before.get(), "Before: no 25 ns subtraction", "l");
        const std::string label = "After: T=" + std::to_string(r.threshold);
        legend.AddEntry(after.get(), label.c_str(), "l");
        legend.AddEntry(smoothBefore.get(), "Before, smoothed for scoring", "l");
        legend.AddEntry(smoothAfter.get(), "After, smoothed for scoring", "l");
        legend.Draw(); canvas.Modified(); canvas.Update(); detail::writeOrThrow(canvas);
        if (!pdfPrefix.empty()) canvas.SaveAs((pdfPrefix + "_distribution.pdf").c_str());
        canvas.Clear();
    }
    {
        TCanvas canvas(("canvas_scan_map" + suffix).c_str(), "All threshold distributions", 1000, 700);
        canvas.SetLeftMargin(0.12); canvas.SetBottomMargin(0.12); canvas.SetRightMargin(0.16);
        canvas.SetLogz(); map->SetMinimum(0.5); map->Draw("COLZ");
        TLine line(r.threshold, r.xminNs, r.threshold, r.xmaxNs);
        line.SetLineColor(kRed + 1); line.SetLineStyle(2); line.SetLineWidth(2); line.Draw();
        canvas.Modified(); canvas.Update(); detail::writeOrThrow(canvas);
        if (!pdfPrefix.empty()) canvas.SaveAs((pdfPrefix + "_map.pdf").c_str());
        canvas.Clear();
    }
    return r;
}
} // namespace h2g_toa
#endif // H2GCROC_TOA_SCAN_NO_ROOT
#endif // H2GCROC_TOA_THRESHOLD_SCAN_HXX
