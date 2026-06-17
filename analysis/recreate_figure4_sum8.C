// recreate_figure4_sum8.C
//
// Recreates the timing-resolution-vs-position analysis inspired by Fig. 4 of
// Betancourt et al., JINST 12 (2017) P11023, using the SHiP optical-simulation
// TTree "sipm_hits".
//
// Main definitions used here:
//   END left  = global_id 0..7   (analog sum of 8 SiPMs)
//   END right = global_id 8..15  (analog sum of 8 SiPMs)
//   TOP       = global_id 16..85
//   TOP8      = analog sum of the N nearest top SiPMs (N=8 by default)
//   TOP-all   = analog sum of every top SiPM, stored as a diagnostic
//
// The timing pickoff is a leading-edge crossing of a normalized single-PE
// response (difference of exponentials), with an absolute threshold expressed
// in PE-equivalent amplitude. Defaults reproduce the provisional electronics
// model already used in the project: tau_r=0.5 ns, tau_f=5 ns, threshold=4 PE.
//
// Outputs:
//   * TGraphErrors, TH1D distributions, TF1 fits, TTrees and TCanvas objects
//     in one ROOT file.
//   * PNG/PDF copies of the canvases in <output_stem>_plots/.
//   * CSV summary next to the ROOT file.
//
// Run with defaults:
//   root -l -b -q 'recreate_figure4_sum8.C+()'
//
// Interactive loading (the macro still writes all canvases to the ROOT file):
//   root -l
//   .L recreate_figure4_sum8.C+
//   recreate_figure4_sum8();
//
// The four default data directories are taken from the orchestrator/project:
//   EJ-204 END-only : /home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614
//   EJ-230 END-only : /home/reriosto/SHiP/t0minidaq/endonly_mylar_230
//   EJ-204 END+TOP  : /home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000
//   EJ-230 END+TOP  : /home/reriosto/SHiP/t0minidaq/results_ej230/data

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TMultiGraph.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLatex.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace fig4sum8 {

constexpr double kSprRiseNs = 0.5;
constexpr double kSprFallNs = 5.0;
constexpr double kBarHalfLengthMm = 700.0;
constexpr double kReferenceVeffCmNs = 15.5;

struct DatasetConfig {
    std::string key;
    std::string material;
    std::string topology;
    std::string inputDir;
    bool expectTop = false;
};

struct EventHits {
    std::vector<double> endLeft;
    std::vector<double> endRight;
    std::vector<double> topNearest;
    std::vector<double> topAll;
};

struct FitResult {
    double meanNs = std::numeric_limits<double>::quiet_NaN();
    double meanErrNs = std::numeric_limits<double>::quiet_NaN();
    double sigmaPs = std::numeric_limits<double>::quiet_NaN();
    double sigmaErrPs = std::numeric_limits<double>::quiet_NaN();
    double rmsPs = std::numeric_limits<double>::quiet_NaN();
    double chi2Ndf = std::numeric_limits<double>::quiet_NaN();
    int n = 0;
    bool usedFit = false;
};

struct SpeedResult {
    double slopeNsPerCm = std::numeric_limits<double>::quiet_NaN();
    double slopeErrNsPerCm = std::numeric_limits<double>::quiet_NaN();
    double velocityCmNs = std::numeric_limits<double>::quiet_NaN();
    double velocityErrCmNs = std::numeric_limits<double>::quiet_NaN();
    double chi2Ndf = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
};

struct PositionResult {
    double xMm = std::numeric_limits<double>::quiet_NaN();
    double distanceFromLeftCm = std::numeric_limits<double>::quiet_NaN();
    int nEvents = 0;
    int nEndLeft = 0;
    int nEndRight = 0;
    int nEndBoth = 0;
    int nTopNearest = 0;
    int nTopAll = 0;
    std::string topIds;

    FitResult endLeft;
    FitResult endRight;
    FitResult endMean;
    FitResult endWeightedMean;
    FitResult endDelta;
    FitResult topNearest;
    FitResult topAll;
};

struct DatasetResult {
    DatasetConfig cfg;
    std::vector<PositionResult> points;
    SpeedResult speedLeft;
    SpeedResult speedRight;
    SpeedResult speedDelta;
    double speedConsensusCmNs = std::numeric_limits<double>::quiet_NaN();
    double speedConsensusErrCmNs = std::numeric_limits<double>::quiet_NaN();
};

bool IsFinite(double x) { return std::isfinite(x); }

std::string Sanitize(const std::string& input) {
    std::string out;
    out.reserve(input.size());
    for (char c : input) {
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') || c == '_') {
            out.push_back(c);
        } else {
            out.push_back('_');
        }
    }
    return out;
}

bool EndsWith(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string JoinPath(const std::string& a, const std::string& b) {
    if (a.empty()) return b;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

std::string StemWithoutRoot(const std::string& path) {
    return EndsWith(path, ".root") ? path.substr(0, path.size() - 5) : path;
}

std::vector<std::string> FindRootFiles(const std::string& directory) {
    std::vector<std::string> files;
    if (gSystem->AccessPathName(directory.c_str())) return files;

    TSystemDirectory dir(directory.c_str(), directory.c_str());
    std::unique_ptr<TList> list(dir.GetListOfFiles());
    if (!list) return files;

    TIter next(list.get());
    while (auto* obj = next()) {
        auto* file = dynamic_cast<TSystemFile*>(obj);
        if (!file || file->IsDirectory()) continue;
        const std::string name = file->GetName();
        if (!EndsWith(name, ".root")) continue;
        if (name.find("photon_hits") == std::string::npos) continue;
        files.push_back(JoinPath(directory, name));
    }
    std::sort(files.begin(), files.end());
    return files;
}

double Mean(const std::vector<double>& values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    return std::accumulate(values.begin(), values.end(), 0.0) /
           static_cast<double>(values.size());
}

double StdDev(const std::vector<double>& values) {
    if (values.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    const double mean = Mean(values);
    double sum = 0.0;
    for (double v : values) sum += (v - mean) * (v - mean);
    return std::sqrt(sum / static_cast<double>(values.size() - 1));
}

double Median(std::vector<double> values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    const std::size_t mid = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    double med = values[mid];
    if (values.size() % 2 == 0) {
        const auto lowerMax = std::max_element(values.begin(), values.begin() + mid);
        med = 0.5 * (med + *lowerMax);
    }
    return med;
}

double RobustSigma(const std::vector<double>& values) {
    if (values.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    const double med = Median(values);
    std::vector<double> deviations;
    deviations.reserve(values.size());
    for (double v : values) deviations.push_back(std::abs(v - med));
    return 1.4826 * Median(std::move(deviations));
}

double SprPeakTime() {
    return kSprRiseNs * kSprFallNs / (kSprFallNs - kSprRiseNs) *
           std::log(kSprFallNs / kSprRiseNs);
}

double SprNorm() {
    const double peak = SprPeakTime();
    return 1.0 /
           (std::exp(-peak / kSprFallNs) - std::exp(-peak / kSprRiseNs));
}

double Pulse(double slowState, double fastState, double dt) {
    return SprNorm() *
           (slowState * std::exp(-dt / kSprFallNs) -
            fastState * std::exp(-dt / kSprRiseNs));
}

// Analog sum of all selected SiPM hits. Each detected PE contributes one
// normalized single-PE pulse. The returned time is the first threshold crossing.
double LeadingEdgeTime(std::vector<double> arrivals, double thresholdPe) {
    if (arrivals.empty() || thresholdPe <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::sort(arrivals.begin(), arrivals.end());

    double slowState = 0.0;
    double fastState = 0.0;
    std::size_t index = 0;

    while (index < arrivals.size()) {
        const double current = arrivals[index];
        std::size_t nextIndex = index;
        while (nextIndex < arrivals.size() && arrivals[nextIndex] == current) {
            slowState += 1.0;
            fastState += 1.0;
            ++nextIndex;
        }

        const double interval = nextIndex < arrivals.size()
                                    ? arrivals[nextIndex] - current
                                    : std::numeric_limits<double>::infinity();
        const double derivative0 =
            fastState / kSprRiseNs - slowState / kSprFallNs;

        if (derivative0 > 0.0) {
            const double peakDt =
                std::log((fastState * kSprFallNs) /
                         (slowState * kSprRiseNs)) /
                (1.0 / kSprRiseNs - 1.0 / kSprFallNs);
            const double reach = std::min(peakDt, interval);
            if (reach >= 0.0 &&
                Pulse(slowState, fastState, reach) >= thresholdPe) {
                double low = 0.0;
                double high = reach;
                for (int iteration = 0; iteration < 60; ++iteration) {
                    const double middle = 0.5 * (low + high);
                    if (Pulse(slowState, fastState, middle) >= thresholdPe)
                        high = middle;
                    else
                        low = middle;
                }
                return current + high;
            }
        }

        if (nextIndex >= arrivals.size()) break;
        slowState *= std::exp(-interval / kSprFallNs);
        fastState *= std::exp(-interval / kSprRiseNs);
        index = nextIndex;
    }

    return std::numeric_limits<double>::quiet_NaN();
}

double TopPositionMm(int globalId) {
    const int index = globalId - 16;
    if (index < 0 || index >= 70) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (index < 35) return -692.0 + 20.0 * index;
    return 12.0 + 20.0 * (index - 35);
}

std::vector<int> NearestTopIds(double xMm, int count) {
    count = std::max(1, std::min(count, 70));
    std::vector<std::pair<double, int>> ranked;
    ranked.reserve(70);
    for (int id = 16; id <= 85; ++id) {
        ranked.emplace_back(std::abs(TopPositionMm(id) - xMm), id);
    }
    std::sort(ranked.begin(), ranked.end(),
              [](const auto& a, const auto& b) {
                  if (a.first != b.first) return a.first < b.first;
                  return a.second < b.second;
              });
    std::vector<int> ids;
    ids.reserve(count);
    for (int i = 0; i < count; ++i) ids.push_back(ranked[i].second);
    std::sort(ids.begin(), ids.end());
    return ids;
}

std::string IdListString(const std::vector<int>& ids) {
    std::ostringstream out;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (i) out << ";";
        out << ids[i];
    }
    return out.str();
}

FitResult FitDistribution(const std::vector<double>& values,
                          const std::string& name,
                          const std::string& title,
                          TDirectory* histDir) {
    FitResult result;
    result.n = static_cast<int>(values.size());
    result.rmsPs = values.size() > 1 ? 1000.0 * StdDev(values)
                                     : std::numeric_limits<double>::quiet_NaN();
    if (values.empty()) return result;

    const double center0 = Median(values);
    double robustSigma = RobustSigma(values);
    if (!IsFinite(robustSigma) || robustSigma < 1e-5) robustSigma = 1e-4;
    const double lo = center0 - 8.0 * robustSigma;
    const double hi = center0 + 8.0 * robustSigma;

    auto hist = std::make_unique<TH1D>(name.c_str(), title.c_str(), 120, lo, hi);
    hist->SetDirectory(nullptr);
    for (double v : values) hist->Fill(v);

    if (values.size() >= 20) {
        double mean = center0;
        double sigma = robustSigma;
        TF1 gaus((name + "_gaus").c_str(), "gaus", lo, hi);
        int status = 1;
        for (int iteration = 0; iteration < 4; ++iteration) {
            const double fitLo = std::max(lo, mean - 2.0 * sigma);
            const double fitHi = std::min(hi, mean + 2.0 * sigma);
            gaus.SetRange(fitLo, fitHi);
            gaus.SetParameters(hist->GetMaximum(), mean, sigma);
            gaus.SetParLimits(2, 1e-7, hi - lo);
            TFitResultPtr fit = hist->Fit(&gaus, "QRSN");
            status = static_cast<int>(fit);
            if (status != 0 || gaus.GetParameter(2) <= 0.0) break;
            mean = gaus.GetParameter(1);
            sigma = std::abs(gaus.GetParameter(2));
        }

        if (status == 0 && sigma > 0.0) {
            result.meanNs = mean;
            result.meanErrNs = gaus.GetParError(1);
            result.sigmaPs = 1000.0 * sigma;
            result.sigmaErrPs = 1000.0 * gaus.GetParError(2);
            result.chi2Ndf = gaus.GetNDF() > 0
                                 ? gaus.GetChisquare() / gaus.GetNDF()
                                 : std::numeric_limits<double>::quiet_NaN();
            result.usedFit = true;
        }
    }

    if (!result.usedFit) {
        result.meanNs = Mean(values);
        const double sigmaNs = StdDev(values);
        result.meanErrNs = IsFinite(sigmaNs)
                               ? sigmaNs / std::sqrt(static_cast<double>(values.size()))
                               : std::numeric_limits<double>::quiet_NaN();
        result.sigmaPs = IsFinite(sigmaNs) ? 1000.0 * sigmaNs
                                           : std::numeric_limits<double>::quiet_NaN();
        result.sigmaErrPs = values.size() > 1 && IsFinite(result.sigmaPs)
                                ? result.sigmaPs /
                                      std::sqrt(2.0 * (values.size() - 1.0))
                                : std::numeric_limits<double>::quiet_NaN();
        result.chi2Ndf = -1.0;
    }

    if (histDir) {
        histDir->cd();
        hist->Write();
    }
    return result;
}

void SetGraphStyle(TGraphErrors* graph, int color, int marker) {
    if (!graph) return;
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetLineWidth(2);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(0.9);
}

TGraphErrors* MakeGraph(const std::vector<PositionResult>& points,
                        const std::string& name,
                        const std::string& title,
                        const std::function<bool(const PositionResult&)>& accept,
                        const std::function<double(const PositionResult&)>& xValue,
                        const std::function<double(const PositionResult&)>& yValue,
                        const std::function<double(const PositionResult&)>& yError) {
    auto* graph = new TGraphErrors();
    graph->SetName(name.c_str());
    graph->SetTitle(title.c_str());
    int index = 0;
    for (const auto& p : points) {
        if (!accept(p)) continue;
        const double x = xValue(p);
        const double y = yValue(p);
        const double ey = yError(p);
        if (!IsFinite(x) || !IsFinite(y)) continue;
        graph->SetPoint(index, x, y);
        graph->SetPointError(index, 0.0, IsFinite(ey) ? ey : 0.0);
        ++index;
    }
    return graph;
}

SpeedResult FitSpeed(TGraphErrors* graph,
                     const std::string& fitName,
                     double numerator,
                     TDirectory* outputDir) {
    SpeedResult result;
    if (!graph || graph->GetN() < 2) return result;

    double xmin = 0.0, xmax = 0.0;
    double x = 0.0, y = 0.0;
    graph->GetPoint(0, xmin, y);
    xmax = xmin;
    for (int i = 1; i < graph->GetN(); ++i) {
        graph->GetPoint(i, x, y);
        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
    }
    if (!(xmax > xmin)) return result;

    // The fit is deliberately attached to the graph so that it is also visible
    // when the stored TCanvas is reopened from the output ROOT file.
    auto* fit = new TF1(fitName.c_str(), "pol1", xmin, xmax);
    TFitResultPtr fitPtr = graph->Fit(fit, "QRS");
    if (static_cast<int>(fitPtr) != 0) {
        // ROOT may already have attached the TF1 to the graph; keep it alive.
        return result;
    }

    const double slope = fit->GetParameter(1);
    const double slopeErr = fit->GetParError(1);
    if (std::abs(slope) < 1e-12) return result;

    result.slopeNsPerCm = slope;
    result.slopeErrNsPerCm = slopeErr;
    result.velocityCmNs = numerator / std::abs(slope);
    result.velocityErrCmNs = numerator * slopeErr / (slope * slope);
    result.chi2Ndf = fit->GetNDF() > 0
                         ? fit->GetChisquare() / fit->GetNDF()
                         : std::numeric_limits<double>::quiet_NaN();
    result.valid = true;

    fit->SetLineWidth(2);
    if (outputDir) {
        outputDir->cd();
        fit->Write();
    }
    return result;
}

std::pair<double, double> WeightedVelocity(const std::vector<SpeedResult>& speeds) {
    double sumW = 0.0;
    double sumWV = 0.0;
    std::vector<double> validValues;
    for (const auto& s : speeds) {
        if (!s.valid || !IsFinite(s.velocityCmNs)) continue;
        validValues.push_back(s.velocityCmNs);
        if (IsFinite(s.velocityErrCmNs) && s.velocityErrCmNs > 0.0) {
            const double w = 1.0 / (s.velocityErrCmNs * s.velocityErrCmNs);
            sumW += w;
            sumWV += w * s.velocityCmNs;
        }
    }
    if (sumW > 0.0) return {sumWV / sumW, std::sqrt(1.0 / sumW)};
    if (validValues.empty()) {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }
    const double mean = Mean(validValues);
    const double err = validValues.size() > 1
                           ? StdDev(validValues) /
                                 std::sqrt(static_cast<double>(validValues.size()))
                           : std::numeric_limits<double>::quiet_NaN();
    return {mean, err};
}

PositionResult ProcessFile(const std::string& fileName,
                           const DatasetConfig& dataset,
                           double thresholdPe,
                           int topGroupSize,
                           int expectedEvents,
                           TDirectory* histDir) {
    PositionResult result;

    TFile input(fileName.c_str(), "READ");
    if (input.IsZombie()) {
        std::cerr << "[WARN] Cannot open " << fileName << "\n";
        return result;
    }
    auto* tree = dynamic_cast<TTree*>(input.Get("sipm_hits"));
    if (!tree) {
        std::cerr << "[WARN] TTree sipm_hits not found in " << fileName << "\n";
        return result;
    }

    TTreeReader reader(tree);
    TTreeReaderValue<int> eventId(reader, "event_id");
    TTreeReaderValue<int> globalId(reader, "global_id");
    TTreeReaderValue<double> timeNs(reader, "time_ns");
    TTreeReaderValue<double> gunX(reader, "gun_x_mm");

    std::vector<EventHits> events(std::max(0, expectedEvents));
    std::vector<int> nearestTopIds;
    std::set<int> nearestTopSet;
    bool initialized = false;
    int maxEventId = -1;

    while (reader.Next()) {
        if (!initialized) {
            result.xMm = *gunX;
            result.distanceFromLeftCm =
                (result.xMm + kBarHalfLengthMm) / 10.0;
            nearestTopIds = NearestTopIds(result.xMm, topGroupSize);
            nearestTopSet.insert(nearestTopIds.begin(), nearestTopIds.end());
            result.topIds = IdListString(nearestTopIds);
            initialized = true;
        }

        if (*eventId < 0) continue;
        maxEventId = std::max(maxEventId, *eventId);
        if (static_cast<std::size_t>(*eventId) >= events.size()) {
            events.resize(static_cast<std::size_t>(*eventId) + 1);
        }

        auto& event = events[*eventId];
        if (*globalId >= 0 && *globalId <= 7) {
            event.endLeft.push_back(*timeNs);
        } else if (*globalId >= 8 && *globalId <= 15) {
            event.endRight.push_back(*timeNs);
        } else if (*globalId >= 16 && *globalId <= 85) {
            event.topAll.push_back(*timeNs);
            if (nearestTopSet.count(*globalId)) {
                event.topNearest.push_back(*timeNs);
            }
        }
    }

    if (!initialized) {
        std::cerr << "[WARN] Empty tree in " << fileName << "\n";
        return result;
    }

    result.nEvents = std::max(expectedEvents, maxEventId + 1);
    if (static_cast<int>(events.size()) < result.nEvents) events.resize(result.nEvents);

    std::vector<double> tLeft(result.nEvents,
                              std::numeric_limits<double>::quiet_NaN());
    std::vector<double> tRight(result.nEvents,
                               std::numeric_limits<double>::quiet_NaN());
    std::vector<double> tTopNearest(result.nEvents,
                                    std::numeric_limits<double>::quiet_NaN());
    std::vector<double> tTopAll(result.nEvents,
                                std::numeric_limits<double>::quiet_NaN());

    std::vector<double> leftValues, rightValues, topNearestValues, topAllValues;
    leftValues.reserve(result.nEvents);
    rightValues.reserve(result.nEvents);
    topNearestValues.reserve(result.nEvents);
    topAllValues.reserve(result.nEvents);

    for (int ev = 0; ev < result.nEvents; ++ev) {
        tLeft[ev] = LeadingEdgeTime(std::move(events[ev].endLeft), thresholdPe);
        tRight[ev] = LeadingEdgeTime(std::move(events[ev].endRight), thresholdPe);
        tTopNearest[ev] =
            LeadingEdgeTime(std::move(events[ev].topNearest), thresholdPe);
        tTopAll[ev] = LeadingEdgeTime(std::move(events[ev].topAll), thresholdPe);

        if (IsFinite(tLeft[ev])) leftValues.push_back(tLeft[ev]);
        if (IsFinite(tRight[ev])) rightValues.push_back(tRight[ev]);
        if (IsFinite(tTopNearest[ev])) topNearestValues.push_back(tTopNearest[ev]);
        if (IsFinite(tTopAll[ev])) topAllValues.push_back(tTopAll[ev]);
    }

    result.nEndLeft = static_cast<int>(leftValues.size());
    result.nEndRight = static_cast<int>(rightValues.size());
    result.nTopNearest = static_cast<int>(topNearestValues.size());
    result.nTopAll = static_cast<int>(topAllValues.size());

    const std::string xTag = "x" + std::to_string(static_cast<int>(std::llround(result.xMm)));
    result.endLeft = FitDistribution(
        leftValues, "h_end_left_" + xTag,
        dataset.material + " " + dataset.topology + " END-left SUM8;time [ns];events/bin",
        histDir);
    result.endRight = FitDistribution(
        rightValues, "h_end_right_" + xTag,
        dataset.material + " " + dataset.topology + " END-right SUM8;time [ns];events/bin",
        histDir);
    result.topNearest = FitDistribution(
        topNearestValues, "h_top_nearest_" + xTag,
        dataset.material + " " + dataset.topology + " TOP nearest SiPM sum;time [ns];events/bin",
        histDir);
    result.topAll = FitDistribution(
        topAllValues, "h_top_all_" + xTag,
        dataset.material + " " + dataset.topology + " TOP-all;time [ns];events/bin",
        histDir);

    std::vector<double> meanValues, weightedValues, deltaValues;
    meanValues.reserve(result.nEvents);
    weightedValues.reserve(result.nEvents);
    deltaValues.reserve(result.nEvents);

    const double sigmaLeftNs = result.endLeft.sigmaPs / 1000.0;
    const double sigmaRightNs = result.endRight.sigmaPs / 1000.0;
    const double wLeft = IsFinite(sigmaLeftNs) && sigmaLeftNs > 0.0
                             ? 1.0 / (sigmaLeftNs * sigmaLeftNs)
                             : 0.0;
    const double wRight = IsFinite(sigmaRightNs) && sigmaRightNs > 0.0
                              ? 1.0 / (sigmaRightNs * sigmaRightNs)
                              : 0.0;

    for (int ev = 0; ev < result.nEvents; ++ev) {
        if (!IsFinite(tLeft[ev]) || !IsFinite(tRight[ev])) continue;
        ++result.nEndBoth;
        meanValues.push_back(0.5 * (tLeft[ev] + tRight[ev]));
        deltaValues.push_back(tLeft[ev] - tRight[ev]);
        if (wLeft + wRight > 0.0) {
            weightedValues.push_back(
                (wLeft * tLeft[ev] + wRight * tRight[ev]) / (wLeft + wRight));
        }
    }

    result.endMean = FitDistribution(
        meanValues, "h_end_mean_" + xTag,
        dataset.material + " " + dataset.topology + " END mean time;mean time [ns];events/bin",
        histDir);
    result.endWeightedMean = FitDistribution(
        weightedValues, "h_end_weighted_" + xTag,
        dataset.material + " " + dataset.topology + " END weighted mean;weighted time [ns];events/bin",
        histDir);
    result.endDelta = FitDistribution(
        deltaValues, "h_end_delta_" + xTag,
        dataset.material + " " + dataset.topology + " END #Delta t; t_{L}-t_{R} [ns];events/bin",
        histDir);

    std::cout << std::fixed << std::setprecision(2)
              << "[" << dataset.key << "] x=" << result.xMm << " mm"
              << "  N=" << result.nEvents
              << "  eff(L/R/TOP8)="
              << (result.nEvents > 0 ? 100.0 * result.nEndLeft / result.nEvents : 0.0)
              << "/"
              << (result.nEvents > 0 ? 100.0 * result.nEndRight / result.nEvents : 0.0)
              << "/"
              << (result.nEvents > 0 ? 100.0 * result.nTopNearest / result.nEvents : 0.0)
              << "%  sigma(L/R/mean/TOP8)="
              << result.endLeft.sigmaPs << "/" << result.endRight.sigmaPs << "/"
              << result.endMean.sigmaPs << "/" << result.topNearest.sigmaPs
              << " ps\n";

    return result;
}

void WritePointTree(const DatasetResult& dataset, TDirectory* directory) {
    if (!directory) return;
    directory->cd();
    auto tree = std::make_unique<TTree>("summary", "Timing summary by interaction position");

    double x_mm = 0.0, distance_left_cm = 0.0;
    int n_events = 0, n_end_left = 0, n_end_right = 0, n_end_both = 0;
    int n_top_nearest = 0, n_top_all = 0;
    std::string top_ids;

    double mean_left_ns = 0.0, mean_left_err_ns = 0.0, sigma_left_ps = 0.0,
           sigma_left_err_ps = 0.0;
    double mean_right_ns = 0.0, mean_right_err_ns = 0.0, sigma_right_ps = 0.0,
           sigma_right_err_ps = 0.0;
    double mean_mean_ns = 0.0, sigma_mean_ps = 0.0, sigma_mean_err_ps = 0.0;
    double mean_weighted_ns = 0.0, sigma_weighted_ps = 0.0,
           sigma_weighted_err_ps = 0.0;
    double mean_delta_ns = 0.0, mean_delta_err_ns = 0.0, sigma_delta_ps = 0.0;
    double mean_top8_ns = 0.0, sigma_top8_ps = 0.0, sigma_top8_err_ps = 0.0;
    double mean_topall_ns = 0.0, sigma_topall_ps = 0.0,
           sigma_topall_err_ps = 0.0;

    tree->Branch("x_mm", &x_mm);
    tree->Branch("distance_from_left_cm", &distance_left_cm);
    tree->Branch("n_events", &n_events);
    tree->Branch("n_end_left", &n_end_left);
    tree->Branch("n_end_right", &n_end_right);
    tree->Branch("n_end_both", &n_end_both);
    tree->Branch("n_top_nearest", &n_top_nearest);
    tree->Branch("n_top_all", &n_top_all);
    tree->Branch("top_ids", &top_ids);

    tree->Branch("mean_left_ns", &mean_left_ns);
    tree->Branch("mean_left_err_ns", &mean_left_err_ns);
    tree->Branch("sigma_left_ps", &sigma_left_ps);
    tree->Branch("sigma_left_err_ps", &sigma_left_err_ps);
    tree->Branch("mean_right_ns", &mean_right_ns);
    tree->Branch("mean_right_err_ns", &mean_right_err_ns);
    tree->Branch("sigma_right_ps", &sigma_right_ps);
    tree->Branch("sigma_right_err_ps", &sigma_right_err_ps);
    tree->Branch("mean_mean_ns", &mean_mean_ns);
    tree->Branch("sigma_mean_ps", &sigma_mean_ps);
    tree->Branch("sigma_mean_err_ps", &sigma_mean_err_ps);
    tree->Branch("mean_weighted_ns", &mean_weighted_ns);
    tree->Branch("sigma_weighted_ps", &sigma_weighted_ps);
    tree->Branch("sigma_weighted_err_ps", &sigma_weighted_err_ps);
    tree->Branch("mean_delta_ns", &mean_delta_ns);
    tree->Branch("mean_delta_err_ns", &mean_delta_err_ns);
    tree->Branch("sigma_delta_ps", &sigma_delta_ps);
    tree->Branch("mean_top8_ns", &mean_top8_ns);
    tree->Branch("sigma_top8_ps", &sigma_top8_ps);
    tree->Branch("sigma_top8_err_ps", &sigma_top8_err_ps);
    tree->Branch("mean_topall_ns", &mean_topall_ns);
    tree->Branch("sigma_topall_ps", &sigma_topall_ps);
    tree->Branch("sigma_topall_err_ps", &sigma_topall_err_ps);

    for (const auto& p : dataset.points) {
        x_mm = p.xMm;
        distance_left_cm = p.distanceFromLeftCm;
        n_events = p.nEvents;
        n_end_left = p.nEndLeft;
        n_end_right = p.nEndRight;
        n_end_both = p.nEndBoth;
        n_top_nearest = p.nTopNearest;
        n_top_all = p.nTopAll;
        top_ids = p.topIds;

        mean_left_ns = p.endLeft.meanNs;
        mean_left_err_ns = p.endLeft.meanErrNs;
        sigma_left_ps = p.endLeft.sigmaPs;
        sigma_left_err_ps = p.endLeft.sigmaErrPs;
        mean_right_ns = p.endRight.meanNs;
        mean_right_err_ns = p.endRight.meanErrNs;
        sigma_right_ps = p.endRight.sigmaPs;
        sigma_right_err_ps = p.endRight.sigmaErrPs;
        mean_mean_ns = p.endMean.meanNs;
        sigma_mean_ps = p.endMean.sigmaPs;
        sigma_mean_err_ps = p.endMean.sigmaErrPs;
        mean_weighted_ns = p.endWeightedMean.meanNs;
        sigma_weighted_ps = p.endWeightedMean.sigmaPs;
        sigma_weighted_err_ps = p.endWeightedMean.sigmaErrPs;
        mean_delta_ns = p.endDelta.meanNs;
        mean_delta_err_ns = p.endDelta.meanErrNs;
        sigma_delta_ps = p.endDelta.sigmaPs;
        mean_top8_ns = p.topNearest.meanNs;
        sigma_top8_ps = p.topNearest.sigmaPs;
        sigma_top8_err_ps = p.topNearest.sigmaErrPs;
        mean_topall_ns = p.topAll.meanNs;
        sigma_topall_ps = p.topAll.sigmaPs;
        sigma_topall_err_ps = p.topAll.sigmaErrPs;
        tree->Fill();
    }
    tree->Write();
}

void SaveCanvas(TCanvas* canvas, TDirectory* directory,
                const std::string& plotDir) {
    if (!canvas || !directory) return;
    directory->cd();
    canvas->Write();
    canvas->SaveAs(JoinPath(plotDir, std::string(canvas->GetName()) + ".png").c_str());
    canvas->SaveAs(JoinPath(plotDir, std::string(canvas->GetName()) + ".pdf").c_str());
}

void CreateDatasetObjects(DatasetResult& dataset,
                          TDirectory* directory,
                          const std::string& plotDir,
                          int topGroupSize) {
    if (!directory) return;
    directory->cd();

    const std::string key = Sanitize(dataset.cfg.key);

    auto* gSigmaLeft = MakeGraph(
        dataset.points, "g_sigma_end_left_" + key,
        dataset.cfg.material + " " + dataset.cfg.topology +
            ";distance from left end [cm];time resolution [ps]",
        [](const PositionResult& p) { return p.endLeft.n >= 20; },
        [](const PositionResult& p) { return p.distanceFromLeftCm; },
        [](const PositionResult& p) { return p.endLeft.sigmaPs; },
        [](const PositionResult& p) { return p.endLeft.sigmaErrPs; });
    auto* gSigmaRight = MakeGraph(
        dataset.points, "g_sigma_end_right_" + key, "",
        [](const PositionResult& p) { return p.endRight.n >= 20; },
        [](const PositionResult& p) { return p.distanceFromLeftCm; },
        [](const PositionResult& p) { return p.endRight.sigmaPs; },
        [](const PositionResult& p) { return p.endRight.sigmaErrPs; });
    auto* gSigmaMean = MakeGraph(
        dataset.points, "g_sigma_end_mean_" + key, "",
        [](const PositionResult& p) { return p.endMean.n >= 20; },
        [](const PositionResult& p) { return p.distanceFromLeftCm; },
        [](const PositionResult& p) { return p.endMean.sigmaPs; },
        [](const PositionResult& p) { return p.endMean.sigmaErrPs; });
    auto* gSigmaWeighted = MakeGraph(
        dataset.points, "g_sigma_end_weighted_" + key, "",
        [](const PositionResult& p) { return p.endWeightedMean.n >= 20; },
        [](const PositionResult& p) { return p.distanceFromLeftCm; },
        [](const PositionResult& p) { return p.endWeightedMean.sigmaPs; },
        [](const PositionResult& p) { return p.endWeightedMean.sigmaErrPs; });
    auto* gSigmaTop8 = MakeGraph(
        dataset.points, "g_sigma_top_nearest_" + key,
        dataset.cfg.material + " " + dataset.cfg.topology +
            ";x_{gun} [cm];time resolution [ps]",
        [](const PositionResult& p) { return p.topNearest.n >= 20; },
        [](const PositionResult& p) { return p.xMm / 10.0; },
        [](const PositionResult& p) { return p.topNearest.sigmaPs; },
        [](const PositionResult& p) { return p.topNearest.sigmaErrPs; });
    auto* gSigmaTopAll = MakeGraph(
        dataset.points, "g_sigma_top_all_" + key, "",
        [](const PositionResult& p) { return p.topAll.n >= 20; },
        [](const PositionResult& p) { return p.xMm / 10.0; },
        [](const PositionResult& p) { return p.topAll.sigmaPs; },
        [](const PositionResult& p) { return p.topAll.sigmaErrPs; });

    auto* gMeanLeft = MakeGraph(
        dataset.points, "g_mean_end_left_" + key,
        dataset.cfg.material + " " + dataset.cfg.topology +
            ";x_{gun} [cm];mean trigger time [ns]",
        [](const PositionResult& p) { return p.endLeft.n >= 20; },
        [](const PositionResult& p) { return p.xMm / 10.0; },
        [](const PositionResult& p) { return p.endLeft.meanNs; },
        [](const PositionResult& p) { return p.endLeft.meanErrNs; });
    auto* gMeanRight = MakeGraph(
        dataset.points, "g_mean_end_right_" + key, "",
        [](const PositionResult& p) { return p.endRight.n >= 20; },
        [](const PositionResult& p) { return p.xMm / 10.0; },
        [](const PositionResult& p) { return p.endRight.meanNs; },
        [](const PositionResult& p) { return p.endRight.meanErrNs; });
    auto* gMeanDelta = MakeGraph(
        dataset.points, "g_mean_end_delta_" + key,
        dataset.cfg.material + " " + dataset.cfg.topology +
            ";x_{gun} [cm];mean(t_{L}-t_{R}) [ns]",
        [](const PositionResult& p) { return p.endDelta.n >= 20; },
        [](const PositionResult& p) { return p.xMm / 10.0; },
        [](const PositionResult& p) { return p.endDelta.meanNs; },
        [](const PositionResult& p) { return p.endDelta.meanErrNs; });

    auto* gEffLeft = MakeGraph(
        dataset.points, "g_eff_end_left_" + key,
        dataset.cfg.material + " " + dataset.cfg.topology +
            ";x_{gun} [cm];trigger efficiency [%]",
        [](const PositionResult& p) { return p.nEvents > 0; },
        [](const PositionResult& p) { return p.xMm / 10.0; },
        [](const PositionResult& p) {
            return 100.0 * p.nEndLeft / static_cast<double>(p.nEvents);
        },
        [](const PositionResult&) { return 0.0; });
    auto* gEffRight = MakeGraph(
        dataset.points, "g_eff_end_right_" + key, "",
        [](const PositionResult& p) { return p.nEvents > 0; },
        [](const PositionResult& p) { return p.xMm / 10.0; },
        [](const PositionResult& p) {
            return 100.0 * p.nEndRight / static_cast<double>(p.nEvents);
        },
        [](const PositionResult&) { return 0.0; });
    auto* gEffTop8 = MakeGraph(
        dataset.points, "g_eff_top_nearest_" + key, "",
        [](const PositionResult& p) { return p.nEvents > 0; },
        [](const PositionResult& p) { return p.xMm / 10.0; },
        [](const PositionResult& p) {
            return 100.0 * p.nTopNearest / static_cast<double>(p.nEvents);
        },
        [](const PositionResult&) { return 0.0; });

    SetGraphStyle(gSigmaLeft, kRed + 1, 20);
    SetGraphStyle(gSigmaRight, kBlue + 1, 21);
    SetGraphStyle(gSigmaMean, kGreen + 2, 22);
    SetGraphStyle(gSigmaWeighted, kMagenta + 1, 23);
    SetGraphStyle(gSigmaTop8, kOrange + 7, 20);
    SetGraphStyle(gSigmaTopAll, kCyan + 2, 24);
    SetGraphStyle(gMeanLeft, kRed + 1, 20);
    SetGraphStyle(gMeanRight, kBlue + 1, 21);
    SetGraphStyle(gMeanDelta, kBlack, 22);
    SetGraphStyle(gEffLeft, kRed + 1, 20);
    SetGraphStyle(gEffRight, kBlue + 1, 21);
    SetGraphStyle(gEffTop8, kOrange + 7, 22);

    for (auto* graph : {gSigmaLeft, gSigmaRight, gSigmaMean, gSigmaWeighted,
                        gSigmaTop8, gSigmaTopAll, gMeanLeft, gMeanRight,
                        gMeanDelta, gEffLeft, gEffRight, gEffTop8}) {
        directory->cd();
        graph->Write();
    }

    dataset.speedLeft = FitSpeed(gMeanLeft, "fit_speed_left_" + key, 1.0, directory);
    dataset.speedRight = FitSpeed(gMeanRight, "fit_speed_right_" + key, 1.0, directory);
    dataset.speedDelta = FitSpeed(gMeanDelta, "fit_speed_delta_" + key, 2.0, directory);

    // Rewrite the mean-time graphs after fitting so their attached TF1 objects
    // are present in the standalone graph keys, not only inside the canvases.
    directory->cd();
    gMeanLeft->Write(gMeanLeft->GetName(), TObject::kOverwrite);
    gMeanRight->Write(gMeanRight->GetName(), TObject::kOverwrite);
    gMeanDelta->Write(gMeanDelta->GetName(), TObject::kOverwrite);

    const auto consensus = WeightedVelocity(
        {dataset.speedLeft, dataset.speedRight, dataset.speedDelta});
    dataset.speedConsensusCmNs = consensus.first;
    dataset.speedConsensusErrCmNs = consensus.second;

    auto* cFig4 = new TCanvas(("c_fig4_" + key).c_str(),
                              ("Figure 4 recreation " + dataset.cfg.key).c_str(),
                              1100, 720);
    cFig4->SetGrid();
    auto* mgFig4 = new TMultiGraph();
    mgFig4->SetName(("mg_fig4_" + key).c_str());
    mgFig4->SetTitle((dataset.cfg.material + " " + dataset.cfg.topology +
                      " — END analog SUM8;distance from left end [cm];time resolution [ps]")
                         .c_str());
    if (gSigmaLeft->GetN()) mgFig4->Add(gSigmaLeft, "LP");
    if (gSigmaRight->GetN()) mgFig4->Add(gSigmaRight, "LP");
    if (gSigmaMean->GetN()) mgFig4->Add(gSigmaMean, "LP");
    if (gSigmaWeighted->GetN()) mgFig4->Add(gSigmaWeighted, "LP");
    mgFig4->Draw("A");
    mgFig4->SetMinimum(0.0);
    auto* legFig4 = new TLegend(0.57, 0.67, 0.89, 0.89);
    legFig4->SetBorderSize(0);
    legFig4->AddEntry(gSigmaLeft, "8xSiPM END-left", "lp");
    legFig4->AddEntry(gSigmaRight, "8xSiPM END-right", "lp");
    legFig4->AddEntry(gSigmaMean, "mean (t_{L}+t_{R})/2", "lp");
    legFig4->AddEntry(gSigmaWeighted, "inverse-variance weighted mean", "lp");
    legFig4->Draw();
    SaveCanvas(cFig4, directory, plotDir);

    if (gSigmaTop8->GetN() || gSigmaTopAll->GetN()) {
        auto* cTop = new TCanvas(("c_top_" + key).c_str(),
                                 ("TOP timing " + dataset.cfg.key).c_str(), 1100, 720);
        cTop->SetGrid();
        auto* mgTop = new TMultiGraph();
        mgTop->SetName(("mg_top_" + key).c_str());
        mgTop->SetTitle((dataset.cfg.material + " " + dataset.cfg.topology +
                         " — TOP timing;x_{gun} [cm];time resolution [ps]")
                            .c_str());
        if (gSigmaTop8->GetN()) mgTop->Add(gSigmaTop8, "LP");
        if (gSigmaTopAll->GetN()) mgTop->Add(gSigmaTopAll, "LP");
        mgTop->Draw("A");
        mgTop->SetMinimum(0.0);
        auto* legTop = new TLegend(0.55, 0.73, 0.89, 0.89);
        legTop->SetBorderSize(0);
        std::ostringstream label;
        label << "nearest " << topGroupSize << " TOP SiPMs";
        legTop->AddEntry(gSigmaTop8, label.str().c_str(), "lp");
        legTop->AddEntry(gSigmaTopAll, "all 70 TOP SiPMs (diagnostic)", "lp");
        legTop->Draw();
        SaveCanvas(cTop, directory, plotDir);
    }

    auto* cSpeed = new TCanvas(("c_speed_" + key).c_str(),
                               ("Effective propagation speed " + dataset.cfg.key).c_str(),
                               1100, 900);
    cSpeed->Divide(1, 2);
    cSpeed->cd(1);
    gPad->SetGrid();
    auto* mgSpeed = new TMultiGraph();
    mgSpeed->SetName(("mg_speed_" + key).c_str());
    mgSpeed->SetTitle((dataset.cfg.material + " " + dataset.cfg.topology +
                       ";x_{gun} [cm];mean trigger time [ns]")
                          .c_str());
    if (gMeanLeft->GetN()) mgSpeed->Add(gMeanLeft, "P");
    if (gMeanRight->GetN()) mgSpeed->Add(gMeanRight, "P");
    mgSpeed->Draw("A");
    auto* legSpeed = new TLegend(0.56, 0.68, 0.89, 0.89);
    legSpeed->SetBorderSize(0);
    legSpeed->AddEntry(gMeanLeft, "END-left SUM8", "lp");
    legSpeed->AddEntry(gMeanRight, "END-right SUM8", "lp");
    legSpeed->Draw();

    TLatex text;
    text.SetNDC(true);
    text.SetTextSize(0.035);
    double yText = 0.62;
    if (dataset.speedLeft.valid) {
        text.DrawLatex(0.56, yText,
                       Form("v_{L}=%.2f #pm %.2f cm/ns",
                            dataset.speedLeft.velocityCmNs,
                            dataset.speedLeft.velocityErrCmNs));
        yText -= 0.05;
    }
    if (dataset.speedRight.valid) {
        text.DrawLatex(0.56, yText,
                       Form("v_{R}=%.2f #pm %.2f cm/ns",
                            dataset.speedRight.velocityCmNs,
                            dataset.speedRight.velocityErrCmNs));
        yText -= 0.05;
    }

    cSpeed->cd(2);
    gPad->SetGrid();
    gMeanDelta->Draw("AP");
    if (dataset.speedDelta.valid) {
        text.DrawLatex(0.15, 0.84,
                       Form("v_{#Delta t}=2/|d#LT#Delta t#GT/dx|=%.2f #pm %.2f cm/ns",
                            dataset.speedDelta.velocityCmNs,
                            dataset.speedDelta.velocityErrCmNs));
    }
    if (IsFinite(dataset.speedConsensusCmNs)) {
        text.DrawLatex(0.15, 0.77,
                       Form("diagnostic combined v_{eff}=%.2f #pm %.2f cm/ns",
                            dataset.speedConsensusCmNs,
                            dataset.speedConsensusErrCmNs));
        text.DrawLatex(0.15, 0.70,
                       Form("reference: %.1f cm/ns; deviation = %.1f%%",
                            kReferenceVeffCmNs,
                            100.0 * (dataset.speedConsensusCmNs -
                                     kReferenceVeffCmNs) /
                                kReferenceVeffCmNs));
    }
    SaveCanvas(cSpeed, directory, plotDir);

    auto* cEff = new TCanvas(("c_efficiency_" + key).c_str(),
                             ("Trigger efficiency " + dataset.cfg.key).c_str(),
                             1100, 720);
    cEff->SetGrid();
    auto* mgEff = new TMultiGraph();
    mgEff->SetName(("mg_efficiency_" + key).c_str());
    mgEff->SetTitle((dataset.cfg.material + " " + dataset.cfg.topology +
                     ";x_{gun} [cm];trigger efficiency [%]")
                        .c_str());
    if (gEffLeft->GetN()) mgEff->Add(gEffLeft, "LP");
    if (gEffRight->GetN()) mgEff->Add(gEffRight, "LP");
    if (gEffTop8->GetN()) mgEff->Add(gEffTop8, "LP");
    mgEff->Draw("A");
    mgEff->SetMinimum(0.0);
    mgEff->SetMaximum(105.0);
    auto* legEff = new TLegend(0.60, 0.72, 0.89, 0.89);
    legEff->SetBorderSize(0);
    legEff->AddEntry(gEffLeft, "END-left SUM8", "lp");
    legEff->AddEntry(gEffRight, "END-right SUM8", "lp");
    if (gEffTop8->GetN()) legEff->AddEntry(gEffTop8, "TOP nearest group", "lp");
    legEff->Draw();
    SaveCanvas(cEff, directory, plotDir);

    directory->cd();
    TParameter<double>("v_eff_left_cm_ns", dataset.speedLeft.velocityCmNs).Write();
    TParameter<double>("v_eff_left_err_cm_ns", dataset.speedLeft.velocityErrCmNs).Write();
    TParameter<double>("v_eff_right_cm_ns", dataset.speedRight.velocityCmNs).Write();
    TParameter<double>("v_eff_right_err_cm_ns", dataset.speedRight.velocityErrCmNs).Write();
    TParameter<double>("v_eff_delta_cm_ns", dataset.speedDelta.velocityCmNs).Write();
    TParameter<double>("v_eff_delta_err_cm_ns", dataset.speedDelta.velocityErrCmNs).Write();
    TParameter<double>("v_eff_consensus_cm_ns", dataset.speedConsensusCmNs).Write();
    TParameter<double>("v_eff_consensus_err_cm_ns", dataset.speedConsensusErrCmNs).Write();
    TParameter<double>("v_eff_reference_cm_ns", kReferenceVeffCmNs).Write();

    WritePointTree(dataset, directory);
}

void WriteCsv(const std::string& csvName,
              const std::vector<DatasetResult>& datasets) {
    std::ofstream csv(csvName);
    csv << "dataset,material,topology,x_mm,distance_from_left_cm,n_events,"
           "eff_end_left,eff_end_right,eff_end_both,eff_top_nearest,"
           "mean_left_ns,sigma_left_ps,sigma_left_err_ps,"
           "mean_right_ns,sigma_right_ps,sigma_right_err_ps,"
           "mean_end_ns,sigma_end_mean_ps,sigma_end_weighted_ps,"
           "mean_delta_ns,sigma_delta_ps,mean_top_nearest_ns,sigma_top_nearest_ps,"
           "mean_top_all_ns,sigma_top_all_ps,top_ids,"
           "v_left_cm_ns,v_right_cm_ns,v_delta_cm_ns,v_consensus_cm_ns\n";
    csv << std::setprecision(10);
    for (const auto& ds : datasets) {
        for (const auto& p : ds.points) {
            const auto eff = [&](int n) {
                return p.nEvents > 0 ? static_cast<double>(n) / p.nEvents
                                     : std::numeric_limits<double>::quiet_NaN();
            };
            csv << ds.cfg.key << "," << ds.cfg.material << "," << ds.cfg.topology
                << "," << p.xMm << "," << p.distanceFromLeftCm << ","
                << p.nEvents << "," << eff(p.nEndLeft) << ","
                << eff(p.nEndRight) << "," << eff(p.nEndBoth) << ","
                << eff(p.nTopNearest) << "," << p.endLeft.meanNs << ","
                << p.endLeft.sigmaPs << "," << p.endLeft.sigmaErrPs << ","
                << p.endRight.meanNs << "," << p.endRight.sigmaPs << ","
                << p.endRight.sigmaErrPs << "," << p.endMean.meanNs << ","
                << p.endMean.sigmaPs << "," << p.endWeightedMean.sigmaPs << ","
                << p.endDelta.meanNs << "," << p.endDelta.sigmaPs << ","
                << p.topNearest.meanNs << "," << p.topNearest.sigmaPs << ","
                << p.topAll.meanNs << "," << p.topAll.sigmaPs << ","
                << '"' << p.topIds << '"' << ","
                << ds.speedLeft.velocityCmNs << ","
                << ds.speedRight.velocityCmNs << ","
                << ds.speedDelta.velocityCmNs << ","
                << ds.speedConsensusCmNs << "\n";
        }
    }
}

void WriteVelocityTree(TFile& output, const std::vector<DatasetResult>& datasets) {
    output.cd();
    auto tree = std::make_unique<TTree>("velocity_summary",
                                        "Effective optical propagation speed summary");
    std::string dataset, material, topology;
    double v_left = 0.0, ev_left = 0.0, v_right = 0.0, ev_right = 0.0;
    double v_delta = 0.0, ev_delta = 0.0, v_consensus = 0.0,
           ev_consensus = 0.0, reference = kReferenceVeffCmNs;

    tree->Branch("dataset", &dataset);
    tree->Branch("material", &material);
    tree->Branch("topology", &topology);
    tree->Branch("v_left_cm_ns", &v_left);
    tree->Branch("v_left_err_cm_ns", &ev_left);
    tree->Branch("v_right_cm_ns", &v_right);
    tree->Branch("v_right_err_cm_ns", &ev_right);
    tree->Branch("v_delta_cm_ns", &v_delta);
    tree->Branch("v_delta_err_cm_ns", &ev_delta);
    tree->Branch("v_consensus_cm_ns", &v_consensus);
    tree->Branch("v_consensus_err_cm_ns", &ev_consensus);
    tree->Branch("reference_cm_ns", &reference);

    for (const auto& ds : datasets) {
        dataset = ds.cfg.key;
        material = ds.cfg.material;
        topology = ds.cfg.topology;
        v_left = ds.speedLeft.velocityCmNs;
        ev_left = ds.speedLeft.velocityErrCmNs;
        v_right = ds.speedRight.velocityCmNs;
        ev_right = ds.speedRight.velocityErrCmNs;
        v_delta = ds.speedDelta.velocityCmNs;
        ev_delta = ds.speedDelta.velocityErrCmNs;
        v_consensus = ds.speedConsensusCmNs;
        ev_consensus = ds.speedConsensusErrCmNs;
        tree->Fill();
    }
    tree->Write();
}

}  // namespace fig4sum8

void recreate_figure4_sum8(
    const char* outputRoot =
        "/home/reriosto/SHiP/orchestrator/outputs/figure4_sum8_validation.root",
    const char* ej204EndOnlyDir =
        "/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614",
    const char* ej230EndOnlyDir =
        "/home/reriosto/SHiP/t0minidaq/endonly_mylar_230",
    const char* ej204EndTopDir =
        "/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000",
    const char* ej230EndTopDir =
        "/home/reriosto/SHiP/t0minidaq/results_ej230/data",
    double thresholdPe = 4.0,
    int topGroupSize = 8,
    int expectedEvents = 2000) {
    using namespace fig4sum8;

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gROOT->SetBatch(kTRUE);

    const std::string outputName = outputRoot;
    const std::string outputStem = StemWithoutRoot(outputName);
    const std::string plotDir = outputStem + "_plots";
    const std::string outputParent = gSystem->DirName(outputName.c_str());
    gSystem->mkdir(outputParent.c_str(), true);
    gSystem->mkdir(plotDir.c_str(), true);

    std::vector<DatasetConfig> configs = {
        {"ej204_endonly", "EJ-204", "END-only", ej204EndOnlyDir, false},
        {"ej230_endonly", "EJ-230", "END-only", ej230EndOnlyDir, false},
        {"ej204_endtop", "EJ-204", "END+TOP", ej204EndTopDir, true},
        {"ej230_endtop", "EJ-230", "END+TOP", ej230EndTopDir, true},
    };

    TFile output(outputName.c_str(), "RECREATE");
    if (output.IsZombie()) {
        std::cerr << "[FATAL] Cannot create output ROOT file: " << outputName << "\n";
        return;
    }

    output.cd();
    TNamed methodology(
        "methodology",
        "END-left=global_id 0..7 analog SUM8; END-right=8..15 analog SUM8; "
        "TOP=nearest configurable group (default 8) and all 70 diagnostic; "
        "single-PE response=normalized exp(-t/5ns)-exp(-t/0.5ns); "
        "leading-edge threshold=4 PE by default; Gaussian core=4 iterations in +/-2 sigma; "
        "v_eff from slopes of <t_L>, <t_R>, and <t_L-t_R> versus x; "
        "the combined value is diagnostic because the three estimators are correlated.");
    methodology.Write();
    TParameter<double>("threshold_pe", thresholdPe).Write();
    TParameter<int>("top_group_size", topGroupSize).Write();
    TParameter<int>("expected_events_per_position", expectedEvents).Write();
    TParameter<double>("bar_half_length_mm", kBarHalfLengthMm).Write();

    std::vector<DatasetResult> datasets;

    for (const auto& cfg : configs) {
        std::cout << "\n============================================================\n"
                  << "Dataset: " << cfg.key << "\n"
                  << "Input:   " << cfg.inputDir << "\n"
                  << "============================================================\n";

        const auto files = FindRootFiles(cfg.inputDir);
        if (files.empty()) {
            std::cerr << "[WARN] No photon_hits*.root files found; dataset skipped.\n";
            continue;
        }

        DatasetResult dataset;
        dataset.cfg = cfg;

        output.cd();
        auto* dsDir = output.mkdir(Sanitize(cfg.key).c_str());
        auto* histDir = dsDir->mkdir("time_distributions");

        for (const auto& file : files) {
            auto point = ProcessFile(file, cfg, thresholdPe, topGroupSize,
                                     expectedEvents, histDir);
            if (IsFinite(point.xMm)) dataset.points.push_back(std::move(point));
        }

        std::sort(dataset.points.begin(), dataset.points.end(),
                  [](const PositionResult& a, const PositionResult& b) {
                      return a.xMm < b.xMm;
                  });

        if (dataset.points.empty()) {
            std::cerr << "[WARN] No valid positions for " << cfg.key << "\n";
            continue;
        }

        CreateDatasetObjects(dataset, dsDir, plotDir, topGroupSize);
        datasets.push_back(std::move(dataset));
    }

    WriteVelocityTree(output, datasets);
    output.Write();
    output.Close();

    WriteCsv(outputStem + ".csv", datasets);

    std::cout << "\n============================================================\n"
              << "Analysis finished\n"
              << "ROOT: " << outputName << "\n"
              << "CSV:  " << outputStem << ".csv\n"
              << "Plots:" << plotDir << "/\n"
              << "============================================================\n";

    for (const auto& ds : datasets) {
        std::cout << std::fixed << std::setprecision(3)
                  << ds.cfg.key << ": v_eff = " << ds.speedConsensusCmNs
                  << " +/- " << ds.speedConsensusErrCmNs << " cm/ns"
                  << " (reference " << kReferenceVeffCmNs << " cm/ns)\n";
    }
}
