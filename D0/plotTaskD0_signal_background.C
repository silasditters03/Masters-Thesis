// Modular TaskD0 plotting macro for reconstructed-data and MC analysis outputs.
//
// Auto mode (default) detects which histogram families were filled:
//   Data-mode: original reconstructed-candidate histograms (hMass, hMassFinerBinning,
//              hMassVsPhi, etc.) -> 1D X projections, including an inclusive mass
//              histogram in raw counts for the mass fitter.
//   MC-mode:   14 cutQA signal/background comparisons; original MC mass spectra
//              (separate truth categories and their inclusive sum).
//
// Both data and MC inputs use the SAME macro; no MC truth is required for data.
// The 14 cutQA observables use the new 28-hist task when available, or the older
// 112-hist Pre distributions as a fallback. Normalization affects the DISPLAYED
// MC curves only; all saved histograms and the inclusive mass remain raw counts.
//
// Examples (run from the D0 directory):
// root -l -b -q 'plotTaskD0_signal_background.C("work-reco-qa/AnalysisResults.root","D0PlotData.root")'
// root -l -b -q 'plotTaskD0_signal_background.C("work-cutqa-all200/AnalysisResults.root","D0PlotMC.root",true)'
// Optional fifth argument: "auto" (default), "data", "mc", or "both".

#include <TCanvas.h>
#include <TCollection.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TObject.h>
#include <TStyle.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace d0sigbkg {

// Handle TDirectories, TLists, and names containing literal slashes.
TObject* find(TObject* parent, const std::string& path)
{
  if (!parent || path.empty()) return nullptr;
  TObject* exact = nullptr;
  if (auto* dir = dynamic_cast<TDirectory*>(parent)) {
    exact = dir->Get(path.c_str());
  } else if (auto* coll = dynamic_cast<TCollection*>(parent)) {
    exact = coll->FindObject(path.c_str());
  }
  if (exact) return exact;

  const auto slash = path.find('/');
  if (slash == std::string::npos) return nullptr;
  const std::string head = path.substr(0, slash);
  const std::string tail = path.substr(slash + 1);
  TObject* child = nullptr;
  if (auto* dir = dynamic_cast<TDirectory*>(parent)) {
    child = dir->Get(head.c_str());
  } else if (auto* coll = dynamic_cast<TCollection*>(parent)) {
    child = coll->FindObject(head.c_str());
  }
  return find(child, tail);
}

// Project any 1D/2D/3D histogram onto X (the observable or invariant-mass axis).
// Include underflow/overflow along the integrated Y/Z dimensions.
std::unique_ptr<TH1> projectX(TObject* registry, const std::string& source,
                              const std::string& name)
{
  auto* hist = dynamic_cast<TH1*>(find(registry, source));
  if (!hist) return nullptr;
  TH1* result = nullptr;
  if (auto* h3 = dynamic_cast<TH3*>(hist)) {
    result = h3->ProjectionX(name.c_str(), 0, h3->GetNbinsY() + 1,
                             0, h3->GetNbinsZ() + 1, "e");
  } else if (auto* h2 = dynamic_cast<TH2*>(hist)) {
    result = h2->ProjectionX(name.c_str(), 0, h2->GetNbinsY() + 1, "e");
  } else if (hist->GetDimension() == 1) {
    result = static_cast<TH1*>(hist->Clone(name.c_str()));
  }
  if (result) result->SetDirectory(nullptr);
  return std::unique_ptr<TH1>(result);
}

bool populated(const TH1* hist)
{
  return hist && hist->GetEntries() > 0;
}

struct Single {
  const char* source;
  const char* title;
  const char* axis;
};

// Same 14 reconstructed-observable projections as the older data macro.
const Single reconstructed[] = {
  {"hMass", "Reconstructed candidates: invariant mass", "m(K#pi) (GeV/c^{2})"},
  {"hMassFinerBinning", "Reconstructed candidates: invariant mass, narrow range", "m(K#pi) (GeV/c^{2})"},
  {"hMassVsPhi", "Reconstructed candidates: invariant mass, integrated over p_{T} and #phi", "m(K#pi) (GeV/c^{2})"},
  {"hPtProng0", "Reconstructed candidates: positive-prong p_{T}", "p_{T}(prong 0) (GeV/c)"},
  {"hPtProng1", "Reconstructed candidates: negative-prong p_{T}", "p_{T}(prong 1) (GeV/c)"},
  {"hd0Prong0", "Reconstructed candidates: positive-prong impact parameter", "Signed d_{0}(prong 0) (cm)"},
  {"hd0Prong1", "Reconstructed candidates: negative-prong impact parameter", "Signed d_{0}(prong 1) (cm)"},
  {"hd0d0", "Reconstructed candidates: impact-parameter product", "d_{0,0}d_{0,1} (cm^{2})"},
  {"hCTS", "Reconstructed candidates: helicity cosine, D^{0} hypothesis", "cos #theta*"},
  {"hCPA", "Reconstructed candidates: cosine of pointing angle", "cos pointing angle"},
  {"hCPAXYFinerBinning", "Reconstructed candidates: cosine of pointing angle XY", "cos pointing angle XY"},
  {"hNormalisedDecLengthxy", "Reconstructed candidates: normalized decay length XY", "L_{XY}/#sigma_{L_{XY}}"},
  {"hDecLength", "Reconstructed candidates: decay length", "L (cm)"},
  {"hDecLengthxy", "Reconstructed candidates: decay length XY", "L_{XY} (cm)"}
};

struct Cut {
  const char* key;
  const char* source;
  const char* title;
  const char* axis;
};

// 15 selector cuts use 14 distinct observables (both 3D decay-length cuts use decLen).
const Cut cuts[] = {
  {"mass_window", "m", "Invariant-mass difference", "|m(K#pi) - m(D^{0})| (GeV/c^{2})"},
  {"candidate_DCA_XY", "DCA", "D^{0} candidate transverse impact parameter", "|d_{0,XY}(D^{0})| (cm)"},
  {"cos_theta_star", "cosThetaStar", "Absolute helicity cosine", "|cos #theta*|"},
  {"kaon_pt", "ptK", "Kaon transverse momentum", "p_{T}(K) (GeV/c)"},
  {"pion_pt", "ptPi", "Pion transverse momentum", "p_{T}(#pi) (GeV/c)"},
  {"kaon_DCA_XY", "d0K", "Kaon track DCA XY to primary vertex", "|DCA_{XY}(K)| (cm)"},
  {"pion_DCA_XY", "d0pi", "Pion track DCA XY to primary vertex", "|DCA_{XY}(#pi)| (cm)"},
  {"dca_product", "d0d0", "Signed daughter impact-parameter product", "d_{0,0} d_{0,1} (cm^{2})"},
  {"pointing_angle", "cpa", "Cosine of pointing angle", "cos pointing angle"},
  {"pointing_angle_XY", "cpaXY", "Cosine of pointing angle XY", "cos pointing angle XY"},
  {"norm_decay_length_XY", "normDecLenXY", "Normalised transverse decay length", "L_{XY}/#sigma_{L_{XY}}"},
  {"decay_length", "decLen", "Three-dimensional decay length", "L (cm)"},
  {"decay_length_XY", "decLenXY", "Transverse decay length", "L_{XY} (cm)"},
  {"norm_daughter_IP", "normDauIPXYMin", "Minimum normalised daughter impact parameter", "min(|d_{0,0}/#sigma_{0}|, |d_{0,1}/#sigma_{1}|)"}
};

void saveHistogram(TH1* histogram, TDirectory* directory)
{
  directory->cd();
  histogram->Write();
}

// Both 28-hist and legacy 112-hist MC cutQA layouts are supported.
std::unique_ptr<TH1> combinedQA(TObject* registry, const Cut& cut,
                                const char* truth, const std::string& name)
{
  const std::string direct = std::string("cutQA/") + cut.source + "_" + truth;
  if (find(registry, direct)) return projectX(registry, direct, name);

  const std::string prefix = std::string("cutQA/") + cut.source + "_";
  auto d0 = projectX(registry, prefix + "D0_Pre_" + truth, name + "_D0_tmp");
  auto d0bar = projectX(registry, prefix + "D0bar_Pre_" + truth, name + "_D0bar_tmp");
  if (!d0 || !d0bar) return nullptr;
  d0->SetName(name.c_str());
  d0->Add(d0bar.get());
  return d0;
}

// Source spectra are original MC hMass{Sig,Bkg}{D0,D0bar} TH3 histograms.
std::unique_ptr<TH1> selectedMCMass(TObject* registry, const char* truth,
                                    const std::string& name)
{
  std::unique_ptr<TH1> result;
  for (const char* hypothesis : {"D0", "D0bar"}) {
    const std::string source = std::string("hMass") + truth + hypothesis;
    auto projected = projectX(registry, source, name + "_" + hypothesis + "_tmp");
    if (!projected) return nullptr;
    if (!result) {
      result = std::move(projected);
      result->SetName(name.c_str());
    } else {
      result->Add(projected.get());
    }
  }
  return result;
}

// An inclusive histogram is deliberately kept unnormalized, for fitting.
void saveInclusive(const TH1* source, TDirectory* histDir, TDirectory* plotDir,
                   const char* sourceDescription)
{
  if (!populated(source)) return;
  std::unique_ptr<TH1> hist(static_cast<TH1*>(source->Clone("invariant_mass_inclusive")));
  hist->SetDirectory(nullptr);
  hist->SetTitle("Selected D^{0} + #bar{D}^{0} candidates; m(K#pi) (GeV/c^{2}); Candidates / bin");
  saveHistogram(hist.get(), histDir);

  TCanvas canvas("c_invariant_mass_inclusive", "Inclusive invariant mass for fitting", 1000, 750);
  canvas.SetLeftMargin(0.14);
  canvas.SetBottomMargin(0.12);
  hist->SetLineColor(kBlack);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.65);
  hist->SetMinimum(0.);
  hist->Draw("E1");
  plotDir->cd();
  canvas.Write("invariant_mass_inclusive");
  std::cout << "Inclusive invariant mass saved (raw counts; source: "
            << sourceDescription << ")\n";
}

// Data processing: no truth labels, only reconstructed-candidate quantities.
// Return the best available raw inclusive mass spectrum for the fitter.
std::unique_ptr<TH1> drawData(TObject* registry, TDirectory* histBase,
                              TDirectory* plotBase)
{
  auto* histDir = histBase->mkdir("reco");
  auto* plotDir = plotBase->mkdir("reco");
  int count = 0;
  for (const auto& spec : reconstructed) {
    const std::string name = std::string(spec.source) + "_1D";
    auto histogram = projectX(registry, spec.source, name);
    if (!populated(histogram.get())) continue;

    histogram->SetTitle(spec.title);
    histogram->GetXaxis()->SetTitle(spec.axis);
    histogram->GetYaxis()->SetTitle("Entries");
    saveHistogram(histogram.get(), histDir);

    TCanvas canvas(("c_reco_" + name).c_str(), spec.title, 1000, 750);
    canvas.SetLeftMargin(0.14);
    canvas.SetBottomMargin(0.12);
    histogram->SetLineColor(kBlack);
    histogram->SetLineWidth(2);
    histogram->SetMinimum(0.);
    histogram->Draw("E1");
    plotDir->cd();
    canvas.Write(name.c_str());
    ++count;
  }
  std::cout << "Reconstructed-data plots: " << count << "/14\n";

  // Prefer narrow mass axis (1.5848-2.1848 GeV/c2, 120 bins in the current task)
  // over the wide hMass axis (0-5 GeV/c2, 500 bins).
  auto mass = projectX(registry, "hMassFinerBinning", "data_mass_for_fitter");
  if (populated(mass.get())) return mass;
  return projectX(registry, "hMass", "data_mass_for_fitter");
}

int drawCutQA(TObject* registry, TDirectory* histDir,
              TDirectory* plotDir, bool normalizeShapes)
{
  int count = 0;
  for (const auto& cut : cuts) {
    const std::string key(cut.key);
    auto signal = combinedQA(registry, cut, "Sig", key + "_signal");
    auto background = combinedQA(registry, cut, "Bkg", key + "_background");
    if (!signal || !background ||
        (!populated(signal.get()) && !populated(background.get()))) {
      std::cerr << "Skipping cutQA " << key << " (missing or empty)\n";
      continue;
    }
    signal->GetXaxis()->SetTitle(cut.axis);
    background->GetXaxis()->SetTitle(cut.axis);
    saveHistogram(signal.get(), histDir);
    saveHistogram(background.get(), histDir);

    std::unique_ptr<TH1> drawSig(static_cast<TH1*>(signal->Clone((key + "_draw_sig").c_str())));
    std::unique_ptr<TH1> drawBkg(static_cast<TH1*>(background->Clone((key + "_draw_bkg").c_str())));
    drawSig->SetDirectory(nullptr);
    drawBkg->SetDirectory(nullptr);
    if (normalizeShapes) {
      const double sigSum = drawSig->Integral(1, drawSig->GetNbinsX());
      const double bkgSum = drawBkg->Integral(1, drawBkg->GetNbinsX());
      if (sigSum > 0.) drawSig->Scale(1. / sigSum);
      if (bkgSum > 0.) drawBkg->Scale(1. / bkgSum);
    }

    TCanvas canvas(("c_" + key).c_str(), cut.title, 1000, 750);
    canvas.SetLogy();
    canvas.SetLeftMargin(0.14);
    canvas.SetBottomMargin(0.12);
    drawSig->SetTitle((std::string("D^{0} + #bar{D}^{0}: ") + cut.title).c_str());
    drawSig->GetXaxis()->SetTitle(cut.axis);
    drawSig->GetYaxis()->SetTitle(normalizeShapes ? "Fraction of candidates / bin" : "Entries");
    drawSig->GetYaxis()->SetTitleOffset(1.45);
    drawSig->SetLineColor(kRed + 1);
    drawBkg->SetLineColor(kBlue + 1);
    drawSig->SetLineWidth(2);
    drawBkg->SetLineWidth(2);
    const double top = std::max(drawSig->GetMaximum(), drawBkg->GetMaximum());
    const double floor = normalizeShapes ? 1.e-8 : 0.1;
    drawSig->SetMinimum(std::min(floor, std::max(1.e-12, top * 0.001)));
    drawSig->SetMaximum(std::max(1.35 * top, 1.e-10));
    drawSig->Draw("HIST");
    drawBkg->Draw("HIST SAME");
    TLegend legend(0.57, 0.76, 0.89, 0.89);
    legend.SetBorderSize(0);
    legend.AddEntry(drawSig.get(), "MC signal", "l");
    legend.AddEntry(drawBkg.get(), "MC background", "l");
    legend.Draw();

    plotDir->cd();
    canvas.Write(key.c_str());
    ++count;
  }
  std::cout << "MC cut-QA comparisons: " << count << "/14\n";
  return count;
}

std::unique_ptr<TH1> drawMCMass(TObject* registry, TDirectory* histDir,
                                TDirectory* plotDir, bool normalizeShapes)
{
  auto signal = selectedMCMass(registry, "Sig", "invariant_mass_signal");
  auto background = selectedMCMass(registry, "Bkg", "invariant_mass_background");
  if (!signal || !background ||
      (!populated(signal.get()) && !populated(background.get()))) {
    std::cout << "MC truth-separated mass histograms missing/empty; skipping MC mass plot.\n";
    return nullptr;
  }
  saveHistogram(signal.get(), histDir);
  saveHistogram(background.get(), histDir);

  std::unique_ptr<TH1> drawSig(static_cast<TH1*>(signal->Clone("draw_mc_mass_sig")));
  std::unique_ptr<TH1> drawBkg(static_cast<TH1*>(background->Clone("draw_mc_mass_bkg")));
  drawSig->SetDirectory(nullptr);
  drawBkg->SetDirectory(nullptr);
  if (normalizeShapes) {
    const double s = drawSig->Integral(1, drawSig->GetNbinsX());
    const double b = drawBkg->Integral(1, drawBkg->GetNbinsX());
    if (s > 0.) drawSig->Scale(1. / s);
    if (b > 0.) drawBkg->Scale(1. / b);
  }

  TCanvas canvas("c_invariant_mass_selected", "Selected MC mass: signal vs background", 1000, 750);
  canvas.SetLogy();
  canvas.SetLeftMargin(0.14);
  canvas.SetBottomMargin(0.12);
  drawSig->SetTitle("Selected-candidate invariant mass: MC signal and background");
  drawSig->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
  drawSig->GetYaxis()->SetTitle(normalizeShapes ? "Fraction of candidates / bin" : "Entries");
  drawSig->GetYaxis()->SetTitleOffset(1.45);
  drawSig->SetLineColor(kRed + 1);
  drawBkg->SetLineColor(kBlue + 1);
  drawSig->SetLineWidth(2);
  drawBkg->SetLineWidth(2);
  const double top = std::max(drawSig->GetMaximum(), drawBkg->GetMaximum());
  const double floor = normalizeShapes ? 1.e-8 : 0.1;
  drawSig->SetMinimum(std::min(floor, std::max(1.e-12, top * 0.001)));
  drawSig->SetMaximum(std::max(1.35 * top, 1.e-10));
  drawSig->Draw("HIST");
  drawBkg->Draw("HIST SAME");
  TLegend legend(0.57, 0.76, 0.89, 0.89);
  legend.SetBorderSize(0);
  legend.AddEntry(drawSig.get(), "MC signal", "l");
  legend.AddEntry(drawBkg.get(), "MC background", "l");
  legend.Draw();
  plotDir->cd();
  canvas.Write("invariant_mass_selected");

  auto inclusive = std::unique_ptr<TH1>(
      static_cast<TH1*>(signal->Clone("mc_mass_for_fitter")));
  inclusive->SetDirectory(nullptr);
  inclusive->Add(background.get());
  return inclusive;
}

} // namespace d0sigbkg

// Compatible with previous calls using 2, 3, or 4 arguments.
// mode=auto: plot each available histogram family, independently.
// mode=data: reconstructed histograms only.
// mode=mc: MC cutQA and MC mass histograms only.
// mode=both: run both modules, skipping absent histogram families.
void plotTaskD0_signal_background(
  const char* input = "AnalysisResults.root",
  const char* output = "D0SignalBackground.root",
  bool normalizeShapes = false,
  const char* taskPath = "hf-task-d0",
  const char* mode = "auto")
{
  const std::string choice(mode);
  if (choice != "auto" && choice != "data" && choice != "mc" && choice != "both") {
    throw std::runtime_error("Mode must be auto, data, mc, or both");
  }
  if (std::string(input) == output) {
    throw std::runtime_error("Input and output ROOT filenames must differ");
  }
  std::unique_ptr<TFile> in(TFile::Open(input, "READ"));
  if (!in || in->IsZombie()) {
    throw std::runtime_error(std::string("Cannot open input: ") + input);
  }
  TObject* task = d0sigbkg::find(in.get(), taskPath);
  if (!task) {
    throw std::runtime_error(std::string("Missing task directory: ") + taskPath);
  }
  TObject* registry = d0sigbkg::find(task, "registry");
  if (!registry) registry = task;

  // Detection depends on filled histograms, not only definitions (both sets
  // of histograms can be registered by the task even if one process is disabled).
  auto dataProbe = d0sigbkg::projectX(registry, "hMassFinerBinning", "probe_reco_mass");
  if (!d0sigbkg::populated(dataProbe.get())) {
    dataProbe = d0sigbkg::projectX(registry, "hMass", "probe_reco_mass_fallback");
  }
  const bool haveData = d0sigbkg::populated(dataProbe.get());
  auto mcQAProbeBkg = d0sigbkg::combinedQA(registry, d0sigbkg::cuts[0], "Bkg", "probe_mc_qa_bkg");
  auto mcQAProbeSig = d0sigbkg::combinedQA(registry, d0sigbkg::cuts[0], "Sig", "probe_mc_qa_sig");
  const bool haveQA = d0sigbkg::populated(mcQAProbeBkg.get()) ||
                      d0sigbkg::populated(mcQAProbeSig.get());
  // The task registers MC histograms even during data processing, so checking
  // existence alone would incorrectly enable the MC module for data output.
  bool haveMCMass = false;
  for (const char* key : {"hMassSigD0", "hMassBkgD0", "hMassSigD0bar", "hMassBkgD0bar"}) {
    auto probe = d0sigbkg::projectX(registry, key, std::string("probe_") + key);
    if (d0sigbkg::populated(probe.get())) {
      haveMCMass = true;
      break;
    }
  }

  const bool enableData = choice == "data" || choice == "both" || (choice == "auto" && haveData);
  const bool enableMC = choice == "mc" || choice == "both" || (choice == "auto" && (haveQA || haveMCMass));

  std::cout << "Input mode (auto detection): reco=" << (haveData ? "filled" : "empty")
            << ", cutQA=" << (haveQA ? "filled" : "empty") << '\n';

  TFile out(output, "RECREATE");
  if (out.IsZombie()) {
    throw std::runtime_error(std::string("Cannot create output: ") + output);
  }
  TDirectory* histDir = out.mkdir("histograms");
  TDirectory* plotDir = out.mkdir("plots");
  gStyle->SetOptStat(0);

  std::unique_ptr<TH1> dataMass;
  if (enableData) {
    dataMass = d0sigbkg::drawData(registry, histDir, plotDir);
    if (!d0sigbkg::populated(dataMass.get())) {
      std::cerr << "No reconstructed mass histogram filled; data module has no mass to fit.\n";
    }
  }

  std::unique_ptr<TH1> mcMass;
  if (enableMC) {
    if (haveQA) {
      d0sigbkg::drawCutQA(registry, histDir, plotDir, normalizeShapes);
    } else {
      std::cout << "No filled MC cutQA histograms; skipping cutQA overlays.\n";
    }
    mcMass = d0sigbkg::drawMCMass(registry, histDir, plotDir, normalizeShapes);
  }

  // Prefer truth-agnostic reconstructed mass whenever available; this is the
  // spectrum to use with a conventional invariant-mass fitter.
  if (d0sigbkg::populated(dataMass.get())) {
    d0sigbkg::saveInclusive(dataMass.get(), histDir, plotDir,
                            "hMassFinerBinning (or hMass) from data processing");
  } else if (d0sigbkg::populated(mcMass.get())) {
    d0sigbkg::saveInclusive(mcMass.get(), histDir, plotDir,
                            "MC signal + background after selection");
  } else {
    std::cerr << "No inclusive invariant mass could be constructed.\n";
  }
  out.Write();
  out.Close();
  std::cout << "Output: " << output << '\n';
}
