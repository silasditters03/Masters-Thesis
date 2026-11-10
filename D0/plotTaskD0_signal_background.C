// 14 selector-observable overlays: true D0/D0bar signal vs background.
// Supports BOTH source layouts:
//   New 28-hist task: cutQA/<observable>_Sig / _Bkg
//   Previous 112-hist task: cutQA/<observable>_D0_Pre_* and _D0bar_Pre_*
// The new layout is preferred when available. Histograms are projected over pT.
//
// ROOT usage:
// root -l -b -q 'plotTaskD0_signal_background.C("work-cutqa-all200/AnalysisResults.root","D0SignalBackground.root")'
// Normalize displayed curves (RAW histograms are still saved):
// root -l -b -q 'plotTaskD0_signal_background.C("work-cutqa-all200/AnalysisResults.root","D0SignalBackground_shapes.root",true)'

#include <TCanvas.h>
#include <TCollection.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TObject.h>
#include <TStyle.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace d0sigbkg {

// Works for both ROOT directories and histogram TLists. Also handles
// names stored literally with a slash, such as "cutQA/cpa_D0_Pre_Sig".
TObject* find(TObject* parent, const std::string& path)
{
  if (!parent || path.empty()) return nullptr;

  TObject* exact = nullptr;
  if (auto* dir = dynamic_cast<TDirectory*>(parent)) {
    exact = dir->Get(path.c_str());
  } else if (auto* list = dynamic_cast<TCollection*>(parent)) {
    exact = list->FindObject(path.c_str());
  }
  if (exact) return exact;

  const auto slash = path.find('/');
  if (slash == std::string::npos) return nullptr;

  const auto first = path.substr(0, slash);
  const auto rest = path.substr(slash + 1);
  TObject* child = nullptr;
  if (auto* dir = dynamic_cast<TDirectory*>(parent)) {
    child = dir->Get(first.c_str());
  } else if (auto* list = dynamic_cast<TCollection*>(parent)) {
    child = list->FindObject(first.c_str());
  }
  return find(child, rest);
}

struct Cut {
  const char* key;
  const char* source;
  const char* title;
  const char* xAxis;
};

// 15 selector cuts correspond to 14 different plotted observables:
// min and max 3D decay length use the same observable (decLen).
const Cut cuts[] = {
  {"mass_window",        "m",             "Invariant-mass difference",               "|m(K#pi) - m(D^{0})| (GeV/c^{2})"},
  {"candidate_DCA_XY",   "DCA",           "D^{0} candidate transverse impact parameter", "|d_{0,XY}(D^{0})| (cm)"},
  {"cos_theta_star",     "cosThetaStar",  "Absolute helicity cosine",                "|cos #theta*|"},
  {"kaon_pt",            "ptK",           "Kaon transverse momentum",                "p_{T}(K) (GeV/c)"},
  {"pion_pt",            "ptPi",          "Pion transverse momentum",                "p_{T}(#pi) (GeV/c)"},
  {"kaon_DCA_XY",        "d0K",           "Kaon track DCA XY to primary vertex",    "|DCA_{XY}(K)| (cm)"},
  {"pion_DCA_XY",        "d0pi",          "Pion track DCA XY to primary vertex",    "|DCA_{XY}(#pi)| (cm)"},
  {"dca_product",        "d0d0",          "Signed daughter impact-parameter product", "d_{0,0} d_{0,1} (cm^{2})"},
  {"pointing_angle",     "cpa",           "Cosine of pointing angle",              "cos pointing angle"},
  {"pointing_angle_XY",  "cpaXY",         "Cosine of pointing angle XY",           "cos pointing angle XY"},
  {"norm_decay_length_XY","normDecLenXY", "Normalised transverse decay length",     "L_{XY}/#sigma_{L_{XY}}"},
  {"decay_length",       "decLen",        "Three-dimensional decay length",          "L (cm)"},
  {"decay_length_XY",    "decLenXY",      "Transverse decay length",                "L_{XY} (cm)"},
  {"norm_daughter_IP",   "normDauIPXYMin", "Minimum normalised daughter impact parameter", "min(|d_{0,0}/#sigma_{0}|, |d_{0,1}/#sigma_{1}|)"}
};

// Read one hypothesis' Pre histogram and project it onto the observable axis.
// pT underflow/overflow are included, matching the previous plotting macro.
std::unique_ptr<TH1> projectPre(TObject* registry, const char* source,
                                const char* hypothesis, const char* truth,
                                const std::string& outputName)
{
  const std::string path = std::string("cutQA/") + source + "_" + hypothesis + "_Pre_" + truth;
  auto* h2 = dynamic_cast<TH2*>(find(registry, path));
  if (!h2) {
    std::cerr << "Missing TH2: " << path << '\n';
    return nullptr;
  }
  auto* projection = h2->ProjectionX(outputName.c_str(), 0, h2->GetNbinsY() + 1, "e");
  if (!projection) return nullptr;
  projection->SetDirectory(nullptr);
  return std::unique_ptr<TH1>(projection);
}

// Prefer the direct combined histogram used by the 28-hist task. For
// legacy analysis outputs, combine the two hypotheses from the Pre stage.
std::unique_ptr<TH1> combined(TObject* registry, const Cut& cut,
                              const char* truth, const std::string& outputName)
{
  const std::string combinedPath = std::string("cutQA/") + cut.source + "_" + truth;
  if (auto* direct = dynamic_cast<TH2*>(find(registry, combinedPath))) {
    auto* projection = direct->ProjectionX(outputName.c_str(), 0, direct->GetNbinsY() + 1, "e");
    if (!projection) return nullptr;
    projection->SetDirectory(nullptr);
    return std::unique_ptr<TH1>(projection);
  }

  // Backward compatibility with existing results from the previous 112-hist code.
  auto d0 = projectPre(registry, cut.source, "D0", truth, outputName + "_D0_temp");
  auto d0bar = projectPre(registry, cut.source, "D0bar", truth, outputName + "_D0bar_temp");
  if (!d0 || !d0bar) return nullptr;
  d0->SetName(outputName.c_str());
  d0->Add(d0bar.get());
  return d0;
}

} // namespace d0sigbkg

// normalizeShapes=false (default): raw candidate counts in both curves.
// normalizeShapes=true: unit-area curves, to compare distribution shapes.
// In either case, histograms saved in the output file keep RAW counts.
void plotTaskD0_signal_background(
  const char* input = "AnalysisResults.root",
  const char* output = "D0SignalBackground.root",
  bool normalizeShapes = false,
  const char* taskPath = "hf-task-d0")
{
  if (std::string(input) == output) {
    throw std::runtime_error("Input and output files must be different");
  }

  std::unique_ptr<TFile> in(TFile::Open(input, "READ"));
  if (!in || in->IsZombie()) {
    throw std::runtime_error(std::string("Cannot open input: ") + input);
  }

  TObject* task = d0sigbkg::find(in.get(), taskPath);
  if (!task) {
    throw std::runtime_error(std::string("Missing ROOT task: ") + taskPath);
  }
  TObject* registry = d0sigbkg::find(task, "registry");
  if (!registry) registry = task;

  TFile out(output, "RECREATE");
  if (out.IsZombie()) {
    throw std::runtime_error(std::string("Cannot create output: ") + output);
  }
  TDirectory* histDir = out.mkdir("histograms");
  TDirectory* plotDir = out.mkdir("plots");
  gStyle->SetOptStat(0);

  int produced = 0;
  for (const auto& cut : d0sigbkg::cuts) {
    const std::string base = cut.key;
    auto signal = d0sigbkg::combined(registry, cut, "Sig", base + "_signal");
    auto background = d0sigbkg::combined(registry, cut, "Bkg", base + "_background");
    if (!signal || !background) {
      std::cerr << "Skipping " << base << " (missing input histogram)\n";
      continue;
    }

    // Preserve raw counts for quantitative comparisons in the output ROOT file.
    histDir->cd();
    signal->Write();
    background->Write();

    // Separate draw copies so normalizing the display never changes saved counts.
    std::unique_ptr<TH1> drawSig(static_cast<TH1*>(signal->Clone((base + "_draw_signal").c_str())));
    std::unique_ptr<TH1> drawBkg(static_cast<TH1*>(background->Clone((base + "_draw_background").c_str())));
    drawSig->SetDirectory(nullptr);
    drawBkg->SetDirectory(nullptr);

    if (normalizeShapes) {
      const double sumSig = drawSig->Integral(1, drawSig->GetNbinsX());
      const double sumBkg = drawBkg->Integral(1, drawBkg->GetNbinsX());
      if (sumSig > 0) drawSig->Scale(1.0 / sumSig);
      if (sumBkg > 0) drawBkg->Scale(1.0 / sumBkg);
    }

    TCanvas canvas(("c_" + base).c_str(), cut.title, 1000, 750);
    canvas.SetLogy();
    canvas.SetLeftMargin(0.14);
    canvas.SetBottomMargin(0.12);

    drawSig->SetTitle((std::string("D^{0} + #bar{D}^{0}: ") + cut.title).c_str());
    drawSig->GetXaxis()->SetTitle(cut.xAxis);
    drawSig->GetYaxis()->SetTitle(normalizeShapes ? "Fraction of candidates / bin" : "Entries");
    drawSig->GetYaxis()->SetTitleOffset(1.45);
    drawSig->SetLineColor(kRed + 1);
    drawBkg->SetLineColor(kBlue + 1);
    drawSig->SetLineWidth(2);
    drawBkg->SetLineWidth(2);
    drawSig->SetMinimum(0.0);
    drawSig->SetMaximum(std::max(1.e-12, 1.25 * std::max(drawSig->GetMaximum(), drawBkg->GetMaximum())));

    drawSig->Draw("HIST");
    drawBkg->Draw("HIST SAME");
    TLegend legend(0.56, 0.76, 0.89, 0.89);
    legend.SetBorderSize(0);
    legend.AddEntry(drawSig.get(), "MC signal", "l");
    legend.AddEntry(drawBkg.get(), "MC background", "l");
    legend.Draw();

    plotDir->cd();
    canvas.Write(base.c_str());
    ++produced;

    std::cout << base
              << ": signal=" << signal->GetEntries()
              << ", background=" << background->GetEntries()
              << " (raw entries; D0 + D0bar)\n";
  }

  out.Write();
  out.Close();
  std::cout << "Saved " << produced << "/14 comparison canvases to " << output << '\n';
}
