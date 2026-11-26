#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>

#include <algorithm>
#include <iostream>
#include <memory>

void plotGenPairPt()
{
  const char* input =
    "/Users/ml2024/SilasThesis/Masters-Thesis/D0Dstar/"
    "work-full/AnalysisResults.root";

  const char* output =
    "/Users/ml2024/SilasThesis/Masters-Thesis/D0Dstar/"
    "GenPairPlots.root";

  std::unique_ptr<TFile> file(TFile::Open(input));
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open: " << input << '\n';
    return;
  }

  const char* names[] = {
    "hPtD0DstarGen",
    "hPtD0InSelectedGenPairs",
    "hPtDstarInSelectedGenPairs"
  };

  const char* labels[] = {
    "D^{0}D^{*} pair",
    "D^{0} in selected pairs",
    "D^{*} in selected pairs"
  };

  const int colours[] = {kBlack, kBlue + 1, kRed + 1};
  TH1* histograms[3] = {};
  double maximum = 0.;

  for (int i = 0; i < 3; ++i) {
    TH1* source = nullptr;
    file->GetObject(
      Form("hf-task-d0-dstar/%s", names[i]), source);

    if (!source) {
      std::cerr << "Missing histogram: " << names[i] << '\n';
      return;
    }

    histograms[i] =
      static_cast<TH1*>(source->Clone(Form("ptOverlay_%d", i)));
    histograms[i]->SetDirectory(nullptr);
    histograms[i]->SetStats(false);
    histograms[i]->SetLineColor(colours[i]);
    histograms[i]->SetLineWidth(2);

    maximum = std::max(maximum, histograms[i]->GetMaximum());
  }

  auto* canvas = new TCanvas("cGenPairPt", "Generated pT spectra", 950, 700);
  canvas->SetLeftMargin(0.13);
  canvas->SetBottomMargin(0.13);
  canvas->SetTicks(1, 1);

  histograms[0]->SetTitle(
    "Selected generated D^{0}-D^{*} pairs;p_{T} (GeV/c);Entries / bin");
  histograms[0]->SetMinimum(0.);
  histograms[0]->SetMaximum(maximum > 0. ? 1.25 * maximum : 1.);

  histograms[0]->Draw("HIST");
  histograms[1]->Draw("HIST SAME");
  histograms[2]->Draw("HIST SAME");

  auto* legend = new TLegend(0.57, 0.69, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  for (int i = 0; i < 3; ++i) {
    legend->AddEntry(histograms[i], labels[i], "l");
  }

  legend->Draw();

  // Delta phi with statistical error bars.
  TH1* sourceDeltaPhi = nullptr;
  file->GetObject(
    "hf-task-d0-dstar/hDeltaPhiGenPairs", sourceDeltaPhi);

  if (!sourceDeltaPhi) {
    std::cerr << "Missing histogram: hDeltaPhiGenPairs\n";
    return;
  }

  auto* hDeltaPhi =
    static_cast<TH1*>(sourceDeltaPhi->Clone("deltaPhiPlot"));
  hDeltaPhi->SetDirectory(nullptr);

  // Keep the original binning. Set to 4 for 16 bins instead of 64.
  const int deltaPhiRebin = 1;
  if (deltaPhiRebin > 1) {
    hDeltaPhi->Rebin(deltaPhiRebin);
  }

  hDeltaPhi->SetStats(false);
  hDeltaPhi->SetTitle(
    "Selected generated D^{0}-D^{*} pairs;"
    "#Delta#varphi = #varphi(D^{0}) - #varphi(D*) (rad);"
    "Pairs / bin");
  hDeltaPhi->SetLineColor(kBlue + 1);
  hDeltaPhi->SetMarkerColor(kBlue + 1);
  hDeltaPhi->SetMarkerStyle(20);
  hDeltaPhi->SetMarkerSize(0.8);
  hDeltaPhi->SetMinimum(0.);

  double deltaPhiMaximum = 0.;
  for (int bin = 1; bin <= hDeltaPhi->GetNbinsX(); ++bin) {
    deltaPhiMaximum = std::max(
      deltaPhiMaximum,
      hDeltaPhi->GetBinContent(bin) + hDeltaPhi->GetBinError(bin));
  }
  hDeltaPhi->SetMaximum(
    deltaPhiMaximum > 0. ? 1.2 * deltaPhiMaximum : 1.);

  auto* canvasDeltaPhi =
    new TCanvas("cGenPairDeltaPhi", "Generated pair delta phi", 950, 700);
  canvasDeltaPhi->SetLeftMargin(0.13);
  canvasDeltaPhi->SetBottomMargin(0.13);
  canvasDeltaPhi->SetTicks(1, 1);

  hDeltaPhi->Draw("E1");

  TFile outputFile(output, "UPDATE");
  if (outputFile.IsZombie()) {
    std::cerr << "Cannot open output file: " << output << '\n';
    return;
  }

  outputFile.cd();

  // Save the canvas with all three spectra overlaid.
  canvas->Write("cGenPairPt", TObject::kOverwrite);

  // Also save the individual histograms for later use.
  for (int i = 0; i < 3; ++i) {
    histograms[i]->Write(names[i], TObject::kOverwrite);
  }

  canvasDeltaPhi->Write("cGenPairDeltaPhi", TObject::kOverwrite);
  hDeltaPhi->Write("hDeltaPhiGenPairs", TObject::kOverwrite);
  outputFile.Close();

  std::cout << "Saved: " << output << '\n';
}