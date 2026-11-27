#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <iostream>
#include <algorithm>

void projectDstarSelections(const char* input)
{
  TFile file(input, "READ");
  if (file.IsZombie()) {
    std::cerr << "Cannot open: " << input << '\n';
    return;
  }

  TH2* source = nullptr;
  file.GetObject(
    "hf-candidate-selector-dstar-to-d0-pi/QA/hSelections",
    source);

  if (!source) {
    std::cerr << "Cannot find the Dstar selector QA/hSelections\n";
    return;
  }

  // Sum all pT bins, including underflow and overflow.
  auto* projection = source->ProjectionX(
    "hDstarSelectionProjection", 0, source->GetNbinsY() + 1, "e");

  projection->SetDirectory(nullptr);
  projection->SetStats(false);
  projection->SetTitle("Dstar selector cut flow;Selection;Candidates");
  projection->SetMinimum(0.5);
  projection->SetMaximum(2. * std::max(1., projection->GetMaximum()));
  projection->SetMarkerStyle(20);

  std::cout << "\nDstar selector counts:\n";
  for (int bin = 1; bin <= projection->GetNbinsX(); ++bin) {
    std::cout << projection->GetXaxis()->GetBinLabel(bin)
              << ": " << projection->GetBinContent(bin) << '\n';
  }

  auto* canvas = new TCanvas(
    "cDstarSelectionProjection", "Dstar selector cut flow", 1200, 750);
  canvas->SetBottomMargin(0.35);
  canvas->SetLogy();
  projection->GetXaxis()->LabelsOption("v");
  projection->GetXaxis()->SetLabelSize(0.025);
  projection->Draw("E1");

  TFile output(
    "/Users/ml2024/SilasThesis/Masters-Thesis/D0Dstar/GenPairPlots.root",
    "UPDATE");

  if (output.IsZombie()) {
    std::cerr << "Cannot open GenPairPlots.root for writing\n";
    return;
  }

  projection->Write("hDstarSelectionProjection", TObject::kOverwrite);
  canvas->Write("cDstarSelectionProjection", TObject::kOverwrite);
}