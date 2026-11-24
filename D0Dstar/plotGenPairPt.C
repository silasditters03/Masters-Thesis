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
    "work-tcc-angular-all/AnalysisResults.root";

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

  outputFile.Close();

  std::cout << "Saved: " << output << '\n';
}