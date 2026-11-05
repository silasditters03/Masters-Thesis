#include <TCanvas.h>
#include <TCollection.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TStyle.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace d0plot {

// Read from ROOT directories or serialized histogram collections.
TObject* get(TObject* parent, const std::string& path)
{
  if (!parent || path.empty()) return nullptr;

  const auto slash = path.find('/');
  const auto name = path.substr(0, slash);
  TObject* child = nullptr;

  if (auto* directory = dynamic_cast<TDirectory*>(parent)) {
    child = directory->Get(name.c_str());
  } else if (auto* collection = dynamic_cast<TCollection*>(parent)) {
    child = collection->FindObject(name.c_str());
  }

  return slash == std::string::npos
           ? child
           : get(child, path.substr(slash + 1));
}

// Project the observable on X.
// Integrate all other axes, including their underflow and overflow.
std::unique_ptr<TH1> project(
  TObject* task,
  const char* source,
  const std::string& name)
{
  auto* histogram = dynamic_cast<TH1*>(get(task, source));
  if (!histogram) return nullptr;

  TH1* result = nullptr;

  if (auto* h3 = dynamic_cast<TH3*>(histogram)) {
    result = h3->ProjectionX(
      name.c_str(),
      0, h3->GetNbinsY() + 1,
      0, h3->GetNbinsZ() + 1,
      "e");
  } else if (auto* h2 = dynamic_cast<TH2*>(histogram)) {
    result = h2->ProjectionX(
      name.c_str(),
      0, h2->GetNbinsY() + 1,
      "e");
  } else if (histogram->GetDimension() == 1) {
    result = static_cast<TH1*>(
      histogram->Clone(name.c_str()));
  }

  if (!result) return nullptr;

  result->SetDirectory(nullptr);
  return std::unique_ptr<TH1>(result);
}

bool populated(const TH1* histogram)
{
  return histogram &&
         histogram->Integral(
           0, histogram->GetNbinsX() + 1) > 0.;
}

// Definitions are separate from projection and drawing code.
// Add further available observables to the tables below.
struct Single {
  const char* source;
  const char* title;
  const char* axis;
};

struct Pair {
  const char* key;
  const char* title;
  const char* axis;
  const char* signal;
  const char* background;
};

class Writer {
  TObject* task;
  std::string tag;
  TDirectory* histograms;
  TDirectory* plots;
  TCanvas canvas;

  void store(TH1* histogram, const char* axis)
  {
    histogram->GetXaxis()->SetTitle(axis);
    histogram->GetYaxis()->SetTitle("Entries");

    histograms->cd();
    histogram->Write();

    std::cout
      << tag << "/" << histogram->GetName()
      << ": visible=" << histogram->Integral()
      << ", underflow=" << histogram->GetBinContent(0)
      << ", overflow="
      << histogram->GetBinContent(
           histogram->GetNbinsX() + 1)
      << '\n';
  }

  void saveCanvas(const std::string& key)
  {
    canvas.SetName((tag + "_" + key).c_str());
    canvas.SetTitle(key.c_str());
    canvas.Modified();
    canvas.Update();

    plots->cd();
    canvas.Write(key.c_str());
    ++count;
  }

public:
  int count = 0;

  Writer(
    TObject* source,
    TFile& output,
    const std::string& category)
    : task(source),
      tag(category),
      histograms(
        output.GetDirectory("histograms")->mkdir(
          category.c_str())),
      plots(
        output.GetDirectory("plots")->mkdir(
          category.c_str())),
      canvas(
        (category + "_canvas").c_str(),
        "TaskD0",
        1000, 750)
  {
    canvas.SetLeftMargin(0.14);
  }

  void single(const Single& spec)
  {
    const std::string key =
      std::string(spec.source) + "_1D";

    auto histogram = project(task, spec.source, key);

    if (!populated(histogram.get())) {
      std::cout
        << "Skip " << tag << ": " << spec.source
        << " missing/empty\n";
      return;
    }

    histogram->SetTitle(spec.title);
    store(histogram.get(), spec.axis);

    canvas.Clear();

    histogram->SetLineColor(kBlack);
    histogram->SetLineWidth(2);
    histogram->SetMinimum(0.);
    histogram->Draw("E1");

    saveCanvas(key);
  }

  void pair(const Pair& spec)
  {
    auto signal = project(
      task, spec.signal,
      std::string(spec.key) + "_signal");

    auto background = project(
      task, spec.background,
      std::string(spec.key) + "_background");

    if (!signal || !background ||
        (!populated(signal.get()) &&
         !populated(background.get()))) {
      std::cout
        << "Skip MC pair: " << spec.key
        << " missing/empty\n";
      return;
    }

    signal->SetTitle(spec.title);
    background->SetTitle(spec.title);

    store(signal.get(), spec.axis);
    store(background.get(), spec.axis);

    canvas.Clear();

    // Raw counts: no normalization.
    signal->SetLineColor(kBlue + 1);
    background->SetLineColor(kRed + 1);
    signal->SetLineWidth(2);
    background->SetLineWidth(2);

    signal->SetMinimum(0.);
    signal->SetMaximum(
      std::max(
        1.,
        1.35 * std::max(
          signal->GetMaximum(),
          background->GetMaximum())));

    signal->Draw("HIST");
    background->Draw("HIST SAME");

    TLegend legend(0.60, 0.76, 0.88, 0.88);
    legend.SetBorderSize(0);

    legend.AddEntry(
      signal.get(),
      populated(signal.get())
        ? "MC signal"
        : "MC signal (empty)",
      "l");

    legend.AddEntry(
      background.get(),
      populated(background.get())
        ? "MC background"
        : "MC background (empty)",
      "l");

    legend.Draw();
    saveCanvas(spec.key);
  }
};

// Histograms filled by TaskD0's reconstructed/data-processing branch.
// These can also exist in simulation if that branch was enabled.
void reconstructed(Writer& writer)
{
  const std::vector<Single> specifications{
    {
      "hMass",
      "Reconstructed candidates: invariant mass",
      "m(K#pi) (GeV/c^{2})"
    },
    {
      "hMassFinerBinning",
      "Reconstructed candidates: invariant mass, narrow range",
      "m(K#pi) (GeV/c^{2})"
    },
    {
      "hMassVsPhi",
      "Reconstructed candidates: mass, integrated over p_{T} and #phi",
      "m(K#pi) (GeV/c^{2})"
    },
    {
      "hPtProng0",
      "Reconstructed candidates: positive-prong p_{T}",
      "p_{T}(prong 0) (GeV/c)"
    },
    {
      "hPtProng1",
      "Reconstructed candidates: negative-prong p_{T}",
      "p_{T}(prong 1) (GeV/c)"
    },
    {
      "hd0Prong0",
      "Reconstructed candidates: positive-prong impact parameter",
      "Signed candidate d_{0}(prong 0) (cm)"
    },
    {
      "hd0Prong1",
      "Reconstructed candidates: negative-prong impact parameter",
      "Signed candidate d_{0}(prong 1) (cm)"
    },
    {
      "hd0d0",
      "Reconstructed candidates: impact-parameter product",
      "d_{0,0}d_{0,1} (cm^{2})"
    },
    {
      "hCTS",
      "Reconstructed candidates: helicity cosine, D^{0} hypothesis",
      "cos #theta*"
    },
    {
      "hCPA",
      "Reconstructed candidates: cosine of pointing angle",
      "cos pointing angle"
    },
    {
      "hCPAXYFinerBinning",
      "Reconstructed candidates: cosine of pointing angle XY",
      "cos pointing angle XY"
    },
    {
      "hNormalisedDecLengthxy",
      "Reconstructed candidates: normalized decay length XY",
      "L_{XY}/#sigma_{L_{XY}}"
    },
    {
      "hDecLength",
      "Reconstructed candidates: decay length (min/max cuts)",
      "L (cm)"
    },
    {
      "hDecLengthxy",
      "Reconstructed candidates: decay length XY",
      "L_{XY} (cm)"
    }
  };

  for (const auto& spec : specifications) {
    writer.single(spec);
  }
}

// Truth-separated MC histograms.
// Non-mass pairs in the current TaskD0 describe the D0 hypothesis.
void simulation(Writer& writer)
{
  for (const auto* source : {
         "hMassSigBkgD0",
         "hMassSigBkgD0bar",
         "hMassReflBkgD0",
         "hMassReflBkgD0bar"}) {
    writer.single({
      source,
      source,
      "m(K#pi) (GeV/c^{2})"
    });
  }

  const std::vector<Pair> specifications{
    {
      "mass_D0", "MC: D^{0} mass",
      "m(K#pi) (GeV/c^{2})",
      "hMassSigD0", "hMassBkgD0"
    },
    {
      "mass_D0bar", "MC: #bar{D}^{0} mass",
      "m(K#pi) (GeV/c^{2})",
      "hMassSigD0bar", "hMassBkgD0bar"
    },
    {
      "cos_theta_star", "MC D^{0}: helicity cosine",
      "cos #theta*",
      "hCTSSig", "hCTSBkg"
    },
    {
      "pt_K", "MC D^{0}: kaon p_{T}",
      "p_{T}(K) (GeV/c)",
      "hPtProng1Sig", "hPtProng1Bkg"
    },
    {
      "pt_Pi", "MC D^{0}: pion p_{T}",
      "p_{T}(#pi) (GeV/c)",
      "hPtProng0Sig", "hPtProng0Bkg"
    },
    {
      "d0K_candidate",
      "MC D^{0}: kaon candidate-prong impact parameter",
      "Signed candidate d_{0}(K) (cm)",
      "hd0Prong1Sig", "hd0Prong1Bkg"
    },
    {
      "d0pi_candidate",
      "MC D^{0}: pion candidate-prong impact parameter",
      "Signed candidate d_{0}(#pi) (cm)",
      "hd0Prong0Sig", "hd0Prong0Bkg"
    },
    {
      "d0d0", "MC D^{0}: impact-parameter product",
      "d_{0,0}d_{0,1} (cm^{2})",
      "hd0d0Sig", "hd0d0Bkg"
    },
    {
      "cos_pointing_angle",
      "MC D^{0}: cosine of pointing angle",
      "cos pointing angle",
      "hCPASig", "hCPABkg"
    },
    {
      "cos_pointing_angle_xy",
      "MC D^{0}: cosine of pointing angle XY",
      "cos pointing angle XY",
      "hCPAxySig", "hCPAxyBkg"
    },
    {
      "norm_decay_length_xy",
      "MC D^{0}: normalized decay length XY",
      "L_{XY}/#sigma_{L_{XY}}",
      "hNormalisedDecLengthXYSig",
      "hNormalisedDecLengthXYBkg"
    },
    {
      "decay_length",
      "MC D^{0}: decay length (min/max cuts)",
      "L (cm)",
      "hDecLengthSig", "hDecLengthBkg"
    },
    {
      "decay_length_xy",
      "MC D^{0}: decay length XY",
      "L_{XY} (cm)",
      "hDecLengthXYSig", "hDecLengthXYBkg"
    }
  };

  for (const auto& spec : specifications) {
    writer.pair(spec);
  }
}

} // namespace d0plot

void plotTaskD0(
  const char* input = "AnalysisResults.root",
  const char* output = "TaskD0Plots.root",
  const char* mode = "data",
  const char* taskPath = "hf-task-d0")
{
  const std::string choice(mode);

  if (choice != "data" &&
      choice != "mc" &&
      choice != "both") {
    throw std::runtime_error(
      "Mode must be data, mc, or both");
  }

  if (std::string(input) == output) {
    throw std::runtime_error(
      "Use different input and output files");
  }

  std::unique_ptr<TFile> file(
    TFile::Open(input, "READ"));

  if (!file || file->IsZombie()) {
    throw std::runtime_error(
      "Cannot open input ROOT file");
  }

  auto* task = d0plot::get(file.get(), taskPath);

  if (!task) {
    throw std::runtime_error(
      std::string("Missing task path: ") + taskPath);
  }

  if (auto* registry = d0plot::get(task, "registry")) {
    task = registry;
  }

  TFile result(output, "RECREATE");

  if (result.IsZombie()) {
    throw std::runtime_error(
      "Cannot create output ROOT file");
  }

  result.mkdir("histograms");
  result.mkdir("plots");

  gStyle->SetOptStat(0);

  if (choice == "data" || choice == "both") {
    d0plot::Writer writer(task, result, "reco");
    d0plot::reconstructed(writer);

    std::cout
      << "Saved " << writer.count
      << " reconstructed-candidate plots\n";
  }

  if (choice == "mc" || choice == "both") {
    d0plot::Writer writer(task, result, "mc");
    d0plot::simulation(writer);

    std::cout
      << "Saved " << writer.count
      << " MC plots\n";
  }

  std::cout << "Output: " << output << '\n';
}