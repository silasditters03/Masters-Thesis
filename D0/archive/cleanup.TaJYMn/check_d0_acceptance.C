
    
  
// Read-only diagnostic. Prints only presence/absence status, never values/counts.
#include <TFile.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TString.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <vector>
void check_d0_acceptance(const char* filename = "data/002/AO2D.root")
{
  TFile file(filename);
  if (file.IsZombie()) {
    std::cout << "CHECK_FAILED: cannot open input\n";
    return;
  }
  bool inRapidity = false, inRapidityNotBackground = false, failed = false;
  bool checked = false;
  TIter keys(file.GetListOfKeys());
  while (auto* key = static_cast<TKey*>(keys())) {
    if (!TString(key->GetName()).BeginsWith("DF_")) continue;
    auto* tree = file.Get<TTree>(Form("%s/O2mcparticle_001", key->GetName()));
    if (!tree) { failed = true; continue; }
    const char* names[] = {"fPdgCode", "fIndexSlice_Daughters", "fPx", "fPy", "fPz", "fFlags"};
    bool schemaOk = true;
    for (auto name : names) if (!tree->GetLeaf(name)) schemaOk = false;
    if (!schemaOk) { failed = true; continue; }
    checked = true;
    auto n = tree->GetEntries();
    if (!n) continue;
    tree->SetEstimate(n + 1);
    if (tree->Draw("fPdgCode", "", "goff") != n) { failed = true; continue; }
    const std::vector<int> codes(tree->GetV1(), tree->GetV1() + n);
    auto* daughters = tree->GetLeaf("fIndexSlice_Daughters");
    for (Long64_t i = 0; i < n; ++i) {
      if (std::abs(codes[i]) != 421) continue;
      if (tree->GetEntry(i) <= 0) { failed = true; continue; }
      auto a = daughters->GetValueLong64(0), b = daughters->GetValueLong64(1);
      if (a < 0) continue;
      if (b < a || b >= n) { failed = true; continue; }
      if (b != a + 1) continue;
      int sign = codes[i] > 0 ? 1 : -1;
      if (!((codes[a] == 211 * sign && codes[b] == -321 * sign) ||
            (codes[b] == 211 * sign && codes[a] == -321 * sign))) continue;
      const double px = tree->GetLeaf("fPx")->GetValue();
      const double py = tree->GetLeaf("fPy")->GetValue();
      const double pz = tree->GetLeaf("fPz")->GetValue();
      // Same fixed D0 mass and default generated-rapidity bound as the task.
      constexpr double massD0 = 1.86484;
      const double y = std::asinh(pz / std::sqrt(massD0 * massD0 + px * px + py * py));
      if (std::abs(y) > 0.5 || !std::isfinite(y)) continue;
      inRapidity = true;
      // Framework/DataTypes.h: MCParticleFlags::FromBackgroundEvent == 0x2.
      const auto flags = static_cast<unsigned int>(tree->GetLeaf("fFlags")->GetValueLong64());
      if ((flags & 0x2u) == 0) inRapidityNotBackground = true;
    }
  }
  auto status = [&](bool present) { return present ? "YES" : failed || !checked ? "CHECK_FAILED" : "NO"; };
  std::cout << "D0_KPI_WITHIN_DEFAULT_RAPIDITY: " << status(inRapidity) << '\n';
  std::cout << "D0_KPI_WITHIN_RAPIDITY_NOT_BACKGROUND: " << status(inRapidityNotBackground) << '\n';
}
