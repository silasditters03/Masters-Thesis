#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <TH1F.h>
#include <TFile.h>
#include <cmath>

using namespace o2;
using namespace o2::framework;

struct MyDeltaPhiTask {
    TH1F* hDeltaPhi = nullptr;
    TFile* outFile = nullptr;

    void init(InitContext& ic) {
        outFile = new TFile("DeltaPhiResults.root", "RECREATE");
        hDeltaPhi = new TH1F("hDeltaPhi", "Delta-Phi Distribution;#Delta#phi;Counts", 100, -M_PI, M_PI);
    }

    void process(aod::Tracks const& tracks) {
        for (auto it1 = tracks.begin(); it1 != tracks.end(); ++it1) {
            auto it2 = it1;
            ++it2; 
            
            for (; it2 != tracks.end(); ++it2) {
                // Use dot operator for O2 table iterators
                float dphi = it1.phi() - it2.phi();
                
                while (dphi > M_PI) dphi -= 2 * M_PI;
                while (dphi < -M_PI) dphi += 2 * M_PI;
                
                hDeltaPhi->Fill(dphi);
            }
        }
    }

    void endOfStream(EndOfStreamContext& ec) {
        if (outFile) {
            outFile->cd();
            hDeltaPhi->Write();
            outFile->Close();
        }
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg) {
    return WorkflowSpec{
        adaptAnalysisTask<MyDeltaPhiTask>(cfg)
    };
}