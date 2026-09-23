#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h" // Provides the workflow hooks
#include <TH1F.h>
#include <TFile.h>
#include <cmath>

using namespace o2;
using namespace o2::framework;

struct MyDeltaPhiTask {
    TH1F* hDeltaPhi = nullptr;
    TFile* outFile = nullptr;

    // 1. Initialize the file and histogram directly
    void init(InitContext& ic) {
        outFile = new TFile("DeltaPhiResults.root", "RECREATE");
        hDeltaPhi = new TH1F("hDeltaPhi", "Delta-Phi Distribution;#Delta#phi;Counts", 100, -M_PI, M_PI);
    }

    // 2. Process the tracks using C++ iterators instead of array subscripts
    void process(aod::Tracks const& tracks) {
        // Loop over the track table safely
        for (auto it1 = tracks.begin(); it1 != tracks.end(); ++it1) {
            auto it2 = it1;
            ++it2; // Start inner loop one element ahead to avoid self-pairing
            
            for (; it2 != tracks.end(); ++it2) {
                float dphi = it1->phi() - it2->phi();
                
                // Normalize angle to [-pi, pi]
                while (dphi > M_PI) dphi -= 2 * M_PI;
                while (dphi < -M_PI) dphi += 2 * M_PI;
                
                hDeltaPhi->Fill(dphi);
            }
        }
    }

    // 3. Save and close the file when the data stream ends
    void endOfStream(EndOfStreamContext& ec) {
        if (outFile) {
            outFile->cd();
            hDeltaPhi->Write();
            outFile->Close();
        }
    }
};

// 4. Hook the task into the O2 workflow explicitly (replaces WORKFLOW_APP)
WorkflowSpec defineDataProcessing(ConfigContext const& cfg) {
    return WorkflowSpec{
        adaptAnalysisTask<MyDeltaPhiTask>(cfg)
    };
}