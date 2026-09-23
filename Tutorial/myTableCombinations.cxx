#include "Framework/AnalysisTask.h"
#include <TH1F.h>
#include <cmath>

using namespace o2;
using namespace o2::framework;

struct MyDeltaPhiTask {
    // ALICE O2 macro to automatically register and save the histogram to AnalysisResults.root
    HISTO1D(hDeltaPhi, "hDeltaPhi", "Delta-Phi Distribution;#Delta#phi;Counts", 100, -M_PI, M_PI);

    // The framework automatically feeds track data from the AO2D file into this function
    void process(aod::Tracks const& tracks) {
        
        // Loop over all track combinations in the event
        for (int i = 0; i < tracks.size(); ++i) {
            for (int j = i + 1; j < tracks.size(); ++j) {
                
                // Calculate delta-phi
                float dphi = tracks[i].phi() - tracks[j].phi();
                
                // Normalize the angle to the [-pi, pi] range
                while (dphi > M_PI) dphi -= 2 * M_PI;
                while (dphi < -M_PI) dphi += 2 * M_PI;
                
                hDeltaPhi->Fill(dphi);
            }
        }
    }
};

// Generates the main() function and hooks your task into the O2 data pipeline
WORKFLOW_APP(
    o2::framework::AnalysisTask<MyDeltaPhiTask>
);