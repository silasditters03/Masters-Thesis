#include "Framework/AnalysisTask.h"
#include "Framework/Core.h"
#include "Framework/Logger.h" // Required for the LOGF print statements

using namespace o2;
using namespace o2::framework;

// The task structure
struct TracksPerCollision {
  
  // The framework automatically feeds this function one collision at a time,
  // along with ONLY the tracks that belong to it.
  void process(aod::Collision const& collision, aod::Tracks const& tracks) {
    
    // Print info about the collision itself
    LOGF(info, "Collision time: %f", collision.collisionTime());
    LOGF(info, "Number of tracks in this collision: %d", tracks.size());
    
    // Loop over the tracks for this specific collision
    for (auto& track : tracks) {
      // We add a quick pT filter just so we don't spam your terminal with millions of lines
      if (track.pt() > 5.0f) {
         LOGF(info, "  -> Found High-pT Track: pT = %f GeV/c", track.pt());
      }
    }
  }
};