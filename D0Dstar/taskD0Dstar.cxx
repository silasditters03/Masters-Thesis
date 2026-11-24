// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskD0Dstar.cxx
/// \brief D0 analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Minjung Kim <minjung.kim@cern.ch>, CERN
/// \author Silas Ditters

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsUpcHf.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGUD/Core/UPCHelpers.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>

#include <Rtypes.h>

#include <algorithm> // std::min
#include <array>
#include <cstdint>
#include <numeric>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>
#include <Math/Vector4D.h>
#include <set>
#include <utility>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::hf_evsel;
using namespace o2::analysis::hf_upc;

/// D0 analysis task
namespace
{
enum CandTypeSel {
  SigD0 = 0,      // Signal D0
  SigD0bar,       // Signal D0bar
  ReflectedD0,    // Reflected D0
  ReflectedD0bar, // Reflected D0bar
  PureSigD0,      // Signal D0 exclude Reflected D0bar
  PureSigD0bar    // Signal D0bar exclude Reflected D0
};
} // namespace
struct HfTaskD0Dstar {
  enum RecoEventStatus {
    NoRecoCollision = 0,
    NoSelectedRecoCollision,
    SelectedRecoCollision
  };

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Final selection flag for D* candidates"};
  Configurable<bool> selectionFlagHfD0ToPiK{"selectionFlagHfD0ToPiK", true, "D0 decay-type flag for D* candidates"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<double> yTccRecoMax{"yTccRecoMax", 0.8, "Maximum absolute reconstructed D0-D* pair rapidity; negative disables cut"};
  Configurable<double> yTccGenMax{"yTccGenMax", 0.8, "Maximum absolute generated D0-Dstar pair rapidity; negative disables cut"};
  Configurable<double> etaGenMesonMax{"etaGenMesonMax", 0.8, "Maximum absolute eta of each generated D0 and Dstar; negative disables cut"};
  Configurable<int> selectionFlagHf{"selectionFlagHf", 1, "Selection Flag for HF flagged candidates"};
  Configurable<int> selectionTopol{"selectionTopol", 1, "Selection Flag for topologically selected candidates"};
  Configurable<int> selectionCand{"selectionCand", 1, "Selection Flag for conj. topol. selected candidates"};
  Configurable<int> selectionPid{"selectionPid", 1, "Selection Flag for reco PID candidates"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};
  Configurable<int> centEstimator{"centEstimator", 0, "Centrality estimation (None: 0, FT0C: 2, FT0M: 3)"};
  Configurable<bool> fillEventStatus{"fillEventStatus", false, "Use generated centrality in hSparseAcc and append reco-event status (requires storeCentrality)"};
  Configurable<int> occEstimator{"occEstimator", 0, "Occupancy estimation (None: 0, ITS: 1, FT0C: 2)"};
  Configurable<bool> storeCentrality{"storeCentrality", false, "Flag to store centrality information"};
  Configurable<bool> storeOccupancyAndIR{"storeOccupancyAndIR", false, "Flag to store occupancy information and interaction rate"};
  Configurable<bool> storeTrackQuality{"storeTrackQuality", false, "Flag to store track quality information"};
  Configurable<bool> storeZdcEnergy{"storeZdcEnergy", false, "Flag to store ZDC energy info"};
  Configurable<bool> storeZdcTime{"storeZdcTime", false, "Flag to store ZDC time info"};
  Configurable<bool> storeEtaProngs{"storeEtaProngs", false, "Flag to store eta of the prongs"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> irSource{"irSource", "ZNC hadronic", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  HfEventSelection hfEvSel;         // event selection and monitoring
  HfUpcGapThresholds upcThresholds; // UPC gap determination thresholds
  ctpRateFetcher mRateFetcher;

  SliceCache cache;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  using D0Candidates = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
  using D0CandidatesMc = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>;
  using D0CandidatesKF = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngKF>>;
  using D0CandidatesMcKF = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngKF, aod::HfCand2ProngMcRec>>;

  using D0CandidatesMl = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
  using D0CandidatesMlMc = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0, aod::HfCand2ProngMcRec>>;
  using D0CandidatesMlKF = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0, aod::HfCand2ProngKF>>;
  using D0CandidatesMlMcKF = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0, aod::HfCand2ProngKF, aod::HfCand2ProngMcRec>>;

  // Reconstructed D* candidates with their D0 daughters and selection flags.
  using DstarCandidates = soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi>;
  // Reconstructed candidates with MC matching information.
  using DstarCandidatesMc = soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfCandDstarMcRec>;
  // Reconstructed candidates with ML scores.
  using DstarCandidatesMl = soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi>;
  // Reconstructed candidates with ML scores and MC matching information.
  using DstarCandidatesMlMc = soa::Join<aod::HfD0FromDstar, aod::HfCandDstars, aod::HfSelDstarToD0Pi, aod::HfMlDstarToD0Pi, aod::HfCandDstarMcRec>;
  // Generated MC particles with D* decay matching information.
  using DstarParticlesMc = soa::Join<aod::McParticles,aod::HfCandDstarMcGen>;
  
  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>;
  using CollisionsCent = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms, aod::CentFT0Cs>;
  using CollisionsWithMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::PVMults>;
  using CollisionsWithMcLabelsCent = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::PVMults, aod::CentFT0Ms, aod::CentFT0Cs>;
  using McCollisionsWithCentrality = soa::Join<aod::McCollisions, aod::McCentFT0Cs, aod::McCentFT0Ms>;
  using TracksSelQuality = soa::Join<aod::TracksExtra, aod::TracksWMc>;
  using TracksWPid = soa::Join<o2::aod::FullTracks, aod::TracksDCA, o2::aod::TrackSelection, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  // using TracksWithExtra = o2::soa::Join<o2::aod::FullTracks, o2::aod::TrackSelection>;
  Filter filterD0Flag = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_2prong::DecayType::D0ToPiK))) != static_cast<uint8_t>(0);
  Preslice<TracksWPid> perCol = aod::track::collisionId;
  Preslice<aod::HfCand2Prong> candD0PerCollision = aod::hf_cand::collisionId; //D0 candidate per collision

  PresliceUnsorted<CollisionsWithMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<CollisionsWithMcLabelsCent> colPerMcCollisionCent = aod::mccollisionlabel::mcCollisionId;

  Partition<D0Candidates> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<D0CandidatesKF> selectedD0CandidatesKF = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<D0CandidatesMc> selectedD0CandidatesMc = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;
  Partition<D0CandidatesMcKF> selectedD0CandidatesMcKF = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;

  Partition<D0CandidatesMl> selectedD0CandidatesMl = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<D0CandidatesMlKF> selectedD0CandidatesMlKF = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<D0CandidatesMlMc> selectedD0CandidatesMlMc = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;
  Partition<D0CandidatesMlMcKF> selectedD0CandidatesMlMcKF = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;

  // Reconstructed data: candidates passing the final D* selection.
  Partition<DstarCandidates> selectedDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi;
  Partition<DstarCandidatesMl> selectedDstarCandidatesMl = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi;
  // Reconstructed MC: retain candidates passing the D0 decay-type flag,
  // allowing the later analysis to study each selection stage.
  Partition<DstarCandidatesMc> selectedDstarCandidatesMc = aod::hf_sel_candidate_dstar::isRecoD0Flag == selectionFlagHfD0ToPiK;
  Partition<DstarCandidatesMlMc> selectedDstarCandidatesMlMc = aod::hf_sel_candidate_dstar::isRecoD0Flag == selectionFlagHfD0ToPiK;

  // ThnSparse for ML outputScores and Vars
  ConfigurableAxis thnConfigAxisBkgScore{"thnConfigAxisBkgScore", {50, 0, 1}, "Bkg score bins"};
  ConfigurableAxis thnConfigAxisNonPromptScore{"thnConfigAxisNonPromptScore", {50, 0, 1}, "Non-prompt score bins"};
  ConfigurableAxis thnConfigAxisPromptScore{"thnConfigAxisPromptScore", {50, 0, 1}, "Prompt score bins"};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {120, 1.5848, 2.1848}, "Cand. inv-mass bins"};
  ConfigurableAxis thnConfigAxisPtB{"thnConfigAxisPtB", {1000, 0, 100}, "Cand. beauty mother pTB bins"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {500, 0, 50}, "Cand. pT bins"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {20, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {3, -0.5, 2.5}, "Cand. origin type"};
  ConfigurableAxis thnConfigAxisCandType{"thnConfigAxisCandType", {6, -0.5, 5.5}, "D0 type"};
  ConfigurableAxis thnConfigAxisGenPtD{"thnConfigAxisGenPtD", {500, 0, 50}, "Gen Pt D"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {1000, 0, 100}, "Gen Pt B"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {110, 0., 110.}, ""};
  ConfigurableAxis thnConfigAxisOccupancy{"thnConfigAxisOccupancy", {14, 0, 14000}, "axis for centrality"};
  ConfigurableAxis thnConfigAxisMinItsNCls{"thnConfigAxisMinItsNCls", {5, 3, 8}, "axis for minimum ITS NCls of candidate prongs"};
  ConfigurableAxis thnConfigAxisMinTpcNCrossedRows{"thnConfigAxisMinTpcNCrossedRows", {10, 70, 180}, "axis for minimum TPC NCls crossed rows of candidate prongs"};
  ConfigurableAxis thnConfigAxisIR{"thnConfigAxisIR", {5000, 0, 500}, "Interaction rate (kHz)"};
  ConfigurableAxis thnConfigAxisGapType{"thnConfigAxisGapType", {7, -1.5, 5.5}, "axis for UPC gap type (see TrueGap enum in o2::aod::sgselector)"};
  ConfigurableAxis thnConfigAxisFT0{"thnConfigAxisFT0", {1001, -1.5, 999.5}, "axis for FT0 amplitude (a.u.)"};
  ConfigurableAxis thnConfigAxisFV0A{"thnConfigAxisFV0A", {2001, -1.5, 1999.5}, "axis for FV0-A amplitude (a.u.)"};
  ConfigurableAxis thnConfigAxisFDD{"thnConfigAxisFDD", {200, 0., 4000.}, "axis for FDD amplitude (a.u.)"};
  ConfigurableAxis thnConfigAxisZN{"thnConfigAxisZN", {510, -1.5, 49.5}, "axis for ZN energy (a.u.)"};
  ConfigurableAxis thnConfigAxisTimeZN{"thnConfigAxisTimeZN", {700, -35., 35.}, "axis for ZN energy (a.u.)"};
  ConfigurableAxis thnConfigAxisEtaProng{"thnConfigAxisEtaProng", {200, -1., 1.}, "axis for eta of D0 prongs"};

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSig", "2-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSigPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSigNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecBg", "2-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hYGenPrompt", "MC particles (matched, prompt);#it{y}^{gen.};entries", {HistType::kTH1F, {{300, -1.5, 1.5}}}},
     {"hYGenNonPrompt", "MC particles (matched, non-prompt);#it{y}^{gen.};entries", {HistType::kTH1F, {{300, -1.5, 1.5}}}},
     {"hPtGenSig", "2-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hCPARecSig", "2-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hCPARecBg", "2-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hEtaRecSig", "2-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hEtaRecBg", "2-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hPtGenVsPtRecSig", "2-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH2F, {{360, 0., 36.}, {360, 0., 36.}}}},
     {"hYGenVsYRecSig", "2-prong candidates (matched);#it{y}^{gen.} ;#it{y}^{rec.} ;entries", {HistType::kTH2F, {{300, -1.5, 1.5}, {300, -1.5, 1.5}}}},
     {"hPtProng0Sig", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hPtProng1Sig", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hDecLengthSig", "2-prong candidates (matched);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hDecLengthXYSig", "2-prong candidates (matched);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthSig", "2-prong candidates (matched);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthXYSig", "2-prong candidates (matched);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hd0Prong0Sig", "2-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0Prong1Sig", "2-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0d0Sig", "2-prong candidates (matched);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}}},
     {"hCTSSig", "2-prong candidates (matched);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCtSig", "2-prong candidates (matched);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}}},
     {"hCPASig", "2-prong candidates (matched);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCPAxySig", "2-prong candidates (matched);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hPtProng0Bkg", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hPtProng1Bkg", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hDecLengthBkg", "2-prong candidates (checked);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hDecLengthXYBkg", "2-prong candidates (checked);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthBkg", "2-prong candidates (checked);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthXYBkg", "2-prong candidates (checked);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hd0Prong0Bkg", "2-prong candidates (checked);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0Prong1Bkg", "2-prong candidates (checked);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0d0Bkg", "2-prong candidates (checked);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}}},
     {"hCTSBkg", "2-prong candidates (checked);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCtBkg", "2-prong candidates (checked);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}}},
     {"hCPABkg", "2-prong candidates (checked);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCPAxyBkg", "2-prong candidates (checked);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hPtVsYRecSig_RecoPID", "2-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoPID", "2-prong candidates (RecoPID - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoPID", "2-prong candidates (RecoPID - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoCand", "2-prong candidates (RecoCand - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoCand", "2-prong candidates (RecoCand - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoCand", "2-prong candidates (RecoCand - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoTopol", "2-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoTopol", "2-prong candidates (RecoTopol - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoTopol", "2-prong candidates (RecoTopol - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoHFFlag", "2-prong candidates (RecoHFFlag - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoHFFlag", "2-prong candidates (RecoHFFlag - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoHFFlag", "2-prong candidates (RecoHFFlag - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGen", "2-prong candidates (matched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGenPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGenNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hMassVsPtGenVsPtRecSig", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{gen.} (GeV/#it{c});#it{p}_{T}^{rec.} (GeV/#it{c})", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {150, 0., 30.}}}},
     {"hMassSigD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassBkgD0", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassReflBkgD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigBkgD0", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassBkgD0bar", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassReflBkgD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigBkgD0bar", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}}}};

  
  void init(InitContext& ic)
  {
    // LOGP(info, "Registry: DCAFitterN={} Cent={}", ic.options().get<bool>("processDataWithDCAFitterN"), ic.options().get<bool>("processDataWithDCAFitterNCent"));
    // LOGP(info, "Members: DCAFitterN={} Cent={}", static_cast<bool>(doprocessDataWithDCAFitterN), static_cast<bool>(doprocessDataWithDCAFitterNCent));

    // // Keep your existing initialization below.
    // LOGP(info, "Active switches: DCAFitterN={}, DCAFitterNCent={}",
    //  static_cast<bool>(doprocessDataWithDCAFitterN),
    //  static_cast<bool>(doprocessDataWithDCAFitterNCent));
    std::array<bool, 14> doprocess{doprocessDataWithDCAFitterN, doprocessDataWithDCAFitterNCent, doprocessDataWithKFParticle, doprocessMcWithDCAFitterN, doprocessMcWithDCAFitterNCent, doprocessMcWithKFParticle, doprocessDataWithDCAFitterNMl, doprocessDataWithDCAFitterNMlCent, doprocessDataWithKFParticleMl, doprocessMcWithDCAFitterNMl, doprocessMcWithDCAFitterNMlCent, doprocessMcWithKFParticleMl, doprocessDataWithDCAFitterNWithUpc, doprocessDataWithDCAFitterNMlWithUpc};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) == 0) {
      LOGP(fatal, "At least one process function should be enabled at a time.");
    }
    if ((doprocessDataWithDCAFitterN || doprocessDataWithDCAFitterNCent || doprocessMcWithDCAFitterN || doprocessMcWithDCAFitterNCent || doprocessDataWithDCAFitterNMl || doprocessDataWithDCAFitterNMlCent || doprocessMcWithDCAFitterNMl || doprocessMcWithDCAFitterNMlCent) && (doprocessDataWithKFParticle || doprocessMcWithKFParticle || doprocessDataWithKFParticleMl || doprocessMcWithKFParticleMl)) {
      LOGP(fatal, "DCAFitterN and KFParticle can not be enabled at a time.");
    }
    if ((storeCentrality || storeOccupancyAndIR) && !(doprocessDataWithDCAFitterNCent || doprocessMcWithDCAFitterNCent || doprocessDataWithDCAFitterNMlCent || doprocessMcWithDCAFitterNMlCent || doprocessDataWithDCAFitterNWithUpc || doprocessDataWithDCAFitterNMlWithUpc)) {
      LOGP(fatal, "Can't enable the storeCentrality and storeOccupancy without cent process or UPC process");
    }
    if ((storeZdcEnergy || storeZdcTime) && !(doprocessDataWithDCAFitterNWithUpc || doprocessDataWithDCAFitterNMlWithUpc)) {
      LOGP(fatal, "Can't enable the storeZdcEnergy and storeZdcTime without UPC process");
    }
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassVsPhi", "2-prong candidates vs phi;inv. mass (#pi K) (GeV/#it{c}^{2});phi (rad);entries", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {32, 0, o2::constants::math::TwoPI}}});
    registry.add("hDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecLengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0ErrProng0", "2-prong candidates;prong 0 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0ErrProng1", "2-prong candidates;prong 1 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTS", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassFinerBinning", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthFinerBinning", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{400, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxyFinerBinning", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{400, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0FinerBinning", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1FinerBinning", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0FinerBinning", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -0.1, 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTSFinerBinning", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtFinerBinning", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{500, -0., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAFinerBinning", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAXYFinerBinning", "2-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0VsPtSig", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1VsPtSig", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0VsPtSig", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2}) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -0.1, 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTSVsPtSig", "2-prong candidates;cos #it{#theta}* (D^{0}) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAVsPtSig", "2-prong candidates;cosine of pointing angle vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAXYVsPtSig", "2-prong candidates;cosine of pointing angle xy vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecLengthxyVsPtSig", "2-prong candidates;decay length xy (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthVsPtSig", "2-prong candidates;decay length (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxyVsPtSig", "2-prong candidates;decay length xy (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    //New plots Silas
    registry.add("hCollisionsData", "Input collisions in processData;Counter;Collisions", {HistType::kTH1D, {{1, -0.5, 0.5}}});
    registry.add("hCollisionsMcReco", "Input reconstructed collisions in processMc;Counter;Collisions", {HistType::kTH1D, {{1, -0.5, 0.5}}});
    registry.add("hCollisionsMcGen", "Input generated collisions in processMc;Counter;Collisions", {HistType::kTH1D, {{1, -0.5, 0.5}}});
    registry.add("hMassTcc", "Selected D0-D* pairs; #it{M}(D^{0}D*) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{600, 3.6, 4.2}}});
    registry.add("hPtTcc", "Selected D0-D* pairs; #it{p}_{T}(D^{0}D*) (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 100.}}});
    registry.add("hMassVsPtTcc", "Selected D0-D* pairs; #it{M}(D^{0}D*) (GeV/#it{c}^{2});" "#it{p}_{T}(D^{0}D*) (GeV/#it{c})", {HistType::kTH2F, {{600, 3.6, 4.2}, {200, 0., 100.}}});
    registry.add("hTccInputs", "Input counts;counter;count", {HistType::kTH1D, {{3, -0.5, 2.5}}});
    registry.add("hTccDstarFlow", "D* selection for pairing;stage;count", {HistType::kTH1D, {{4, -0.5, 3.5}}});
    registry.add("hTccPairFlow", "D0-D* pairing selections;stage;count", {HistType::kTH1D, {{5, -0.5, 4.5}}});
    registry.add("hMcPairTruthStatus", "Selected reconstructed pairs;constituent matching category;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}});
    registry.add("hMassTccBothMatched", "Pairs with both constituents truth matched;#it{M}(D^{0}D*) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{600, 3.6, 4.2}}});
    registry.add("hMcGenCounts", "Generated particle counts;category;entries", {HistType::kTH1F, {{3, -0.5, 2.5}}});
    registry.add("hMassGenD0DstarParent", "Generated parents decaying directly to D0-D*;parent mass (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{1000, 0., 10.}}});
    registry.add("hPtGenD0DstarParent", "Generated parents decaying directly to D0-D*;parent #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 100.}}});
    registry.add("hTccPairFlowGen", "Generated D0-Dstar pair cut flow;Stage;Pairs", {HistType::kTH1D, {{7, 0.5, 7.5}}});

    auto* genPairFlowAxis = registry.get<TH1>(HIST("hTccPairFlowGen"))->GetXaxis();
    genPairFlowAxis->SetBinLabel(1, "Same collision");
    genPairFlowAxis->SetBinLabel(2, "Both charm combinations");
    genPairFlowAxis->SetBinLabel(3, "Independent D0");
    genPairFlowAxis->SetBinLabel(4, "Decay channels");
    genPairFlowAxis->SetBinLabel(5, "Distinct daughters");
    genPairFlowAxis->SetBinLabel(6, "Rapidity");
    genPairFlowAxis->SetBinLabel(7, "Meson eta");

    registry.add("hMassD0DstarGen", "Selected generated D0-Dstar pairs;#it{M}(D^{0}D*) (GeV/#it{c}^{2});Pairs", {HistType::kTH1D, {{1200, 3.6, 15.6}}});
    registry.add("hPtD0DstarGen", "Selected generated D0-Dstar pairs;#it{p}_{T} (GeV/#it{c});Pairs", {HistType::kTH1D, {{200, 0., 100.}}});
    registry.add("hMassVsPtD0DstarGen", "Selected generated D0-Dstar pairs;#it{M}(D^{0}D*) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2D, {{1200, 3.6, 15.6}, {200, 0., 100.}}});
    registry.add("hPtD0InSelectedGenPairs", "Independent D0 in selected generated pairs;#it{p}_{T}(D^{0}) (GeV/#it{c});Pairs", {HistType::kTH1D, {{200, 0., 100.}}}); 
    registry.add("hPtDstarInSelectedGenPairs", "Dstar in selected generated pairs;#it{p}_{T}(D*) (GeV/#it{c});Pairs", {HistType::kTH1D, {{200, 0., 100.}}});
    registry.add("hDeltaEtaGenPairs", "Selected generated D0-Dstar pairs;#Delta#eta = #eta(D^{0}) - #eta(D*);Pairs", {HistType::kTH1D, {{200, -5., 5.}}});
    registry.add("hDeltaPhiGenPairs", "Selected generated D0-Dstar pairs;#Delta#varphi = #varphi(D^{0}) - #varphi(D*) (rad);Pairs", {HistType::kTH1D, {{64, -0.5 * o2::constants::math::PI, 1.5 * o2::constants::math::PI}}});
    registry.add("hDeltaPhiVsDeltaEtaGenPairs", "Selected generated D0-Dstar pairs;#Delta#varphi (rad); #Delta#eta;", {HistType::kTH2D, {{64, -0.5 * o2::constants::math::PI, 1.5 * o2::constants::math::PI}, {200, -5., 5.}}});

    registry.add("hGenPairUniqueness", "Selected generated pairs: repeated row-pair check;Category;Pairs",{HistType::kTH1D, {{2, 0.5, 2.5}}});

    auto* uniquenessAxis = registry.get<TH1>(HIST("hGenPairUniqueness"))->GetXaxis();
    uniquenessAxis->SetBinLabel(1, "First occurrence");
    uniquenessAxis->SetBinLabel(2, "Repeated pair");
    
    const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (#pi K) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPtB{thnConfigAxisPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisY{thnConfigAxisY, "y"};
    const AxisSpec thnAxisOrigin{thnConfigAxisOrigin, "Origin"};
    const AxisSpec thnAxisCandType{thnConfigAxisCandType, "D0 type"};
    const AxisSpec thnAxisGenPtD{thnConfigAxisGenPtD, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisGenPtB{thnConfigAxisGenPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisNumPvContr{thnConfigAxisNumPvContr, "Number of PV contributors"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality"};
    const AxisSpec thnAxisOccupancy{thnConfigAxisOccupancy, "Occupancy"};
    const AxisSpec thnAxisMinItsNCls{thnConfigAxisMinItsNCls, "Minimum ITS cluster found"};
    const AxisSpec thnAxisMinTpcNCrossedRows{thnConfigAxisMinTpcNCrossedRows, "Minimum TPC crossed rows"};
    const AxisSpec thnAxisIR{thnConfigAxisIR, "Interaction rate"};
    const AxisSpec thnAxisGapType{thnConfigAxisGapType, "Gap type"};
    const AxisSpec thnAxisFT0A{thnConfigAxisFT0, "FT0-A amplitude"};
    const AxisSpec thnAxisFT0C{thnConfigAxisFT0, "FT0-C amplitude"};
    const AxisSpec thnAxisFV0A{thnConfigAxisFV0A, "FV0-A amplitude"};
    const AxisSpec thnAxisFDDA{thnConfigAxisFDD, "FDD-A amplitude"};
    const AxisSpec thnAxisFDDC{thnConfigAxisFDD, "FDD-C amplitude"};
    const AxisSpec thnAxisEnergyZNA{thnConfigAxisZN, "ZNA energy"};
    const AxisSpec thnAxisEnergyZNC{thnConfigAxisZN, "ZNC energy"};
    const AxisSpec thnAxisTimeZNA{thnConfigAxisTimeZN, "ZNA Time"};
    const AxisSpec thnAxisTimeZNC{thnConfigAxisTimeZN, "ZNC Time"};
    const AxisSpec thnAxisEtaProngs{thnConfigAxisEtaProng, "Eta of prongs"};

    if (doprocessMcWithDCAFitterN || doprocessMcWithDCAFitterNCent || doprocessMcWithKFParticle || doprocessMcWithDCAFitterNMl || doprocessMcWithDCAFitterNMlCent || doprocessMcWithKFParticleMl) {
      std::vector<AxisSpec> axesAcc = {thnAxisGenPtD, thnAxisGenPtB, thnAxisY, thnAxisOrigin, thnAxisNumPvContr};

      if (storeCentrality) {
        axesAcc.emplace_back(thnConfigAxisCent, fillEventStatus ? "Generated centrality (%)" : "Centrality");
      }
      // interaction rate only store in Data and MC Reco. Level
      if (storeOccupancyAndIR) {
        axesAcc.push_back(thnAxisOccupancy);
      }

      if (fillEventStatus) {
        if (!storeCentrality) {
          LOGP(fatal, "fillEventStatus requires storeCentrality=true.");
        }
        if (centEstimator != CentralityEstimator::FT0C && centEstimator != CentralityEstimator::FT0M) {
          LOGP(fatal, "centEstimator must be FT0C (2) or FT0M (3).");
        }
        axesAcc.emplace_back(3, -0.5, 2.5, "Reco event status (0: no reco, 1: none selected, 2: selected)");
      }
      registry.add("hSparseAcc", "Thn for generated D0 from charm and beauty", HistType::kTHnSparseD, axesAcc);
      registry.get<THnSparse>(HIST("hSparseAcc"))->Sumw2();
    }

    std::vector<AxisSpec> axes = {thnAxisMass, thnAxisPt, thnAxisY, thnAxisCandType};
    if (doprocessMcWithDCAFitterN || doprocessMcWithDCAFitterNCent || doprocessMcWithKFParticle || doprocessMcWithDCAFitterNMl || doprocessMcWithDCAFitterNMlCent || doprocessMcWithKFParticleMl) {
      axes.push_back(thnAxisPtB);
      axes.push_back(thnAxisOrigin);
    }
    if (storeCentrality) {
      axes.push_back(thnAxisCent);
      axes.push_back(thnAxisNumPvContr);
    }
    if (storeOccupancyAndIR) {
      axes.push_back(thnAxisOccupancy);
      axes.push_back(thnAxisIR);
    }
    if (storeTrackQuality) {
      axes.push_back(thnAxisMinItsNCls);
      axes.push_back(thnAxisMinTpcNCrossedRows);
    }
    if (applyMl) {
      const AxisSpec thnAxisBkgScore{thnConfigAxisBkgScore, "BDT score bkg."};
      const AxisSpec thnAxisNonPromptScore{thnConfigAxisNonPromptScore, "BDT score non-prompt."};
      const AxisSpec thnAxisPromptScore{thnConfigAxisPromptScore, "BDT score prompt."};

      // Insert ML scores after pt (position 2) to match taskDplus structure: [mass, pt, mlScores, ...]
      if (doprocessDataWithDCAFitterNMlWithUpc) {
        axes.insert(axes.begin() + 2, thnAxisNonPromptScore);
        axes.insert(axes.begin() + 2, thnAxisPromptScore);
        axes.insert(axes.begin() + 2, thnAxisBkgScore);
      } else {
        axes.insert(axes.begin(), thnAxisPromptScore);
        axes.insert(axes.begin(), thnAxisNonPromptScore);
        axes.insert(axes.begin(), thnAxisBkgScore);
      }
    }
    if (doprocessDataWithDCAFitterNMlWithUpc || doprocessDataWithDCAFitterNWithUpc) {
      axes.push_back(thnAxisGapType);
      axes.push_back(thnAxisFT0A);
      axes.push_back(thnAxisFT0C);
      axes.push_back(thnAxisFV0A);
      axes.push_back(thnAxisFDDA);
      axes.push_back(thnAxisFDDC);
      axes.push_back(thnAxisNumPvContr);
      if (storeZdcEnergy) {
        axes.push_back(thnAxisEnergyZNA);
        axes.push_back(thnAxisEnergyZNC);
      }
      if (storeZdcTime) {
        axes.push_back(thnAxisTimeZNA);
        axes.push_back(thnAxisTimeZNC);
      }
      if (storeEtaProngs) {
        axes.push_back(thnAxisEtaProngs);
        axes.push_back(thnAxisEtaProngs);
      }
    }

    if (applyMl) {
      registry.add("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type", "Thn for D0 candidates", HistType::kTHnSparseD, axes);
      registry.get<THnSparse>(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"))->Sumw2();
    } else {
      registry.add("hMassVsPtVsPtBVsYVsOriginVsD0Type", "Thn for D0 candidates", HistType::kTHnSparseD, axes);
      registry.get<THnSparse>(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"))->Sumw2();
    }

    registry.add("Data/fitInfo/ampFT0A_vs_ampFT0C", "FT0-A vs FT0-C amplitude;FT0-A amplitude (a.u.);FT0-C amplitude (a.u.)", {HistType::kTH2F, {{2500, 0., 250}, {2500, 0., 250}}});
    registry.add("Data/zdc/energyZNA_vs_energyZNC", "ZNA vs ZNC common energy;E_{ZNA}^{common} (a.u.);E_{ZNC}^{common} (a.u.)", {HistType::kTH2F, {{1000, 0., 20000}, {1000, 0., 20000}}});
    registry.add("Data/zdc/timeZNA_vs_timeZNC", "ZNA vs ZNC time;ZNA Time;ZNC time", {HistType::kTH2F, {{200, -10., 10.}, {200, -10., 10.}}});
    registry.add("Data/hUpcGapAfterSelection", "UPC gap type after selection;Gap type;Counts", {HistType::kTH1F, {{7, -1.5, 5.5}}});
    registry.add("Data/hGapVsEtaTrack0", "UPC gap vs Eta;Gap type;Eta", {HistType::kTH2F, {{7, -1.5, 5.5}, {50, -1., 1.}}});
    registry.add("Data/hGapVsEtaTrack1", "UPC gap vs Eta;Gap type;Eta", {HistType::kTH2F, {{7, -1.5, 5.5}, {50, -1., 1.}}});
    registry.add("QAtracks/hEtaTrackVsGap", "", {HistType::kTH2F, {{7, -1.5, 5.5}, {100, -1.8, 1.8}}});
    registry.add("QAtracks/hPtTrackVsGap", "", {HistType::kTH2F, {{7, -1.5, 5.5}, {100, 0, 50}}});

    registry.add("Data/hTPCnSigProng0Pion_GapA", "Gap A Prong 0;P (GeV/c) ;TPC nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTPCnSigProng1Kaon_GapA", "Gap A Prong 1;P (GeV/c) ;TPC nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTPCnSigProng0Kaon_GapA", "Gap A Prong 0;P (GeV/c) ;TPC nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTPCnSigProng1Pion_GapA", "Gap A Prong 1;P (GeV/c) ;TPC nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTPCnSigProng0Pion_GapC", "Gap C Prong 0;P (GeV/c) ;TPC nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTPCnSigProng1Kaon_GapC", "Gap C Prong 1;P (GeV/c) ;TPC nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTPCnSigProng0Kaon_GapC", "Gap C Prong 0;P (GeV/c) ;TPC nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTPCnSigProng1Pion_GapC", "Gap C Prong 1;P (GeV/c) ;TPC nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});

    registry.add("Data/hTOFnSigProng0Pion_GapA", "Gap A Prong 0;P (GeV/c) ;TOF nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTOFnSigProng1Kaon_GapA", "Gap A Prong 1;P (GeV/c) ;TOF nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTOFnSigProng0Kaon_GapA", "Gap A Prong 0;P (GeV/c) ;TOF nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTOFnSigProng1Pion_GapA", "Gap A Prong 1;P (GeV/c) ;TOF nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTOFnSigProng0Pion_GapC", "Gap C Prong 0;P (GeV/c) ;TOF nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTOFnSigProng1Kaon_GapC", "Gap C Prong 1;P (GeV/c) ;TOF nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTOFnSigProng0Kaon_GapC", "Gap C Prong 0;P (GeV/c) ;TOF nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});
    registry.add("Data/hTOFnSigProng1Pion_GapC", "Gap C Prong 1;P (GeV/c) ;TOF nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {120, -6., 6.}}});

    registry.add("Data/hTpcTofnSigProng0Pion_GapA", "Gap A Prong 0;P (GeV/c) ;TpcTof nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {50, 0., 10.}}});
    registry.add("Data/hTpcTofnSigProng1Kaon_GapA", "Gap A Prong 1;P (GeV/c) ;TpcTof nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {50, 0., 10.}}});
    registry.add("Data/hTpcTofnSigProng0Kaon_GapA", "Gap A Prong 0;P (GeV/c) ;TpcTof nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {50, 0., 10.}}});
    registry.add("Data/hTpcTofnSigProng1Pion_GapA", "Gap A Prong 1;P (GeV/c) ;TpcTof nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {50, 0., 10.}}});
    registry.add("Data/hTpcTofnSigProng0Pion_GapC", "Gap C Prong 0;P (GeV/c) ;TpcTof nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {50, 0., 10.}}});
    registry.add("Data/hTpcTofnSigProng1Kaon_GapC", "Gap C Prong 1;P (GeV/c) ;TpcTof nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {50, 0., 10.}}});
    registry.add("Data/hTpcTofnSigProng0Kaon_GapC", "Gap C Prong 0;P (GeV/c) ;TpcTof nSigma Pion", {HistType::kTH2F, {{100, 0, 50}, {50, 0., 10.}}});
    registry.add("Data/hTpcTofnSigProng1Pion_GapC", "Gap C Prong 1;P (GeV/c) ;TpcTof nSigma Kaon", {HistType::kTH2F, {{100, 0, 50}, {50, 0., 10.}}});
    registry.add("Data/hGapVsRap", "UPC gap vs Eta;Gap type;Eta", {HistType::kTH2F, {{7, -1.5, 5.5}, {50, -1., 1.}}});
    // QA histograms for Event level info
    registry.add("QAevents/hPVcontrVsGap", " ;Gap ; N PV contributors", {HistType::kTH2F, {{7, -1.5, 5.5}, {800, 0., 200.}}});
    registry.add("QAevents/ampFT0AVsC_GapA", "FT0-A vs FT0-C amplitude;FT0-A amplitude (a.u.);FT0-C amplitude (a.u.) Gap A", {HistType::kTH2F, {{2500, 0., 250}, {2500, 0., 250}}});
    registry.add("QAevents/ampFT0AVsC_GapC", "FT0-A vs FT0-C amplitude;FT0-A amplitude (a.u.);FT0-C amplitude (a.u.) Gap C", {HistType::kTH2F, {{2500, 0., 250}, {2500, 0., 250}}});
    registry.add("QAevents/energyZNAvsC_GapA", "ZNA vs ZNC common energy;E_{ZNA}^{common} (a.u.);E_{ZNC}^{common} (a.u.) Gap A", {HistType::kTH2F, {{1000, 0., 20000}, {1000, 0., 20000}}});
    registry.add("QAevents/energyZNAvsC_GapC", "ZNA vs ZNC common energy;E_{ZNA}^{common} (a.u.);E_{ZNC}^{common} (a.u.) Gap A", {HistType::kTH2F, {{1000, 0., 20000}, {1000, 0., 20000}}});
    registry.add("QAevents/timeZNAvsC_GapA", "ZNA vs ZNC time;ZNA Time;ZNC time Gap A", {HistType::kTH2F, {{200, -10., 10.}, {200, -10., 10.}}});
    registry.add("QAevents/timeZNAvsC_GapC", "ZNA vs ZNC time;ZNA Time;ZNC time Gap C", {HistType::kTH2F, {{200, -10., 10.}, {200, -10., 10.}}});
    // QA histograms for tracks
    // A side gap
    registry.add("QAtracks/hTPCnSigmaPi_GapA", "Gap A;P (GeV/c) ;TPC nSigma Pi", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTPCnSigmaKa_GapA", "Gap A;P (GeV/c) ;TPC nSigma Ka", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTPCnSigmaPr_GapA", "Gap A;P (GeV/c) ;TPC nSigma Pr", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTOFnSigmaPi_GapA", "Gap A;P (GeV/c) ;TOF nSigma Pi", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTOFnSigmaKa_GapA", "Gap A;P (GeV/c) ;TOF nSigma Ka", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTOFnSigmaPr_GapA", "Gap A;P (GeV/c) ;TOF nSigma Pr", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTPCTOFnSigmaPi_GapA", "Gap A;P (GeV/c) ;TPCTOF nSigma Pi", {HistType::kTH2F, {{360, 0, 36}, {120, 0., 6.}}});
    registry.add("QAtracks/hTPCTOFnSigmaKa_GapA", "Gap A;P (GeV/c) ;TPCTOF nSigma Ka", {HistType::kTH2F, {{360, 0, 36}, {120, 0., 6.}}});
    registry.add("QAtracks/hTPCTOFnSigmaPr_GapA", "Gap A;P (GeV/c) ;TPCTOF nSigma Pr", {HistType::kTH2F, {{360, 0, 36}, {120, 0., 6.}}});

    registry.add("QAtracks/hTPCNCls_GapA", "Gap A; TPC Cls", {HistType::kTH1F, {{160, 0., 160.}}});
    registry.add("QAtracks/hTPCChi2_GapA", "Gap A; TPC chi2", {HistType::kTH1F, {{24, 0., 6.}}});
    registry.add("QAtracks/hITSNCls_GapA", "Gap A; TPC Cls", {HistType::kTH1F, {{8, -1., 7.}}});
    registry.add("QAtracks/hDCAxy_GapA", "Gap A; DCA xy", {HistType::kTH1F, {{400, -2, 2.}}});
    registry.add("QAtracks/hDCAz_GapA", "Gap A; DCA z", {HistType::kTH1F, {{400, -4, 4.}}});

    // C side gap
    registry.add("QAtracks/hTPCnSigmaPi_GapC", "Gap C;P (GeV/c) ;TPC nSigma Pi", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTPCnSigmaKa_GapC", "Gap C;P (GeV/c) ;TPC nSigma Ka", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTPCnSigmaPr_GapC", "Gap C;P (GeV/c) ;TPC nSigma Pr", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTOFnSigmaPi_GapC", "Gap C;P (GeV/c) ;TOF nSigma Pi", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTOFnSigmaKa_GapC", "Gap C;P (GeV/c) ;TOF nSigma Ka", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTOFnSigmaPr_GapC", "Gap C;P (GeV/c) ;TOF nSigma Pr", {HistType::kTH2F, {{360, 0, 36}, {120, -6., 6.}}});
    registry.add("QAtracks/hTPCTOFnSigmaPi_GapC", "Gap C;P (GeV/c) ;TPCTOF nSigma Pi", {HistType::kTH2F, {{360, 0, 36}, {120, 0, 6.}}});
    registry.add("QAtracks/hTPCTOFnSigmaKa_GapC", "Gap C;P (GeV/c) ;TPCTOF nSigma Ka", {HistType::kTH2F, {{360, 0, 36}, {120, 0., 6.}}});
    registry.add("QAtracks/hTPCTOFnSigmaPr_GapC", "Gap C;P (GeV/c) ;TPCTOF nSigma Pr", {HistType::kTH2F, {{360, 0, 36}, {120, 0., 6.}}});

    registry.add("QAtracks/hTPCNCls_GapC", "Gap C; TPC Cls", {HistType::kTH1F, {{160, 0., 160.}}});
    registry.add("QAtracks/hTPCChi2_GapC", "Gap C; TPC chi2", {HistType::kTH1F, {{24, 0., 6.}}});
    registry.add("QAtracks/hITSNCls_GapC", "Gap C; TPC Cls", {HistType::kTH1F, {{8, -1., 7.}}});
    registry.add("QAtracks/hDCAxy_GapC", "Gap C; DCA xy", {HistType::kTH1F, {{400, -2, 2.}}});
    registry.add("QAtracks/hDCAz_GapC", "Gap C; DCA z", {HistType::kTH1F, {{400, -4, 4.}}});

    if (fillEventStatus && (doprocessMcWithDCAFitterN || doprocessMcWithDCAFitterNCent || doprocessMcWithKFParticle || doprocessMcWithDCAFitterNMl || doprocessMcWithDCAFitterNMlCent || doprocessMcWithKFParticleMl)) {
      hfEvSel.init(registry);
    } else {
      hfEvSel.addHistograms(registry);
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  //Template function: Allows processData to work with different configurations
  template <int ReconstructionType, bool ApplyMl, typename D0CandType, typename DstarCandType, typename CollType, typename BCsType>
  void processData(D0CandType const& d0Candidates,
                   DstarCandType const& dstarCandidates,
                   CollType const& collisions,
                   aod::TracksWExtra const&,
                   BCsType const&)
  {
    registry.fill(HIST("hTccInputs"), 0.);
    registry.fill(HIST("hTccInputs"), 1., static_cast<double>(d0Candidates.size()));
    registry.fill(HIST("hTccInputs"), 2., static_cast<double>(dstarCandidates.size()));
    registry.fill(HIST("hCollisionsData"), 0., static_cast<double>(collisions.size()));
    //Dstar Loop, first because its rarer
    for (const auto& dstar : dstarCandidates) {

      registry.fill(HIST("hTccDstarFlow"), 0);
      //Check if candidate is flagged as correct Dstar, if not: skip
      if (!dstar.isSelDstarToD0Pi()) {
        continue;
      }
      registry.fill(HIST("hTccDstarFlow"), 1);

      //Check if there is a soft pion; if not: skip
      const auto softPiSign = dstar.signSoftPi();
      if (softPiSign == 0){
        continue;
      }
      registry.fill(HIST("hTccDstarFlow"), 2);

      //Get D0 from Dstar's collision
      const auto d0s = d0Candidates.sliceBy(candD0PerCollision, dstar.collisionId());

      if (d0s.size() > 0) {
        registry.fill(HIST("hTccDstarFlow"), 3);
      }

      //Loop over D0
      for (const auto& d0 : d0s) {
        registry.fill(HIST("hTccPairFlow"), 0);
        //Check if D0 satisfies D0 -> K Pi flag
        if (!(d0.hfflag() & (1 << aod::hf_cand_2prong::DecayType::D0ToPiK))) {
          continue;
        }
        registry.fill(HIST("hTccPairFlow"), 1);

        //Checks if Tetraquark has right charge
        const bool isTccPlus = softPiSign > 0 && d0.isSelD0() >= selectionFlagD0;
        const bool isTccMinus = softPiSign < 0 && d0.isSelD0bar() >= selectionFlagD0bar;
        if (!isTccPlus && !isTccMinus) {
          continue;
        }
        registry.fill(HIST("hTccPairFlow"), 2);

        //Checks if no tracks overlap
        if (d0.prong0Id() == dstar.prong0Id() ||
          d0.prong0Id() == dstar.prong1Id() ||
          d0.prong0Id() == dstar.prongPiId() ||
          d0.prong1Id() == dstar.prong0Id() ||
          d0.prong1Id() == dstar.prong1Id() ||
          d0.prong1Id() == dstar.prongPiId()) {
            continue;
          }
          registry.fill(HIST("hTccPairFlow"), 3);
        
        //Get D0 mass
        double massD0;

        //Ternary expression, we know if we have D0 or D0bar from charge of Tcc
        if constexpr (ReconstructionType == aod::hf_cand::VertexerType::KfParticle){
          massD0 = isTccPlus ? d0.kfGeoMassD0() : d0.kfGeoMassD0bar();
        }
        else {
          massD0 = isTccPlus ? HfHelper::invMassD0ToPiK(d0) : HfHelper::invMassD0barToKPi(d0);
        }

        //Get D* mass
        const double massDstar = isTccPlus ? dstar.invMassDstar() : dstar.invMassAntiDstar();

        //Get three momenta
        const auto pD0 = d0.pVector();
        const auto pDstar = dstar.pVector();

        //Get four-momenta
        const ROOT::Math::PxPyPzMVector fourVectorD0(pD0[0], pD0[1], pD0[2], massD0);
        const ROOT::Math::PxPyPzMVector fourVectorDstar(pDstar[0], pDstar[1], pDstar[2], massDstar);

        const auto fourVectorTcc = fourVectorD0 + fourVectorDstar; //Tcc 4 momentum

        const double massTcc = fourVectorTcc.M(); //Tcc mass from 4-vector
        const double ptTcc = fourVectorTcc.Pt(); //Tcc pT from 4-vector
        const double yTcc = fourVectorTcc.Rapidity(); //Tcc rapidity from 4-vector

        //Rapidity cut for TCC
        if (yTccRecoMax >= 0. && std::abs(yTcc) > yTccRecoMax){
          continue;
        }
        registry.fill(HIST("hTccPairFlow"), 4);

        //Fills histograms
        registry.fill(HIST("hMassTcc"), massTcc);
        registry.fill(HIST("hPtTcc"), ptTcc);
        registry.fill(HIST("hMassVsPtTcc"), massTcc, ptTcc);

      }


    }
  }


  template <bool FillMl, typename CollType, typename CandType, typename BCsType>
  void runAnalysisPerCollisionDataWithUpc(CollType const& collisions,
                                          CandType const& candidates,
                                          BCsType const& bcs,
                                          aod::FT0s const& ft0s,
                                          aod::FV0As const& fv0as,
                                          aod::FDDs const& fdds,
                                          TracksWPid const& tracks)
  {
    for (const auto& collision : collisions) {
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMaskWithUpc<true, CentralityEstimator::None, BCsType>(collision, centrality, ccdb, registry, bcs);
      if (rejectionMask != 0) {
        continue;
      }
      const auto& bc = collision.template bc_as<BCsType>();

      // Determine gap type using SGSelector with BC range checking
      const auto gapResult = hf_upc::determineGapType(collision, bcs, upcThresholds);
      const int gap = gapResult.value;
      const auto numPvContributors = collision.numContrib();

      // Use the BC with FIT activity if available from SGSelector
      auto bcForUPC = bc;
      if (gapResult.bc) {
        bcForUPC = *(gapResult.bc);
      }
      // Get FIT information from the UPC BC
      upchelpers::FITInfo fitInfo{};
      udhelpers::getFITinfo(fitInfo, bcForUPC, bcs, ft0s, fv0as, fdds);

      // Get ZDC energies if available (extract once and reuse)
      const bool hasZdc = bcForUPC.has_zdc();
      float zdcEnergyZNA = -1.f;
      float zdcEnergyZNC = -1.f;
      float zdcTimeZNA = -999.f;
      float zdcTimeZNC = -999.f;

      if (hasZdc) {
        const auto& zdc = bcForUPC.zdc();
        zdcEnergyZNA = zdc.energyCommonZNA();
        zdcEnergyZNC = zdc.energyCommonZNC();
        zdcTimeZNA = zdc.timeZNA();
        zdcTimeZNC = zdc.timeZNC();
        registry.fill(HIST("Data/zdc/energyZNA_vs_energyZNC"), zdcEnergyZNA, zdcEnergyZNC);
        registry.fill(HIST("Data/zdc/timeZNA_vs_timeZNC"), zdcTimeZNA, zdcTimeZNC);
      }
      registry.fill(HIST("Data/fitInfo/ampFT0A_vs_ampFT0C"), fitInfo.ampFT0A, fitInfo.ampFT0C);
      registry.fill(HIST("Data/hUpcGapAfterSelection"), gap);
      registry.fill(HIST("QAevents/hPVcontrVsGap"), gap, collision.numContrib());

      if (gap == 0) {
        registry.fill(HIST("QAevents/ampFT0AVsC_GapA"), fitInfo.ampFT0A, fitInfo.ampFT0C);
        registry.fill(HIST("QAevents/energyZNAvsC_GapA"), zdcEnergyZNA, zdcEnergyZNC);
        registry.fill(HIST("QAevents/timeZNAvsC_GapA"), zdcTimeZNA, zdcTimeZNC);
      }

      if (gap == 1) {
        registry.fill(HIST("QAevents/ampFT0AVsC_GapC"), fitInfo.ampFT0A, fitInfo.ampFT0C);
        registry.fill(HIST("QAevents/energyZNAvsC_GapC"), zdcEnergyZNA, zdcEnergyZNC);
        registry.fill(HIST("QAevents/timeZNAvsC_GapC"), zdcTimeZNA, zdcTimeZNC);
      }

      const auto thisCollId = collision.globalIndex();
      const auto& groupedD0Candidates = candidates.sliceBy(candD0PerCollision, thisCollId);

      // Calculate occupancy and interaction rate if needed
      float occ{-1.f};
      float ir{-1.f};
      if (storeOccupancyAndIR && occEstimator != OccupancyEstimator::None) {
        occ = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
        ir = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource, true) * 1.e-3; // kHz
      }

      auto tracksample = tracks.sliceBy(perCol, collision.globalIndex());
      if (collision.numContrib() > 1) {
        for (const auto& track : tracksample) {
          if (!track.isGlobalTrack()) {
            continue;
          }
          registry.fill(HIST("QAtracks/hEtaTrackVsGap"), gap, track.eta());
          registry.fill(HIST("QAtracks/hPtTrackVsGap"), gap, track.pt());

          if (gap == 0) { // QA for Gap A
            registry.fill(HIST("QAtracks/hTPCnSigmaPi_GapA"), track.p(), track.tpcNSigmaPi());
            registry.fill(HIST("QAtracks/hTPCnSigmaKa_GapA"), track.p(), track.tpcNSigmaKa());
            registry.fill(HIST("QAtracks/hTPCnSigmaPr_GapA"), track.p(), track.tpcNSigmaPr());
            registry.fill(HIST("QAtracks/hTOFnSigmaPi_GapA"), track.p(), track.tofNSigmaPi());
            registry.fill(HIST("QAtracks/hTOFnSigmaKa_GapA"), track.p(), track.tofNSigmaKa());
            registry.fill(HIST("QAtracks/hTOFnSigmaPr_GapA"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("QAtracks/hTPCTOFnSigmaPi_GapA"), track.p(), track.tpcTofNSigmaPi());
            registry.fill(HIST("QAtracks/hTPCTOFnSigmaKa_GapA"), track.p(), track.tpcTofNSigmaKa());
            registry.fill(HIST("QAtracks/hTPCTOFnSigmaPr_GapA"), track.p(), track.tpcTofNSigmaPr());
            registry.fill(HIST("QAtracks/hTPCNCls_GapA"), track.tpcNClsFound());
            registry.fill(HIST("QAtracks/hTPCChi2_GapA"), track.tpcChi2NCl());
            registry.fill(HIST("QAtracks/hITSNCls_GapA"), track.itsNCls());
            registry.fill(HIST("QAtracks/hDCAxy_GapA"), track.dcaXY());
            registry.fill(HIST("QAtracks/hDCAz_GapA"), track.dcaZ());
          }
          if (gap == 1) { // QA for Gap C
            registry.fill(HIST("QAtracks/hTPCnSigmaPi_GapC"), track.p(), track.tpcNSigmaPi());
            registry.fill(HIST("QAtracks/hTPCnSigmaKa_GapC"), track.p(), track.tpcNSigmaKa());
            registry.fill(HIST("QAtracks/hTPCnSigmaPr_GapC"), track.p(), track.tpcNSigmaPr());
            registry.fill(HIST("QAtracks/hTOFnSigmaPi_GapC"), track.p(), track.tofNSigmaPi());
            registry.fill(HIST("QAtracks/hTOFnSigmaKa_GapC"), track.p(), track.tofNSigmaKa());
            registry.fill(HIST("QAtracks/hTOFnSigmaPr_GapC"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("QAtracks/hTPCTOFnSigmaPi_GapC"), track.p(), track.tpcTofNSigmaPi());
            registry.fill(HIST("QAtracks/hTPCTOFnSigmaKa_GapC"), track.p(), track.tpcTofNSigmaKa());
            registry.fill(HIST("QAtracks/hTPCTOFnSigmaPr_GapC"), track.p(), track.tpcTofNSigmaPr());
            registry.fill(HIST("QAtracks/hTPCNCls_GapC"), track.tpcNClsFound());
            registry.fill(HIST("QAtracks/hTPCChi2_GapC"), track.tpcChi2NCl());
            registry.fill(HIST("QAtracks/hITSNCls_GapC"), track.itsNCls());
            registry.fill(HIST("QAtracks/hDCAxy_GapC"), track.dcaXY());
            registry.fill(HIST("QAtracks/hDCAz_GapC"), track.dcaZ());
          }
        }
      }

      for (const auto& candidate : groupedD0Candidates) {
        if (yCandRecoMax >= 0. && std::abs(HfHelper::yD0(candidate)) > yCandRecoMax) {
          continue;
        }

        const float massD0 = HfHelper::invMassD0ToPiK(candidate);
        const float massD0bar = HfHelper::invMassD0barToKPi(candidate);
        const auto ptCandidate = candidate.pt();

        auto track0 = candidate.template prong0_as<TracksWPid>();
        auto track1 = candidate.template prong1_as<TracksWPid>();

        registry.fill(HIST("Data/hGapVsEtaTrack0"), gap, track0.eta());
        registry.fill(HIST("Data/hGapVsEtaTrack1"), gap, track1.eta());
        registry.fill(HIST("Data/hGapVsRap"), gap, HfHelper::yD0(candidate));

        if (gap == 0 && candidate.isSelD0() >= selectionFlagD0) { // A side // D0 --> K-Pi+
          registry.fill(HIST("Data/hTPCnSigProng0Pion_GapA"), track0.p(), track0.tpcNSigmaPi());
          registry.fill(HIST("Data/hTPCnSigProng1Kaon_GapA"), track1.p(), track1.tpcNSigmaKa());
          registry.fill(HIST("Data/hTOFnSigProng0Pion_GapA"), track0.p(), track0.tofNSigmaPi());
          registry.fill(HIST("Data/hTOFnSigProng1Kaon_GapA"), track1.p(), track1.tofNSigmaKa());
          registry.fill(HIST("Data/hTpcTofnSigProng0Pion_GapA"), track0.p(), track0.tpcTofNSigmaPi());
          registry.fill(HIST("Data/hTpcTofnSigProng1Kaon_GapA"), track1.p(), track1.tpcTofNSigmaKa());
        }

        if (gap == 0 && candidate.isSelD0bar() >= selectionFlagD0) { // A side // D0-bar --> K+Pi-
          registry.fill(HIST("Data/hTPCnSigProng0Kaon_GapA"), track0.p(), track0.tpcNSigmaKa());
          registry.fill(HIST("Data/hTPCnSigProng1Pion_GapA"), track1.p(), track1.tpcNSigmaPi());
          registry.fill(HIST("Data/hTOFnSigProng0Kaon_GapA"), track0.p(), track0.tofNSigmaKa());
          registry.fill(HIST("Data/hTOFnSigProng1Pion_GapA"), track1.p(), track1.tofNSigmaPi());
          registry.fill(HIST("Data/hTpcTofnSigProng0Kaon_GapA"), track0.p(), track0.tpcTofNSigmaKa());
          registry.fill(HIST("Data/hTpcTofnSigProng1Pion_GapA"), track1.p(), track1.tpcTofNSigmaPi());
        }

        if (gap == 1 && candidate.isSelD0() >= selectionFlagD0) { // C side // D0 --> K-Pi+
          registry.fill(HIST("Data/hTPCnSigProng0Pion_GapC"), track0.p(), track0.tpcNSigmaPi());
          registry.fill(HIST("Data/hTPCnSigProng1Kaon_GapC"), track1.p(), track1.tpcNSigmaKa());
          registry.fill(HIST("Data/hTOFnSigProng0Pion_GapC"), track0.p(), track0.tofNSigmaPi());
          registry.fill(HIST("Data/hTOFnSigProng1Kaon_GapC"), track1.p(), track1.tofNSigmaKa());
          registry.fill(HIST("Data/hTpcTofnSigProng0Pion_GapC"), track0.p(), track0.tpcTofNSigmaPi());
          registry.fill(HIST("Data/hTpcTofnSigProng1Kaon_GapC"), track1.p(), track1.tpcTofNSigmaKa());
        }

        if (gap == 1 && candidate.isSelD0bar() >= selectionFlagD0) { // C side // D0-bar --> K+Pi-
          registry.fill(HIST("Data/hTPCnSigProng0Kaon_GapC"), track0.p(), track0.tpcNSigmaKa());
          registry.fill(HIST("Data/hTPCnSigProng1Pion_GapC"), track1.p(), track1.tpcNSigmaPi());
          registry.fill(HIST("Data/hTOFnSigProng0Kaon_GapC"), track0.p(), track0.tofNSigmaKa());
          registry.fill(HIST("Data/hTOFnSigProng1Pion_GapC"), track1.p(), track1.tofNSigmaPi());
          registry.fill(HIST("Data/hTpcTofnSigProng0Kaon_GapC"), track0.p(), track0.tpcTofNSigmaKa());
          registry.fill(HIST("Data/hTpcTofnSigProng1Pion_GapC"), track1.p(), track1.tpcTofNSigmaPi());
        }

        if (candidate.isSelD0() >= selectionFlagD0) {
          registry.fill(HIST("hMass"), massD0, ptCandidate);
          registry.fill(HIST("hMassFinerBinning"), massD0, ptCandidate);
          registry.fill(HIST("hMassVsPhi"), massD0, ptCandidate, candidate.phi());
        }
        if (candidate.isSelD0bar() >= selectionFlagD0bar) {
          registry.fill(HIST("hMass"), massD0bar, ptCandidate);
          registry.fill(HIST("hMassFinerBinning"), massD0bar, ptCandidate);
          registry.fill(HIST("hMassVsPhi"), massD0bar, ptCandidate, candidate.phi());
        }

        // Fill THnSparse with structure matching histogram axes: [mass, pt, (mlScores if FillMl), rapidity, d0Type, (cent if storeCentrality), (occ, ir if storeOccupancyAndIR), gapType, FT0A, FT0C, FV0A, FDDA, FDDC, ZNA, ZNC]
        auto fillTHnData = [&](float mass, int d0Type) {
          // Pre-calculate vector size to avoid reallocations
          constexpr int NAxesBase = 11;                       // mass, pt, rapidity, d0Type, gapType, FT0A, FT0C, FV0A, FDDA, FDDC, ZNA, ZNC, nPVcontr
          constexpr int NAxesMl = FillMl ? 3 : 0;             // 3 ML scores if FillMl
          int const nAxesCent = storeCentrality ? 1 : 0;      // centrality if storeCentrality
          int const nAxesOccIR = storeOccupancyAndIR ? 2 : 0; // occupancy and IR if storeOccupancyAndIR
          int const nAxesZdcEnergy = storeZdcEnergy ? 2 : 0;  // ZDC energy if storeZdcEnergy
          int const nAxesZdcTime = storeZdcTime ? 2 : 0;      // ZDC time if storeZdctime
          int const nAxesEtaProngs = storeEtaProngs ? 2 : 0;  // Eta distributions of prong
          int const nAxesTotal = NAxesBase + NAxesMl + nAxesCent + nAxesOccIR + nAxesZdcEnergy + nAxesZdcTime + nAxesEtaProngs;

          std::vector<double> valuesToFill;
          valuesToFill.reserve(nAxesTotal);

          // Fill values in order matching histogram axes
          valuesToFill.push_back(static_cast<double>(mass));
          valuesToFill.push_back(static_cast<double>(ptCandidate));
          if constexpr (FillMl) {
            auto const& mlScores = candidate.mlProbD0();
            if (mlScores.size() == NAxesMl) {
              valuesToFill.push_back(mlScores[0]);
              valuesToFill.push_back(mlScores[1]);
              valuesToFill.push_back(mlScores[2]);
            }
          }
          valuesToFill.push_back(static_cast<double>(HfHelper::yD0(candidate)));
          valuesToFill.push_back(static_cast<double>(d0Type));
          if (storeCentrality) {
            valuesToFill.push_back(centrality);
          }
          if (storeOccupancyAndIR) {
            valuesToFill.push_back(occ);
            valuesToFill.push_back(ir);
          }
          valuesToFill.push_back(static_cast<double>(gap));
          valuesToFill.push_back(static_cast<double>(fitInfo.ampFT0A));
          valuesToFill.push_back(static_cast<double>(fitInfo.ampFT0C));
          valuesToFill.push_back(static_cast<double>(fitInfo.ampFV0A));
          valuesToFill.push_back(static_cast<double>(fitInfo.ampFDDA));
          valuesToFill.push_back(static_cast<double>(fitInfo.ampFDDC));
          valuesToFill.push_back(static_cast<double>(numPvContributors));
          if (storeZdcEnergy) {
            valuesToFill.push_back(static_cast<double>(zdcEnergyZNA));
            valuesToFill.push_back(static_cast<double>(zdcEnergyZNC));
          }
          if (storeZdcTime) {
            valuesToFill.push_back(static_cast<double>(zdcTimeZNA));
            valuesToFill.push_back(static_cast<double>(zdcTimeZNC));
          }
          if (storeEtaProngs) {
            valuesToFill.push_back(static_cast<double>(track0.eta()));
            valuesToFill.push_back(static_cast<double>(track1.eta()));
          }
          if constexpr (FillMl) {
            registry.get<THnSparse>(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type"))->Fill(valuesToFill.data());
          } else {
            registry.get<THnSparse>(HIST("hMassVsPtVsPtBVsYVsOriginVsD0Type"))->Fill(valuesToFill.data());
          }
        };

        if (candidate.isSelD0() >= selectionFlagD0) {
          fillTHnData(massD0, SigD0);
          fillTHnData(massD0, candidate.isSelD0bar() ? ReflectedD0 : PureSigD0);
        }
        if (candidate.isSelD0bar() >= selectionFlagD0bar) {
          fillTHnData(massD0bar, SigD0bar);
          fillTHnData(massD0bar, candidate.isSelD0() ? ReflectedD0bar : PureSigD0bar);
        }
      }
    }
  }

  // DCAFitterN, without centrality.
  void processDataWithDCAFitterN(D0Candidates const&, DstarCandidates const& dstarCandidates, Collisions const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter, false>(selectedD0Candidates, dstarCandidates, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processDataWithDCAFitterN, "Process reconstructed D0-D* pairs with DCAFitterN", true);

  // DCAFitterN, with centrality tables.
  void processDataWithDCAFitterNCent(D0Candidates const&, DstarCandidates const& dstarCandidates, CollisionsCent const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter, false>(selectedD0Candidates, dstarCandidates, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processDataWithDCAFitterNCent, "Process reconstructed D0-D* pairs with DCAFitterN and centrality", false);

  // KFParticle D0 reconstruction, without centrality.
  void processDataWithKFParticle(D0CandidatesKF const&, DstarCandidates const& dstarCandidates, Collisions const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::KfParticle, false>(selectedD0CandidatesKF, dstarCandidates, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processDataWithKFParticle, "Process reconstructed D0-D* pairs with KFParticle D0 reconstruction", false);

  // DCAFitterN with D0 ML information, without centrality.
  void processDataWithDCAFitterNMl(D0CandidatesMl const&, DstarCandidates const& dstarCandidates, Collisions const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter, true>(selectedD0CandidatesMl, dstarCandidates, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processDataWithDCAFitterNMl, "Process reconstructed D0-D* pairs with DCAFitterN and D0 ML information", false);

  // DCAFitterN with D0 ML information and centrality.
  void processDataWithDCAFitterNMlCent(D0CandidatesMl const&, DstarCandidates const& dstarCandidates, CollisionsCent const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter, true>(selectedD0CandidatesMl, dstarCandidates, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processDataWithDCAFitterNMlCent, "Process reconstructed D0-D* pairs with DCAFitterN, D0 ML information and centrality", false);

  // KFParticle D0 reconstruction with D0 ML information.
  void processDataWithKFParticleMl(D0CandidatesMlKF const&, DstarCandidates const& dstarCandidates, Collisions const& collisions, aod::TracksWExtra const& tracks, aod::BcFullInfos const& bcs)
  {
    processData<aod::hf_cand::VertexerType::KfParticle, true>(selectedD0CandidatesMlKF, dstarCandidates, collisions, tracks, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processDataWithKFParticleMl, "Process reconstructed D0-D* pairs with KFParticle D0 reconstruction and D0 ML information", false);
  // TODO: add processKFParticleMlCent

  
  //New ProcessMC function for Tcc
  template <int ReconstructionType, bool ApplyML, typename D0CandType, typename DstarCandType, typename CollType, typename BCsType>
  void processMc(D0CandType const& d0Candidates,
                 DstarCandType const& dstarCandidates,
                 aod::McParticles const& mcParticles,
                 TracksSelQuality const& tracks,
                 CollType const& collisions,
                 McCollisionsWithCentrality const& mcCollisions,
                 BCsType const& bcs)
  {
    //Reconstructed level
    registry.fill(HIST("hCollisionsMcReco"), 0., static_cast<double>(collisions.size()));
    registry.fill(HIST("hCollisionsMcGen"), 0., static_cast<double>(mcCollisions.size()));
    for (const auto& dstar : dstarCandidates) {

      //Check Dstar -> D0 pi candidate flag
      if (!dstar.isSelDstarToD0Pi()) {
        continue;
      }

      //Check if there is a proper soft pion
      const auto softPiSign = dstar.signSoftPi();
      if (softPiSign == 0) {
        continue;
      }

      //Get D0 candidates in D*'s collision
      const auto d0s = d0Candidates.sliceBy(candD0PerCollision, dstar.collisionId());

      for (const auto& d0 : d0s) {
        
        //Check if D0 -> K Pi flag is satisfied
        if (!(d0.hfflag() & (1 << aod::hf_cand_2prong::DecayType::D0ToPiK))) {
          continue;
        }

        const bool isTccPlus = softPiSign > 0 && d0.isSelD0() >= selectionFlagD0;
        const bool isTccMinus = softPiSign < 0 && d0.isSelD0bar() >= selectionFlagD0bar;

        //Check if Tcc has right charge
        if (!isTccPlus && !isTccMinus) {
          continue;
        }

        // Reject shared tracks.
        if (d0.prong0Id() == dstar.prong0Id() ||
            d0.prong0Id() == dstar.prong1Id() ||
            d0.prong0Id() == dstar.prongPiId() ||
            d0.prong1Id() == dstar.prong0Id() ||
            d0.prong1Id() == dstar.prong1Id() ||
            d0.prong1Id() == dstar.prongPiId()) {
          continue;
        }

        //Reconstructed constituent masses
        double massD0;

        if constexpr (ReconstructionType == aod::hf_cand::VertexerType::KfParticle) {
          massD0 = isTccPlus ? d0.kfGeoMassD0() : d0.kfGeoMassD0bar();
        }
        else{
          massD0 = isTccPlus ? HfHelper::invMassD0ToPiK(d0) : HfHelper::invMassD0barToKPi(d0);
        }

        const double massDstar = isTccPlus ?  dstar.invMassDstar() : dstar.invMassAntiDstar();
        
        const auto pD0 = d0.pVector();
        const auto pDstar = dstar.pVector();

        const ROOT::Math::PxPyPzMVector fourVectorD0(pD0[0], pD0[1], pD0[2], massD0);
        const ROOT::Math::PxPyPzMVector fourVectorDstar(pDstar[0], pDstar[1], pDstar[2], massDstar);

        const auto fourVectorTcc = fourVectorD0 + fourVectorDstar;

        const double massTcc = fourVectorTcc.M();
        const double ptTcc = fourVectorTcc.Pt();
        const double yTcc = fourVectorTcc.Rapidity();

        if (yTccRecoMax >= 0. && std::abs(yTcc) > yTccRecoMax) {
          continue;
        }

        //Inclusive reconstructed pairs
        registry.fill(HIST("hMassTcc"), massTcc);
        registry.fill(HIST("hPtTcc"), ptTcc);
        registry.fill(HIST("hMassVsPtTcc"), massTcc, ptTcc);

        //MC-specific information starts here
        constexpr int channelD0 = static_cast<int>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
        constexpr int channelDstar = static_cast<int>(o2::hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi);

        const int expectedSign = isTccPlus ? 1 : -1;

        //Signed matching also checks selected particle/antiparticle hypthesis
        const bool matchedD0 = static_cast<int>(d0.flagMcMatchRec()) == expectedSign * channelD0;
        const bool matchedDstar = static_cast<int>(dstar.flagMcMatchRec()) == expectedSign * channelDstar;

        const int truthStatus = static_cast<int>(matchedD0) + 2 * static_cast<int>(matchedDstar);
        registry.fill(HIST("hMcPairTruthStatus"), truthStatus);

        if (matchedD0 && matchedDstar) {
          registry.fill(HIST("hMassTccBothMatched"), massTcc);

          //Get on daugther track from independant D0, and soft pion from D*
          const auto d0Track = d0.template prong0_as<TracksSelQuality>();
          const auto softPiTrack = dstar.template prongPi_as<TracksSelQuality>();

          if (!d0Track.has_mcParticle() || !softPiTrack.has_mcParticle()) {
            continue;
          }

          //Follow the daugther tracks MC ancestry to their meson parents
          //"True" also accepts the charge-conjugate particles.
          const int indexMcD0 = RecoDecay::getMother(mcParticles, d0Track.template mcParticle_as<aod::McParticles>(), o2::constants::physics::Pdg::kD0, true);
          const int indexMcDstar = RecoDecay::getMother(mcParticles, softPiTrack.template mcParticle_as<aod::McParticles>(), o2::constants::physics::Pdg::kDStar, true);

          if (indexMcD0 < 0 || indexMcDstar < 0) {
            continue;
          }

          //getMother returns global indices; rawIteratorAt uses local row indices
          const auto mcD0 = mcParticles.rawIteratorAt(indexMcD0 - mcParticles.offset());
          const auto mcDstar = mcParticles.rawIteratorAt(indexMcDstar - mcParticles.offset());

          //Check generated flavour against reconstructed hypothesis
          if (mcD0.pdgCode() != expectedSign * o2::constants::physics::Pdg::kD0 ||
              mcDstar.pdgCode() != expectedSign * o2::constants::physics::Pdg::kDStar){
                continue;
              }

          //Both candidates must have same generated collision
          if (mcD0.mcCollisionId() != mcDstar.mcCollisionId()) {
            continue;
          }

          //Check if they share the same mother
          //IndexCommonMother, default = -1, means "no found"
          int indexCommonMother = -1;
          for (const auto& motherD0 : mcD0.mothers_as<aod::McParticles>()) {
            for (const auto& motherDstar : mcDstar.mothers_as<aod::McParticles>()) {
              if (motherD0.globalIndex() == motherDstar.globalIndex()) {
                indexCommonMother = motherD0.globalIndex();
                break;
              }
            }
            if (indexCommonMother >= 0) {
              break;
            }
          }
          if (indexCommonMother < 0) {
            continue;
          }

          //Common mother of D0 and D*
          const auto commonMother = mcParticles.rawIteratorAt(indexCommonMother - mcParticles.offset());

          // Temporary diagnostic: identify what actually produced the pair.
          LOGP(info, "Matched D0-Dstar pair: common mother PDG={}, pT={}", commonMother.pdgCode(), commonMother.pt());
        }
      }
    }

    //Generated level
    for (const auto& particle : mcParticles) {
      const int absPdg = std::abs(particle.pdgCode());

      if (absPdg == o2::constants::physics::Pdg::kD0) {
        registry.fill(HIST("hMcGenCounts"), 0.);
      }

      if (absPdg == o2::constants::physics::Pdg::kDStar) {
        registry.fill(HIST("hMcGenCounts"), 1.);
      }

      if (!particle.has_daughters()) {
        continue;
      }

      //Check if MC particle has parent -> D0 + D* decay
      const bool hasD0DstarDecay = RecoDecay::isMatchedMCGen(mcParticles, particle, absPdg, std::array<int,2>{o2::constants::physics::Pdg::kD0, o2::constants::physics::Pdg::kDStar}, true, nullptr, 1);

      if (!hasD0DstarDecay) {
        continue;
      }

      registry.fill(HIST("hMcGenCounts"), 2.);

      //Get parents 4-momentum
      const ROOT::Math::PxPyPzEVector fourVectorParent(particle.px(), particle.py(), particle.pz(), particle.e());

      registry.fill(HIST("hMassGenD0DstarParent"), fourVectorParent.M());
      registry.fill(HIST("hPtGenD0DstarParent"), fourVectorParent.Pt());

       // Temporary diagnostic to establish the parent's identity.
      LOGP(info, "Generated D0-Dstar parent: PDG={}, mass={}, pT={}", particle.pdgCode(), fourVectorParent.M(), fourVectorParent.Pt());  
    }
     
    // Local to this processMc call, shared across all generated-pair iterations.
    std::set<std::pair<int64_t, int64_t>> seenSelectedGenPairs;
    //Generated D0 D* pairs without common parent

    std::unordered_map<int64_t, std::vector<int64_t>> genD0IdsByCollision;

    for (const auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != o2::constants::physics::Pdg::kD0) {
        continue;
      }

      if (particle.mcCollisionId() < 0) {
        continue;
      }

      genD0IdsByCollision[particle.mcCollisionId()].push_back(particle.globalIndex());
    }

    //Loop over generated D*
    for (const auto& genDstar : mcParticles) {
      if (std::abs(genDstar.pdgCode()) != o2::constants::physics::Pdg::kDStar) {
        continue;
      }

      const auto collisionIt = genD0IdsByCollision.find(genDstar.mcCollisionId());

      if (collisionIt == genD0IdsByCollision.end()) {
        continue;
      }

      const int expectedSign = genDstar.pdgCode() > 0 ? 1 : -1;

      //Check Dstar -> D0 pi and get its daugthers
      std::vector<int> dstarDaughterIds;

      const bool matchesDstarDecay = RecoDecay::isMatchedMCGen(mcParticles, genDstar, o2::constants::physics::Pdg::kDStar,
                    std::array<int, 2>{o2::constants::physics::Pdg::kD0, +kPiPlus}, true, nullptr, 1, &dstarDaughterIds);

      int internalD0Id = -1;
      int softPiId = -1;

      if (matchesDstarDecay) {
        for (const int daughterId : dstarDaughterIds) {
          const auto daughter = mcParticles.rawIteratorAt(daughterId - mcParticles.offset());

          if (std::abs(daughter.pdgCode()) == o2::constants::physics::Pdg::kD0) {
            internalD0Id = daughterId;
          } else if (std::abs(daughter.pdgCode()) == kPiPlus) {
            softPiId = daughterId;
          }
        }
      }

      //Check the Dstar's internal D0 -> pi K decay
      std::vector<int> internalD0DaughterIds;
      bool matchesInternalD0Decay = false;

      if (matchesDstarDecay && internalD0Id >= 0 && softPiId >= 0) {
        const auto internalD0 = mcParticles.rawIteratorAt(internalD0Id - mcParticles.offset());

        matchesInternalD0Decay = RecoDecay::isMatchedMCGen(mcParticles, internalD0, o2::constants::physics::Pdg::kD0,
                                    std::array<int, 2>{+kPiPlus, -kKPlus}, true, nullptr, 1, &internalD0DaughterIds);
        }

        //Pair with every generated D0 in the same MC collision
        for (const int64_t genD0Id : collisionIt->second) {
          const auto genD0 = mcParticles.rawIteratorAt(genD0Id - mcParticles.offset());

          //Stage 1: same collision mesons
          registry.fill(HIST("hTccPairFlowGen"), 1.);

          // Stage 2: accept both D0 and anti-D0 with either Dstar charge.
          registry.fill(HIST("hTccPairFlowGen"), 2.);

          //Stage 3: exclude Dstar's own direct D0 daugther
          bool isDstarDaughter = false;

          for (const auto motherId : genD0.mothersIds()) {
            if (motherId == genDstar.globalIndex()) {
              isDstarDaughter = true;
              break;
            }
          }

          if (isDstarDaughter) {
            continue;
          }

          registry.fill(HIST("hTccPairFlowGen"), 3.);

          //Stage 4: require full decay chain
          if (!matchesDstarDecay || !matchesInternalD0Decay) {
            continue;
          }

          std::vector<int> independentD0DaughterIds;

          const bool matchesIndependentD0Decay = RecoDecay::isMatchedMCGen(mcParticles, genD0, o2::constants::physics::Pdg::kD0, std::array<int, 2>{+kPiPlus, -kKPlus}, true, nullptr, 1, &independentD0DaughterIds);

          if (!matchesIndependentD0Decay) {
          continue;
          }

          registry.fill(HIST("hTccPairFlowGen"), 4.);

          //Stage 5: check no shared-track selection
          std::array<int, 5> finalDaughterIds{
            independentD0DaughterIds[0],
            independentD0DaughterIds[1],
            internalD0DaughterIds[0],
            internalD0DaughterIds[1],
            softPiId};

          std::sort(finalDaughterIds.begin(), finalDaughterIds.end());

          if (std::adjacent_find(finalDaughterIds.begin(), finalDaughterIds.end()) != finalDaughterIds.end()) {
            continue;
          }

          registry.fill(HIST("hTccPairFlowGen"), 5.);

          //construct the pair from generated momenta and energies
          const ROOT::Math::PxPyPzEVector fourVectorGenD0(
            genD0.px(), genD0.py(), genD0.pz(), genD0.e());

          const ROOT::Math::PxPyPzEVector fourVectorGenDstar(
            genDstar.px(), genDstar.py(), genDstar.pz(), genDstar.e());

          const auto fourVectorGenPair = fourVectorGenD0 + fourVectorGenDstar;

          //Stage 6: generated pair rapidity selection
          if (yTccGenMax >= 0. && std::abs(fourVectorGenPair.Rapidity()) > yTccGenMax) {
            continue;
          }

          registry.fill(HIST("hTccPairFlowGen"), 6.);

          // Stage 7: both generated mesons inside the eta acceptance.
          const double etaD0 = genD0.eta();
          const double etaDstar = genDstar.eta();

          if (!std::isfinite(etaD0) || !std::isfinite(etaDstar)) {
            continue;
          }

          if (etaGenMesonMax >= 0. && (std::abs(etaD0) >= etaGenMesonMax || std::abs(etaDstar) >= etaGenMesonMax)) {
            continue;
          }

          registry.fill(HIST("hTccPairFlowGen"), 7.);

          // Diagnose repeated fills of the exact same generated meson pair.
          const std::pair<int64_t, int64_t> pairKey{genD0.globalIndex(), genDstar.globalIndex()};

          const bool firstOccurrence = seenSelectedGenPairs.insert(pairKey).second;

          registry.fill(HIST("hGenPairUniqueness"), firstOccurrence ? 1. : 2.);
          

          const double deltaEtaGen = genD0.eta() - genDstar.eta();
          const double deltaPhiGen = RecoDecay::constrainAngle(static_cast<double>(genD0.phi()) - genDstar.phi(), -0.5 * o2::constants::math::PI);

          //D0 and D* momenta and pair observables
          registry.fill(HIST("hPtD0InSelectedGenPairs"), genD0.pt());
          registry.fill(HIST("hPtDstarInSelectedGenPairs"), genDstar.pt());
          registry.fill(HIST("hMassD0DstarGen"), fourVectorGenPair.M());
          registry.fill(HIST("hPtD0DstarGen"), fourVectorGenPair.Pt());
          registry.fill(HIST("hMassVsPtD0DstarGen"), fourVectorGenPair.M(), fourVectorGenPair.Pt());
          registry.fill(HIST("hDeltaEtaGenPairs"), deltaEtaGen);
          registry.fill(HIST("hDeltaPhiGenPairs"), deltaPhiGen);
          registry.fill(HIST("hDeltaPhiVsDeltaEtaGenPairs"), deltaPhiGen, deltaEtaGen);
      }
    }
  }


    void processMcWithDCAFitterN(D0CandidatesMc const&,
                              DstarCandidatesMc const& dstarCandidates,
                              aod::McParticles const& mcParticles,
                              TracksSelQuality const& tracks,
                              CollisionsWithMcLabels const& collisions,
                              McCollisionsWithCentrality const& mcCollisions,
                              aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, false>(selectedD0CandidatesMc, dstarCandidates, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processMcWithDCAFitterN, "Process MC D0-Dstar pairs with DCAFitterN", false);

  void processMcWithDCAFitterNCent(D0CandidatesMc const&,
                                  DstarCandidatesMc const& dstarCandidates,
                                  aod::McParticles const& mcParticles,
                                  TracksSelQuality const& tracks,
                                  CollisionsWithMcLabelsCent const& collisions,
                                  McCollisionsWithCentrality const& mcCollisions,
                                  aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, false>(selectedD0CandidatesMc, dstarCandidates, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processMcWithDCAFitterNCent, "Process MC D0-Dstar pairs with DCAFitterN and centrality", false);

  void processMcWithKFParticle(D0CandidatesMcKF const&,
                              DstarCandidatesMc const& dstarCandidates,
                              aod::McParticles const& mcParticles,
                              TracksSelQuality const& tracks,
                              CollisionsWithMcLabels const& collisions,
                              McCollisionsWithCentrality const& mcCollisions,
                              aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, false>(selectedD0CandidatesMcKF, dstarCandidates, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processMcWithKFParticle, "Process MC D0-Dstar pairs with KFParticle D0 reconstruction", false);
  // TODO: add processMcWithKFParticleCent

  void processMcWithDCAFitterNMl(D0CandidatesMlMc const&,
                                DstarCandidatesMc const& dstarCandidates,
                                aod::McParticles const& mcParticles,
                                TracksSelQuality const& tracks,
                                CollisionsWithMcLabels const& collisions,
                                McCollisionsWithCentrality const& mcCollisions,
                                aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, true>(selectedD0CandidatesMlMc, dstarCandidates, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processMcWithDCAFitterNMl, "Process MC D0-Dstar pairs with DCAFitterN and D0 ML information", false);

  void processMcWithDCAFitterNMlCent(D0CandidatesMlMc const&,
                                    DstarCandidatesMc const& dstarCandidates,
                                    aod::McParticles const& mcParticles,
                                    TracksSelQuality const& tracks,
                                    CollisionsWithMcLabelsCent const& collisions,
                                    McCollisionsWithCentrality const& mcCollisions,
                                    aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter, true>(selectedD0CandidatesMlMc, dstarCandidates, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processMcWithDCAFitterNMlCent, "Process MC D0-Dstar pairs with DCAFitterN, D0 ML information and centrality", false);

  void processMcWithKFParticleMl(D0CandidatesMlMcKF const&,
                                DstarCandidatesMc const& dstarCandidates,
                                aod::McParticles const& mcParticles,
                                TracksSelQuality const& tracks,
                                CollisionsWithMcLabels const& collisions,
                                McCollisionsWithCentrality const& mcCollisions,
                                aod::BcFullInfos const& bcs)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle, true>(selectedD0CandidatesMlMcKF, dstarCandidates, mcParticles, tracks, collisions, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processMcWithKFParticleMl, "Process MC D0-Dstar pairs with KFParticle D0 reconstruction and D0 ML information", false);
  // TODO: add the processMcWithKFParticleMlCent

  void processDataWithDCAFitterNWithUpc(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults> const& collisions,
                                        aod::BcFullInfos const& bcs,
                                        D0Candidates const&,
                                        aod::TracksWExtra const&,
                                        aod::FT0s const& ft0s,
                                        aod::FV0As const& fv0as,
                                        aod::FDDs const& fdds,
                                        TracksWPid const& tracks,
                                        aod::Zdcs const& /*zdcs*/
  )
  {
    runAnalysisPerCollisionDataWithUpc<false>(collisions, selectedD0Candidates, bcs, ft0s, fv0as, fdds, tracks);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processDataWithDCAFitterNWithUpc, "Process real data with DCAFitterN w/o ML with UPC", false);

  void processDataWithDCAFitterNMlWithUpc(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults> const& collisions,
                                          aod::BcFullInfos const& bcs,
                                          D0CandidatesMl const&,
                                          aod::TracksWExtra const&,
                                          aod::FT0s const& ft0s,
                                          aod::FV0As const& fv0as,
                                          aod::FDDs const& fdds,
                                          TracksWPid const& tracks,
                                          aod::Zdcs const& /*zdcs*/)
  {
    runAnalysisPerCollisionDataWithUpc<true>(collisions, selectedD0CandidatesMl, bcs, ft0s, fv0as, fdds, tracks);
  }
  PROCESS_SWITCH(HfTaskD0Dstar, processDataWithDCAFitterNMlWithUpc, "Process real data with DCAFitterN and ML with UPC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskD0Dstar>(cfgc, TaskName{"hf-task-d0-dstar"})};
}
