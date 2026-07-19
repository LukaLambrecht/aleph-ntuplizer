// D* meson finder.
//
// Reconstructs D* -> D0 pi -> K pi pi candidates by combining an already-found D0
// candidate (see analyzer_dzerofinder.cxx) with a third track (the "soft"/slow pion).
//
// Translated from the CMSSW / nanoAOD version at:
// https://github.com/LukaLambrecht/HcNano/blob/main/HcNano/plugins/HToDStarMesonProducer.cc
// (see also the companion header HToDStarMesonProducer.h for the constants).
//
// Adaptations made for ALEPH / FCCAnalyses conventions:
// - combinatorics run over ReconstructedParticles rather than raw tracks:
//   ALEPH reco particles already carry the (magnetic-field-corrected) momentum
//   and charge needed for the K/pi mass-hypothesis logic, while the associated
//   edm4hep::TrackState (via rp.tracks_begin, a direct link after
//   RecoTools::getModifiedRecoParticles) provides the parameters for the vertex fit.
// - the CMS-specific track/vertex "reference point" pre-selection cuts were dropped.
//   In CMS miniAOD, pat::PackedCandidate::referencePoint() encodes the primary
//   vertex a candidate was assigned to by particle flow, so comparing reference
//   points is a cheap way to check whether two tracks were assigned to the same
//   vertex. edm4hep::TrackState.referencePoint has no such meaning here (it is
//   documented as "reference point of the track parameters, e.g. the origin at
//   the IP", i.e. a shared coordinate anchor, not a per-candidate vertex
//   assignment), so there is no meaningful equivalent to translate. The actual
//   vertex-compatibility requirement is still fully enforced below via the chi2
//   cut on the fitted vertex, exactly as in the reference; dropping the cheap
//   pre-check only costs some extra combinatorics, not correctness.
// - the reference code's dstarmass constant (1.96847 GeV) is actually the D_s
//   meson mass, not the D*(2010)+/- mass (2.01026 GeV, PDG); this looks like a
//   copy-paste bug in the source repo (the wide +-0.1 GeV window meant it still
//   mostly worked). The correct value is used here.
// - gen-matching (present in the reference) is not translated: it depends on
//   CMS-specific H->D* gen-chain matching machinery that has no equivalent in
//   this repository, and does not apply to ALEPH qqbar simulation anyway.
// - vertex fitting reuses VertexFitterMod(...), defined in analyzer_svfinder.cxx,
//   the same function already used by the secondary vertex / V0 finders.
//   Note: this file does not #include analyzer_svfinder.cxx (it has no
//   top-level include guard of its own, unlike its namespaced counterparts,
//   so including it here in addition to analysis.py's own
//   ROOT.gInterpreter.Declare(...) of it would cause duplicate-definition
//   errors). Instead, analysis.py must Declare this file *after*
//   analyzer_svfinder.cxx, so VertexFitterMod is already available in the
//   (persistent, cross-Declare-call) interpreter state by the time this file
//   is loaded -- exactly the same mechanism that already lets analysis.py's
//   own .Define(...) calls reference namespaces declared in earlier files.
// - also note: VertexFitterMod never populates FCCAnalysesVertex.vertex.ndf
//   (only .vertex.chi2, already normalized to chi2/ndof internally). The
//   existing SV finder code already relies only on .vertex.chi2 for this
//   reason (see check_constraints() in analyzer_svfinder.cxx); this file
//   does the same.
//
// Update: the D0-pair-finding half of this file (the K/pi combinatorics, mass
// hypothesis logic, and D0 vertex fit) has been factored out into
// analyzer_dzerofinder.cxx, since the bare D0 -> K pi decay turned out to be a
// more promising b/c-jet calibration channel on its own (much larger branching
// fraction, at the cost of more combinatorial background -- the D* chain's extra
// pi-slow leg suppressed the yield too much to be useful here). getDStarCandidates
// below now takes an already-found DZeroMesonFinder::DZeroCandidates collection as
// input, exactly the same pattern already used by analyzer_bmesonfinder.cxx (which
// takes an already-found DStarCandidates collection rather than re-deriving D*
// candidates itself) -- this lets analysis.py compute the D0 candidates exactly
// once and share them between the standalone D0 output, the D* finder, and (via the
// D* finder) the B finder, regardless of which combination of do_dzero_candidates /
// do_dstar_candidates / do_bmeson_candidates flags is switched on.

#include <ROOT/RVec.hxx>
#include <cmath>

#include "FCCAnalyses/VertexingUtils.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"

#include <iostream>
#include <algorithm>

// safe to include: analyzer_dzerofinder.cxx has its own top-level include guard
// (DZeroMesonFinder_H), so this is a no-op if analysis.py has already separately
// Declared it (which it must do anyway, after analyzer_svfinder.cxx, whenever
// do_dzero_candidates/do_dstar_candidates/do_bmeson_candidates need it).
#include "analyzer_dzerofinder.cxx"


#ifndef DStarMesonFinder_H
#define DStarMesonFinder_H

namespace DStarMesonFinder{

// particle masses [GeV]
const double m_pi = 0.13957039;
const double m_DStar = 2.01026; // D*(2010)+/-

// output collection: one entry per D* candidate found in the event
struct DStarCandidates{
    ROOT::VecOps::RVec<float> mass;
    ROOT::VecOps::RVec<float> pt;
    ROOT::VecOps::RVec<float> eta;
    ROOT::VecOps::RVec<float> phi;
    ROOT::VecOps::RVec<float> d0_mass;
    ROOT::VecOps::RVec<float> d0_pt;
    ROOT::VecOps::RVec<float> d0_eta;
    ROOT::VecOps::RVec<float> d0_phi;
    ROOT::VecOps::RVec<float> d0_massDiff;
    ROOT::VecOps::RVec<float> pi1_pt; // soft pion (from D* decay)
    ROOT::VecOps::RVec<float> pi1_eta;
    ROOT::VecOps::RVec<float> pi1_phi;
    ROOT::VecOps::RVec<int>   pi1_charge;
    ROOT::VecOps::RVec<float> k_pt; // kaon (from D0 decay)
    ROOT::VecOps::RVec<float> k_eta;
    ROOT::VecOps::RVec<float> k_phi;
    ROOT::VecOps::RVec<int>   k_charge;
    ROOT::VecOps::RVec<float> pi2_pt; // pion (from D0 decay)
    ROOT::VecOps::RVec<float> pi2_eta;
    ROOT::VecOps::RVec<float> pi2_phi;
    ROOT::VecOps::RVec<int>   pi2_charge;
    ROOT::VecOps::RVec<float> tr1tr2_deltaR;
    ROOT::VecOps::RVec<float> tr3d0_deltaR;
    ROOT::VecOps::RVec<float> d0vtx_chi2Normalized;
    ROOT::VecOps::RVec<float> dstarvtx_chi2Normalized;

    // distance between the D0 vertex (2-track K/pi2 fit) and the primary vertex, and
    // the same for the D* vertex (3-track fit incl. the soft pion). Both are provided
    // since it's not obvious a priori which is the better b-vs-c-jet discriminator:
    // - d0vtx is a clean fit of the true D0 decay point (the K/pi2 genuinely originate
    //   there), so its distance to the PV is (D0 flight length) for a primary/c-jet D0,
    //   or (B flight length + D0 flight length) for a b-jet cascade D0 -- always offset
    //   by at least the D0's own flight length, even for c-jets.
    // - dstarvtx additionally includes the slow pion, which genuinely originates from
    //   the D* decay point itself (essentially the D* *production* point too, since the
    //   D* does not fly any measurable distance) -- i.e. at the primary vertex for a
    //   prompt/c-jet D*, or at the B decay vertex for a b-jet cascade D*. But since the
    //   K/pi2 do NOT originate from that same point (they come from the displaced D0
    //   decay point further downstream), forcing all 3 tracks into one common vertex
    //   fit is itself an approximation -- the fitted position ends up some resolution-
    //   weighted compromise between the D* and D0 decay points, not a clean fit of
    //   either. Which of the two ends up more discriminating is an empirical question.
    ROOT::VecOps::RVec<float> d0vtx_dxy;  // transverse (x-y) distance
    ROOT::VecOps::RVec<float> d0vtx_dxyz; // full 3D distance
    ROOT::VecOps::RVec<float> dstarvtx_dxy;
    ROOT::VecOps::RVec<float> dstarvtx_dxyz;

    // true for genuine opposite-charge K/pi2 candidates (signal-like); false for
    // same-charge candidates, which are kept alongside as a combinatorial background
    // estimate (see getDStarCandidates) and should generally be excluded from anything
    // downstream that assumes a genuine D0 (e.g. the B finder only uses candidates
    // with isOppositeSign == true).
    ROOT::VecOps::RVec<bool> isOppositeSign;

    // indices into ReconstructedParticles for the 3 decay products, one entry per
    // candidate. For internal cross-referencing only (e.g. by BMesonFinder, to reuse
    // an already-found D* candidate's tracks without re-deriving them); not meant to
    // be exposed as their own ntuple branches, since a bare RP index has no standalone
    // meaning without the RP collection itself also being saved.
    ROOT::VecOps::RVec<int> k_idx;
    ROOT::VecOps::RVec<int> pi2_idx;
    ROOT::VecOps::RVec<int> pi1_idx;
};

// note: the isGoodTrack track-quality helper and the K/pi pair-finding combinatorics
// that used to live here have moved to analyzer_dzerofinder.cxx (see DZeroMesonFinder::
// isGoodTrack / getSelectedIndices / getDZeroCandidates); this file now only adds the
// third (soft pion) track to an already-found D0 candidate.

DStarCandidates getDStarCandidates(
        const DZeroMesonFinder::DZeroCandidates& d0Candidates,
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recoParticles,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& trackDatas,
        const TVector3& primaryVertex,
        double minTrackPt = 0.3,
        double minSoftPionPt = 0.5,
        double dstarMassWindow = 0.1,
        double maxSoftPionDeltaR = 0.1,
        double maxVertexChi2Normalized = 5.,
        unsigned int maxCandidates = 30){
    /*
    Find D* -> D0 pi -> K pi pi candidates, by combining each entry of d0Candidates
    (see analyzer_dzerofinder.cxx::getDZeroCandidates) with a third track (the soft
    pion). minTrackPt is used only to re-derive the same selected-track set that
    d0Candidates was built from (for the soft-pion loop below); it should match the
    minTrackPt value used to compute d0Candidates, or the soft-pion candidate pool
    will be inconsistent with the D0 pair selection.
    Note: cut values above are carried over from the CMS reference implementation
    as defaults, since they were tuned for the CMS silicon tracker; ALEPH's
    TPC-based tracking has different resolution, so these will likely need
    retuning. This only affects efficiency/purity, not the correctness of the
    translation.
    */

    DStarCandidates result;

    // preselect reco particles with a good associated track (same selection that
    // d0Candidates was built from -- needed here again to draw the soft-pion
    // candidate from the same quality-selected pool)
    std::vector<int> selectedIndices = DZeroMesonFinder::getSelectedIndices(
        recoParticles, trackStates, trackDatas, minTrackPt);

    // loop over already-found D0 candidates
    for(size_t d0idx = 0; d0idx < d0Candidates.mass.size(); d0idx++){

        int kIdx = d0Candidates.k_idx[d0idx];
        int pi2Idx = d0Candidates.pi_idx[d0idx];

        // re-derive the D0 4-vector from the K/pi momenta directly (rather than via
        // SetPtEtaPhiM from d0Candidates' stored float pt/eta/phi/mass), so this
        // reproduces the exact same double-precision value getDZeroCandidates used
        // internally, bit-for-bit -- avoids a float32 round-trip precision loss that
        // would otherwise leak into dstarP4 = d0P4 + pi1P4 below.
        const auto& rpK = recoParticles[kIdx];
        const auto& rpPi2 = recoParticles[pi2Idx];
        TLorentzVector kP4_forD0; kP4_forD0.SetXYZM(rpK.momentum.x, rpK.momentum.y, rpK.momentum.z, DZeroMesonFinder::m_K);
        TLorentzVector pi2P4_forD0; pi2P4_forD0.SetXYZM(rpPi2.momentum.x, rpPi2.momentum.y, rpPi2.momentum.z, DZeroMesonFinder::m_pi);
        TLorentzVector d0P4 = kP4_forD0 + pi2P4_forD0;

        // fit a D0 vertex from the pair of tracks (re-derived here since d0Candidates
        // only stores summary observables, not the intermediate fit result -- same
        // tracks, so this reproduces the exact same fit as in getDZeroCandidates)
        ROOT::VecOps::RVec<edm4hep::TrackState> pairTracks;
        pairTracks.push_back(trackStates[recoParticles[kIdx].tracks_begin]);
        pairTracks.push_back(trackStates[recoParticles[pi2Idx].tracks_begin]);

        // loop over third track (soft pion candidate)
        for(size_t kk = 0; kk < selectedIndices.size(); kk++){
            int k = selectedIndices[kk];
            if(k == kIdx || k == pi2Idx) continue;
            const auto& rp3 = recoParticles[k];

            TLorentzVector pi1P4;
            pi1P4.SetXYZM(rp3.momentum.x, rp3.momentum.y, rp3.momentum.z, m_pi);

            // pi candidate must have a minimum pt
            if(pi1P4.Pt() < minSoftPionPt) continue;

            // candidate must point approximately in the same direction as the D0
            double dR3D0 = pi1P4.DeltaR(d0P4);
            if(dR3D0 > maxSoftPionDeltaR) continue;

            // make invariant mass (under the pi mass hypothesis for the third track)
            TLorentzVector dstarP4 = d0P4 + pi1P4;
            double dstarMass = dstarP4.M();

            // check if mass is close enough to the D* mass
            if(std::abs(dstarMass - m_DStar) > dstarMassWindow) continue;

            // fit a D* vertex from the triplet of tracks
            ROOT::VecOps::RVec<edm4hep::TrackState> tripletTracks = pairTracks;
            tripletTracks.push_back(trackStates[rp3.tracks_begin]);
            FCCAnalyses::VertexingUtils::FCCAnalysesVertex dstarvtx = VertexFitterMod(2, tripletTracks);
            double dstarvtxChi2 = dstarvtx.vertex.chi2;
            if(dstarvtxChi2 < 0. || dstarvtxChi2 > maxVertexChi2Normalized) continue;

            // distance between the D* vertex and the primary vertex (see struct docs)
            TVector3 dstarvtxPos(dstarvtx.vertex.position.x, dstarvtx.vertex.position.y, dstarvtx.vertex.position.z);
            TVector3 dstarvtxDiff = dstarvtxPos - primaryVertex;
            float dstarvtxDxy = std::sqrt(dstarvtxDiff.X()*dstarvtxDiff.X() + dstarvtxDiff.Y()*dstarvtxDiff.Y());
            float dstarvtxDxyz = dstarvtxDiff.Mag();

            // store the candidate
            result.mass.push_back(dstarP4.M());
            result.pt.push_back(dstarP4.Pt());
            result.eta.push_back(dstarP4.Eta());
            result.phi.push_back(dstarP4.Phi());
            result.d0_mass.push_back(d0Candidates.mass[d0idx]);
            result.d0_pt.push_back(d0Candidates.pt[d0idx]);
            result.d0_eta.push_back(d0Candidates.eta[d0idx]);
            result.d0_phi.push_back(d0Candidates.phi[d0idx]);
            // uses the (double-precision, just-reconstructed) d0P4.M() rather than
            // d0Candidates.mass[d0idx] (already float-rounded), to avoid a spurious
            // extra rounding step in the subtraction
            result.d0_massDiff.push_back(dstarP4.M() - d0P4.M());
            result.pi1_pt.push_back(pi1P4.Pt());
            result.pi1_eta.push_back(pi1P4.Eta());
            result.pi1_phi.push_back(pi1P4.Phi());
            result.pi1_charge.push_back((int)rp3.charge);
            result.k_pt.push_back(d0Candidates.k_pt[d0idx]);
            result.k_eta.push_back(d0Candidates.k_eta[d0idx]);
            result.k_phi.push_back(d0Candidates.k_phi[d0idx]);
            result.k_charge.push_back(d0Candidates.k_charge[d0idx]);
            result.pi2_pt.push_back(d0Candidates.pi_pt[d0idx]);
            result.pi2_eta.push_back(d0Candidates.pi_eta[d0idx]);
            result.pi2_phi.push_back(d0Candidates.pi_phi[d0idx]);
            result.pi2_charge.push_back(d0Candidates.pi_charge[d0idx]);
            result.tr1tr2_deltaR.push_back(d0Candidates.tr1tr2_deltaR[d0idx]);
            result.tr3d0_deltaR.push_back(dR3D0);
            result.d0vtx_chi2Normalized.push_back(d0Candidates.vtx_chi2Normalized[d0idx]);
            result.dstarvtx_chi2Normalized.push_back(dstarvtxChi2);
            result.d0vtx_dxy.push_back(d0Candidates.vtx_dxy[d0idx]);
            result.d0vtx_dxyz.push_back(d0Candidates.vtx_dxyz[d0idx]);
            result.dstarvtx_dxy.push_back(dstarvtxDxy);
            result.dstarvtx_dxyz.push_back(dstarvtxDxyz);
            result.isOppositeSign.push_back(d0Candidates.isOppositeSign[d0idx]);
            result.k_idx.push_back(kIdx);
            result.pi2_idx.push_back(pi2Idx);
            result.pi1_idx.push_back(k);

            if(result.mass.size() >= maxCandidates) return result;
        } // end loop over third track
    } // end loop over D0 candidates

    return result;
}

// "poor man's" gen-level matching: check whether there is a gen-level D*(2010)+/-
// meson (PDG code 413) within dR of each reconstructed D* candidate.
// note: this is a purely geometric, charge-blind match on the D* candidate as a
// whole (not on its individual decay products), since mother-daughter links are
// not populated in the gen-particle collection for this ALEPH EDM4HEP conversion
// (checked directly: MCParticles.parents_begin/end and daughters_begin/end are
// zero-length for every particle, and the _MCParticles_parents/_daughters link
// collections are empty). This means it cannot distinguish a genuine D* decay
// from an accidental nearby gen D* in a busy event, and cannot be used to tune
// selection cuts against individual decay products -- only to get a rough,
// first-order handle on how much of the reconstructed peak is real.
ROOT::VecOps::RVec<bool> matchToGenDStar(
        const ROOT::VecOps::RVec<float>& cand_pt,
        const ROOT::VecOps::RVec<float>& cand_eta,
        const ROOT::VecOps::RVec<float>& cand_phi,
        const ROOT::VecOps::RVec<float>& cand_mass,
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles,
        double dRThreshold = 0.1){

    ROOT::VecOps::RVec<bool> result;

    // collect 4-vectors of gen-level D*(2010)+/- mesons
    std::vector<TLorentzVector> genDStars;
    for(const auto& gp : genParticles){
        if(std::abs(gp.PDG) != 413) continue;
        TLorentzVector p4;
        p4.SetXYZM(gp.momentum.x, gp.momentum.y, gp.momentum.z, gp.mass);
        genDStars.push_back(p4);
    }

    // match each candidate to the closest gen D* (if any within threshold)
    for(size_t i = 0; i < cand_pt.size(); i++){
        TLorentzVector cand;
        cand.SetPtEtaPhiM(cand_pt[i], cand_eta[i], cand_phi[i], cand_mass[i]);
        bool matched = false;
        for(const auto& genD : genDStars){
            if(cand.DeltaR(genD) < dRThreshold){ matched = true; break; }
        }
        result.push_back(matched);
    }
    return result;
}

// same as above, but a dummy for running on data (no gen particles available)
ROOT::VecOps::RVec<bool> matchToGenDStarDummy(const ROOT::VecOps::RVec<float>& cand_pt){
    ROOT::VecOps::RVec<bool> result(cand_pt.size(), false);
    return result;
}

}

#endif
