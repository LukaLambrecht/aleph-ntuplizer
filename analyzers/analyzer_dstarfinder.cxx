// D* meson finder.
//
// Reconstructs D* -> D0 pi -> K pi pi candidates from triplets of tracks
// (via their associated reconstructed particles).
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

#include <ROOT/RVec.hxx>
#include <cmath>

#include "FCCAnalyses/VertexingUtils.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"

#include <iostream>
#include <algorithm>


#ifndef DStarMesonFinder_H
#define DStarMesonFinder_H

namespace DStarMesonFinder{

// particle masses [GeV]
const double m_pi = 0.13957039;
const double m_K = 0.493677;
const double m_D0 = 1.86484;
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

// helper: same track quality selection as TrackTools::getSelectedTracks
// (covariance sanity checks + chi2/ndof), applied here directly on the
// reco-particle's associated track, plus a minimum pt requirement.
// kept self-contained rather than factored out, to avoid coupling this file
// to the exact interface of TrackTools::getSelectedTracks.
bool isGoodTrack(
        const edm4hep::TrackState& trst,
        const edm4hep::TrackData& trkobj,
        const TVector3& p,
        double minPt){
    const auto& c = trst.covMatrix;
    if (c[0] <= 0 || c[2] <= 0 || c[9] <= 0) return false;
    if (!std::isfinite(c[0]) || !std::isfinite(c[2]) || !std::isfinite(c[9])) return false;
    if (trkobj.ndf == 0) return false;
    if (trkobj.chi2 / trkobj.ndf > 10.) return false;
    if (p.Pt() < minPt) return false;
    return true;
}

DStarCandidates getDStarCandidates(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recoParticles,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& trackDatas,
        double minTrackPt = 0.3,
        double minKaonPt = 1.0,
        double minSoftPionPt = 0.5,
        double maxTrackPairDeltaR = 0.4,
        double d0MassWindow = 0.035,
        double dstarMassWindow = 0.1,
        double maxSoftPionDeltaR = 0.1,
        double maxVertexChi2Normalized = 5.,
        unsigned int maxCandidates = 30){
    /*
    Find D* -> D0 pi -> K pi pi candidates.
    Note: cut values above are carried over from the CMS reference implementation
    as defaults, since they were tuned for the CMS silicon tracker; ALEPH's
    TPC-based tracking has different resolution, so these will likely need
    retuning. This only affects efficiency/purity, not the correctness of the
    translation.
    */

    DStarCandidates result;

    // preselect reco particles with a good associated track
    std::vector<int> selectedIndices;
    for(size_t idx = 0; idx < recoParticles.size(); idx++){
        const auto& rp = recoParticles[idx];
        size_t trackIdx = rp.tracks_begin;
        if(trackIdx >= trackStates.size()) continue;
        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);
        if(!isGoodTrack(trackStates[trackIdx], trackDatas[trackIdx], p, minTrackPt)) continue;
        selectedIndices.push_back((int)idx);
    }

    // loop over pairs of tracks (K/pi2 candidates, i.e. the D0 decay products)
    for(size_t ii = 0; ii < selectedIndices.size(); ii++){
        for(size_t jj = ii+1; jj < selectedIndices.size(); jj++){
            int i = selectedIndices[ii];
            int j = selectedIndices[jj];
            const auto& rp1 = recoParticles[i];
            const auto& rp2 = recoParticles[j];

            // same-charge pairs are kept (rather than required to be opposite-charge)
            // to additionally provide a same-sign combinatorial background estimate,
            // following the same approach as the reference implementation (see file
            // docstring); isOppositeSign is stored per-candidate below so downstream
            // code (e.g. the B finder) can select on it rather than have it hard-cut here.
            bool isOppositeSign = (rp1.charge * rp2.charge) < 0;

            TLorentzVector p1; p1.SetXYZM(rp1.momentum.x, rp1.momentum.y, rp1.momentum.z, rp1.mass);
            TLorentzVector p2; p2.SetXYZM(rp2.momentum.x, rp2.momentum.y, rp2.momentum.z, rp2.mass);

            // candidates must point approximately in the same direction
            double dR12 = p1.DeltaR(p2);
            if(dR12 > maxTrackPairDeltaR) continue;

            // find which track plays the "positive" and "negative" role in the K/pi
            // mass-hypothesis test below. For opposite-charge pairs this follows the
            // actual physical charge; for same-charge pairs (kept only for the
            // combinatorial background estimate, see above) there is no physical
            // meaning to this assignment either way, so it is fixed deterministically
            // by track index instead of the reference implementation's random choice
            // (for reproducibility; the choice is arbitrary in this case regardless).
            int posIdx = isOppositeSign ? ((rp1.charge > 0) ? i : j) : i;
            int negIdx = isOppositeSign ? ((rp1.charge > 0) ? j : i) : j;
            const auto& rpPos = recoParticles[posIdx];
            const auto& rpNeg = recoParticles[negIdx];

            // make invariant mass under both K-pi hypotheses
            // (although the D0 decays preferentially to K- pi+, the original
            // particle could be an anti-D0, which decays preferentially to K+ pi-)
            TLorentzVector piPlusP4;  piPlusP4.SetXYZM(rpPos.momentum.x, rpPos.momentum.y, rpPos.momentum.z, m_pi);
            TLorentzVector KMinusP4;  KMinusP4.SetXYZM(rpNeg.momentum.x, rpNeg.momentum.y, rpNeg.momentum.z, m_K);
            TLorentzVector KPlusP4;   KPlusP4.SetXYZM(rpPos.momentum.x, rpPos.momentum.y, rpPos.momentum.z, m_K);
            TLorentzVector piMinusP4; piMinusP4.SetXYZM(rpNeg.momentum.x, rpNeg.momentum.y, rpNeg.momentum.z, m_pi);

            TLorentzVector d0P4 = piPlusP4 + KMinusP4;
            TLorentzVector d0barP4 = piMinusP4 + KPlusP4;
            double d0Mass = d0P4.M();
            double d0barMass = d0barP4.M();

            // invariant mass must be close to the D0 resonance mass
            TLorentzVector pi2P4, kP4;
            int pi2Idx = -1, kIdx = -1;
            bool matched = false;
            if( (std::abs(d0Mass - m_D0) < d0MassWindow)
                && (std::abs(d0Mass - m_D0) < std::abs(d0barMass - m_D0)) ){
                pi2P4 = piPlusP4; kP4 = KMinusP4;
                pi2Idx = posIdx; kIdx = negIdx;
                matched = true;
            } else if( (std::abs(d0barMass - m_D0) < d0MassWindow)
                       && (std::abs(d0barMass - m_D0) < std::abs(d0Mass - m_D0)) ){
                pi2P4 = piMinusP4; kP4 = KPlusP4;
                pi2Idx = negIdx; kIdx = posIdx;
                d0P4 = d0barP4; d0Mass = d0barMass;
                matched = true;
            }
            if(!matched) continue;

            // K candidate must have a minimum pt
            if(kP4.Pt() < minKaonPt) continue;

            // fit a D0 vertex from the pair of tracks
            ROOT::VecOps::RVec<edm4hep::TrackState> pairTracks;
            pairTracks.push_back(trackStates[rp1.tracks_begin]);
            pairTracks.push_back(trackStates[rp2.tracks_begin]);
            FCCAnalyses::VertexingUtils::FCCAnalysesVertex d0vtx = VertexFitterMod(2, pairTracks);
            double d0vtxChi2 = d0vtx.vertex.chi2;
            if(d0vtxChi2 < 0. || d0vtxChi2 > maxVertexChi2Normalized) continue;

            // loop over third track (soft pion candidate)
            for(size_t kk = 0; kk < selectedIndices.size(); kk++){
                int k = selectedIndices[kk];
                if(k == i || k == j) continue;
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

                // store the candidate
                result.mass.push_back(dstarP4.M());
                result.pt.push_back(dstarP4.Pt());
                result.eta.push_back(dstarP4.Eta());
                result.phi.push_back(dstarP4.Phi());
                result.d0_mass.push_back(d0P4.M());
                result.d0_pt.push_back(d0P4.Pt());
                result.d0_eta.push_back(d0P4.Eta());
                result.d0_phi.push_back(d0P4.Phi());
                result.d0_massDiff.push_back(dstarP4.M() - d0P4.M());
                result.pi1_pt.push_back(pi1P4.Pt());
                result.pi1_eta.push_back(pi1P4.Eta());
                result.pi1_phi.push_back(pi1P4.Phi());
                result.pi1_charge.push_back((int)rp3.charge);
                result.k_pt.push_back(kP4.Pt());
                result.k_eta.push_back(kP4.Eta());
                result.k_phi.push_back(kP4.Phi());
                result.k_charge.push_back((int)recoParticles[kIdx].charge);
                result.pi2_pt.push_back(pi2P4.Pt());
                result.pi2_eta.push_back(pi2P4.Eta());
                result.pi2_phi.push_back(pi2P4.Phi());
                result.pi2_charge.push_back((int)recoParticles[pi2Idx].charge);
                result.tr1tr2_deltaR.push_back(dR12);
                result.tr3d0_deltaR.push_back(dR3D0);
                result.d0vtx_chi2Normalized.push_back(d0vtxChi2);
                result.dstarvtx_chi2Normalized.push_back(dstarvtxChi2);
                result.isOppositeSign.push_back(isOppositeSign);
                result.k_idx.push_back(kIdx);
                result.pi2_idx.push_back(pi2Idx);
                result.pi1_idx.push_back(k);

                if(result.mass.size() >= maxCandidates) return result;
            } // end loop over third track
        }
    } // end loop over first and second track

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
