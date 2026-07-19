// D0 meson finder.
//
// Reconstructs D0 -> K pi candidates directly (no intermediate D* resonance),
// from pairs of tracks (via their associated reconstructed particles).
//
// Motivation: the D* -> D0 pi -> K pi pi chain (analyzer_dstarfinder.cxx) has a
// fragmentation-fraction * branching-fraction product that turned out too small
// to be useful as a b/c-jet calibration handle on this dataset. The bare D0 -> K pi
// decay has a much larger branching fraction (no extra pi-slow leg to multiply in),
// at the cost of more combinatorial background (no D*-D0 mass-difference constraint
// to reject fakes). This file exists to let that trade-off be evaluated empirically.
//
// This is a straight extraction of the D0-pair-finding half of
// analyzer_dstarfinder.cxx (same track preselection, same K/pi mass-hypothesis
// logic, same same-sign background-estimate approach, same vertex-fit machinery),
// factored out so it can be used standalone here AND reused by the D* finder
// (which now takes an already-found DZeroCandidates collection as input rather
// than re-deriving D0 pairs itself -- mirroring how analyzer_bmesonfinder.cxx
// already takes an already-found DStarCandidates collection as input). See
// analyzer_dstarfinder.cxx for the original translation notes (CMSSW reference,
// ALEPH-specific adaptations, vertex-fit unit/ndf caveats); those all still apply
// here unchanged, since the underlying per-pair logic is untouched.

#include <ROOT/RVec.hxx>
#include <cmath>

#include "FCCAnalyses/VertexingUtils.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"

#include <iostream>
#include <algorithm>


#ifndef DZeroMesonFinder_H
#define DZeroMesonFinder_H

namespace DZeroMesonFinder{

// particle masses [GeV]
const double m_pi = 0.13957039;
const double m_K = 0.493677;
const double m_D0 = 1.86484;

// output collection: one entry per D0 candidate found in the event
struct DZeroCandidates{
    ROOT::VecOps::RVec<float> mass;
    ROOT::VecOps::RVec<float> pt;
    ROOT::VecOps::RVec<float> eta;
    ROOT::VecOps::RVec<float> phi;
    ROOT::VecOps::RVec<float> k_pt;
    ROOT::VecOps::RVec<float> k_eta;
    ROOT::VecOps::RVec<float> k_phi;
    ROOT::VecOps::RVec<int>   k_charge;
    ROOT::VecOps::RVec<float> pi_pt;
    ROOT::VecOps::RVec<float> pi_eta;
    ROOT::VecOps::RVec<float> pi_phi;
    ROOT::VecOps::RVec<int>   pi_charge;
    ROOT::VecOps::RVec<float> tr1tr2_deltaR;
    ROOT::VecOps::RVec<float> vtx_chi2Normalized;

    // distance between the fitted D0 (K/pi) vertex and the primary vertex -- for a
    // genuine D0, this is (D0 flight length) in a primary/c-jet production, or
    // (B flight length + D0 flight length) for a b-jet cascade D0.
    ROOT::VecOps::RVec<float> vtx_dxy;  // transverse (x-y) distance
    ROOT::VecOps::RVec<float> vtx_dxyz; // full 3D distance

    // true for genuine opposite-charge K/pi candidates (signal-like); false for
    // same-charge candidates, which are kept alongside as a combinatorial background
    // estimate (see getDZeroCandidates) and should generally be excluded from
    // anything downstream that assumes a genuine D0 (e.g. the D* finder only uses
    // candidates with isOppositeSign == true).
    ROOT::VecOps::RVec<bool> isOppositeSign;

    // indices into ReconstructedParticles for the 2 decay products, one entry per
    // candidate. For internal cross-referencing only (e.g. by DStarMesonFinder, to
    // reuse an already-found D0 candidate's tracks without re-deriving them); not
    // meant to be exposed as their own ntuple branches, since a bare RP index has no
    // standalone meaning without the RP collection itself also being saved.
    ROOT::VecOps::RVec<int> k_idx;
    ROOT::VecOps::RVec<int> pi_idx;
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

// helper: indices of reco particles passing isGoodTrack above. Exposed (not static
// to this file) so DStarMesonFinder can re-derive the same selected-track set for
// its own soft-pion loop, consistent with the set used to find the D0 candidates
// it receives as input.
std::vector<int> getSelectedIndices(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recoParticles,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& trackDatas,
        double minTrackPt){
    std::vector<int> selectedIndices;
    for(size_t idx = 0; idx < recoParticles.size(); idx++){
        const auto& rp = recoParticles[idx];
        size_t trackIdx = rp.tracks_begin;
        if(trackIdx >= trackStates.size()) continue;
        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);
        if(!isGoodTrack(trackStates[trackIdx], trackDatas[trackIdx], p, minTrackPt)) continue;
        selectedIndices.push_back((int)idx);
    }
    return selectedIndices;
}

DZeroCandidates getDZeroCandidates(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recoParticles,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& trackDatas,
        const TVector3& primaryVertex,
        double minTrackPt = 0.3,
        double minKaonPt = 1.0,
        double maxTrackPairDeltaR = 0.4,
        double d0MassWindow = 0.035,
        double maxVertexChi2Normalized = 5.,
        // deliberately generous: this caps D0 candidates per event, but the D* finder
        // also consumes this same collection (see analyzer_dstarfinder.cxx), where an
        // overly tight cap here would silently drop otherwise-valid D* candidates
        // built from D0 pairs found later in the loop. In practice, realistic ALEPH
        // events are nowhere near this many opposite/same-sign pairs passing the mass
        // window and vertex chi2 cut, so this is not expected to bite.
        unsigned int maxCandidates = 500){
    /*
    Find D0 -> K pi candidates.
    Note: cut values above are carried over from the CMS reference implementation
    that analyzer_dstarfinder.cxx was translated from (see there), since they were
    tuned for the CMS silicon tracker; ALEPH's TPC-based tracking has different
    resolution, so these will likely need retuning. This only affects
    efficiency/purity, not the correctness of the translation.
    */

    DZeroCandidates result;

    // preselect reco particles with a good associated track
    std::vector<int> selectedIndices = getSelectedIndices(recoParticles, trackStates, trackDatas, minTrackPt);

    // loop over pairs of tracks (K/pi candidates)
    for(size_t ii = 0; ii < selectedIndices.size(); ii++){
        for(size_t jj = ii+1; jj < selectedIndices.size(); jj++){
            int i = selectedIndices[ii];
            int j = selectedIndices[jj];
            const auto& rp1 = recoParticles[i];
            const auto& rp2 = recoParticles[j];

            // same-charge pairs are kept (rather than required to be opposite-charge)
            // to additionally provide a same-sign combinatorial background estimate,
            // following the same approach as analyzer_dstarfinder.cxx (see there);
            // isOppositeSign is stored per-candidate below so downstream code (e.g.
            // the D* finder) can select on it rather than have it hard-cut here.
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
            // by track index instead (for reproducibility; the choice is arbitrary in
            // this case regardless).
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
            TLorentzVector piP4, kP4;
            int piIdx = -1, kIdx = -1;
            bool matched = false;
            if( (std::abs(d0Mass - m_D0) < d0MassWindow)
                && (std::abs(d0Mass - m_D0) < std::abs(d0barMass - m_D0)) ){
                piP4 = piPlusP4; kP4 = KMinusP4;
                piIdx = posIdx; kIdx = negIdx;
                matched = true;
            } else if( (std::abs(d0barMass - m_D0) < d0MassWindow)
                       && (std::abs(d0barMass - m_D0) < std::abs(d0Mass - m_D0)) ){
                piP4 = piMinusP4; kP4 = KPlusP4;
                piIdx = negIdx; kIdx = posIdx;
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

            // distance between the D0 vertex and the primary vertex (see struct docs)
            TVector3 d0vtxPos(d0vtx.vertex.position.x, d0vtx.vertex.position.y, d0vtx.vertex.position.z);
            TVector3 d0vtxDiff = d0vtxPos - primaryVertex;
            float d0vtxDxy = std::sqrt(d0vtxDiff.X()*d0vtxDiff.X() + d0vtxDiff.Y()*d0vtxDiff.Y());
            float d0vtxDxyz = d0vtxDiff.Mag();

            // store the candidate
            result.mass.push_back(d0P4.M());
            result.pt.push_back(d0P4.Pt());
            result.eta.push_back(d0P4.Eta());
            result.phi.push_back(d0P4.Phi());
            result.k_pt.push_back(kP4.Pt());
            result.k_eta.push_back(kP4.Eta());
            result.k_phi.push_back(kP4.Phi());
            result.k_charge.push_back((int)recoParticles[kIdx].charge);
            result.pi_pt.push_back(piP4.Pt());
            result.pi_eta.push_back(piP4.Eta());
            result.pi_phi.push_back(piP4.Phi());
            result.pi_charge.push_back((int)recoParticles[piIdx].charge);
            result.tr1tr2_deltaR.push_back(dR12);
            result.vtx_chi2Normalized.push_back(d0vtxChi2);
            result.vtx_dxy.push_back(d0vtxDxy);
            result.vtx_dxyz.push_back(d0vtxDxyz);
            result.isOppositeSign.push_back(isOppositeSign);
            result.k_idx.push_back(kIdx);
            result.pi_idx.push_back(piIdx);

            if(result.mass.size() >= maxCandidates) return result;
        }
    } // end loop over first and second track

    return result;
}

// "poor man's" gen-level matching: check whether there is a gen-level D0 (PDG code
// 421) within dR of each reconstructed D0 candidate. Same caveats as
// DStarMesonFinder::matchToGenDStar (see analyzer_dstarfinder.cxx): purely
// geometric, charge-blind, whole-candidate match, since mother-daughter links are
// not populated in the gen-particle collection for this ALEPH EDM4HEP conversion.
ROOT::VecOps::RVec<bool> matchToGenDZero(
        const ROOT::VecOps::RVec<float>& cand_pt,
        const ROOT::VecOps::RVec<float>& cand_eta,
        const ROOT::VecOps::RVec<float>& cand_phi,
        const ROOT::VecOps::RVec<float>& cand_mass,
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles,
        double dRThreshold = 0.1){

    ROOT::VecOps::RVec<bool> result;

    // collect 4-vectors of gen-level D0 mesons
    std::vector<TLorentzVector> genDZeros;
    for(const auto& gp : genParticles){
        if(std::abs(gp.PDG) != 421) continue;
        TLorentzVector p4;
        p4.SetXYZM(gp.momentum.x, gp.momentum.y, gp.momentum.z, gp.mass);
        genDZeros.push_back(p4);
    }

    // match each candidate to the closest gen D0 (if any within threshold)
    for(size_t i = 0; i < cand_pt.size(); i++){
        TLorentzVector cand;
        cand.SetPtEtaPhiM(cand_pt[i], cand_eta[i], cand_phi[i], cand_mass[i]);
        bool matched = false;
        for(const auto& genD : genDZeros){
            if(cand.DeltaR(genD) < dRThreshold){ matched = true; break; }
        }
        result.push_back(matched);
    }
    return result;
}

// same as above, but a dummy for running on data (no gen particles available)
ROOT::VecOps::RVec<bool> matchToGenDZeroDummy(const ROOT::VecOps::RVec<float>& cand_pt){
    ROOT::VecOps::RVec<bool> result(cand_pt.size(), false);
    return result;
}

}

#endif
