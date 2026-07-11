// Semileptonic B meson finder.
//
// Reconstructs the "tag side" of B -> D* l nu decays (D* -> D0 pi -> K pi pi, plus a
// nearby lepton with the standard charge correlation), for use as a clean,
// exclusively-tagged b-jet to calibrate an inclusive b-tagger against (the classic LEP
// double-tag / hemisphere technique).
//
// Note: since the neutrino is unmeasured, this does *not* reconstruct the B mass --
// only the visible (D* + lepton) system. The purpose here is a clean tag, not a mass
// measurement.
//
// This deliberately reuses the D* finder (analyzer_dstarfinder.cxx) rather than
// reimplementing D0/D* reconstruction: analysis.py computes the D* candidates once
// (whenever D* finding and/or B finding is switched on) and passes them in here, so
// this file never repeats that combinatorial search itself.
//
// Note on the lepton charge convention: the standard "right-sign" combination is
// lepton charge opposite the D* candidate's kaon charge (equivalently, same sign as
// the soft pion). This is stored as a flag (isRightSign) rather than applied as a hard
// cut, so wrong-sign candidates remain available as the standard background/mixing
// control sample. The convention itself is recalled from memory, not pulled from a
// live reference -- it should be checked empirically (e.g. gen-matched candidates
// should be predominantly right-sign) before being trusted for an actual analysis.

#include <ROOT/RVec.hxx>
#include <cmath>

#include "FCCAnalyses/VertexingUtils.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/ParticleIDData.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>

// safe to include directly: analyzer_dstarfinder.cxx has its own top-level include
// guard (DStarMesonFinder_H), unlike analyzer_svfinder.cxx. This gives compile-time
// access to DStarMesonFinder::DStarCandidates and its mass constants.
#include "analyzer_dstarfinder.cxx"


#ifndef BMesonFinder_H
#define BMesonFinder_H

namespace BMesonFinder{

// lepton masses [GeV]
const double m_e = 0.000511;
const double m_mu = 0.105658;

// output collection: one entry per (D* candidate, lepton) pair found in the event,
// after keeping only the best D* per distinct lepton (see getBCandidates)
struct BCandidates{
    ROOT::VecOps::RVec<int>   dstar_idx; // index into the DStarCandidates arrays this came from
    ROOT::VecOps::RVec<int>   lepton_idx; // index into ReconstructedParticles for the lepton
    ROOT::VecOps::RVec<float> lepton_pt;
    ROOT::VecOps::RVec<float> lepton_eta;
    ROOT::VecOps::RVec<float> lepton_phi;
    ROOT::VecOps::RVec<int>   lepton_charge;
    ROOT::VecOps::RVec<int>   lepton_type; // 1 = electron, 2 = muon (matches ParticleID.type)
    ROOT::VecOps::RVec<float> mass; // visible (D* + lepton) mass -- NOT the true B mass
    ROOT::VecOps::RVec<float> pt;
    ROOT::VecOps::RVec<float> eta;
    ROOT::VecOps::RVec<float> phi;
    ROOT::VecOps::RVec<float> dstarlepton_deltaR;
    ROOT::VecOps::RVec<bool>  isRightSign;
    ROOT::VecOps::RVec<float> bvtx_chi2Normalized;
};

// internal helper: one raw (D*, lepton) pair before the per-lepton best-candidate
// selection (see getBCandidates)
struct RawBCandidate{
    int dstar_idx;
    int lepton_idx;
    float lepton_pt, lepton_eta, lepton_phi;
    int lepton_charge, lepton_type;
    float mass, pt, eta, phi;
    float dstarlepton_deltaR;
    bool isRightSign;
    float bvtx_chi2Normalized;
};

BCandidates getBCandidates(
        const DStarMesonFinder::DStarCandidates& dstarCandidates,
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& recoParticles,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& trackDatas,
        const ROOT::VecOps::RVec<edm4hep::ParticleIDData>& particleIDs,
        double minLeptonPt = 1.0,
        double maxDstarLeptonDeltaR = 0.6,
        double maxVertexChi2Normalized = 5.,
        unsigned int maxCandidates = 30){
    /*
    Find (D*, lepton) pairs consistent with a semileptonic B -> D* l nu decay.
    Note: cut values above are starting defaults (no reference implementation to
    translate this time, unlike the D* finder); meant to be retuned.
    */

    // collect all passing (D*, lepton) pairs first; the same lepton can otherwise
    // end up attached to several different D* candidates (real or fake) found in
    // the same event, which are not independent observations of anything -- see
    // the per-lepton best-candidate selection below.
    std::vector<RawBCandidate> rawCandidates;

    // loop over already-found D* candidates
    for(size_t d = 0; d < dstarCandidates.mass.size(); d++){
        // only genuine opposite-charge D* candidates make sense as the D* leg of a
        // semileptonic B decay; same-charge ones are kept in DStarCandidates purely
        // as a combinatorial background estimate for the D* peak itself (see
        // analyzer_dstarfinder.cxx) and are not meaningful here.
        if(!dstarCandidates.isOppositeSign[d]) continue;

        TLorentzVector dstarP4;
        dstarP4.SetPtEtaPhiM(dstarCandidates.pt[d], dstarCandidates.eta[d], dstarCandidates.phi[d], dstarCandidates.mass[d]);
        int kIdx = dstarCandidates.k_idx[d];
        int pi2Idx = dstarCandidates.pi2_idx[d];
        int pi1Idx = dstarCandidates.pi1_idx[d];
        int kCharge = dstarCandidates.k_charge[d];

        // loop over reco particles to find a lepton candidate
        for(size_t idx = 0; idx < recoParticles.size(); idx++){
            // skip the D* candidate's own 3 tracks
            if((int)idx == kIdx || (int)idx == pi2Idx || (int)idx == pi1Idx) continue;

            // require electron or muon PID type
            int pidType = particleIDs[idx].type;
            if(pidType != 1 && pidType != 2) continue;

            const auto& rp = recoParticles[idx];
            size_t trackIdx = rp.tracks_begin;
            if(trackIdx >= trackStates.size()) continue;

            double leptonMass = (pidType == 1) ? m_e : m_mu;
            TLorentzVector leptonP4;
            leptonP4.SetXYZM(rp.momentum.x, rp.momentum.y, rp.momentum.z, leptonMass);
            if(leptonP4.Pt() < minLeptonPt) continue;

            // lepton must point approximately in the same direction as the D*
            double dR = leptonP4.DeltaR(dstarP4);
            if(dR > maxDstarLeptonDeltaR) continue;

            // fit a vertex from the D* candidate's 3 tracks plus the lepton's track
            ROOT::VecOps::RVec<edm4hep::TrackState> quadTracks;
            quadTracks.push_back(trackStates[recoParticles[kIdx].tracks_begin]);
            quadTracks.push_back(trackStates[recoParticles[pi2Idx].tracks_begin]);
            quadTracks.push_back(trackStates[recoParticles[pi1Idx].tracks_begin]);
            quadTracks.push_back(trackStates[trackIdx]);
            FCCAnalyses::VertexingUtils::FCCAnalysesVertex bvtx = VertexFitterMod(2, quadTracks);
            double bvtxChi2 = bvtx.vertex.chi2;
            if(bvtxChi2 < 0. || bvtxChi2 > maxVertexChi2Normalized) continue;

            // store the raw candidate
            TLorentzVector visibleP4 = dstarP4 + leptonP4;
            bool isRightSign = (rp.charge * kCharge) < 0;

            RawBCandidate cand;
            cand.dstar_idx = (int)d;
            cand.lepton_idx = (int)idx;
            cand.lepton_pt = leptonP4.Pt();
            cand.lepton_eta = leptonP4.Eta();
            cand.lepton_phi = leptonP4.Phi();
            cand.lepton_charge = (int)rp.charge;
            cand.lepton_type = pidType;
            cand.mass = visibleP4.M();
            cand.pt = visibleP4.Pt();
            cand.eta = visibleP4.Eta();
            cand.phi = visibleP4.Phi();
            cand.dstarlepton_deltaR = dR;
            cand.isRightSign = isRightSign;
            cand.bvtx_chi2Normalized = bvtxChi2;
            rawCandidates.push_back(cand);
        } // end loop over lepton candidates
    } // end loop over D* candidates

    // per-lepton best-candidate selection: a given lepton can be geometrically
    // compatible with several D* candidates found in the same event (most of which
    // would then be spurious for that lepton); keep only the one with the best
    // (lowest) 4-track vertex chi2, since that directly tests whether the 4 tracks
    // are actually consistent with a common vertex. Still allows up to one B
    // candidate per lepton found (e.g. one per b-quark in a Z -> bb event).
    std::unordered_map<int, size_t> bestForLepton; // lepton_idx -> index into rawCandidates
    for(size_t i = 0; i < rawCandidates.size(); i++){
        int lidx = rawCandidates[i].lepton_idx;
        auto it = bestForLepton.find(lidx);
        if(it == bestForLepton.end() ||
           rawCandidates[i].bvtx_chi2Normalized < rawCandidates[it->second].bvtx_chi2Normalized){
            bestForLepton[lidx] = i;
        }
    }

    BCandidates result;
    for(const auto& entry : bestForLepton){
        const RawBCandidate& cand = rawCandidates[entry.second];
        result.dstar_idx.push_back(cand.dstar_idx);
        result.lepton_idx.push_back(cand.lepton_idx);
        result.lepton_pt.push_back(cand.lepton_pt);
        result.lepton_eta.push_back(cand.lepton_eta);
        result.lepton_phi.push_back(cand.lepton_phi);
        result.lepton_charge.push_back(cand.lepton_charge);
        result.lepton_type.push_back(cand.lepton_type);
        result.mass.push_back(cand.mass);
        result.pt.push_back(cand.pt);
        result.eta.push_back(cand.eta);
        result.phi.push_back(cand.phi);
        result.dstarlepton_deltaR.push_back(cand.dstarlepton_deltaR);
        result.isRightSign.push_back(cand.isRightSign);
        result.bvtx_chi2Normalized.push_back(cand.bvtx_chi2Normalized);
        if(result.mass.size() >= maxCandidates) break;
    }

    return result;
}

// "poor man's" gen-level match (see analyzer_dstarfinder.cxx for the same approach and
// its caveats): checks whether there is a gen-level B0/B+ (PDG code 511/521) within dR
// of the *D* candidate's* direction (not the visible D*+lepton system -- the unmeasured
// neutrino makes the visible system's direction a worse proxy for the true B direction
// than the D* alone).
ROOT::VecOps::RVec<bool> matchToGenB(
        const ROOT::VecOps::RVec<int>& dstar_idx,
        const DStarMesonFinder::DStarCandidates& dstarCandidates,
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles,
        double dRThreshold = 0.1){

    ROOT::VecOps::RVec<bool> result;

    // collect 4-vectors of gen-level B0/B+ mesons
    std::vector<TLorentzVector> genBs;
    for(const auto& gp : genParticles){
        int pdg = std::abs(gp.PDG);
        if(pdg != 511 && pdg != 521) continue;
        TLorentzVector p4;
        p4.SetXYZM(gp.momentum.x, gp.momentum.y, gp.momentum.z, gp.mass);
        genBs.push_back(p4);
    }

    for(size_t i = 0; i < dstar_idx.size(); i++){
        int d = dstar_idx[i];
        TLorentzVector dstarP4;
        dstarP4.SetPtEtaPhiM(dstarCandidates.pt[d], dstarCandidates.eta[d], dstarCandidates.phi[d], dstarCandidates.mass[d]);
        bool matched = false;
        for(const auto& genB : genBs){
            if(dstarP4.DeltaR(genB) < dRThreshold){ matched = true; break; }
        }
        result.push_back(matched);
    }
    return result;
}

// same as above, but a dummy for running on data (no gen particles available)
ROOT::VecOps::RVec<bool> matchToGenBDummy(const ROOT::VecOps::RVec<int>& dstar_idx){
    ROOT::VecOps::RVec<bool> result(dstar_idx.size(), false);
    return result;
}

}

#endif
