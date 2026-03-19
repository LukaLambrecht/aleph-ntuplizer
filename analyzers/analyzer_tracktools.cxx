// Tools for dealing with collections of tracks

#include <ROOT/RVec.hxx>

#include "FCCAnalyses/JetConstituentsUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"

#include <iostream>
#include <algorithm>


#ifndef TrackTools_H
#define TrackTools_H

namespace TrackTools{

// get a collection of track states with some track parameters modified
ROOT::VecOps::RVec<edm4hep::TrackState>
getModifiedTrackStates(const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks) {
    ROOT::VecOps::RVec<edm4hep::TrackState> out;
    out.reserve(tracks.size());
    for (const auto &track : tracks) {
        edm4hep::TrackState newtrack = track;

        // flip sign of D0 and omega
        newtrack.D0 = -newtrack.D0;
        newtrack.omega = -newtrack.omega;

        // propagate sign flip of D0 and omega to covariance matrix
        newtrack.covMatrix[1]  = -newtrack.covMatrix[1]; // cov(D0, phi)
        newtrack.covMatrix[4]  = -newtrack.covMatrix[4]; // cov(phi, omega)
        newtrack.covMatrix[6]  = -newtrack.covMatrix[6]; // cov(D0, z0)
        newtrack.covMatrix[8]  = -newtrack.covMatrix[8]; // cov(omega, z0)
        newtrack.covMatrix[10] = -newtrack.covMatrix[10]; // cov(D0, tanLambda)
        newtrack.covMatrix[12] = -newtrack.covMatrix[12]; // cov(omega, tanLambda)

        // express D0, Z0 and omega in mm instead of cm
        // note: disabled again because it is hard to sync with earlier version:
        //       the vertex reconstruction gives slightly different results,
        //       because it depends on both the values and the covariances (must be in consistent units),
        //       but in the earlier version these covariances were not yet modified;
        //       so just to be consistent, keep everything in cm here,
        //       and modify in the code where needed (similar to what was done in earlier version);
        //       can update evertying to a more convenient approach later.
        //newtrack.D0 *= 10;
        //newtrack.Z0 *= 10;
        //newtrack.omega *= 0.1;

        // note: should probably also progagate this to covariance matrix,
        // but not yet done in order to stay in sync with earlier version.

        out.push_back(std::move(newtrack));
    }
    return out;
}

// get the set complement for a collection of tracks
ROOT::VecOps::RVec<edm4hep::TrackState>
get_complement(ROOT::VecOps::RVec<edm4hep::TrackState> allTracks,
                   ROOT::VecOps::RVec<edm4hep::TrackState> selTracks){

    ROOT::VecOps::RVec<edm4hep::TrackState> result;
    for (auto &track : allTracks) {
        bool isSel = false;
        for (auto &selTrack : selTracks) {
            if( track.D0 == selTrack.D0 && track.phi == selTrack.phi && track.omega == selTrack.omega &&
                track.Z0 == selTrack.Z0 && track.tanLambda == selTrack.tanLambda &&
                track.time == selTrack.time && track.referencePoint.x == selTrack.referencePoint.x &&
                track.referencePoint.y == selTrack.referencePoint.y &&
                track.referencePoint.z == selTrack.referencePoint.z ){
                isSel = true;
                break;
            }
        }
        if (!isSel) result.push_back(track);
    }
    return result;
    
}

// same as above but for Track objects instead of TrackStates
ROOT::VecOps::RVec<edm4hep::TrackData>
get_complement(ROOT::VecOps::RVec<edm4hep::TrackData> allTracks,
                   ROOT::VecOps::RVec<edm4hep::TrackState> allTrackStates,
                   ROOT::VecOps::RVec<edm4hep::TrackData> selTracks,
                   ROOT::VecOps::RVec<edm4hep::TrackState> selTrackStates){

    // safety checks
    if( allTracks.size() != allTrackStates.size()) {
        std::string msg = "Vector sizes do not match (in get_complement):";
        msg += " allTracks has size " + std::to_string(allTracks.size());
        msg += " while allTrackStates has size " + std::to_string(allTrackStates.size());
        throw std::runtime_error(msg);
    }
    if( selTracks.size() != selTrackStates.size()) {
        std::string msg = "Vector sizes do not match (in get_complement):";
        msg += " selTracks has size " + std::to_string(selTracks.size());
        msg += " while selTrackStates has size " + std::to_string(selTrackStates.size());
        throw std::runtime_error(msg);
    }

    ROOT::VecOps::RVec<edm4hep::TrackData> result;
    for (int idx = 0; idx < allTracks.size(); idx++) {
        const edm4hep::TrackData track = allTracks.at(idx);
        const edm4hep::TrackState trackState = allTrackStates.at(idx);
        bool isSel = false;
        for (auto &selTrackState : selTrackStates) {
            if( trackState.D0 == selTrackState.D0 && trackState.phi == selTrackState.phi && trackState.omega == selTrackState.omega &&
                trackState.Z0 == selTrackState.Z0 && trackState.tanLambda == selTrackState.tanLambda &&
                trackState.time == selTrackState.time && trackState.referencePoint.x == selTrackState.referencePoint.x &&
                trackState.referencePoint.y == selTrackState.referencePoint.y &&
                trackState.referencePoint.z == selTrackState.referencePoint.z ){
                isSel = true;
                break;
            }
        }
        if (!isSel) result.push_back(track);
    }
    return result;
    
}

// helper function to find tracks passing some baseline selection
// note: this function uses both the track states and the track objects (the latter for chi2);
//       it is assumed that both collections are corresponding trivially, i.e. one-to-one by same index.
ROOT::VecOps::RVec<edm4hep::TrackState> getSelectedTracks(
        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& trackDatas,
        const float D0max, const float Z0max){

        // safety check
        if( tracks.size() != trackDatas.size() ){
            throw std::runtime_error("Vector sizes do not match (in getSelectedTracks plain).");
        }

        // do filtering
        ROOT::VecOps::RVec<edm4hep::TrackState> selectedTracks;
        for (int idx = 0; idx < tracks.size(); idx++) {
            const edm4hep::TrackState trk = tracks.at(idx);
            const edm4hep::TrackData trkobj = trackDatas.at(idx);
            const auto& c = trk.covMatrix;
            if (c[0] <= 0 || c[2] <= 0 || c[9] <= 0) continue;
            if (c[0] < 1e-12 || c[2] < 1e-12 || c[9] <= 1e-12) continue;
            if (!std::isfinite(c[0]) || !std::isfinite(c[2]) || !std::isfinite(c[9])) continue;
            if (D0max > 0 && std::abs(trk.D0)>D0max) continue;
            if (Z0max > 0 && std::abs(trk.Z0)>Z0max) continue;
            if (trkobj.ndf==0) continue;
            if (trkobj.chi2 / trkobj.ndf > 10.) continue;
            selectedTracks.push_back(trk);
        }
        return selectedTracks;
}

// same as above but with tracks per jet
ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> getSelectedTracks(
        const ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>>& tracksPerJet,
        const ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackData>>& trackDatasPerJet,
        const float D0max, const float Z0max){

        // safety check
        if( tracksPerJet.size() != trackDatasPerJet.size() ){
            throw std::runtime_error("Vector sizes do not match (in getSelectedTracks per jet).");
        }

        // initialization
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> selectedTracks;

        // loop over jets
        for (int idx = 0; idx < tracksPerJet.size(); idx++) {
            ROOT::VecOps::RVec<edm4hep::TrackState> tracks = tracksPerJet.at(idx);
            ROOT::VecOps::RVec<edm4hep::TrackData> trackDatas = trackDatasPerJet.at(idx);
            selectedTracks.push_back( getSelectedTracks(tracks, trackDatas, D0max, Z0max) );
        }
        return selectedTracks;
}

// helper function to distribute tracks without reco particle over tracks per jet
ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>>
append_losttracks_to_tracksperjet(
        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracksWithoutReco,
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> tracksWithRecoPerJet,
        const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets){

        // safety checks
        if (jets.size()==0){ return tracksWithRecoPerJet; }
        if (tracksWithRecoPerJet.size() != jets.size()) {
            throw std::runtime_error("Vector sizes do not match (in append_losttracks_to_tracksperjet for TrackStates).");
        }

        // make a vector for each jet
        ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector> jetP4s;
        for (const auto& jet : jets) {
            jetP4s.emplace_back(jet.pt(), jet.eta(), jet.phi(), jet.m());
        }

        // loop over tracks to append
        for(const edm4hep::TrackState& track : tracksWithoutReco){
            float phi_track = track.phi;
            float eta_track = std::asinh(track.tanLambda);
            ROOT::Math::PtEtaPhiMVector trackP4(
                1.0,       // dummy pT
                eta_track,
                phi_track,
                0.0
            );

            // find index of closest jet
            float minDR = 99.;
            int bestJetIdx = 0;
            for (size_t i = 0; i < jetP4s.size(); ++i) {
                ROOT::Math::PtEtaPhiMVector jetP4 = jetP4s.at(i);
                const float dR = ROOT::Math::VectorUtil::DeltaR(trackP4, jetP4);
                if (dR < minDR){ minDR = dR; bestJetIdx = i; }
            }

            // append
            if (bestJetIdx >= 0){ tracksWithRecoPerJet.at(bestJetIdx).push_back(track); }
        } // end of loop over tracks

        return tracksWithRecoPerJet;
}

// same as above but for tracks rather than track states...
ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackData>>
append_losttracks_to_tracksperjet(
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracksWithoutReco,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStatesWithoutReco,
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackData>> tracksWithRecoPerJet,
        const ROOT::VecOps::RVec<fastjet::PseudoJet>& jets){

        // safety checks
        if (jets.size()==0){ return tracksWithRecoPerJet; }
        if (tracksWithRecoPerJet.size() != jets.size()) {
            std::string msg = "Vector sizes do not match (in append_losttracks_to_tracksperjet for TrackData):";
            msg += " tracksWithRecoPerJet has size " + std::to_string(tracksWithRecoPerJet.size());
            msg += " while jets has size " + std::to_string(jets.size());
            throw std::runtime_error(msg);
        }
        if (tracksWithoutReco.size() != trackStatesWithoutReco.size()) {
            std::string msg = "Vector sizes do not match (in append_losttracks_to_tracksperjet for TrackData):";
            msg += " tracksWithoutReco has size " + std::to_string(tracksWithoutReco.size());
            msg += " while trackStatesWithoutReco has size " + std::to_string(trackStatesWithoutReco.size());
            throw std::runtime_error(msg);
        }

        // make a vector for each jet
        ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector> jetP4s;
        for (const auto& jet : jets) {
            jetP4s.emplace_back(jet.pt(), jet.eta(), jet.phi(), jet.m());
        }

        // loop over tracks to append
        for(int idx = 0; idx < tracksWithoutReco.size(); idx++){
            const edm4hep::TrackData track = tracksWithoutReco.at(idx);
            const edm4hep::TrackState trackState = trackStatesWithoutReco.at(idx);
            float phi_track = trackState.phi;
            float eta_track = std::asinh(trackState.tanLambda);
            ROOT::Math::PtEtaPhiMVector trackP4(
                1.0,       // dummy pT
                eta_track,
                phi_track,
                0.0
            );

            // find index of closest jet
            float minDR = 99.;
            int bestJetIdx = 0;
            for (size_t i = 0; i < jetP4s.size(); ++i) {
                ROOT::Math::PtEtaPhiMVector jetP4 = jetP4s.at(i);
                const float dR = ROOT::Math::VectorUtil::DeltaR(trackP4, jetP4);
                if (dR < minDR){ minDR = dR; bestJetIdx = i; }
            }

            // append
            if (bestJetIdx >= 0){ tracksWithRecoPerJet.at(bestJetIdx).push_back(track); }
        } // end of loop over tracks

        return tracksWithRecoPerJet;
}

// helper function to find tracks not compatible with primary vertex.
ROOT::VecOps::RVec<edm4hep::TrackState> getSecondaryTracks(
        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& primaryTracks){

        // skip tracks compatible with primary vertex
        ROOT::VecOps::RVec<edm4hep::TrackState> secondaryTracks;
        secondaryTracks = FCCAnalyses::VertexFitterSimple::get_NonPrimaryTracks(tracks, primaryTracks);
        return secondaryTracks;
}

// same as above but with tracks per jet
ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> getSecondaryTracks(
        const ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>>& tracksPerJet,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& primaryTracks){

        // initialization
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> secondaryTracks;

        // loop over jets
        for ( const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks : tracksPerJet) {
            secondaryTracks.push_back( getSecondaryTracks(tracks, primaryTracks) );
        }
        return secondaryTracks;
}

// helper function to count the number of tracks per jet
ROOT::VecOps::RVec<int> countTracks(
        const ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>>& tracks){
        ROOT::VecOps::RVec<int> nTracks;
        for( const ROOT::VecOps::RVec<edm4hep::TrackState>& this_tracks : tracks) {
            nTracks.push_back( this_tracks.size() );
        }
        return nTracks;
}

}

#endif
