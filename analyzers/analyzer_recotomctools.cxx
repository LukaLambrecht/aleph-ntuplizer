// Tools for dealing with gen particles

#include <ROOT/RVec.hxx>

#include "FCCAnalyses/JetConstituentsUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"

#include <iostream>
#include <algorithm>


#ifndef RecoToMCTools_H
#define RecoToMCTools_H

namespace RecoToMCTools{

// Helper function to make track to gen particle mapping
ROOT::VecOps::RVec<int> makeTrackToMCMapping(
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<podio::ObjectID>& tracktomc_links,
        const ROOT::VecOps::RVec<podio::ObjectID>& mctotrack_links){

    // initialize
    ROOT::VecOps::RVec<int> trackToMCMap;
    for(int idx = 0; idx < tracks.size(); idx++){ trackToMCMap.push_back(9999); }

    // check
    if( tracktomc_links.size() != mctotrack_links.size() ){
        std::string msg = "These vectors should have the same size";
        throw std::runtime_error(msg);
    }

    // loop over index pairs
    for(int idx = 0; idx < tracktomc_links.size(); idx++){
        int track_idx = mctotrack_links.at(idx).index;
        int mc_idx = tracktomc_links.at(idx).index;
        trackToMCMap.at(track_idx) = mc_idx;
    }
    return trackToMCMap;
}

// Same as above, but the other way round
ROOT::VecOps::RVec<int> makeMCToTrackMapping(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles,
        const ROOT::VecOps::RVec<podio::ObjectID>& tracktomc_links,
        const ROOT::VecOps::RVec<podio::ObjectID>& mctotrack_links){

    // initialize
    ROOT::VecOps::RVec<int> MCToTrackMap;
    for(int idx = 0; idx < genParticles.size(); idx++){ MCToTrackMap.push_back(9999); }

    // check
    if( tracktomc_links.size() != mctotrack_links.size() ){
        std::string msg = "These vectors should have the same size";
        throw std::runtime_error(msg);
    }

    // loop over index pairs
    for(int idx = 0; idx < tracktomc_links.size(); idx++){
        int track_idx = mctotrack_links.at(idx).index;
        int mc_idx = tracktomc_links.at(idx).index;
        MCToTrackMap.at(mc_idx) = track_idx;
    }
    return MCToTrackMap;
}

// Get MC particle index for each reconstructed particle
// Note: for ALEPH data, there exists a direct link between tracks and MC particles
//       but neutral reco particles do not seem to have something like that,
//       so these are given a dummy value.
ROOT::VecOps::RVec<int> get_genParticleIndex(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles,
        const ROOT::VecOps::RVec<int>& trackToMCMap){

    // initializations
    ROOT::VecOps::RVec<int> result;

    // loop over reco particles
    for(auto & rp: rps){

        // find track
        size_t trackIndex = rp.tracks_begin;
        if(trackIndex < tracks.size()){
            // find gen particle index
            size_t genIndex = trackToMCMap.at(trackIndex);
            if(genIndex < genParticles.size()){
                result.push_back(genIndex);
            }
            else{ result.push_back(-1); }
        }
        else{ result.push_back(-1); }
    }
    return result;
}


// Get MC particle PDG ID for each reconstructed particle
// Note: for ALEPH data, there exists a direct link between tracks and MC particles
//       but neutral reco particles do not seem to have something like that,
//       so these are given a dummy value.
ROOT::VecOps::RVec<int> get_pdgid(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles,
        const ROOT::VecOps::RVec<int>& trackToMCMap){

    // initializations
    ROOT::VecOps::RVec<int> result;

    // get gen particle indices for each reco particle
    ROOT::VecOps::RVec<int> genParticleIndices = get_genParticleIndex(
        rps,
        tracks,
        genParticles,
        trackToMCMap
    );

    // loop over reco particles
    for(unsigned int reco_idx = 0; reco_idx < rps.size(); reco_idx++){
        int gen_idx = genParticleIndices.at(reco_idx);
        if(gen_idx >= 0 && gen_idx < genParticles.size()){
            const edm4hep::MCParticleData genParticle = genParticles.at(gen_idx);
            result.push_back(genParticle.PDG);
        }
        else{ result.push_back(0); }
    }
    return result;
}

// Same as above, but per jet
ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> get_pdgid(
        const ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>>& jcs,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles,
        const ROOT::VecOps::RVec<int>& trackToMCMap){

    // loop over jets
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> out;
    for (const auto &jc : jcs){
        out.emplace_back( get_pdgid(jc, tracks, genParticles, trackToMCMap) );
    }
    return out;
}

// Same as above, but dummy for data
ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> get_pdgidDummy(
        const ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>>& jcs){

    // loop over jets and jet constituents
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> out;
    for (const auto &jc : jcs){
        ROOT::VecOps::RVec<int> temp;
        for(int idx = 0; idx < jc.size(); idx++){ temp.push_back(0); }
        out.push_back(temp);
    }
    return out;
}

}

#endif
