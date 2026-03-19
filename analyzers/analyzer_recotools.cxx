// Tools for dealing with collections of reco particles

#include <ROOT/RVec.hxx>

#include "FCCAnalyses/JetConstituentsUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"

#include <iostream>
#include <algorithm>


#ifndef RecoTools_H
#define RecoTools_H

namespace RecoTools{

// helper function
size_t getTrackIndex(
        const edm4hep::ReconstructedParticleData& rp,
        const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){
    /*
    Get index of track in track collection corrsponding to a reco particle.
    Note: depending on the edm4hep version, this is NOT simply RecoParticle.tracks_begin...
    */
    if(rp.tracks_begin < reco2track_links.size()) {
      const auto &oid = reco2track_links.at(rp.tracks_begin);
      size_t trackIndex = oid.index;
      return trackIndex;
    }
    return 9999;
}

// get a collection of reco particles with some parameters modified
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
getModifiedRecoParticles(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<podio::ObjectID>& reco2track_links){

    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> out;
    out.reserve(rps.size());
    for (const auto &rp : rps) {
        edm4hep::ReconstructedParticleData newrp = rp;

        // make the link between reco particle and track direct,
        // (as in newer versions of EDM4HEP (?))
        // instead of having to go through the intermediate collection of linking indices every time.
        newrp.tracks_begin = getTrackIndex(rp, reco2track_links);
        newrp.tracks_end = newrp.tracks_begin + 1;

        out.push_back(std::move(newrp));
    }
    return out;
}

}

#endif
