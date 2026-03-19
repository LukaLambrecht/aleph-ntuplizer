// Tools for dealing with links between reconstructed particles and tracks

#include <ROOT/RVec.hxx>

#include "FCCAnalyses/JetConstituentsUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"

#include <iostream>
#include <algorithm>

#ifndef RecoToTrackTools_H
#define RecoToTrackTools_H

namespace RecoToTrackTools{

  
  /*
  General getters
  */

  // note: no such function in FCCAnalyses ReconstructedParticle2Track.
  //       should this be added?
  ROOT::VecOps::RVec<edm4hep::TrackState> getRP2TRK_trackState(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
    const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks){
    /*
    Get the full TrackState object for each input particle.
    Note: contrary to most other functions of this template,
          no dummy value is returned for neutral particles;
          hence the output vector contains all the tracks of the set of input particles,
          but the length might be different and there is no one-to-one correspondence!
          this function is useful e.g. for getting all tracks of a given jet
          (or other collection of particles) but not for extracting per-particle properties.
    */

    // initializations
    ROOT::VecOps::RVec<edm4hep::TrackState> result;

    // loop over reco particles
    for (auto & p: rps) {
        // find track
        size_t trackIndex = p.tracks_begin;
        if(trackIndex < tracks.size()){
            edm4hep::TrackState tr = tracks.at(trackIndex);
            result.push_back(tr);
        }
    }
    return result;
  }


  // note: no such function in FCCAnalyses ReconstructedParticle2Track.
  //       should this be added?
  ROOT::VecOps::RVec<edm4hep::TrackData> getRP2TRK_track(
    const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
    const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks){
    /*
    Get the full Track object for each input particle.
    Note: contrary to most other functions of this template,
          no dummy value is returned for neutral particles;
          hence the output vector contains all the tracks of the set of input particles,
          but the length might be different and there is no one-to-one correspondence!
          this function is useful e.g. for getting all tracks of a given jet
          (or other collection of particles) but not for extracting per-particle properties.
    */

    // initializations
    ROOT::VecOps::RVec<edm4hep::TrackData> result;

    // loop over reco particles
    for (auto & p: rps) {
        // find track
        size_t trackIndex = p.tracks_begin;
        if(trackIndex < tracks.size()){
            edm4hep::TrackData tr = tracks.at(trackIndex);
            result.push_back(tr);
        }
    }
    return result;
  }

  /*
  Get chi2 and related variables
  */

  ROOT::VecOps::RVec<float>
  getRP2TRK_chi2(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
                   const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks){
    // get chi2 of track fit
    
    // initializations
    ROOT::VecOps::RVec<float> result;

    // loop over reco particles
    for (auto & p: in) {
        bool valid = false;
        // find track
        size_t trackIndex = p.tracks_begin;
        if(trackIndex < tracks.size()){
            edm4hep::TrackData tr = tracks.at(trackIndex);
            float chi2 = tr.chi2;
            result.push_back(chi2);
            valid = true;
        }
        if(!valid){ result.push_back( -1 ); }
    }
    return result;
  }

  ROOT::VecOps::RVec<int>
  getRP2TRK_ndof(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
                   const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks){
    // get number of degrees of freedom of track fit

    // initializations
    ROOT::VecOps::RVec<int> result;

    // loop over reco particles
    for (auto & p: in) {
        bool valid = false;
        // find track
        size_t trackIndex = p.tracks_begin;
        if(trackIndex < tracks.size()){
            edm4hep::TrackData tr = tracks.at(trackIndex);
            int ndof = tr.ndf;
            result.push_back(ndof);
            valid = true;
        }
        if(!valid){ result.push_back( -1 ); }
    }
    return result;
  }

  ROOT::VecOps::RVec<float>
  getRP2TRK_chi2Normalized(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
                   const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks){
    // get chi2 of track fit

    // initializations
    ROOT::VecOps::RVec<float> result;

    // loop over reco particles
    for (auto & p: in) {
        bool valid = false;
        // find track
        size_t trackIndex = p.tracks_begin;
        if(trackIndex < tracks.size()){
            edm4hep::TrackData tr = tracks.at(trackIndex);
            float chi2 = tr.chi2;
            int ndof = tr.ndf;
            if(ndof > 0){
                result.push_back(chi2 / (float)ndof);
                valid = true;
            }
        }
        if(!valid){ result.push_back( -1 ); }
    }
    return result;
  }

  /*
  Get number of track hits in each tracker subdetector
  */

  ROOT::VecOps::RVec<int> getRP2TRK_nTrackHits(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<int>& subdetectorHitNumbers,
        const int subdetectorNumber){
    // Get number of track hits for a given set of reco particles
    // for a given subdetector number

    // initializations
    ROOT::VecOps::RVec<float> out;

    // loop over reco particles
    for(auto & p: rps) {
        bool valid = false;
        // find track
        size_t trackIndex = p.tracks_begin;
        if(trackIndex < tracks.size()){
            edm4hep::TrackData tr = tracks.at(trackIndex);
            // find index in hit collection
            size_t hitIdx = tr.subdetectorHitNumbers_begin + subdetectorNumber;
            if(hitIdx < subdetectorHitNumbers.size()){
                int nHits = subdetectorHitNumbers.at(hitIdx);
                out.push_back(nHits);
                valid = true;
            }
        }
        if(!valid){ out.push_back( -1 ); }
    }
    return out;
  }

  ROOT::VecOps::RVec<int> getRP2TRK_nTrackHits_VDET(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<int>& subdetectorHitNumbers){
    // Get number of VDET (vertex detector) hits for a given set of reco particles.
    // Note: assumes inside-out numbering of the subdetectors,
    // such that VDET is at index 0 (to check).
    return getRP2TRK_nTrackHits(rps, tracks, subdetectorHitNumbers, 0);
  }

  ROOT::VecOps::RVec<int> getRP2TRK_nTrackHits_ITC(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<int>& subdetectorHitNumbers){
    // Get number of ITC (inner tracking chamber) hits for a given set of reco particles.
    // Note: assumes inside-out numbering of the subdetectors,
    // such that ITC is at index 1 (to check).
    return getRP2TRK_nTrackHits(rps, tracks, subdetectorHitNumbers, 1);
  }

  ROOT::VecOps::RVec<int> getRP2TRK_nTrackHits_TPC(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const ROOT::VecOps::RVec<int>& subdetectorHitNumbers){
    // Get number of TPC (time projection chamber) hits for a given set of reco particles.
    // Note: assumes inside-out numbering of the subdetectors,
    // such that TPC is at index 2 (to check).
    return getRP2TRK_nTrackHits(rps, tracks, subdetectorHitNumbers, 2);
  }

  /*
  Magnetic field
  */

  ROOT::VecOps::RVec<float> getRP2TRK_Bz(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates) {
    // Get the magnetic field strength from the momentum of a particle
    // and the curvature of its associated track.
    // Based on the formula p[GeV/c] = c[= 0.3 in these units] B[T] r[m]
    // with p the momentum, c the speed of light, B the magnetic field strength,
    // and r the radius of curvature.

    // initializations
    double cSpeed = 2.99792458e8; // speed of light in m/s
    cSpeed *= 1e-9 * 1e-3; // conversion to appropriate units
    // (see formula above, and also accounting for r in mm instead of in m).
    ROOT::VecOps::RVec<float> out;

    // loop over reco particles
    for(auto & p: rps) {
        bool valid = false;
        size_t trackIndex = p.tracks_begin;
        if(trackIndex < trackStates.size()) {
            edm4hep::TrackState trst = trackStates.at(trackIndex);
            double pt = sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y);
            double omega = trst.omega;
            omega *= 0.1;
            // extra numerical factor needed for Aleph data
            // (because of the curvature being expressed in 1/cm instead of 1/mm)
            double Bz = omega / cSpeed * pt * std::copysign(1.0, p.charge);
            out.push_back(Bz);
            valid = true;
        }
        if(!valid){ out.push_back(-9.); }
    }
    return out;
  }

  float Bz(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
           const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates) {
    // Get the magnetic field strength for an event
    // from the momenta of the particles and the curvatures of their associated tracks.
    // Based on the formula p[GeV/c] = c[= 0.3 in these units] B[T] r[m]
    // with p the momentum, c the speed of light, B the magnetic field strength,
    // and r the radius of curvature.

    // initializations
    double cSpeed = 2.99792458e8; // speed of light in m/s
    cSpeed *= 1e-9 * 1e-3; // conversion to appropriate units
    // (see formula above, and also accounting for r in mm instead of in m).
    double Bz = -9; // dummy value if no valid result found

    // loop over reco particles
    for(auto & p: rps) {
      size_t trackIndex = p.tracks_begin;
      if(trackIndex >= trackStates.size()){ continue; }
      edm4hep::TrackState trst = trackStates.at(trackIndex);
      double pt = sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y);
      double omega = trst.omega;
      omega *= 0.1;
      // extra numerical factor needed for Aleph data
      // (because of the curvature being expressed in 1/cm instead of 1/mm)
      Bz = omega / cSpeed * pt * std::copysign(1.0, p.charge);
      // break after the first particle with valid result
      break;
    }
    return Bz;
  }

  /*
  Impact parameters
  */

  ROOT::VecOps::RVec<float> XPtoPar_dxy(
      const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
      const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
      const TLorentzVector& PV,
      const float& Bz) {
    // Get dxy with respect to primary vertex,
    // given the track state parameters calculated with respect to the origin.

    // initializations
    double cSpeed = 2.99792458e8; // speed of light in m/s
    cSpeed *= 1e-9 * 1e-3; // conversion to appropriate units (see Bz calculation)
    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {
      size_t trackIndex = rp.tracks_begin;
      if(trackIndex < trackStates.size()){
        edm4hep::TrackState trst = trackStates.at(trackIndex);

        float D0_wrt0 = trst.D0;
        float Z0_wrt0 = trst.Z0;
        float phi0_wrt0 = trst.phi;

        // note: phi0 is not the position vector azimuth,
        // but the azimuth of the momentum vector at the point of closest approach!
        TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
        TVector3 x = X - PV.Vect();
        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);
        double pt = p.Pt();

        double a = - rp.charge * Bz * cSpeed;
        // extra numerical factor for aleph, see Bz calculation
        a *= -10;
        double r2 = x(0) * x(0) + x(1) * x(1);
        double cross = x(0) * p(1) - x(1) * p(0);
        double D=-9;
        if (pt * pt - 2 * a * cross + a * a * r2 > 0) {
          double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
          if (pt < 10.0) D = (T - pt) / a;
          else D = (-2 * cross + a * r2) / (T + pt);
        }
        out.push_back(D);
      } else { out.push_back(-9.); }
    }
    return out;
  }

  ROOT::VecOps::RVec<float> XPtoPar_dz(
      const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
      const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
      const TLorentzVector& PV,
      const float& Bz) {
    // Get dz with respect to primary vertex,
    // given the track state parameters calculated with respect to the origin.

    // initializations
    double cSpeed = 2.99792458e8; // speed of light in m/s
    cSpeed *= 1e-9 * 1e-3; // conversion to appropriate units
    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {
      size_t trackIndex = rp.tracks_begin;
      if(trackIndex < trackStates.size()){
        edm4hep::TrackState trst = trackStates.at(trackIndex);

        float D0_wrt0 = trst.D0;
        float Z0_wrt0 = trst.Z0;
        float phi0_wrt0 = trst.phi;

        // note: phi0 is not the position vector azimuth,
        // but the azimuth of the momentum vector at the point of closest approach!
        TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
        TVector3 x = X - PV.Vect();
        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);
        double pt = p.Pt();

        double a = - rp.charge * Bz * cSpeed;
        // extra numerical factor for aleph, see Bz calculation
        a *= -10;
        double C = a/(2 * pt);
        double r2 = x(0) * x(0) + x(1) * x(1);
        double cross = x(0) * p(1) - x(1) * p(0);
        double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
        double D;
        if (pt < 10.0) D = (T - pt) / a;
        else D = (-2 * cross + a * r2) / (T + pt);
        double B = C * TMath::Sqrt(TMath::Max(r2 - D * D, 0.0) / (1 + 2 * C * D));
        if ( TMath::Abs(B) > 1.) B = TMath::Sign(1, B);
        double st = TMath::ASin(B) / C;
        double ct = p(2) / pt;
        double z0;
        double dot = x(0) * p(0) + x(1) * p(1);
        if (dot > 0.0) z0 = x(2) - ct * st;
        else z0 = x(2) + ct * st;

        out.push_back(z0);
      } else { out.push_back(-9.); }
    }
    return out;
  }

  ROOT::VecOps::RVec<float> XPtoPar_phi(
      const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& in,
      const ROOT::VecOps::RVec<edm4hep::TrackState>& trackStates,
      const TLorentzVector& PV,
      const float& Bz) {
    // Get phi with respect to primary vertex,
    // given the track state parameters calculated with respect to the origin.

    // initializations
    double cSpeed = 2.99792458e8; // speed of light in m/s
    cSpeed *= 1e-9 * 1e-3; // conversion to appropriate units
    ROOT::VecOps::RVec<float> out;

    for (const auto & rp: in) {
      size_t trackIndex = rp.tracks_begin;
      if(trackIndex < trackStates.size()){
        edm4hep::TrackState trst = trackStates.at(trackIndex);

        float D0_wrt0 = trst.D0;
        float Z0_wrt0 = trst.Z0;
        float phi0_wrt0 = trst.phi;

        // note: phi0 is not the position vector azimuth,
        // but the azimuth of the momentum vector at the point of closest approach!
        TVector3 X( - D0_wrt0 * TMath::Sin(phi0_wrt0) , D0_wrt0 * TMath::Cos(phi0_wrt0) , Z0_wrt0);
        TVector3 x = X - PV.Vect();
        TVector3 p(rp.momentum.x, rp.momentum.y, rp.momentum.z);
        double pt = p.Pt();

        double a = - rp.charge * Bz * cSpeed;
        // extra numerical factor for aleph, see Bz calculation
        a *= -10;
        double r2 = x(0) * x(0) + x(1) * x(1);
        double cross = x(0) * p(1) - x(1) * p(0);
        double T = TMath::Sqrt(pt * pt - 2 * a * cross + a * a * r2);
        double phi0 = TMath::ATan2((p(1) - a * x(0)) / T, (p(0) + a * x(1)) / T);

        out.push_back(phi0);

      } else { out.push_back(-9.); }
    }
    return out;
  }

}

#endif
