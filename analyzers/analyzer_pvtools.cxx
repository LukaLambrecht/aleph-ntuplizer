// Tools for dealing with primary vertex reconstruction

#include <ROOT/RVec.hxx>

#include "FCCAnalyses/JetConstituentsUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"

#include <iostream>
#include <algorithm>


#ifndef PrimaryVertexTools_H
#define PrimaryVertexTools_H

namespace PrimaryVertexTools{

// helper function to find the primary tracks.
// note: copied from here:
//       https://github.com/HEP-FCC/FCCAnalyses/blob/dabdfaf92f716373dc000509766bdd827dc07392/analyzers/dataframe/src/VertexFitterSimple.cc#L259
//       with the following modifications:
//       - give chi2max as an argument instead of hard-coding it to the value of 25.
//       - change the conversion factor from micrometer to millimeter from 1e-6 (bug?) to 1e-3.
//       maybe later propagate back to main FCCAnalyses...
ROOT::VecOps::RVec<edm4hep::TrackState> get_PrimaryTracks(
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
        double chi2max,
        bool BeamSpotConstraint,
        double bsc_sigmax, double bsc_sigmay, double bsc_sigmaz,
        double bsc_x, double bsc_y, double bsc_z) {
  
  // check if sufficient tracks were provided
  ROOT::VecOps::RVec<edm4hep::TrackState> seltracks = tracks;
  if (seltracks.size() <= 1) return seltracks;
  int Ntr = tracks.size();

  // parse track parameters
  TVectorD **trkPar = new TVectorD *[Ntr];
  TMatrixDSym **trkCov = new TMatrixDSym *[Ntr];
  for (Int_t i = 0; i < Ntr; i++) {
    edm4hep::TrackState t = tracks[i];
    TVectorD par = FCCAnalyses::VertexingUtils::get_trackParam(t);
    trkPar[i] = new TVectorD(par);
    TMatrixDSym Cov = FCCAnalyses::VertexingUtils::get_trackCov(t);
    trkCov[i] = new TMatrixDSym(Cov);
  }

  // initialize the vertex fit
  VertexFit theVertexFit(Ntr, trkPar, trkCov);

  // add beamspot constraint
  if (BeamSpotConstraint) {
    float conv_BSC = 1e-3;
    // (conversion factor from micrometer to millimeter;
    // the beamspot positions and widths are expected to be given in micrometers,
    // while the track- and vertex properties are typically in millimeter.)
    TVectorD xv_BS(3);
    xv_BS[0] = bsc_x * conv_BSC;
    xv_BS[1] = bsc_y * conv_BSC;
    xv_BS[2] = bsc_z * conv_BSC;
    TMatrixDSym cov_BS(3);
    cov_BS[0][0] = pow(bsc_sigmax * conv_BSC, 2);
    cov_BS[1][1] = pow(bsc_sigmay * conv_BSC, 2);
    cov_BS[2][2] = pow(bsc_sigmaz * conv_BSC, 2);
    theVertexFit.AddVtxConstraint(xv_BS, cov_BS);
  }

  // run the fit
  TVectorD x = theVertexFit.GetVtx();

  // temp: printouts for testing
  /*std::cout << "---" << std::endl;
  TVectorD trackschi2 = theVertexFit.GetVtxChi2List();
  for (int i = 0; i < theVertexFit.GetNtrk(); i++) {
      float trackchi2 = trackschi2[i];
      std::cout << trackchi2 << std::endl;
  }*/

  // iteratively remove tracks until the chi2 is below threshold
  float chi2_max = 1e30;
  while (chi2_max >= chi2max) {

    // get maximum chi2 of all tracks currently in the fit
    TVectorD tracks_chi2 = theVertexFit.GetVtxChi2List();
    chi2_max = tracks_chi2.Max();

    // remove track with highest chi2
    int n_removed = 0;
    for (int i = 0; i < theVertexFit.GetNtrk(); i++) {
      float track_chi2 = tracks_chi2[i];
      if (track_chi2 >= chi2_max) {
        theVertexFit.RemoveTrk(i);
        seltracks.erase(seltracks.begin() + i);
        n_removed++;
      }
    }

    // run the fit again
    if (n_removed > 0) {
      if (theVertexFit.GetNtrk() > 1) {
        x = theVertexFit.GetVtx();
        TVectorD new_tracks_chi2 = theVertexFit.GetVtxChi2List();
        chi2_max = new_tracks_chi2.Max();
      } else {
        chi2_max = 0; // exit from the loop w/o crashing..
      }
    }
  } // end while

  // memory cleanup :
  for (Int_t i = 0; i < Ntr; i++) {
    delete trkPar[i];
    delete trkCov[i];
  }
  delete[] trkPar;
  delete[] trkCov;

  return seltracks;
}


// helper function to find primary tracks.
// note: need to keep in sync with fitRecoPrimaryVertex (below)...
ROOT::VecOps::RVec<edm4hep::TrackState> getPrimaryTracks(
        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
        double chi2max = 25.,
        double beamspotX = 0., double beamspotY = 0., double beamspotZ = 0.){

        // convert beamspot position units from centimeter to 10 micrometer
        // note: the FCCAnalyses function expect these values in micrometer,
        //       corresponding to track parameters in mm.
        //       since our track parameters are in cm, we use beamspot units of 10 micrometer.
        beamspotX = beamspotX * 1e3;
        beamspotY = beamspotY * 1e3;
        beamspotZ = beamspotZ * 1e3;

        // define beamspot width
        double sigma_beamspotX = 20; // unit: 10 micrometer
        double sigma_beamspotY = 10; // unit: 10 micrometer
        double sigma_beamspotZ = 2000; // unit: 10 micrometer
        bool doBeamSpotConstraint = true;

        // intitialize output
        ROOT::VecOps::RVec<edm4hep::TrackState> primaryTracks;

        // do filtering
        // note: the input is assumed to have already passed the baseline selection;
        //       here we just apply an extra cut on D0 and Z0 to focus on primary tracks
        ROOT::VecOps::RVec<edm4hep::TrackState> tracksToUse;
        for (const edm4hep::TrackState& trk : tracks) {
            if (std::abs(trk.D0)>0.75 || std::abs(trk.Z0)>2) continue;
            tracksToUse.push_back(trk);
        }
        if( tracksToUse.size() < 2 ){ return primaryTracks; }

        // call primary track finder
        primaryTracks = get_PrimaryTracks(
            tracksToUse,
            chi2max,
            doBeamSpotConstraint,
            sigma_beamspotX, sigma_beamspotY, sigma_beamspotZ,
            beamspotX, beamspotY, beamspotZ
        );
        return primaryTracks;
}


// helper function to re-calculate the primary vertex from the collection of tracks.
// note: no track selection is performed, this is assumed to be done beforehand.
// note: need to keep in sync with getPrimaryTracks (above)...
FCCAnalyses::VertexingUtils::FCCAnalysesVertex fitRecoPrimaryVertex(
        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
        double beamspotX = 0, double beamspotY = 0, double beamspotZ = 0){

        // convert beamspot position units from centimeter to micrometer
        // note: the FCCAnalyses function expect these values in micrometer,
        //       corresponding to track parameters in mm (and current beamspot parameters in cm).
        //       since our track parameters are in cm, we use beamspot units of 10 micrometer.
        beamspotX = beamspotX * 1e3;
        beamspotY = beamspotY * 1e3;
        beamspotZ = beamspotZ * 1e3;

        // define beamspot width
        double sigma_beamspotX = 20; // unit: 10 micrometer
        double sigma_beamspotY = 10; // unit: 10 micrometer
        double sigma_beamspotZ = 2000; // unit: 10 micrometer
        bool doBeamSpotConstraint = true;

        // define dummy vertex in case the fit cannot be performed
        edm4hep::VertexData dummyVertex;
        dummyVertex.chi2 = -1;
        dummyVertex.ndf = 0;
        dummyVertex.position = edm4hep::Vector3f(beamspotX, beamspotY, beamspotZ);
        FCCAnalyses::VertexingUtils::FCCAnalysesVertex dummyVertexObject;
        dummyVertexObject.vertex = dummyVertex;
        dummyVertexObject.ntracks = 0;
        dummyVertexObject.mc_ind = -1;
        if( tracks.size() < 2 ){ return dummyVertexObject; }

        // call primary vertex finder from FCCAnalyses
        FCCAnalyses::VertexingUtils::FCCAnalysesVertex vertex;
        vertex = FCCAnalyses::VertexFitterSimple::VertexFitter_Tk(
            1, tracks, doBeamSpotConstraint,
            sigma_beamspotX, sigma_beamspotY, sigma_beamspotZ,
            beamspotX, beamspotY, beamspotZ
        );
        return vertex;
}

}

#endif
