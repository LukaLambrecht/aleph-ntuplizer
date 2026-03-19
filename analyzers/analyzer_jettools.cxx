// Tools for dealing with collections of jets

#include <ROOT/RVec.hxx>

#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"

#include <iostream>
#include <algorithm>

#include "analyzer_recototracktools.cxx"


#ifndef JetTools_H
#define JetTools_H

namespace JetTools{

    using FCCAnalysesJetConstituents = ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>;
    using FCCAnalysesJetConstituentsData = ROOT::VecOps::RVec<float>;

    // note: this function is not present in FCCAnalyses / JetConstituentUtils.
    //       should it be added?
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> build_trackstates_cluster(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackState>& tracks,
        const std::vector<std::vector<int>>& jet_indices){
        /*
        Build the collection of track states (mapping jet -> track states of (charged) reconstructed particles)
        */
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackState>> tracks_perjet;
        // loop over jets
        for (const auto &this_jet_indices : jet_indices){
            // get constituents for this jet
            FCCAnalysesJetConstituents constituents;
            for (const auto &constituent_index : this_jet_indices){
                constituents.push_back(rps.at(constituent_index));
            }
            // get track states for these constituents
            ROOT::VecOps::RVec<edm4hep::TrackState> this_jet_tracks;
            this_jet_tracks = RecoToTrackTools::getRP2TRK_trackState(constituents, tracks);
            tracks_perjet.push_back(this_jet_tracks);
        }
        return tracks_perjet;
    }

    // note: this function is not present in FCCAnalyses / JetConstituentUtils.
    //       should it be added?
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackData>> build_tracks_cluster(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps,
        const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
        const std::vector<std::vector<int>>& jet_indices){
        /*
        Build the collection of tracks (mapping jet -> tracks of (charged) reconstructed particles)
        */
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::TrackData>> tracks_perjet;
        // loop over jets
        for (const auto &this_jet_indices : jet_indices){
            // get constituents for this jet
            FCCAnalysesJetConstituents constituents;
            for (const auto &constituent_index : this_jet_indices){
                constituents.push_back(rps.at(constituent_index));
            }
            // get track states for these constituents
            ROOT::VecOps::RVec<edm4hep::TrackData> this_jet_tracks;
            this_jet_tracks = RecoToTrackTools::getRP2TRK_track(constituents, tracks);
            tracks_perjet.push_back(this_jet_tracks);
        }
        return tracks_perjet;
    }

    // helper functions
    // (copied from FCCAnalyses / JetConstituentsUtils)
    auto cast_constituent = [](const auto &jcs, auto &&meth)
    {
      ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> out;
      for (const auto &jc : jcs)
        out.emplace_back(meth(jc));
      return out;
    };

    auto cast_constituent_2 = [](const auto &jcs, const auto &coll1, auto &&meth)
    {
      ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> out;
      for (const auto &jc : jcs)
      {
        out.emplace_back(meth(jc, coll1));
      }
      return out;
    };
    
    auto cast_constituent_3 = [](const auto &jcs, const auto &coll1, const auto &coll2, auto &&meth)
    {
      ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> out;
      for (const auto &jc : jcs)
      {
        out.emplace_back(meth(jc, coll1, coll2));
      }
      return out;
    };

    auto cast_constituent_4 = [](const auto &jcs, const auto &coll1, const auto &coll2, const auto &coll3, auto &&meth)
    {
      ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> out;
      for (const auto &jc : jcs)
      {
        out.emplace_back(meth(jc, coll1, coll2, coll3));
      }
      return out;
    };

    // basic getters not present in FCCAnalyses / JetConstituentsUtils
    // (should they be added there?)
    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_px(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, FCCAnalyses::ReconstructedParticle::get_px);
    }

    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_py(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, FCCAnalyses::ReconstructedParticle::get_py);
    }

    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_pz(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, FCCAnalyses::ReconstructedParticle::get_pz);
    }

    // this function is not present in FCCAnalyses / JetConstituentsUtils
    // (should it be added there?)
    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_ptrel_log_cluster(const ROOT::VecOps::RVec<fastjet::PseudoJet> &jets,
                                                                   const ROOT::VecOps::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> out;
      for (size_t i = 0; i < jets.size(); ++i)
      {
        auto &jet_csts = out.emplace_back();
        float pt_jet = jets.at(i).pt();
        auto csts = FCCAnalyses::JetConstituentsUtils::get_jet_constituents(jcs, i);
        for (const auto &jc : csts)
        {
          TLorentzVector jcvec;
          jcvec.SetXYZM(jc.momentum.x, jc.momentum.y, jc.momentum.z, jc.mass);
          float val = (pt_jet > 0.) ? jcvec.Pt() / pt_jet : 1.;
          float ptrel_log = float(std::log10(val));
          jet_csts.emplace_back(ptrel_log);
        }
      }
      return out;
    }

    // this function is not present in FCCAnalyses / JetConstituentsUtils
    // (should it be added there?)
    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_ptrel_cluster(const ROOT::VecOps::RVec<fastjet::PseudoJet> &jets,
                                                               const ROOT::VecOps::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> out;
      for (size_t i = 0; i < jets.size(); ++i)
      {
        auto &jet_csts = out.emplace_back();
        double pt_jet = jets.at(i).pt();
        auto csts = FCCAnalyses::JetConstituentsUtils::get_jet_constituents(jcs, i);
        for (const auto &jc : csts)
        {
          TLorentzVector jcvec;
          jcvec.SetXYZM(jc.momentum.x, jc.momentum.y, jc.momentum.z, jc.mass);
          float val = (pt_jet > 0.) ? jcvec.Pt() / pt_jet : 1.;
          float ptrel = val;
          jet_csts.emplace_back(ptrel);
        }
      }
      return out;
    }

    /*
    Getters for track parameters
    */
    
    // get basic track properties
    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_chi2(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents>& jcs,
                                                    const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks)
    {
      return cast_constituent_2(jcs, tracks, RecoToTrackTools::getRP2TRK_chi2);
    }

    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_ndof(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents>& jcs,
                                                    const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks)
    {
      return cast_constituent_2(jcs, tracks, RecoToTrackTools::getRP2TRK_ndof);
    }

    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_chi2Normalized(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents>& jcs,
                                                    const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks)
    {
      return cast_constituent_2(jcs, tracks, RecoToTrackTools::getRP2TRK_chi2Normalized);
    }

    // get number of track hits in subdetectors
    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_nTrackHits(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents>& jcs,
                                                    const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
                                                    const ROOT::VecOps::RVec<int>& subdetectorHitNumbers,
                                                    const int subdetectorNumber)
    {
      return cast_constituent_4(jcs, tracks, subdetectorHitNumbers, subdetectorNumber,
                                RecoToTrackTools::getRP2TRK_nTrackHits);
    }

    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_nTrackHits_VDET(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents>& jcs,
                                                    const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
                                                    const ROOT::VecOps::RVec<int>& subdetectorHitNumbers)
    {
      return cast_constituent_3(jcs, tracks, subdetectorHitNumbers,
                                RecoToTrackTools::getRP2TRK_nTrackHits_VDET);
    }

    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_nTrackHits_ITC(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents>& jcs,
                                                    const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
                                                    const ROOT::VecOps::RVec<int>& subdetectorHitNumbers)
    {
      return cast_constituent_3(jcs, tracks, subdetectorHitNumbers,
                                RecoToTrackTools::getRP2TRK_nTrackHits_ITC);
    }

    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_nTrackHits_TPC(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents>& jcs,
                                                    const ROOT::VecOps::RVec<edm4hep::TrackData>& tracks,
                                                    const ROOT::VecOps::RVec<int>& subdetectorHitNumbers)
    {
      return cast_constituent_3(jcs, tracks, subdetectorHitNumbers,
                                RecoToTrackTools::getRP2TRK_nTrackHits_TPC);
    }

    // magnetic field
    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> get_Bz(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents> &jcs,
                                                    const ROOT::VecOps::RVec<edm4hep::TrackState> &trackStates)
    {
      return cast_constituent_2(jcs, trackStates, RecoToTrackTools::getRP2TRK_Bz);
    }

    // impact parameters
    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> XPtoPar_dxy(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents> &jcs,
                                                         const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks,
                                                         const TLorentzVector &PV,
                                                         const float &Bz)
    {

      return cast_constituent_4(jcs, tracks, PV, Bz, RecoToTrackTools::XPtoPar_dxy);
    }

    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> XPtoPar_dz(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents> &jcs,
                                                        const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks,
                                                        const TLorentzVector &PV,
                                                        const float &Bz)
    {

      return cast_constituent_4(jcs, tracks, PV, Bz, RecoToTrackTools::XPtoPar_dz);
    }

    ROOT::VecOps::RVec<FCCAnalysesJetConstituentsData> XPtoPar_phi(const ROOT::VecOps::RVec<FCCAnalysesJetConstituents> &jcs,
                                                         const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks,
                                                         const TLorentzVector &PV,
                                                         const float &Bz)
    {

      return cast_constituent_4(jcs, tracks, PV, Bz, RecoToTrackTools::XPtoPar_phi);
    }

 
}

#endif
