// Tools for collecting the dEdx information for jet constituents

// Orginally copied from here:
// https://github.com/Apranikstar/Aleph/blob/main/src/analyzer.h

#include <ROOT/RVec.hxx>

#include "FCCAnalyses/JetConstituentsUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"

#include "edm4hep/MCParticleData.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/RecDqdx.h"

#include <iostream>
#include <algorithm>


#ifndef dEdxTools_H
#define dEdxTools_H

namespace dEdxTools{

namespace rv = ROOT::VecOps;
using FCCAnalysesJetConstituents = rv::RVec<edm4hep::ReconstructedParticleData>;
using FCCAnalysesJetConstituentsData = rv::RVec<float>;

// Embedded Bethe-Bloch parameters {a, b, c, d, e}
std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> params = {
    {"e", {{"pads", {0.9932659976792287, 1.7094419426691188, 0.07384141776786941, 0.0, -2.0}},
           {"wires", {0.7930482536047483, 2.1221617291314856, 0.048686287341931526, 0.0, -2.0}}}},
    {"mu", {{"pads", {0.8590747, 1.32919203, 0.18728265, 0.01248365, -1.96898002}},
            {"wires", {0.55374855, 1.97190215, 0.26416447, 0.01549576, -1.972309}}}},
    {"pi", {{"pads", {0.536889582, 2.10964282, 0.269949484, 0.00310166058, -3.49991524}},
            {"wires", {0.7922792, 1.30952175, 0.19162276, 0.01523584, -2.4295921}}}},
    {"K", {{"pads", {0.25619823, 3.84049172, 0.53610855, 0.65090594, -2.44656219}},
           {"wires", {0.45920637, 1.82346531, 0.34181439, 0.52306056, -2.21132773}}}},
    {"p", {{"pads", {0.73189858, 1.05891917, 0.256201, 1.34293618, -2.02610177}},
           {"wires", {0.68811606, 1.03354162, 0.24500224, 1.44362792, -2.08656854}}}}
};

// Bethe-Bloch: a * (b + c * log(p) + d * p^e)
double bethe_bloch(double p, const std::vector<double>& par) {
    if (p <= 0.0) return 0.0;
    double a = par[0], b = par[1], c = par[2], d = par[3], e_pow = par[4];
    return a * (b + c * std::log(p) + d * std::pow(p, e_pow));
}

// Single hypothesis/single measurement: signed p-value (-9 if invalid/non-finite), is_wires=true for wires
double signed_p_value(double p, double dedx, double err, const std::string& hypothesis, bool is_wires) {
    std::string sensor = is_wires ? "wires" : "pads";
    
    auto part_it = params.find(hypothesis);
    if (part_it == params.end()) return -9.0;

    auto sens_it = part_it->second.find(sensor);
    if (sens_it == part_it->second.end()) return -9.0;

    const auto& par = sens_it->second;

    if (err <= 0.0) return -9.0;

    double expected = bethe_bloch(p, par);
    double residual = (dedx - expected) / err;
    if (!std::isfinite(residual)) return -9.0;

    double abs_z = std::fabs(residual);
    double sf = 0.5 * std::erfc(abs_z / std::sqrt(2.0));  // norm.sf(|z|)
    return 2.0 * sf * (residual > 0 ? 1.0 : -1.0);
}

// All hypotheses for single measurement: array[5] p-values {e, mu, pi, K, p}
std::array<double, 5> all_hypotheses_pvalues(double p, double dedx, double err, bool is_wires) {
    std::string hypos[5] = {"e", "mu", "pi", "K", "p"};
    std::array<double, 5> results;
    for (int i = 0; i < 5; ++i) {
        results[i] = signed_p_value(p, dedx, err, hypos[i], is_wires);
    }
    return results;
}

struct build_constituents_dEdx{
  struct dEdx_and_PID_result{
    rv::RVec<rv::RVec<edm4hep::RecDqdxData>> dedx_constituents;
    rv::RVec<rv::RVec<std::array<double, 5>>> pid_array_constituents;
    };

    dEdx_and_PID_result

    operator()(const rv::RVec<edm4hep::ReconstructedParticleData> &recoParticles,
             const rv::RVec<edm4hep::RecDqdxData> &dEdxCollection,
             const rv::RVec<int> &_dEdxIndicesCollection, 
             const std::vector<std::vector<int>> &jet_indices,
             bool is_wires) const
    { 
        rv::RVec<rv::RVec<edm4hep::RecDqdxData>> dedx_constituents;
        rv::RVec<rv::RVec<std::array<double, 5>>> pid_array_constituents;

        // The links dEdx -> Track and RecoPart -> Track are one-directional, we need a map to store
        // Track.index -> dEdx to not have to loop everytime 
        // in addition, the object itself is stored in <Collection>
        // while the relations (=indices we need for links) are on _<Collection>
        std::unordered_map<int, edm4hep::RecDqdxData> track_index_to_dEdx;
        for (size_t i = 0; i < _dEdxIndicesCollection.size(); ++i) {
          int track_index = _dEdxIndicesCollection[i];
          edm4hep::RecDqdxData dedx_value = dEdxCollection[i];
          track_index_to_dEdx[track_index] = dedx_value;
        }

        // Define dummy for particles with no track or for which no dEdx measurement is found
        edm4hep::RecDqdxData dedx_dummy{};
        dedx_dummy.dQdx.type = -9.0f;
        dedx_dummy.dQdx.value = -9.0f;
        dedx_dummy.dQdx.error = -9.0f;
        std::array<double, 5> pid_array_dummy{{-9.0f, -9.0f, -9.0f, -9.0f, -9.0f}};

        // Now, for each jet loop over the indices of the jet constituents provided by the JetClusteringUtils,
        // retrieve the associated RecoParticle,
        // from there get the link to the Track.
        for (const auto &jet_const_indices : jet_indices) { //loop over jets
          rv::RVec<edm4hep::RecDqdxData> jet_dEdx;
          rv::RVec<std::array<double, 5>> jet_pid_array;

          for (int constituent_index : jet_const_indices) { // loop over jet constituents
            const auto &recoPart = recoParticles[constituent_index];
            TLorentzVector tlv_recoPart; // needed later to get the momentum total
            tlv_recoPart.SetXYZM(recoPart.momentum.x, recoPart.momentum.y, recoPart.momentum.z, recoPart.mass);

            // Initialize flag for good value found (otherwise use dummy)
            bool found = false;

            //loop over tracks associated to the RecoPart (for charged particles should always be exactly one in ALEPH data)
            for (int track_index = recoPart.tracks_begin; track_index < recoPart.tracks_end; ++track_index) {

                  //find the matching dEdx in the map
                  if (track_index_to_dEdx.count(track_index)) {
                    const auto &dEdx = track_index_to_dEdx[track_index];

                    //get the PID hypotheses p-values (= array of five entries fo e, mu, pi, K, p)
                    const auto PID_pvals_array = all_hypotheses_pvalues(tlv_recoPart.P(), dEdx.dQdx.value, dEdx.dQdx.error, is_wires);

                    //check wether the measurement is valid, if not fill default value
                    if (dEdx.dQdx.type == 0) {
                      jet_dEdx.push_back(track_index_to_dEdx[track_index]);
                      jet_pid_array.push_back(PID_pvals_array);
                    }
                    else {
                      jet_dEdx.push_back(dedx_dummy);
                      jet_pid_array.push_back(pid_array_dummy);
                    }

                    found = true;
                    break;
                  }
            }
            // if no track found, i.e. neutral particle, use the dummy 
            if (!found){
              jet_dEdx.push_back(dedx_dummy);
              jet_pid_array.push_back(pid_array_dummy);
            }
          }
          dedx_constituents.push_back(jet_dEdx); 
          pid_array_constituents.push_back(jet_pid_array);
        }
        return {dedx_constituents, pid_array_constituents};
    }
};

//helpers to read the dEdx objects (to check if can reuse existing FCCAna functions instead):

rv::RVec<rv::RVec<float>> get_dEdx_type(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.type);
    }
    values.push_back(inner_values);
  }
  return values;
}

rv::RVec<rv::RVec<float>> get_dEdx_value(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.value);
    }
    values.push_back(inner_values);
  }
  return values;
}

rv::RVec<rv::RVec<float>> get_dEdx_error(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>> &dedx_vec) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : dedx_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d.dQdx.error);
    }
    values.push_back(inner_values);
  }
  return values;
}

rv::RVec<rv::RVec<float>> get_PID_pvalue(const rv::RVec<rv::RVec<std::array<double, 5>>> pid_array_vec, int particle_index) {
  rv::RVec<rv::RVec<float>> values;
  for (const auto &inner_vec : pid_array_vec) {
    rv::RVec<float> inner_values;
    for (const auto &d : inner_vec) {
      inner_values.push_back(d[particle_index]);
    }
    values.push_back(inner_values);
  }
  return values;
}

}

#endif
