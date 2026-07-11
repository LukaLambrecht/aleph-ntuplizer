// Tools for dealing with collections of vertices

#include <ROOT/RVec.hxx>

#include "FCCAnalyses/JetConstituentsUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/ReconstructedParticleData.h"

#include <iostream>
#include <algorithm>


#ifndef VertexTools_H
#define VertexTools_H

namespace VertexTools{

// note: the functions below are not present in FCCAnalyses / VertexingUtils.
//       should they be added there?

// vector of position of all reconstructed SV (in mm)
ROOT::VecOps::RVec<TVector3>
get_position_SV(const ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>& vertices) {
    ROOT::VecOps::RVec<TVector3> result;
    for (FCCAnalyses::VertexingUtils::FCCAnalysesVertex ivtx : vertices) {
        TVector3 xyz(ivtx.vertex.position[0], ivtx.vertex.position[1], ivtx.vertex.position[2]);
        result.push_back(xyz);
    }
    return result;
}

// same as above, but only x-coordinate
ROOT::VecOps::RVec<double>
get_x_SV(const ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>& vertices) {
    ROOT::VecOps::RVec<double> result;
    for (FCCAnalyses::VertexingUtils::FCCAnalysesVertex ivtx : vertices) {
        result.push_back(ivtx.vertex.position[0]);
    }
    return result;
}

// same as above, but relative to PV
ROOT::VecOps::RVec<double>
get_xrel_SV(const ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>& vertices,
            const TVector3& PV) {
    ROOT::VecOps::RVec<double> result;
    for (FCCAnalyses::VertexingUtils::FCCAnalysesVertex ivtx : vertices) {
        result.push_back(ivtx.vertex.position[0] - PV.x()); 
    }
    return result;
}

// same as above, but only y-coordinate
ROOT::VecOps::RVec<double>
get_y_SV(const ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>& vertices) {
    ROOT::VecOps::RVec<double> result;
    for (FCCAnalyses::VertexingUtils::FCCAnalysesVertex ivtx : vertices) {
        result.push_back(ivtx.vertex.position[1]);
    }
    return result;
}

// same as above, but relative to PV
ROOT::VecOps::RVec<double>
get_yrel_SV(const ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>& vertices,
            const TVector3& PV) {
    ROOT::VecOps::RVec<double> result;
    for (FCCAnalyses::VertexingUtils::FCCAnalysesVertex ivtx : vertices) {
        result.push_back(ivtx.vertex.position[1] - PV.y());
    }
    return result;
}

// same as above, but only z-coordinate
ROOT::VecOps::RVec<double>
get_z_SV(const ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>& vertices) {
    ROOT::VecOps::RVec<double> result;
    for (FCCAnalyses::VertexingUtils::FCCAnalysesVertex ivtx : vertices) {
        result.push_back(ivtx.vertex.position[2]);
    }
    return result;
}

// same as above, but relative to PV
ROOT::VecOps::RVec<double>
get_zrel_SV(const ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>& vertices,
            const TVector3& PV) {
    ROOT::VecOps::RVec<double> result;
    for (FCCAnalyses::VertexingUtils::FCCAnalysesVertex ivtx : vertices) {
        result.push_back(ivtx.vertex.position[2] - PV.z());
    }
    return result;
}

// vector of position of all reconstructed SV (for the case of grouping per jet)
ROOT::VecOps::RVec<ROOT::VecOps::RVec<TVector3>>
get_position_SV_jets(ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>> vertices) {
  ROOT::VecOps::RVec<ROOT::VecOps::RVec<TVector3>> result;
  for (ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> this_vertices: vertices) {
    result.push_back(get_position_SV(this_vertices));
  }
  return result;
}

// same as above, but only x-coordinate
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_x_SV_jets(const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>>& vertices) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    for (ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> this_vertices : vertices) {
        result.push_back(get_x_SV(this_vertices));
    }
    return result;
}

// same as above, but relative to PV
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_xrel_SV_jets(const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>>& vertices,
            const TVector3& PV) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    for (ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> this_vertices : vertices) {
        result.push_back(get_xrel_SV(this_vertices, PV));
    }
    return result;
}

// same as above, but only y-coordinate
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_y_SV_jets(const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>>& vertices) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    for (ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> this_vertices : vertices) {
        result.push_back(get_y_SV(this_vertices));
    }
    return result;
}

// same as above, but relative to PV
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_yrel_SV_jets(const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>>& vertices,
            const TVector3& PV) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    for (ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> this_vertices : vertices) {
        result.push_back(get_yrel_SV(this_vertices, PV));
    }
    return result;
}

// same as above, but only z-coordinate
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_z_SV_jets(const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>>& vertices) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    for (ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> this_vertices : vertices) {
        result.push_back(get_z_SV(this_vertices));
    }
    return result;
}

// same as above, but relative to PV
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>
get_zrel_SV_jets(const ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>>& vertices,
            const TVector3& PV) {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    for (ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> this_vertices : vertices) {
        result.push_back(get_zrel_SV(this_vertices, PV));
    }
    return result;
}

// SV relative momentum
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> get_prel_SV_jets(
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>> vertices,
    ROOT::VecOps::RVec<fastjet::PseudoJet> jets) {
  ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
  ROOT::VecOps::RVec<double> i_result;

  // loop over jets
  for (unsigned int i = 0; i < jets.size(); i++) {
    ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> i_vertices = vertices.at(i);
    fastjet::PseudoJet i_jet = jets.at(i);
    // loop over vertices
    for (auto &ivtx : i_vertices) {
        
      // make sum of track momenta
      ROOT::VecOps::RVec<TVector3> p_tracks = ivtx.updated_track_momentum_at_vertex;
      TVector3 p_sum;
      for (TVector3 p_tr : p_tracks){ p_sum += p_tr; }
      double numerator = p_sum.Mag();

      // get jet momentum
      TVector3 jetP(i_jet.px(), i_jet.py(), i_jet.pz());
      double denominator = jetP.Mag();

      // take ratio
      i_result.push_back(numerator / denominator);
    }
    result.push_back(i_result);
    i_result.clear();
  }
  return result;
}

// pointing angle
// (copied from FCCAnalyses but with bugfix in r_vtx_PV)
// SV pointing angle wrt PV
ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> get_pointingangle_SV(
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>> vertices,
    FCCAnalyses::VertexingUtils::FCCAnalysesVertex PV) {
  ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
  ROOT::VecOps::RVec<double> i_result;

  edm4hep::Vector3f r_PV = PV.vertex.position; // in mm

  for (unsigned int i = 0; i < vertices.size(); i++) {
    ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> i_vertices = vertices.at(i);
    for (auto &ivtx : i_vertices) {
      double pointangle = 0.;

      ROOT::VecOps::RVec<TVector3> p_tracks =
          ivtx.updated_track_momentum_at_vertex;
      TVector3 p_sum;
      for (TVector3 p_tr : p_tracks)
        p_sum += p_tr;

      edm4hep::Vector3f r_vtx = ivtx.vertex.position; // in mm

      TVector3 r_vtx_PV(r_vtx[0] - r_PV[0],
                        r_vtx[1] - r_PV[1],
                        r_vtx[2] - r_PV[2]);

      double pDOTr = p_sum.Dot(r_vtx_PV);
      double p_mag = p_sum.Mag();
      double r_mag = r_vtx_PV.Mag();

      pointangle = pDOTr / (p_mag * r_mag);
      i_result.push_back(pointangle);
    }
    result.push_back(i_result);
    i_result.clear();
  }
  return result;
}

}

#endif
