/*
Tools for calculating ALEPH-style impact parameters

These tools include:
- calculating the point of closest approach between a track and a (jet) axis.
- linearizing the track around that point.
- calculating the shortest distance between the linearized track and a (primary) vertex.
*/


#ifndef FCCAnalyses_IPAlephTools_h
#define FCCAnalyses_IPAlephTools_h

#include <cmath>
#include "ROOT/RVec.hxx"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/ReconstructedParticle.h"
#include "podio/ObjectID.h"
#include "fastjet/JetDefinition.hh"


/*
Helper classes
*/

class Helix {
    /*
    Helper class for easily getting the position and tangent to a helix
    for a given path-length parameter.
    */
    public:

    explicit Helix(const edm4hep::TrackState& trackState){
        d0 = trackState.D0;
        z0 = trackState.Z0;
        phi0 = trackState.phi;
        omega = trackState.omega;
        tanL = trackState.tanLambda;
        
        straight = std::abs(omega) < 1e-12;
        R = straight ? 1e12 : 1.0 / omega;

        // Perigee position
        r0.SetXYZ(
            -d0 * std::sin(phi0),
             d0 * std::cos(phi0),
             z0
        );
    }

    // Position on helix at path-length parameter s
    TVector3 position(double s) const {
        if (straight) { return r0 + s * direction(0.0); }
        return TVector3(
            r0.X() + R * (std::sin(phi0 + omega*s) - std::sin(phi0)),
            r0.Y() - R * (std::cos(phi0 + omega*s) - std::cos(phi0)),
            r0.Z() + s * tanL
        );
    }

    // Tangent dr/ds (NOT unit length)
    TVector3 tangent(double s) const {
        if (straight) { return TVector3(std::cos(phi0), std::sin(phi0), tanL); }
        return TVector3(
            std::cos(phi0 + omega*s),
            std::sin(phi0 + omega*s),
            tanL
        );
    }

    // Unit direction at s
    TVector3 direction(double s) const {
        TVector3 t = tangent(s);
        return t.Unit();
    }

    private:

    double d0;
    double z0;
    double phi0;
    double omega;
    double tanL;

    bool straight;
    double R;

    TVector3 r0;
};


/*
Helper functions
*/

double helixToAxisPCA_newton(
    const Helix& helix,
    const TVector3& axisPoint,
    const TVector3& axisDir){
    /*
    Helper function to find the path-length parameter on a helix
    where it is closest to a given axis, using a Newton-Raphson based approach.
    */
    double s_helix = 0.0;
    int max_iter = 20;
    for (int iter = 0; iter < max_iter; ++iter){
        TVector3 r_helix = helix.position(s_helix);
        TVector3 v_helix = helix.tangent(s_helix);

        double s_axis = (r_helix - axisPoint).Dot(axisDir);
        TVector3 r_axis = axisPoint + s_axis * axisDir;

        TVector3 u = r_helix - r_axis;
        TVector3 v_perp = v_helix - axisDir * v_helix.Dot(axisDir);
        double df = 2.0 * u.Dot(v_perp);
        double d2f = 2.0 * v_perp.Mag2(); // approximation

        if (std::abs(df/s_helix) < 1e-6) break;
        if (std::abs(d2f/s_helix) < 1e-6) break;

        s_helix = s_helix - df / d2f;
    }
    return s_helix;
}

double helixToAxisPCA_step(
    const Helix& helix,
    const TVector3& axisPoint,
    const TVector3& axisDir){
    /*
    Helper function to find the path-length parameter on a helix
    where it is closest to a given axis, using a simple iterative step approach.
    */
    double s_helix = 0.0;
    int max_iter = 100;
    for (int iter = 0; iter < max_iter; ++iter){
        
        // Get position and tangent to helix at current point
        TVector3 r_helix = helix.position(s_helix);
        TVector3 t_helix = helix.direction(s_helix);

        // Closest point on axis to current helix point
        double s_axis = (r_helix - axisPoint).Dot(axisDir);
        TVector3 r_axis = axisPoint + s_axis * axisDir;

        // Update step: move along helix towards axis
        double step = (r_axis - r_helix).Dot(t_helix);
        s_helix += step;

        // Convergence criterion
        if( std::abs(step/s_helix) < 1e-5 ) break;
    }
    return s_helix;
}

double helixToAxisPCA(
    const Helix& helix,
    const TVector3& axisPoint,
    const TVector3& axisDir,
    const std::string& method){
    /*
    Helper function to find the path-length parameter on a helix
    where it is closest to a given axis.
    */
    if( method=="newton" ){ return helixToAxisPCA_newton(helix, axisPoint, axisDir); }
    else if( method == "step" ){ return helixToAxisPCA_step(helix, axisPoint, axisDir); }
    // todo: throw an error if none of the above
    return 0.0;
}

TMatrixDSym getCovarianceMatrix(const edm4hep::TrackState& trackState){
    /*
    Helper function to retrieve the covariance matrix from a TrackState.
    Note: this relies crucially on the assumed ordering of the elements
          in the stored covariance array.
          See the LCIO ordering convention and also similar tools in
          ReconstructedParticle2Track.cc.
          Not sure if our ALEPH data actually follows this convention!
    */
    TMatrixDSym C(5);
    const auto& cov = trackState.covMatrix;

    // Diagonal
    C(0,0) = cov[0];    // d0
    C(1,1) = cov[2];    // phi
    C(2,2) = cov[5];    // omega
    C(3,3) = cov[9];    // z0
    C(4,4) = cov[14];   // tanLambda

    // Off-diagonals: symmetrize
    C(1,0) = C(0,1) = cov[1];   // phi-d0
    C(2,0) = C(0,2) = cov[3];   // omega-d0
    C(2,1) = C(1,2) = cov[4];   // omega-phi

    C(3,0) = C(0,3) = cov[6];   // z0-d0
    C(3,1) = C(1,3) = cov[7];   // z0-phi
    C(3,2) = C(2,3) = cov[8];   // z0-omega

    C(4,0) = C(0,4) = cov[10];  // tanλ-d0
    C(4,1) = C(1,4) = cov[11];  // tanλ-phi
    C(4,2) = C(2,4) = cov[12];  // tanλ-omega
    C(4,3) = C(3,4) = cov[13];  // tanλ-z

    return C;
}

TMatrixD getPositionJacobianMatrix(const edm4hep::TrackState& trackState, double s0){
    /*
    Helper function to calculate the Jacobian matrix of the position,
    i.e. the partial derivatives of the track coordinates (x,y,z)
    to the track parameters (d0, phi0, omega, z0, tanLambda),
    at a provided path-length parameter s0.
    */
    TMatrixD J(3,5);

    // Track parameters
    const double d0 = trackState.D0;
    const double phi = trackState.phi;
    const double omega = trackState.omega;
    const double z0 = trackState.Z0;
    const double tanL = trackState.tanLambda;

    // Protect against zero curvature
    bool straight = std::abs(omega) < 1e-12;
    const double R = straight ? 1e12 : 1.0 / omega;
    const double R2 = R * R;

    // Other helpful quantities
    const double Phi = phi + omega * s0;
    const double sinphi = std::sin(phi);
    const double cosphi = std::cos(phi);
    const double sinPhi = std::sin(Phi);
    const double cosPhi = std::cos(Phi);

    // --------------------
    // x derivatives
    // --------------------
    J(0,0) = -sinphi; // dx/dd0
    J(0,1) = -d0 * cosphi + R * (cosPhi - cosphi); // dx/dphi
    J(0,2) = -R2 * (sinPhi - sinphi) + s0 * R * cosPhi; // dx/domega
    J(0,3) = 0.0; // dx/dz0
    J(0,4) = 0.0; // dx/dtanLambda

    // --------------------
    // y derivatives
    // --------------------
    J(1,0) = cosphi; // dy/dd0
    J(1,1) = -d0 * sinphi + R * (sinPhi - sinphi); // dy/dphi
    J(1,2) = R2 * (cosPhi - cosphi) + s0 * R * sinPhi; // dy/domega
    J(1,3) = 0.0; // dy/dz0
    J(1,4) = 0.0; // dy/dtanLambda

    // --------------------
    // z derivatives
    // --------------------
    J(2,0) = 0.0; // dz/dd0
    J(2,1) = 0.0; // dz/dphi
    J(2,2) = 0.0; // dz/domega
    J(2,3) = 1.0; // dz/dz0
    J(2,4) = s0; // dz/dtanLambda

    return J;
}

TMatrixD getPositionCovarianceMatrix(const edm4hep::TrackState& trackState, double s0){
    /*
    Calculate the position covariance matrix at given path-length s0.
    */
    TMatrixDSym C = getCovarianceMatrix(trackState);
    TMatrixD J = getPositionJacobianMatrix(trackState, s0);
    TMatrixD temp(3, 5);
    temp.Mult(J, C);
    TMatrixD res(3, 3);
    res.Mult(temp, J.T());
    return res;
}

double getPositionVariance(const edm4hep::TrackState& trackState,
        double s0, const TVector3& direction){
    /*
    Get position variance in a given direction
    */

    // Get position covariance matrix at a given point
    TMatrixD positionCovariance = getPositionCovarianceMatrix(trackState, s0);

    // Make sure the direction is a unit vector
    TVector3 unitDirection = direction.Unit();

    // Calculate transpose(direction) * covariance * direction
    double sigma2 = 0.0;
    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            sigma2 += unitDirection(i) * positionCovariance(i,j) * unitDirection(j);
        }
    }
    return sigma2;
}

ROOT::VecOps::RVec<TVector3> trackToAxisPCA(
    const edm4hep::TrackState& trackState,
    const TVector3& referencePoint,
    const TVector3& direction){
    /*
    Helper function to calculate the points of closest approach between a track and an axis
    (i.e. a line following a given direction and going through a given point).
    Usually this axis would be a jet axis, i.e. the line following the jet momentum
    and going through the primary vertex.
    This function returns a vector of points, in this order:
     - point on the track closest to the axis.
     - point on the axis closest to the track.
     - point on the linearized track (around the first point) closest to the reference point on the axis.
     - dummy point with additional parameters; these are:
         - the path-length parameter on the track helix where it is closest to the axis.
         - distance of closest approach between the linearized track and the reference point.
         - uncertainty in that distance.
    */

    // Direction must be a unit vector
    TVector3 axisDir = direction.Unit();

    // Make a Helix object from track parameters
    Helix helix = Helix(trackState);

    // Get path-length of POCA
    double s_track = helixToAxisPCA(helix, referencePoint, axisDir, "newton");

    // Final position on helix and axis
    TVector3 trackPCA = helix.position(s_track);
    double s_axis = (trackPCA - referencePoint).Dot(axisDir);
    TVector3 axisPCA = referencePoint + s_axis * axisDir;

    // Linearization of track around POCA to axis
    TVector3 trackDirection = helix.direction(s_track);

    // Calculate POCA between linear extension and referencePoint
    double s_line = (referencePoint - trackPCA).Dot(trackDirection);
    TVector3 linePCA = trackPCA + s_line * trackDirection;

    // Calculate distance of closest approach between linear extension and referencePoint
    double dca = (linePCA - referencePoint).Mag();

    // Estimate uncertainty in that distance
    // Note: we take the covariance matrix at the point on the track helix
    //       where it is closest to the axis, and use this as uncertainty
    //       in the position of the POCA of the linear extension to the referencePoint;
    //       other effects are ignored.
    TVector3 connector = (referencePoint - linePCA).Unit();
    double dcaVar = getPositionVariance(trackState, s_track, connector);
    double dcaUnc = std::sqrt(dcaVar);

    // Make a dummy point with additional parameters to return
    TVector3 params = {s_track, dca, dcaUnc};

    ROOT::VecOps::RVec<TVector3> out = {trackPCA, axisPCA, linePCA, params};
    return out;
}

/*
Interface functions
*/

namespace FCCAnalyses{

namespace IPAlephTools{

namespace rv = ROOT::VecOps;

rv::RVec<rv::RVec<rv::RVec<TVector3>>> getTrackToJetAxisPCA(
        const rv::RVec<fastjet::PseudoJet>& jets,
        const rv::RVec<rv::RVec<edm4hep::ReconstructedParticleData>>& jetConstituents,
        const rv::RVec<edm4hep::TrackState>& tracks,
        const TVector3& primaryVertex){

    // initialize output
    rv::RVec<rv::RVec<rv::RVec<TVector3>>> out;

    // loop over jets and jet constituents
    for (int jetidx = 0; jetidx < jets.size(); ++jetidx){
        TVector3 jetMomentum(jets.at(jetidx).px(), jets.at(jetidx).py(), jets.at(jetidx).pz());
        rv::RVec<rv::RVec<TVector3>> temp;
        for (int jcidx = 0; jcidx < jetConstituents.at(jetidx).size(); ++jcidx){
            edm4hep::ReconstructedParticleData recoParticle = jetConstituents.at(jetidx).at(jcidx);
            // find track for this jet constituent
            size_t trackidx = recoParticle.tracks_begin;
            if(trackidx < tracks.size()){
                edm4hep::TrackState trackState = tracks.at(trackidx);
                // calculate PCA
                rv::RVec<TVector3> pcas = trackToAxisPCA(trackState, primaryVertex, jetMomentum);
                temp.push_back(pcas);
            } else{
                rv::RVec<TVector3> dummy = {
                    TVector3(0, 0, 0), TVector3(0, 0, 0), TVector3(0, 0, 0), TVector3(0, 0, 0)
                };
                temp.push_back(dummy);
            }
        }
        out.push_back(temp);
    }
    return out;
}

rv::RVec<rv::RVec<TVector3>> getTrackPCAToJetAxis(
        const rv::RVec<rv::RVec<rv::RVec<TVector3>>>& pcas){
    rv::RVec<rv::RVec<TVector3>> out;
    for (int i = 0; i < pcas.size(); ++i){
        rv::RVec<TVector3> temp;
        for (int j = 0; j < pcas.at(i).size(); ++j){
            temp.push_back(pcas.at(i).at(j).at(0));
        }
        out.push_back(temp);
    }
    return out;
}

rv::RVec<rv::RVec<TVector3>> getJetAxisPCAToTrack(
        const rv::RVec<rv::RVec<rv::RVec<TVector3>>>& pcas){
    rv::RVec<rv::RVec<TVector3>> out;
    for (int i = 0; i < pcas.size(); ++i){
        rv::RVec<TVector3> temp;
        for (int j = 0; j < pcas.at(i).size(); ++j){
            temp.push_back(pcas.at(i).at(j).at(1));
        }
        out.push_back(temp);
    }
    return out;
}

rv::RVec<rv::RVec<TVector3>> getLinePCAToPrimaryVertex(
        const rv::RVec<rv::RVec<rv::RVec<TVector3>>>& pcas){
    rv::RVec<rv::RVec<TVector3>> out;
    for (int i = 0; i < pcas.size(); ++i){
        rv::RVec<TVector3> temp;
        for (int j = 0; j < pcas.at(i).size(); ++j){
            temp.push_back(pcas.at(i).at(j).at(2));
        }
        out.push_back(temp);
    }
    return out;
}

rv::RVec<rv::RVec<float>> getPCA_x(const rv::RVec<rv::RVec<TVector3>>& pcas){
    rv::RVec<rv::RVec<float>> out;
    for (int i = 0; i < pcas.size(); ++i){
        rv::RVec<float> temp;
        for (int j = 0; j < pcas.at(i).size(); ++j){
            temp.push_back(pcas.at(i).at(j).x());
        }
        out.push_back(temp);
    }
    return out;
}

rv::RVec<rv::RVec<float>> getPCA_y(const rv::RVec<rv::RVec<TVector3>>& pcas){
    rv::RVec<rv::RVec<float>> out;
    for (int i = 0; i < pcas.size(); ++i){
        rv::RVec<float> temp;
        for (int j = 0; j < pcas.at(i).size(); ++j){
            temp.push_back(pcas.at(i).at(j).y());
        }
        out.push_back(temp);
    }
    return out;
}

rv::RVec<rv::RVec<float>> getPCA_z(const rv::RVec<rv::RVec<TVector3>>& pcas){
    rv::RVec<rv::RVec<float>> out;
    for (int i = 0; i < pcas.size(); ++i){
        rv::RVec<float> temp;
        for (int j = 0; j < pcas.at(i).size(); ++j){
            temp.push_back(pcas.at(i).at(j).z());
        }
        out.push_back(temp); 
    }
    return out;
}

rv::RVec<rv::RVec<float>> signedIP3D(
        const rv::RVec<rv::RVec<rv::RVec<TVector3>>>& pcas,
        const rv::RVec<fastjet::PseudoJet>& jets,
        const TVector3& primaryVertex){
    /*
    Calculate signed 3D impact parameter.
    Note: assumed to run on the output of getTrackToJetAxisPCA.
    */
    rv::RVec<rv::RVec<float>> out;
    for (int i = 0; i < pcas.size(); ++i){
        TVector3 jetMomentum(jets.at(i).px(), jets.at(i).py(), jets.at(i).pz());
        rv::RVec<float> temp;
        for (int j = 0; j < pcas.at(i).size(); ++j){
            // get params for this track
            TVector3 axisPCA = pcas.at(i).at(j).at(1);
            TVector3 params = pcas.at(i).at(j).at(3);
            double ip3d = params(1);
            // calculate sign
            double proj = (axisPCA - primaryVertex).Dot(jetMomentum);
            int sign = (proj >= 0) ? 1 : -1;
            // calculate signed impact parameter
            double sip3d = sign * ip3d;
            temp.push_back(sip3d);
        }
        out.push_back(temp);
    }
    return out;
}

rv::RVec<rv::RVec<float>> IP3DUnc(
        const rv::RVec<rv::RVec<rv::RVec<TVector3>>>& pcas,
        const TVector3& primaryVertex){
    /*
    Calculate 3D impact parameter uncertainty.
    Note: assumed to run on the output of getTrackToJetAxisPCA.
    */
    rv::RVec<rv::RVec<float>> out;
    for (int i = 0; i < pcas.size(); ++i){
        rv::RVec<float> temp;
        for (int j = 0; j < pcas.at(i).size(); ++j){
            // get params for this track
            TVector3 params = pcas.at(i).at(j).at(3);
            double ip3dunc = params(2);
            temp.push_back(ip3dunc);
        }
        out.push_back(temp);
    }
    return out;
}

rv::RVec<rv::RVec<float>> signedIP3DSig(
        const rv::RVec<rv::RVec<rv::RVec<TVector3>>>& pcas,
        const rv::RVec<fastjet::PseudoJet>& jets,
        const TVector3& primaryVertex){
    /*
    Calculate signed 3D impact parameter.
    Note: assumed to run on the output of getTrackToJetAxisPCA.
    */
    rv::RVec<rv::RVec<float>> out;
    for (int i = 0; i < pcas.size(); ++i){
        TVector3 jetMomentum(jets.at(i).px(), jets.at(i).py(), jets.at(i).pz());
        rv::RVec<float> temp;
        for (int j = 0; j < pcas.at(i).size(); ++j){
            // get params for this track
            TVector3 axisPCA = pcas.at(i).at(j).at(1);
            TVector3 params = pcas.at(i).at(j).at(3);
            double ip3d = params(1);
            double ip3dunc = params(2);
            // calculate sign
            double proj = (axisPCA - primaryVertex).Dot(jetMomentum);
            int sign = (proj >= 0) ? 1 : -1;
            // safety for 0 uncertainty
            if (ip3dunc < 1e-12){
                ip3d = 0;
                ip3dunc = 1; 
            }
            // calculate signed impact parameter significance
            double sip3d = sign * ip3d / ip3dunc;
            temp.push_back(sip3d);
        }
        out.push_back(temp);
    }
    return out;
}

rv::RVec<rv::RVec<float>> transverseJetDistance(
        const rv::RVec<rv::RVec<rv::RVec<TVector3>>>& pcas){
    /*
    Calculate transverseJetDistance (= distance between track and axis PCAs).
    Note: assumed to run on the output of getTrackToJetAxisPCA.
    */
    rv::RVec<rv::RVec<float>> out;
    for (int i = 0; i < pcas.size(); ++i){
        rv::RVec<float> temp;
        for (int j = 0; j < pcas.at(i).size(); ++j){
            // get points for this track
            TVector3 trackPCA = pcas.at(i).at(j).at(0);
            TVector3 axisPCA = pcas.at(i).at(j).at(1);
            // calculate distance between track and axis PCA
            double dist = (trackPCA - axisPCA).Mag();
            temp.push_back(dist);
        }
        out.push_back(temp);
    }
    return out;
}

rv::RVec<rv::RVec<float>> longitudinalJetDistance(
        const rv::RVec<rv::RVec<rv::RVec<TVector3>>>& pcas,
        const TVector3& primaryVertex){
    /*
    Calculate longitudinal jet distance (= distance between axis PCA and primary vertex).
    Note: assumed to run on the output of getTrackToJetAxisPCA.
    */
    rv::RVec<rv::RVec<float>> out;
    for (int i = 0; i < pcas.size(); ++i){
        rv::RVec<float> temp;
        for (int j = 0; j < pcas.at(i).size(); ++j){
            // get points for this track
            TVector3 axisPCA = pcas.at(i).at(j).at(1);
            // calculate distance between axis PCA and PV
            double dist = (axisPCA - primaryVertex).Mag();
            temp.push_back(dist);
        }
        out.push_back(temp);
    }
    return out;
}

} // end namespace IPAlephTools

} // end namespace FCCAnalyses

#endif
