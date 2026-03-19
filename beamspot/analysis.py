import os
import sys
import ROOT

# load custom analyzer with collections of tracks
analyzer_path = os.path.join(os.path.dirname(__file__), '../analyzers', 'analyzer_tracktools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with primary vertex tools
analyzer_path = os.path.join(os.path.dirname(__file__), '../analyzers', 'analyzer_pvtools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# helper function to re-calculate the primary vertex from the collection of tracks.
ROOT.gInterpreter.Declare("""
    TLorentzVector fitRecoPrimaryVertex(
        ROOT::VecOps::RVec<edm4hep::TrackState> tracks){
        ROOT::VecOps::RVec<edm4hep::TrackState> tracksToUse;
        for (const auto& trk : tracks) {
            const auto& c = trk.covMatrix;
            if (c[0] <= 0 || c[2] <= 0 || c[9] <= 0) continue;
            if (c[0] < 1e-6 || c[2] < 1e-6 || c[9] <= 1e-6) continue;
            if (!std::isfinite(c[0]) || !std::isfinite(c[2]) || !std::isfinite(c[9])) continue;
            if (std::abs(trk.D0)>0.75 || std::abs(trk.Z0)>2) continue;
            tracksToUse.push_back(trk);
        }
        if( tracksToUse.size() < 2 ){ return TLorentzVector(0, 0, 0, 0); }
        ROOT::VecOps::RVec<edm4hep::TrackState> primaryTracks;
        primaryTracks = PrimaryVertexTools::get_PrimaryTracks(tracksToUse,
            25., false, 0, 0, 0, 0, 0, 0);
        if( primaryTracks.size() < 2 ){ return TLorentzVector(0, 0, 0, 0); }
        FCCAnalyses::VertexingUtils::FCCAnalysesVertex fitresult;
        fitresult = FCCAnalyses::VertexFitterSimple::VertexFitter_Tk(1, primaryTracks,
            false, 0, 0, 0, 0, 0, 0);
        edm4hep::VertexData vertex = FCCAnalyses::VertexingUtils::get_VertexData(fitresult);
        TLorentzVector result = {vertex.position.x, vertex.position.y, vertex.position.z, 0.};
        return result;
    }""")


# main analyzer class
class RDFanalysis():

    def analysers(df):

        # initialization
        dfout = df

        # automatically determine detector type based on the presence of some branches
        # (maybe make more robust later)
        det = None
        try:
            df.Alias("_", "ReconstructedParticles")
            det = 'fcc'
        except: pass
        try:
            df.Alias("_", "RecoParticles")
            det = 'aleph'
        except: pass
        if det is None:
            raise Exception('Could not determine detector type.')

        # automatically determine data type based on the presence of some branches
        # (maybe make more robust later)
        dtype = 'data'
        try:
            if det=='fcc': df.Alias("_", "Particle")
            if det=='aleph': df.Alias("_", "MCParticles")
            dtype = 'sim'
        except: pass

        # for Aleph simulation, the collection of MC particles is called "MCParticles",
        # while for FCC simulation it is called "Particle", so define an alias here.
        # and same for "RecoParticles" vs "ReconstructedParticles".
        # and same for "ParticleID" vs "ParticleIDs" (note: not well defined for FCC).
        if det=='aleph':
            if dtype=='sim': dfout = dfout.Alias("Particle", "MCParticles")
            dfout = dfout.Alias("ReconstructedParticles", "RecoParticles")
            dfout = dfout.Alias("ParticleIDs", "ParticleID")

        # for Aleph simulation, the collections EFlowTrack, EFlowTrack_1 and EFlowTrack_2 do not seem to exist,
        # so we need to alias them with other collections.
        # note: the alias for EFlowTrack_2 has not yet been validated, no guarantee that it is correct.
        # note: there are also some modification that need to be made to the Aleph track states:
        #       - flip the sign of D0 and omega (at least so it seems; under investigation)
        #       - unit conversion from cm to mm (this is currently still hard-coded in my branch of FCCAnalyses)
        if det=='aleph':
            dfout = (
                dfout

                .Alias("EFlowTrack", "Tracks")
                # (must be an object of type rv::RVec<edm4hep::TrackData>)
                .Define("EFlowTrack_1", "TrackTools::getModifiedTrackStates(_Tracks_trackStates)")
                # (must be an object of type ROOT::VecOps::RVec<edm4hep::TrackState>)
                .Define("EFlowTrack_2", "1.0 / ReconstructedParticle::get_p(ReconstructedParticles)")
                # (must be an object of type rv::RVec<edm4hep::Quantity>)
            )

        # do the actual analysis
        dfout = (
            dfout

            # get event ID variables
            .Define("runNumber", "EventHeader.runNumber")
            .Define("eventNumber", "EventHeader.eventNumber")

            # recalculate reco primary vertex
            .Define("PrimaryVertexP4", "fitRecoPrimaryVertex(EFlowTrack_1)")

            # store the primary vertex coordinates
            .Define("PV_x", "PrimaryVertexP4.X()")
            .Define("PV_y", "PrimaryVertexP4.Y()")
            .Define("PV_z", "PrimaryVertexP4.Z()")
            
        )
        return dfout

    def output():
        branchList = [
            'runNumber',
            'eventNumber',
            'PV_x',
            'PV_y',
            'PV_z',
        ]
        
        return branchList    
