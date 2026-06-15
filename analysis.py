import os
import sys
import json
import ROOT


### handling of beamspot data ###

# load beamspot json file
beamspotfile = 'data/beamspot.json'
if os.path.exists(beamspotfile):
    print(f'Loading beamspot data from {beamspotfile}...')
    with open(beamspotfile, 'r') as f:
        beamspot_data = json.load(f)
else:
    msg = 'WARNING: trying to load beamspot data, but expected file {beamspotfile} not found.'
    msg += ' Will skip loading beamspot data.'
    print(msg)
    beamspot_data = None

# declare data struct in ROOT and fill it
ROOT.gInterpreter.Declare("""
    #include <map>
    #include <tuple>
    #include <string>

    static std::map<int, std::tuple<double,double,double>> beamspotCoords;
""")
for run, coords in beamspot_data.items():
    ROOT.beamspotCoords.emplace(
        int(run),
        ROOT.std.make_tuple(float(coords['x']), float(coords['y']), float(coords['z']))
    )

# declare retrieval function in ROOT
ROOT.gInterpreter.Declare("""
    TLorentzVector getBeamspotCoords(int run){
        auto it = beamspotCoords.find(run);
        if (it != beamspotCoords.end()){
            std::tuple<double, double, double> coords = it->second;
            return TLorentzVector(std::get<0>(coords), std::get<1>(coords), std::get<2>(coords), 0.0);
        }
        // fallback if run not found
        return TLorentzVector(0.0, 0.0, 0.0, 0.0);
    }
""")

### end of beamspot handling part ###

# load custom analyzer for particle ID retrieval
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_particleid.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with collections of tracks
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_tracktools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with collections of reco particles
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_recotools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with collections of reco particles
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_jettools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with collections of reco particles
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_vertextools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with primary vertex reconstruction
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_pvtools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with collections of secondary vertices
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_svtools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with secondary vertex fitting tools
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_svfinder.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with dE/dx tools
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_dEdx.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with reco to mc linking tools
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_recotomctools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with reco to track linking tools
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_recototracktools.cxx')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# load custom analyzer with Aleph impact parameter definitions
analyzer_path = os.path.join(os.path.dirname(__file__), 'analyzers', 'analyzer_ipalephtools.cc')
ROOT.gInterpreter.Declare(f'#include "{analyzer_path}"')

# helper function for deriving the gen-level event type.
# note: for now, only valid with qqbar simulations,
#       where the event type is between 1 (d dbar) and 5 (b bbar) (see PDG numbering scheme).
# note: the event type is derived simply from the first quark in the list of MCParticles;
#       there is in principle no guarantee for any kind of ordering;
#       we just assume the first quark PDG ID in the MCParticle collection
#       is the one corresponding to the type of quarks produced in the hard scattering.
#       to be checked and refined later.
# note: in the original analyzer that served as source for this one,
#       the event type needed not to be derived, as the simulation was split per quark flavour,
#       so instead the event type was just derived from the file name.
ROOT.gInterpreter.Declare("""
    int get_genEventType(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles) {
        for (const auto& genParticle : genParticles) {
            int pdgid = std::abs(genParticle.PDG);
            if( (pdgid >= 1) && (pdgid <= 6) ){ return pdgid; }
        }
        return -1;
    }""")

# helper function for deriving the reco-level event type.
# for more info on the class bitset encoding this information,
# see here: https://aleph-new.docs.cern.ch/eos/1994-data/
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<unsigned int> get_recoEventType(
        ROOT::VecOps::RVec<unsigned int> classBitsetInt) {
        std::bitset<32> classBitset = std::bitset<32>(classBitsetInt[0]);
        ROOT::VecOps::RVec<unsigned int> res;
        for(unsigned int shift=0; shift<32; shift++){
            if(classBitset[shift]){ res.push_back(shift+1); }
        }
        return res;
    }""")

# helper function for making a dummy RecoParticle-to-Tracks linking collection.
# this is needed for syntax in case the linking collection does not exist,
# and the linking between RecoParticles and Tracks is direct.
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<podio::ObjectID> makeDummyRecoToTracks(
        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& rps){
        ROOT::VecOps::RVec<podio::ObjectID> links;
        links.reserve(rps.size());
        for (size_t i = 0; i < rps.size(); ++i) {
            const auto& rp = rps[i];
            podio::ObjectID oid;
            oid.index = rp.tracks_begin;
            oid.collectionID = 0;
            links.push_back(oid);
        }
        return links;
    }""")

# helper function to retrieve a vertex from the collection of vertices
# in the correct TLorentzVector format.
# note: not clear how to distinguish primary from secondary vertices,
#       this code just relies on the fact that >99% of events (in simulation) have exactly one vertex stored.
# note: what to do in case no vertex is stored for an event?
#       for now just use the dummy of (0, 0, 0) in those cases.
ROOT.gInterpreter.Declare("""
    TLorentzVector getRecoPrimaryVertex(
        const ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex>& vertexCollection){
        edm4hep::VertexData vertex;
        TLorentzVector result = {0., 0., 0., 0.};
        if( vertexCollection.size() == 0 ){ return result; }
        vertex = vertexCollection.at(0).vertex;
        result = {vertex.position.x, vertex.position.y, vertex.position.z, 0.};
        return result;
    }""")

# helper function to count the number of instances per jet
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<int> count(
        const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>& data){
        ROOT::VecOps::RVec<int> counts;
        for( const ROOT::VecOps::RVec<int>& this_data : data) {
            counts.push_back( this_data.size() );
        }
        return counts;
    }""")

# helper function to count the number of instances per jet
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<int> countInstances(
        const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>& data, int targetValue){
        ROOT::VecOps::RVec<int> counts;
        for( const ROOT::VecOps::RVec<int>& this_data : data) {
            int nValues = 0;
            for( int value : this_data ){
                if( value==targetValue ){ nValues++; }
            }
            counts.push_back( nValues );
        }
        return counts;
    }""")

# helper function to get the MC truth-level primary vertex from the first gen-particle.
# note: relies on the assumption that the vertex stored for the first gen-particle is the PV.
ROOT.gInterpreter.Declare("""
    TLorentzVector getMCPV(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles) {
        TLorentzVector result = {0, 0, 0, 0};
        if( genParticles.size() < 1 ){ return result; }
        result = {genParticles[0].vertex.x, genParticles[0].vertex.y, genParticles[0].vertex.z, 0.};
        return result;
    }""")

# helper function to make dummy jet constituent variables
# (e.g. for variables that are available for aleph but not for fcc or the other way round)
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> makeDummyJetConstituentVariable(
        const ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>>& jcs, float dummyValue){
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>> out;
        for (const auto &jc : jcs){
            ROOT::VecOps::RVec<float> temp;
            for (const auto &el : jc){
                temp.emplace_back(dummyValue);
            }
            out.emplace_back(temp);
        }
        return out;
    }""")


# dirty hard-coded settings using global variables
# (try to do more cleanly later, but doesn't seem very trivial with RDFanalysis)
do_secondary_vertices = True
do_v0_candidates = True


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

        # for Aleph simulation, the link between RecoParticles and Tracks is indirect,
        # via an intermediate collection _RecoParticles_tracks;
        # for FCC simulation, we must define it as a transparent mapping for the syntax
        if det=='aleph': dfout = dfout.Alias("Reco2TrackLinks", "_RecoParticles_tracks")
        elif det=='fcc': dfout = dfout.Define("Reco2TrackLinks", "makeDummyRecoToTracks(ReconstructedParticles)")

        # for Aleph simulation, the collection of MC particles is called "MCParticles",
        # while for FCC simulation it is called "Particle", so define an alias here.
        # and same for "RecoParticles" vs "ReconstructedParticles".
        # and same for "ParticleID" vs "ParticleIDs" (note: not well defined for FCC).
        # note: there are also some modification that need to be made to the Aleph reco particles:
        #       - make the linking between reco particles and track states direct
        #         rather than via some intermediate linking indices collection.
        if det=='aleph':
            if dtype=='sim': dfout = dfout.Alias("Particle", "MCParticles")
            dfout = dfout.Define("ReconstructedParticles", "RecoTools::getModifiedRecoParticles(RecoParticles, Reco2TrackLinks)")
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
                # (must be an object of type ROOT::VecOps::RVec<edm4hep::TrackData>)
                .Define("EFlowTrack_1", "TrackTools::getModifiedTrackStates(_Tracks_trackStates)")
                # (must be an object of type ROOT::VecOps::RVec<edm4hep::TrackState>)
                .Define("EFlowTrack_2", "1.0 / ReconstructedParticle::get_p(ReconstructedParticles)")
                # (must be an object of type rv::RVec<edm4hep::Quantity>)
            )

        # add event identifiers
        dfout = (
            dfout
            .Define("runNumber", "EventHeader.runNumber")
            .Define("eventNumber", "EventHeader.eventNumber")
        )

        # do gen-level stuff
        if dtype=='sim':
            dfout = (
                dfout
                
                # store the pdg ID and generator status for all generator particles
                # (mainly for debugging)
                .Define("GenParticle_pdgId", "MCParticle::get_pdg(Particle)")
                .Define("GenParticle_genStatus", "MCParticle::get_genStatus(Particle)")

                # get event type (at generator level)
                .Define("genEventType", "get_genEventType(Particle)")

                # get true primary vertex position
                .Define("GenPrimaryVertexP4", "getMCPV(Particle)")
                .Define("GenPV_x", "GenPrimaryVertexP4.X()")
                .Define("GenPV_y", "GenPrimaryVertexP4.Y()")
                .Define("GenPV_z", "GenPrimaryVertexP4.Z()")
            )
        else:
            # set dummies (maybe later find out how to avoid the need for this)
            dfout = (
                dfout
                .Define("GenParticle_pdgId", "0")
                .Define("GenParticle_genStatus", "0")
                .Define("genEventType", "-1")
                .Define("GenPV_x", "-1")
                .Define("GenPV_y", "-1")
                .Define("GenPV_z", "-1")
            )

        # do reco-level event type
        # (only well-defined for Aleph, put a dummy for FCC)
        if det=='aleph':
            dfout = dfout.Define("recoEventType", "get_recoEventType(ClassBitset)")
        elif det=='fcc':
            dfout = dfout.Define("recoEventType", "-1")

        # get the beamspot position
        if dtype=='data':
            dfout = (
                dfout
                .Define("BeamspotP4", "getBeamspotCoords(EventHeader.runNumber[0])")
                .Define("Beamspot_x", "BeamspotP4.X()")
                .Define("Beamspot_y", "BeamspotP4.Y()")
                .Define("Beamspot_z", "BeamspotP4.Z()")
            )
        else:
            # set to 0 for simulation
            dfout = (
                dfout
                .Define("Beamspot_x", "0.0")
                .Define("Beamspot_y", "0.0")
                .Define("Beamspot_z", "0.0")
            )

        # do track preselection
        dfout = (
            dfout

            # store number of reco particles and tracks per event
            .Define("Event_nRecoParticles", "ReconstructedParticles.size()")
            .Define("Event_nTracks", "EFlowTrack_1.size()")
            .Define("TracksWithReco", "RecoToTrackTools::getRP2TRK_trackState(ReconstructedParticles, EFlowTrack_1)")
            .Define("TracksWithoutReco", "TrackTools::get_complement(EFlowTrack_1, TracksWithReco)")
            .Define("TrackDatasWithReco", "RecoToTrackTools::getRP2TRK_track(ReconstructedParticles, EFlowTrack)")
            .Define("TrackDatasWithoutReco", "TrackTools::get_complement(EFlowTrack, EFlowTrack_1, TrackDatasWithReco, TracksWithReco)")
            .Define("Event_nTracksWithReco", "TracksWithReco.size()")
            .Define("Event_nTracksWithoutReco", "TracksWithoutReco.size()")

            # do track preselection
            .Define("SelectedTracks", "TrackTools::getSelectedTracks(EFlowTrack_1, EFlowTrack, -1, -1)")
        )

        # find the primary vertex
        dfout = (
            dfout

            # get MC primary vertex
            #.Define("PrimaryVertexP4", "FCCAnalyses::MCParticle::get_EventPrimaryVertexP4()(Particle)")

            # alternative: get MC primary vertex from first particle in list of gen particles
            #.Define("PrimaryVertexP4", "getMCPV(Particle)")

            # alternative for running on data: just use a dummy.
            # note: maybe later try to switch to actual reco primary vertex.
            #.Define("PrimaryVertexP4", "TLorentzVector(0.,0.,0.,0.)")

            # alternative for running on data or circumventing other issues with the MC primary vertex:
            # use reco primary vertex.
            # note: not sure how to make this work for FCC sim; there doesn't seem to be an equivalent collection.
            #.Define("PrimaryVertexP4", "getRecoPrimaryVertex(Vertices)")
            
            # alternative: recalculate reco primary vertex
            .Define("PrimaryTracks", "PrimaryVertexTools::getPrimaryTracks(SelectedTracks, 5., Beamspot_x, Beamspot_y, Beamspot_z)")
            .Define("PrimaryVertexObject", "PrimaryVertexTools::fitRecoPrimaryVertex(PrimaryTracks, Beamspot_x, Beamspot_y, Beamspot_z)")
            .Define("PrimaryVertex", "FCCAnalyses::VertexingUtils::get_VertexData(PrimaryVertexObject)")
            .Define("PrimaryVertexP4", "TLorentzVector(PrimaryVertex.position.x, PrimaryVertex.position.y, PrimaryVertex.position.z, 0.)")

            # store the primary vertex coordinates
            # (mainly for debugging)
            .Define("PV_x", "PrimaryVertexP4.X()")
            .Define("PV_y", "PrimaryVertexP4.Y()")
            .Define("PV_z", "PrimaryVertexP4.Z()")

            # make a copy of the primary vertex as a 3-vector rather than a 4-vector
            # (needed in some functions below)
            .Define("PrimaryVertexP3", "TVector3(PV_x, PV_y, PV_z)")

            # store properties of fitted vertex
            # (note: only works for fitted vertex, not for other methods)
            .Define("PV_chi2Normalized", "PrimaryVertex.chi2")
            .Define("PV_ndof", "PrimaryVertex.ndf")
            .Define("Event_nSelectedTracks", "SelectedTracks.size()")
            .Define("Event_nPrimaryTracks", "PrimaryTracks.size()")
        )

        # do jet clustering
        dfout = (
            dfout

            # define the momentum, energy, mass and charge of all reconstructed particles.
            .Define("RP_px",          "ReconstructedParticle::get_px(ReconstructedParticles)")
            .Define("RP_py",          "ReconstructedParticle::get_py(ReconstructedParticles)")
            .Define("RP_pz",          "ReconstructedParticle::get_pz(ReconstructedParticles)")
            .Define("RP_e",           "ReconstructedParticle::get_e(ReconstructedParticles)")
            .Define("RP_m",           "ReconstructedParticle::get_mass(ReconstructedParticles)")
            .Define("RP_q",           "ReconstructedParticle::get_charge(ReconstructedParticles)")
            
            # build "pseudo-jets", meaning each particle is converted to a jet
            # consisting of only that one particle (which can then be clustered in the next step)
            .Define("pseudo_jets",    "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
            
            # run jet clustering
            # note: the JetClustering namespace with its functions is defined here:
            #       FCCAnalyses/addons/FastJet/src/JetClustering.cc.
            # note: the arguments are (in order):
            # - radius: jet radius (typically 1.5)
            # - exclusive: inclusive or exclusive clustering.
            #   "inclusive" = variable number of jets, fixed radius.
            #   "exclusive" = fixed number of jets (typically 2, i.e. hemispheres), variable radius.
            #   in the case of exclusive clustering, the jet radius argument is ignored (?).
            #   use 0 for inclusive clustering, 1 or 2 for exclusive clustering (?)
            #   (see FCCAnalyses/analyzers/dataframe/src/JetClusteringUtils.cc / build_jets)
            # - cut: seems to be a pt threshold for inclusive clustering,
            #   and the number of jets in exclusive clustering (?)
            #   (see FCCAnalyses/analyzers/dataframe/src/JetClusteringUtils.cc / build_jets)
            #   (see https://fastjet.fr/repo/doxygen-3.0.0/classfastjet_1_1ClusterSequence.html)
            # - sorted: jet sorting method; allowed values are:
            #     - 0: sort by pt (default)
            #     - 1: sort by energy
            #   (see FCCAnalyses/analyzers/dataframe/src/JetClusteringUtils.cc / build_jets)
            # - recombination: how jets are summed (leave at default 0).
            # - exponent (of kT algorithm, typically -1, i.e. anti-kT).
            #.Define("FCCAnalysesJets_ee_genkt", "JetClustering::clustering_ee_genkt(1.5, 0, 0, 0, 0, -1)(pseudo_jets)") # inclusive
            #.Define("FCCAnalysesJets_ee_genkt", "JetClustering::clustering_ee_genkt(1.5, 3, 2, 0, 0, -1)(pseudo_jets)") # exclusive

            # alternative jet clustering in sync with CERN team
            # see here: https://github.com/Apranikstar/Aleph/blob/6fe1e52a50a22d3f1e0f978ca1c0e747488bb016/src/stage1.py#L128
            # note: the arguments are (in order):
            # - exclusive
            # - cut
            # - sorted
            # - recombination
            .Define("FCCAnalysesJets_ee_genkt", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")

            # get the jets out of the struct
            .Define("jets_ee_genkt", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_genkt)")

            # get merging parameters out of the struct
            .Define("Event_dmerge1", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt, 1)")
            .Define("Event_dmerge2", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt, 2)")
            .Define("Event_dmerge3", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt, 3)")
            .Define("Event_dmerge4", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt, 4)")

            # define jet-level observables
            .Define("Jets_px", "JetClusteringUtils::get_px(jets_ee_genkt)")
            .Define("Jets_py", "JetClusteringUtils::get_py(jets_ee_genkt)")
            .Define("Jets_pz", "JetClusteringUtils::get_pz(jets_ee_genkt)")
            .Define("Jets_pt", "JetClusteringUtils::get_pt(jets_ee_genkt)")
            .Define("Jets_e", "JetClusteringUtils::get_e(jets_ee_genkt)")
            .Define("Jets_mass", "JetClusteringUtils::get_m(jets_ee_genkt)")
            .Define("Jets_phi", "JetClusteringUtils::get_phi(jets_ee_genkt)")
            .Define("Jets_eta", "JetClusteringUtils::get_eta(jets_ee_genkt)")
            .Define("Jets_theta", "JetClusteringUtils::get_theta(jets_ee_genkt)")
            .Define("Jets_p4", "JetConstituentsUtils::compute_tlv_jets(jets_ee_genkt)")

            # find the reconstructed particles grouped per jet
            # (output struct is a vector of vectors of ReconstructedParticle objects, one vector of ReconstructedParticles for each jet)
            .Define("jetconstituents_ee_genkt", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_genkt)")
            .Define("JetsConstituents", "JetConstituentsUtils::build_constituents_cluster(ReconstructedParticles, jetconstituents_ee_genkt)")
            .Define("Jets_nConstituents", "JetConstituentsUtils::count_consts(JetsConstituents)")

            # find track states grouped per jet, and perform baseline selection
            # (output struct is a vector of vectors of TrackState objects, one vector of TrackStates for each jet)
            # note: the grouping function implicitly only considers tracks connected to a reco particle,
            #       so we have to manually append the tracks without reco particle (based on simple delta-R matching to the closest jet.
            .Define("TracksPerJet", "JetTools::build_trackstates_cluster(ReconstructedParticles, EFlowTrack_1, jetconstituents_ee_genkt)")
            .Redefine("TracksPerJet", "TrackTools::append_losttracks_to_tracksperjet(TracksWithoutReco, TracksPerJet, jets_ee_genkt)")
            .Define("TrackDatasPerJet", "JetTools::build_tracks_cluster(ReconstructedParticles, EFlowTrack, jetconstituents_ee_genkt)")
            .Redefine("TrackDatasPerJet", "TrackTools::append_losttracks_to_tracksperjet(TrackDatasWithoutReco, TracksWithoutReco, TrackDatasPerJet, jets_ee_genkt)")
            .Define("SelectedTracksPerJet", "TrackTools::getSelectedTracks(TracksPerJet, TrackDatasPerJet, -1, -1)")
            .Define("Jets_nTracksPerJet", "TrackTools::countTracks(TracksPerJet)")
            .Define("Jets_nSelectedTracksPerJet", "TrackTools::countTracks(SelectedTracksPerJet)")

            # for debugging: check that all particles and all tracks in the event were clustered
            .Define("Event_nConstituentsSum", "ROOT::VecOps::Sum(Jets_nConstituents)")
            .Define("Event_nTracksPerJetSum", "ROOT::VecOps::Sum(Jets_nTracksPerJet)")
        )

        # some more event level variables
        # (that rely on the above jet clustering info)
        dfout = (
            dfout

            # general properties
            .Define("Event_mass", "JetConstituentsUtils::InvariantMass(Jets_p4[0], Jets_p4[1])")
            .Define("Event_njets", "(int)Jets_p4.size()")
            .Define("Event_Bz", "RecoToTrackTools::Bz(ReconstructedParticles, EFlowTrack_1)")

            # compute the residues for jet-constituents on significant kinematic variables as a check
            # notes:
            # - "tlv_jets" seems to mean: "the lorentz vectors of the jets, calculated directly from the jets"
            # - "sum_tlv_jcs" seems to mean: "the lorentz vectors of the jets, but calculated by summing all constituents"
            .Define("tlv_jets", "JetConstituentsUtils::compute_tlv_jets(jets_ee_genkt)")
            .Define("sum_tlv_jcs", "JetConstituentsUtils::sum_tlv_constituents(JetsConstituents)")
            .Define("Event_de", "JetConstituentsUtils::compute_residue_energy(tlv_jets, sum_tlv_jcs)")
            .Define("Event_dpt", "JetConstituentsUtils::compute_residue_pt(tlv_jets, sum_tlv_jcs)")
            .Define("Event_dphi", "JetConstituentsUtils::compute_residue_phi(tlv_jets, sum_tlv_jcs)")
            .Define("Event_dtheta", "JetConstituentsUtils::compute_residue_theta(tlv_jets, sum_tlv_jcs)")
        )

        # find secondary tracks (for later use in secondary vertex finding)
        dfout = (
            dfout

            # find tracks incompatible with the primary vertex (both per event and per jet).
            # note: primary tracks have already been defined before (when calculating the primary vertex),
            #       but an alternative set is re-calculated here, possibly with a tighter chi2 cut.
            .Define("SecondaryTracks", "TrackTools::getSecondaryTracks(SelectedTracks, PrimaryTracks)")
            .Define("SecondaryTracksPerJet", "TrackTools::getSecondaryTracks(SelectedTracksPerJet, PrimaryTracks)")

            # store counters
            .Define("Event_nSecondaryTracks", "SecondaryTracks.size()")
            .Define("Jets_nSecondaryTracksPerJet", "TrackTools::countTracks(SecondaryTracksPerJet)")
        )

        # find secondary vertices (per jet)
        if do_secondary_vertices:
          dfout = (
            dfout

            # fit secondary vertices per jet.
            # (output struct is a vector of vectors of FCCAnalysesVertex objects, one vector for each jet).
            # the input arguments are (in order):
            # - tracks per jet to use
            # - all tracks
            # - primary vertex
            # - whether to do V0 rejection
            # - chi2 upper bound for vertices
            # - invariant mass upper bound for vertices
            # - track chi2 cut (what does this do exactly?)
            #.Define("SecondaryVertices", "SecondaryVertexFinder::get_SV_jets(SecondaryTracksPerJet, EFlowTrack_1, PrimaryVertexObject, true, 10., 10., 5.)")
            
            # alternative: fit secondary vertices per event (ignoring jets), and associate them to tracks post-rem.
            # note: not sure how to do this best, but for now do simple dR matching.
            .Define("EventSecondaryVertices", "SecondaryVertexFinder::get_SV_event(SecondaryTracks, EFlowTrack_1, PrimaryVertexObject, true, 10., 10., 5.)")
            .Define("SecondaryVertices", "SecondaryVertexTools::distributeSecondaryVerticesOverJets(EventSecondaryVertices, jets_ee_genkt)")

            # calculate properties of secondary vertices to store
            .Define("SecondaryVertices_xrel", "VertexTools::get_xrel_SV_jets(SecondaryVertices, PrimaryVertexP3)")
            .Define("SecondaryVertices_yrel", "VertexTools::get_yrel_SV_jets(SecondaryVertices, PrimaryVertexP3)")
            .Define("SecondaryVertices_zrel", "VertexTools::get_zrel_SV_jets(SecondaryVertices, PrimaryVertexP3)")
            .Define("SecondaryVertices_thetarel", "FCCAnalyses::VertexingUtils::get_relTheta_SV(SecondaryVertices, jets_ee_genkt)")
            .Define("SecondaryVertices_phirel", "FCCAnalyses::VertexingUtils::get_relPhi_SV(SecondaryVertices, jets_ee_genkt)")
            .Define("SecondaryVertices_p", "FCCAnalyses::VertexingUtils::get_pMag_SV(SecondaryVertices)")
            .Define("SecondaryVertices_prel", "VertexTools::get_prel_SV_jets(SecondaryVertices, jets_ee_genkt)")
            .Define("SecondaryVertices_chi2", "FCCAnalyses::VertexingUtils::get_chi2_SV(SecondaryVertices)")
            .Define("SecondaryVertices_chi2Normalized", "FCCAnalyses::VertexingUtils::get_norm_chi2_SV(SecondaryVertices)")
            .Define("SecondaryVertices_ndof", "FCCAnalyses::VertexingUtils::get_nDOF_SV(SecondaryVertices)")
            .Define("SecondaryVertices_nTracks", "FCCAnalyses::VertexingUtils::get_VertexNtrk(SecondaryVertices)")
            .Define("SecondaryVertices_mass", "FCCAnalyses::VertexingUtils::get_invM(SecondaryVertices)")
            .Define("SecondaryVertices_dxy", "FCCAnalyses::VertexingUtils::get_dxy_SV(SecondaryVertices, PrimaryVertexObject)")
            .Define("SecondaryVertices_dxyz", "FCCAnalyses::VertexingUtils::get_d3d_SV(SecondaryVertices, PrimaryVertexObject)")
            .Define("SecondaryVertices_cosPointing", "VertexTools::get_pointingangle_SV(SecondaryVertices, PrimaryVertexObject)")
            .Define("SecondaryVertices_correctedMass", "SecondaryVertexTools::getCorrectedInvMass(SecondaryVertices, PrimaryVertexObject)")

            # store counters
            .Define("Event_nSV", "EventSecondaryVertices.size()")
            .Define("Jets_nSV", "FCCAnalyses::VertexingUtils::get_n_SV_jets(SecondaryVertices)")
          )

        # find V0s (per jet)
        if do_v0_candidates:
          dfout = (
            dfout

            # reconstruct all V0 candidates in the event
            # note: arguments are (in order): tracks, primary vertex, whether to use tight constraints, chi2 threshold
            .Define("EventV0s", "SecondaryVertexFinder::get_V0s(SecondaryTracks, PrimaryVertexObject, false, 10.)")
            .Define("V0Candidates", "SecondaryVertexTools::distributeSecondaryVerticesOverJets(EventV0s.vtx, jets_ee_genkt)")
            .Define("v0_jet_ids", "SecondaryVertexTools::getSecondaryVertexJetIndices(EventV0s.vtx, jets_ee_genkt)")

            # calculate properties of V0 candidates to store
            .Define("V0Candidates_pdgId", "SecondaryVertexTools::distributeOverJets(EventV0s.pdgAbs, v0_jet_ids)")
            .Define("V0Candidates_xrel", "VertexTools::get_xrel_SV_jets(V0Candidates, PrimaryVertexP3)")
            .Define("V0Candidates_yrel", "VertexTools::get_yrel_SV_jets(V0Candidates, PrimaryVertexP3)")
            .Define("V0Candidates_zrel", "VertexTools::get_zrel_SV_jets(V0Candidates, PrimaryVertexP3)")
            .Define("V0Candidates_thetarel", "FCCAnalyses::VertexingUtils::get_relTheta_SV(V0Candidates, jets_ee_genkt)")
            .Define("V0Candidates_phirel", "FCCAnalyses::VertexingUtils::get_relPhi_SV(V0Candidates, jets_ee_genkt)")
            .Define("V0Candidates_p", "FCCAnalyses::VertexingUtils::get_pMag_SV(V0Candidates)")
            .Define("V0Candidates_prel", "VertexTools::get_prel_SV_jets(V0Candidates, jets_ee_genkt)")
            .Define("V0Candidates_chi2", "FCCAnalyses::VertexingUtils::get_chi2_SV(V0Candidates)")
            .Define("V0Candidates_chi2Normalized", "FCCAnalyses::VertexingUtils::get_norm_chi2_SV(V0Candidates)")
            .Define("V0Candidates_ndof", "FCCAnalyses::VertexingUtils::get_nDOF_SV(V0Candidates)")
            .Define("V0Candidates_nTracks", "FCCAnalyses::VertexingUtils::get_VertexNtrk(V0Candidates)")
            .Define("V0Candidates_mass", "SecondaryVertexTools::distributeOverJets(EventV0s.invM, v0_jet_ids)")
            .Define("V0Candidates_dxy", "FCCAnalyses::VertexingUtils::get_dxy_SV(V0Candidates, PrimaryVertexObject)")
            .Define("V0Candidates_dxyz", "FCCAnalyses::VertexingUtils::get_d3d_SV(V0Candidates, PrimaryVertexObject)")
            .Define("V0Candidates_cosPointing", "VertexTools::get_pointingangle_SV(V0Candidates, PrimaryVertexObject)")
            .Define("V0Candidates_correctedMass", "SecondaryVertexTools::getCorrectedInvMass(V0Candidates, PrimaryVertexObject)")

            # store counters
            .Define("Event_nV0Candidates", "EventV0s.vtx.size()")
            .Define("Jets_nV0Candidates", "count(v0_jet_ids)")
            .Define("Jets_nKsCandidates", "countInstances(V0Candidates_pdgId, 310)")
            .Define("Jets_nLambdaCandidates", "countInstances(V0Candidates_pdgId, 3122)")

          )

        # store the PDG ID of every jet constituent
        # (or at least every particle for which it is available, see more details in the helper function)
        if dtype=='sim':
            dfout = (
                dfout
                .Define("TrackToMCMap", "RecoToMCTools::makeTrackToMCMapping(EFlowTrack, _trackMCLink_to, _trackMCLink_from)")
                .Define("JetsConstituents_pdgId", "RecoToMCTools::get_pdgid(JetsConstituents, EFlowTrack, Particle, TrackToMCMap)")
            )

        else:
            # set dummies (maybe later find out how to avoid the need for this)
            dfout = (
                dfout
                .Define("JetsConstituents_pdgId", "RecoToMCTools::get_pdgidDummy(JetsConstituents)")
            )

        # rest of the analysis
        dfout = (
            dfout

            # Extract ParticleID types for all jet constituents
            # ParticleID.type legend: 0:Track, 1:Electron, 2:Muon, 3:Track from V0, 
            #                         4:EM, 5:Ecal hadron/residual, 6:Hcal element, 7:Lcal element
            .Define("JetsConstituents_Types", "getParticleIDTypes(ParticleIDs, jetconstituents_ee_genkt)")
            
            # Convert the ParticleID types to binary flags for particle classification
            # Muon: type 2
            # Electron: type 1
            # Gamma: type 4 (EM)
            # Charged hadron: types 0 (Track), 3 (Track from V0)
            # Neutral hadron: types 5 (Ecal hadron/residual), 6 (Hcal element), 7 (Lcal element)
            .Define("JetsConstituents_isMu", "get_isMu_from_type(JetsConstituents_Types)")
            .Define("JetsConstituents_isEl", "get_isEl_from_type(JetsConstituents_Types)")
            .Define("JetsConstituents_isChargedHad", "get_isChargedHad_from_type(JetsConstituents_Types)")
            .Define("JetsConstituents_isGamma", "get_isGamma_from_type(JetsConstituents_Types)")
            .Define("JetsConstituents_isNeutralHad", "get_isNeutralHad_from_type(JetsConstituents_Types)")

            # counting the types of particles per jet
            .Define("Jets_nMu", "JetConstituentsUtils::count_type(JetsConstituents_isMu)")
            .Define("Jets_nEl", "JetConstituentsUtils::count_type(JetsConstituents_isEl)")
            .Define("Jets_nChargedHad", "JetConstituentsUtils::count_type(JetsConstituents_isChargedHad)")
            .Define("Jets_nPhoton", "JetConstituentsUtils::count_type(JetsConstituents_isGamma)")
            .Define("Jets_nNeutralHad", "JetConstituentsUtils::count_type(JetsConstituents_isNeutralHad)")

            # basic kinematics
            .Define("JetsConstituents_e", "JetConstituentsUtils::get_e(JetsConstituents)")
            .Define("JetsConstituents_pt", "JetConstituentsUtils::get_pt(JetsConstituents)")
            .Define("JetsConstituents_px", "JetTools::get_px(JetsConstituents)")
            .Define("JetsConstituents_py", "JetTools::get_py(JetsConstituents)")
            .Define("JetsConstituents_pz", "JetTools::get_pz(JetsConstituents)")
            .Define("JetsConstituents_theta", "JetConstituentsUtils::get_theta(JetsConstituents)")
            .Define("JetsConstituents_phi", "JetConstituentsUtils::get_phi(JetsConstituents)")
            .Define("JetsConstituents_charge", "JetConstituentsUtils::get_charge(JetsConstituents)")

            # kinematics relative to jet
            .Define("JetsConstituents_erel", "JetConstituentsUtils::get_erel_cluster(jets_ee_genkt, JetsConstituents)")
            .Define("JetsConstituents_erel_log", "JetConstituentsUtils::get_erel_log_cluster(jets_ee_genkt, JetsConstituents)")
            .Define("JetsConstituents_ptrel", "JetTools::get_ptrel_cluster(jets_ee_genkt, JetsConstituents)")
            .Define("JetsConstituents_ptrel_log", "JetTools::get_ptrel_log_cluster(jets_ee_genkt, JetsConstituents)")
            .Define("JetsConstituents_thetarel", "JetConstituentsUtils::get_thetarel_cluster(jets_ee_genkt, JetsConstituents)")
            .Define("JetsConstituents_phirel", "JetConstituentsUtils::get_phirel_cluster(jets_ee_genkt, JetsConstituents)") 
            
            # PID variables
            #.Define("JetsConstituents_dndx", "JetConstituentsUtils::get_dndx(JetsConstituents, EFlowTrack_2, EFlowTrack, JetsConstituents_isChargedHad)")
            #.Define("JetsConstituents_mtof", "JetConstituentsUtils::get_mtof(JetsConstituents, EFlowTrack_L, EFlowTrack, TrackerHits, JetsConstituents_Pids)")
           
            # track properties
            .Define("JetsConstituents_trackChi2", "JetTools::get_chi2(JetsConstituents, EFlowTrack)")
            .Define("JetsConstituents_trackNdof", "JetTools::get_ndof(JetsConstituents, EFlowTrack)")
            .Define("JetsConstituents_trackChi2Normalized", "JetTools::get_chi2Normalized(JetsConstituents, EFlowTrack)")
        )

        # dE/dx values
        dfout = (

            dfout
            .Define("dEdxPadsValue", "dEdxPads.dQdx.value")
            .Define("dEdxPadsError", "dEdxPads.dQdx.error")
            .Define("dEdxWiresValue", "dEdxWires.dQdx.value")
            .Define("dEdxWiresError", "dEdxWires.dQdx.error")

            # basic variables for pads
            .Define("jet_constituents_dEdx_pads_objs", "dEdxTools::build_constituents_dEdx()(ReconstructedParticles, dEdxPads, _dEdxPads_track.index, jetconstituents_ee_genkt, false)")
            .Define("jet_constituents_dEdx_pads_dedx", "jet_constituents_dEdx_pads_objs.dedx_constituents")
            .Define("JetsConstituents_dEdx_pads_type", "dEdxTools::get_dEdx_type(jet_constituents_dEdx_pads_dedx)")
            .Define("JetsConstituents_dEdx_pads_value", "dEdxTools::get_dEdx_value(jet_constituents_dEdx_pads_dedx)")
            .Define("JetsConstituents_dEdx_pads_error", "dEdxTools::get_dEdx_error(jet_constituents_dEdx_pads_dedx)")

            # PID hypotheses compatibility for pads
            .Define("jet_constituents_dEdx_pads_pvals", "jet_constituents_dEdx_pads_objs.pid_array_constituents")
            .Define("JetsConstituents_PID_pval_pads_ele", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_pads_pvals, 0)")
            .Define("JetsConstituents_PID_pval_pads_mu", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_pads_pvals, 1)")
            .Define("JetsConstituents_PID_pval_pads_pi", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_pads_pvals, 2)")
            .Define("JetsConstituents_PID_pval_pads_kaon", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_pads_pvals, 3)")
            .Define("JetsConstituents_PID_pval_pads_proton", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_pads_pvals, 4)")

            # basic variables for wires
            .Define("jet_constituents_dEdx_wires_objs", "dEdxTools::build_constituents_dEdx()(ReconstructedParticles, dEdxWires, _dEdxWires_track.index, jetconstituents_ee_genkt, true)")
            .Define("jet_constituents_dEdx_wires_dedx", "jet_constituents_dEdx_wires_objs.dedx_constituents")
            .Define("JetsConstituents_dEdx_wires_type", "dEdxTools::get_dEdx_type(jet_constituents_dEdx_wires_dedx)")
            .Define("JetsConstituents_dEdx_wires_value", "dEdxTools::get_dEdx_value(jet_constituents_dEdx_wires_dedx)")
            .Define("JetsConstituents_dEdx_wires_error", "dEdxTools::get_dEdx_error(jet_constituents_dEdx_wires_dedx)")

            # PID hypotheses compatibility for wires
            .Define("jet_constituents_dEdx_wires_pvals", "jet_constituents_dEdx_wires_objs.pid_array_constituents")
            .Define("JetsConstituents_PID_pval_wires_ele", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_wires_pvals, 0)")
            .Define("JetsConstituents_PID_pval_wires_mu", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_wires_pvals, 1)")
            .Define("JetsConstituents_PID_pval_wires_pi", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_wires_pvals, 2)")
            .Define("JetsConstituents_PID_pval_wires_kaon", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_wires_pvals, 3)")
            .Define("JetsConstituents_PID_pval_wires_proton", "dEdxTools::get_PID_pvalue(jet_constituents_dEdx_wires_pvals, 4)")
        )

        # number of hits in tracking detectors
        if det=='fcc':
            dfout = (
                dfout
                .Define("JetsConstituents_nTrackHits_VDET", "makeDummyJetConstituentVariable(JetsConstituents, 0)")
                .Define("JetsConstituents_nTrackHits_ITC", "makeDummyJetConstituentVariable(JetsConstituents, 0)")
                .Define("JetsConstituents_nTrackHits_TPC", "makeDummyJetConstituentVariable(JetsConstituents, 0)")
            )

        elif det=='aleph':
            dfout = (
                dfout
                .Define("JetsConstituents_nTrackHits_VDET", "JetTools::get_nTrackHits_VDET(JetsConstituents, EFlowTrack, _Tracks_subdetectorHitNumbers)")
                .Define("JetsConstituents_nTrackHits_ITC", "JetTools::get_nTrackHits_ITC(JetsConstituents, EFlowTrack, _Tracks_subdetectorHitNumbers)")
                .Define("JetsConstituents_nTrackHits_TPC", "JetTools::get_nTrackHits_TPC(JetsConstituents, EFlowTrack, _Tracks_subdetectorHitNumbers)")
            )
            
        dfout = (
            dfout

            # store some track parameters with respect to the nominal origin
            # (mainly for debugging; typically these variables should be re-calculated w.r.t. the primary vertex)
            # note: the parameters have the following meaning:
            #  - d0: transverse impact parameter, i.e. signed transverse distance of closest approach of track to origin
            #  - z0: longitudinal impact parameter, i.e. z-coordinate of the point of closest approach of the track to origin
            #  - phi0: azimuthal angle of track at the point of closest approach of the track to the origin
            #  - omega: track curvature (does not depend on reference point?)
            #  - tan(lambda): pz / pT (related to theta) (does not depend on reference point?)
            # note: the functions below don't do any calculations, they are just getters for the values stored.
            .Define("JetsConstituents_d0_wrt0", "JetConstituentsUtils::get_d0(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_z0_wrt0", "JetConstituentsUtils::get_z0(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_phi0_wrt0", "JetConstituentsUtils::get_phi0(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_omega_wrt0", "JetConstituentsUtils::get_omega(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_tanlambda_wrt0", "JetConstituentsUtils::get_tanLambda(JetsConstituents, EFlowTrack_1)")

            # calculate the magnetic field strength along the z-axis from the curvature of the tracks
            .Define("JetsConstituents_Bz", "JetTools::get_Bz(JetsConstituents, EFlowTrack_1)")
        
            # impact parameters (with respect to the primary vertex) and related quantities (including their uncertainties)
            .Define("JetsConstituents_dxy", "JetTools::XPtoPar_dxy(JetsConstituents, EFlowTrack_1, PrimaryVertexP4, Event_Bz)")
            .Define("JetsConstituents_dz", "JetTools::XPtoPar_dz(JetsConstituents, EFlowTrack_1, PrimaryVertexP4, Event_Bz)")
            .Define("JetsConstituents_phi0", "JetTools::XPtoPar_phi(JetsConstituents, EFlowTrack_1, PrimaryVertexP4, Event_Bz)")
            .Define("JetsConstituents_C", "JetConstituentsUtils::XPtoPar_C(JetsConstituents, EFlowTrack_1, Event_Bz)")
            .Define("JetsConstituents_ct", "JetConstituentsUtils::XPtoPar_ct(JetsConstituents, EFlowTrack_1, Event_Bz)")
            .Define("JetsConstituents_omega_cov", "JetConstituentsUtils::get_omega_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_d0_cov", "JetConstituentsUtils::get_d0_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_z0_cov", "JetConstituentsUtils::get_z0_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_phi0_cov", "JetConstituentsUtils::get_phi0_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_tanlambda_cov", "JetConstituentsUtils::get_tanlambda_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_d0_z0_cov", "JetConstituentsUtils::get_d0_z0_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_phi0_d0_cov", "JetConstituentsUtils::get_phi0_d0_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_phi0_z0_cov", "JetConstituentsUtils::get_phi0_z0_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_tanlambda_phi0_cov", "JetConstituentsUtils::get_tanlambda_phi0_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_tanlambda_d0_cov", "JetConstituentsUtils::get_tanlambda_d0_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_tanlambda_z0_cov", "JetConstituentsUtils::get_tanlambda_z0_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_omega_tanlambda_cov", "JetConstituentsUtils::get_omega_tanlambda_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_omega_phi0_cov", "JetConstituentsUtils::get_omega_phi0_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_omega_d0_cov", "JetConstituentsUtils::get_omega_d0_cov(JetsConstituents, EFlowTrack_1)")
            .Define("JetsConstituents_omega_z0_cov", "JetConstituentsUtils::get_omega_z0_cov(JetsConstituents, EFlowTrack_1)")
            
            # FCC-style signed impact parameters
            .Define("JetsConstituents_Sip2dVal", "JetConstituentsUtils::get_Sip2dVal_clusterV(jets_ee_genkt, JetsConstituents_dxy, JetsConstituents_phi0, 0.0)")
            .Define("JetsConstituents_Sip2dSig", "JetConstituentsUtils::get_Sip2dSig(JetsConstituents_Sip2dVal, JetsConstituents_d0_cov)")
            .Define("JetsConstituents_Sip3dVal", "JetConstituentsUtils::get_Sip3dVal_clusterV(jets_ee_genkt, JetsConstituents_dxy, JetsConstituents_dz, JetsConstituents_phi0, 0.0)")
            .Define("JetsConstituents_Sip3dSig", "JetConstituentsUtils::get_Sip3dSig(JetsConstituents_Sip3dVal, JetsConstituents_d0_cov, JetsConstituents_z0_cov)")
            .Define("JetsConstituents_JetDistVal", "JetConstituentsUtils::get_JetDistVal_clusterV(jets_ee_genkt, JetsConstituents, JetsConstituents_dxy, JetsConstituents_dz, JetsConstituents_phi0, 0.0)")
            .Define("JetsConstituents_JetDistSig", "JetConstituentsUtils::get_JetDistSig(JetsConstituents_JetDistVal, JetsConstituents_d0_cov, JetsConstituents_z0_cov)")

            # ALEPH-style signed impact parameters
            .Define("trackToJetAxisPCAs", "IPAlephTools::getTrackToJetAxisPCA(jets_ee_genkt, JetsConstituents, EFlowTrack_1, PrimaryVertexP3)")
            .Define("JetsConstituents_linearSignedIP3D", "IPAlephTools::signedIP3D(trackToJetAxisPCAs, jets_ee_genkt, PrimaryVertexP3)")
            .Define("JetsConstituents_linearSignedIP3DSig", "IPAlephTools::signedIP3DSig(trackToJetAxisPCAs, jets_ee_genkt, PrimaryVertexP3)")
            .Define("JetsConstituents_transverseJetDistance", "IPAlephTools::transverseJetDistance(trackToJetAxisPCAs)")
            .Define("JetsConstituents_longitudinalJetDistance", "IPAlephTools::longitudinalJetDistance(trackToJetAxisPCAs, PrimaryVertexP3)")

            # extension of the above: store coordinates of PCAs (only for checking and debugging)
            .Define("trackPCAToJetAxis", "IPAlephTools::getTrackPCAToJetAxis(trackToJetAxisPCAs)")
            .Define("jetAxisPCAToTrack", "IPAlephTools::getJetAxisPCAToTrack(trackToJetAxisPCAs)")
            .Define("linePCAToPrimaryVertex", "IPAlephTools::getLinePCAToPrimaryVertex(trackToJetAxisPCAs)")
            .Define("JetsConstituents_trackPCAToJetAxis_x", "IPAlephTools::getPCA_x(trackPCAToJetAxis)")
            .Define("JetsConstituents_trackPCAToJetAxis_y", "IPAlephTools::getPCA_y(trackPCAToJetAxis)")
            .Define("JetsConstituents_trackPCAToJetAxis_z", "IPAlephTools::getPCA_z(trackPCAToJetAxis)")
            .Define("JetsConstituents_jetAxisPCAToTrack_x", "IPAlephTools::getPCA_x(jetAxisPCAToTrack)")
            .Define("JetsConstituents_jetAxisPCAToTrack_y", "IPAlephTools::getPCA_y(jetAxisPCAToTrack)")
            .Define("JetsConstituents_jetAxisPCAToTrack_z", "IPAlephTools::getPCA_z(jetAxisPCAToTrack)")
            .Define("JetsConstituents_linePCAToPrimaryVertex_x", "IPAlephTools::getPCA_x(linePCAToPrimaryVertex)")
            .Define("JetsConstituents_linePCAToPrimaryVertex_y", "IPAlephTools::getPCA_y(linePCAToPrimaryVertex)")
            .Define("JetsConstituents_linePCAToPrimaryVertex_z", "IPAlephTools::getPCA_z(linePCAToPrimaryVertex)")
        )

        return dfout

    def output():

        # define what output to store

        # define output branches
        branchList = []

        # event identifiers
        branchList = [
            'runNumber',
            'eventNumber'
        ]

        # gen-level stuff
        branchList += [
            'genEventType',
            'GenParticle_pdgId',
            'GenParticle_genStatus',
            'GenPV_x',
            'GenPV_y',
            'GenPV_z'
        ]

        # event-level variables
        branchList += [
            'recoEventType',
            'Event_njets',
            'Event_mass',
            'Event_Bz',
            'Event_de',
            'Event_dpt',
            'Event_dphi',
            'Event_dtheta',
            'Event_nRecoParticles',
            'Event_nTracks',
            'Event_nTracksWithReco',
            'Event_nTracksWithoutReco',
            'Event_nSelectedTracks',
            'Event_nPrimaryTracks',
            'Event_nConstituentsSum',
            'Event_nTracksPerJetSum'
        ]

        # primary vertex variables
        branchList += [
            'PV_x',
            'PV_y',
            'PV_z',
            'PV_chi2Normalized',
            'PV_ndof'
        ]
        
        # jet-level variables
        branchList += [
            'Jets_e',
            'Jets_mass',
            'Jets_pt',
            'Jets_px',
            'Jets_py',
            'Jets_pz',
            'Jets_phi',
            'Jets_eta',
            'Jets_theta',
            'Jets_nConstituents',
            'Jets_nTracksPerJet',
            'Jets_nSelectedTracksPerJet',
        ]

        # event-level variables about jet clustering
        branchList += [
            'Event_dmerge1',
            'Event_dmerge2',
            'Event_dmerge3',
            'Event_dmerge4',
        ]

        # secondary track counters
        branchList += [
            'Event_nSecondaryTracks',
            'Jets_nSecondaryTracksPerJet'
        ]

        # secondary vertex variables
        if do_secondary_vertices:
          branchList += [
            'Event_nSV',
            'Jets_nSV',
            'SecondaryVertices_xrel',
            'SecondaryVertices_yrel',
            'SecondaryVertices_zrel',
            'SecondaryVertices_thetarel',
            'SecondaryVertices_phirel',
            'SecondaryVertices_p',
            'SecondaryVertices_prel',
            'SecondaryVertices_chi2',
            'SecondaryVertices_chi2Normalized',
            'SecondaryVertices_ndof',
            'SecondaryVertices_nTracks',
            'SecondaryVertices_mass',
            'SecondaryVertices_dxy',
            'SecondaryVertices_dxyz',
            'SecondaryVertices_cosPointing',
            'SecondaryVertices_correctedMass',
          ]

        # V0-candidate variables
        if do_v0_candidates:
          branchList += [
            'Event_nV0Candidates',
            'Jets_nV0Candidates',
            'Jets_nKsCandidates',
            'Jets_nLambdaCandidates',
            'V0Candidates_pdgId',
            'V0Candidates_xrel',
            'V0Candidates_yrel',
            'V0Candidates_zrel',
            'V0Candidates_thetarel',
            'V0Candidates_phirel',
            'V0Candidates_p',
            'V0Candidates_prel',
            'V0Candidates_chi2',
            'V0Candidates_chi2Normalized',
            'V0Candidates_ndof',
            'V0Candidates_nTracks',
            'V0Candidates_mass',
            'V0Candidates_dxy',
            'V0Candidates_dxyz',
            'V0Candidates_cosPointing',
            'V0Candidates_correctedMass',
          ]

        # jet-constituent-level variables
        branchList += [

            'JetsConstituents_pdgId',

            'JetsConstituents_e', 'JetsConstituents_pt',
            'JetsConstituents_px', 'JetsConstituents_py', 'JetsConstituents_pz',
            'JetsConstituents_theta', 'JetsConstituents_phi',
            'JetsConstituents_charge',
            'JetsConstituents_erel', 'JetsConstituents_erel_log',
            'JetsConstituents_ptrel', 'JetsConstituents_ptrel_log',
            'JetsConstituents_thetarel', 'JetsConstituents_phirel', 
            #'JetsConstituents_dndx',
            #'JetsConstituents_mtof',

            'JetsConstituents_trackChi2',
            'JetsConstituents_trackNdof',
            'JetsConstituents_trackChi2Normalized',
            
            'JetsConstituents_nTrackHits_VDET',
            'JetsConstituents_nTrackHits_ITC',
            'JetsConstituents_nTrackHits_TPC',

            'JetsConstituents_dEdx_pads_type',
            'JetsConstituents_dEdx_pads_value',
            'JetsConstituents_dEdx_pads_error',
            'JetsConstituents_dEdx_wires_type',
            'JetsConstituents_dEdx_wires_value',
            'JetsConstituents_dEdx_wires_error',

            'JetsConstituents_PID_pval_pads_ele',
            'JetsConstituents_PID_pval_pads_mu',
            'JetsConstituents_PID_pval_pads_pi',
            'JetsConstituents_PID_pval_pads_kaon',
            'JetsConstituents_PID_pval_pads_proton',

            'JetsConstituents_PID_pval_wires_ele',
            'JetsConstituents_PID_pval_wires_mu',
            'JetsConstituents_PID_pval_wires_pi',
            'JetsConstituents_PID_pval_wires_kaon',
            'JetsConstituents_PID_pval_wires_proton',

            'JetsConstituents_d0_wrt0',
            'JetsConstituents_z0_wrt0',
            'JetsConstituents_phi0_wrt0',
            'JetsConstituents_omega_wrt0',
            'JetsConstituents_tanlambda_wrt0',
            'JetsConstituents_Bz',
            'JetsConstituents_dxy',
            'JetsConstituents_dz',
            'JetsConstituents_phi0',
            'JetsConstituents_C',
            'JetsConstituents_ct',

            'JetsConstituents_omega_cov',
            'JetsConstituents_d0_cov',
            'JetsConstituents_z0_cov',
            'JetsConstituents_phi0_cov',
            'JetsConstituents_tanlambda_cov',
            'JetsConstituents_d0_z0_cov',
            'JetsConstituents_phi0_d0_cov',
            'JetsConstituents_phi0_z0_cov', 
            'JetsConstituents_tanlambda_phi0_cov',
            'JetsConstituents_tanlambda_d0_cov',
            'JetsConstituents_tanlambda_z0_cov', 
            'JetsConstituents_omega_tanlambda_cov',
            'JetsConstituents_omega_phi0_cov',
            'JetsConstituents_omega_d0_cov',
            'JetsConstituents_omega_z0_cov', 
            
            'JetsConstituents_Sip2dVal',
            'JetsConstituents_Sip2dSig', 
            'JetsConstituents_Sip3dVal',
            'JetsConstituents_Sip3dSig', 
            'JetsConstituents_JetDistVal',
            'JetsConstituents_JetDistSig',

            'JetsConstituents_linearSignedIP3D',
            'JetsConstituents_linearSignedIP3DSig',
            'JetsConstituents_transverseJetDistance',
            'JetsConstituents_longitudinalJetDistance',

            "JetsConstituents_trackPCAToJetAxis_x",
            "JetsConstituents_trackPCAToJetAxis_y",
            "JetsConstituents_trackPCAToJetAxis_z",
            "JetsConstituents_jetAxisPCAToTrack_x",
            "JetsConstituents_jetAxisPCAToTrack_y",
            "JetsConstituents_jetAxisPCAToTrack_z",
            "JetsConstituents_linePCAToPrimaryVertex_x",
            "JetsConstituents_linePCAToPrimaryVertex_y",
            "JetsConstituents_linePCAToPrimaryVertex_z",

            'JetsConstituents_isMu', 
            'JetsConstituents_isEl', 
            'JetsConstituents_isChargedHad',
            'JetsConstituents_isGamma', 
            'JetsConstituents_isNeutralHad',
            'Jets_nMu',
            'Jets_nEl',
            'Jets_nChargedHad',
            'Jets_nPhoton',
            'Jets_nNeutralHad',
        ]

        return branchList    
