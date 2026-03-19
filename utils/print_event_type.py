# Exploratory script to investigate the event type of qq MC samples.

# The provided MC samples are generic qqbar events,
# but we need to find the flavour of the quarks
# in order to define labels for the jets.
# This script explores how to do that.

import ROOT


# helper function for printing some basic properties
# of all gen-level particles in the event
ROOT.gInterpreter.Declare("""
    int printGenParticles(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles) {
        std::cout << "=== begin event ===" << std::endl;
        for (const auto& genParticle : genParticles) {
            std::cout << genParticle.PDG << std::endl;
        }
        std::cout << "=== end event ===" << std::endl;
        std::cout.flush();
        return 0;
    }""")

# helper function for storing the PDG IDs
# of all gen-level particles in the event
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<int> getPdgIds(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles) {
        ROOT::VecOps::RVec<int> res;
        for (const auto& genParticle : genParticles) {
            res.push_back(genParticle.PDG);
        }
        return res;
    }""")

# helper function for storing the PDG ID
# of the first gen-level particle in the event
ROOT.gInterpreter.Declare("""
    int getFirstPdgId(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles) {
        for (const auto& genParticle : genParticles) {
            return genParticle.PDG;
        }
        return -99;
    }""")

# helper function for storing the first quark PDG ID.
# note: there is no formal guarantee for any kind of ordering;
#       we just assume the first quark PDG ID in the MCParticle collection
#       is the one corresponding to the type of quarks produced in the hard scattering.
#       to be refined later.
ROOT.gInterpreter.Declare("""
    int getFirstQuarkPdgId(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles) {
        for (const auto& genParticle : genParticles) {
            int pdgid = std::abs(genParticle.PDG);
            if( (pdgid >= 1) && (pdgid <= 6) ){ return pdgid; }
        }
        return -1;
    }""")

# helper function for storing the first quark generator status
# (just as a first rough check for the assumption outlined above)
ROOT.gInterpreter.Declare("""
    int getFirstQuarkStatus(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles) {
        for (const auto& genParticle : genParticles) {
            int pdgid = std::abs(genParticle.PDG);
            int status = genParticle.generatorStatus;
            if( (pdgid >= 1) && (pdgid <= 6) ){ return status; }
        }
        return -1;
    }""")


# main analyzer class
class RDFanalysis():

    def analysers(df_in):
       
        df_out = (
            df_in
            #.Define("dummy_print", "printGenParticles(MCParticles)")
            #.Define("genPdgIds", "getPdgIds(MCParticles)")
            .Define("firstPdgId", "getFirstPdgId(MCParticles)")
            #.Define("firstQuarkPdgId", "getFirstQuarkPdgId(MCParticles)")
            #.Define("firstQuarkStatus", "getFirstQuarkStatus(MCParticles)")
        )
        return df_out

    def output():
        branchList = ([
          #"genPdgIds",
          "firstPdgId",
          #"firstQuarkPdgId",
          #"firstQuarkStatus"
        ])
        return branchList
