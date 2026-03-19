# Exploratory script to investigate the primary vertex retrieval.

# The standard FCCAnalyses function is here:
# https://github.com/LukaLambrecht/FCCAnalyses/blob/68169d1216b43b3eacf1eb4e8a8298fa9f275417/analyzers/dataframe/src/MCParticle.cc#L102
# but it does not work for Aleph MC (no candidate is found and a dummy is returned).
# This scripts investigates how to fix it.


import ROOT


# helper function for storing the gen status of all electrons
# (hoping one will be the incoming particle)
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<int> getElectronGeneratorStatus(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles) {
        ROOT::VecOps::RVec<int> res;
        for (const auto& genParticle : genParticles) {
            if( genParticle.PDG==11 ){ res.push_back(genParticle.generatorStatus); }
        }
        return res;
    }""")

# helper function for storing the gen status at 10000 + a little bit
# (should correspond to hardest subprocess)
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<int> getSelectedGeneratorStatus(
        const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles) {
        ROOT::VecOps::RVec<int> res;
        for (const auto& genParticle : genParticles) {
            if( genParticle.generatorStatus>=10000 and genParticle.generatorStatus<=10100 ){ res.push_back(genParticle.generatorStatus); }
        }
        return res;
    }""")


# main analyzer class
class RDFanalysis():

    def analysers(df_in):
       
        df_out = (
            df_in
            .Define("electronGeneratorStatus", "getElectronGeneratorStatus(MCParticles)")
            .Define("selectedGeneratorStatus", "getSelectedGeneratorStatus(MCParticles)")
        )
        return df_out

    def output():
        branchList = ([
          "electronGeneratorStatus",
          "selectedGeneratorStatus"
        ])
        return branchList
