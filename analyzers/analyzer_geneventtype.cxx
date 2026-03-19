#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticle.h"
#include <vector>
#include <cmath>


int get_genEventType(const ROOT::VecOps::RVec<edm4hep::MCParticleData>& parts) {
    // note: does not work on our Aleph input files,
    //       since the parents and daughters of MC particles are apparently not stored...

    // initialization
    int primaryBosonIndex = -1;

    // Step 1: find the primary Z or gamma (parent e+e-)
    for (size_t i = 0; i < parts.size(); ++i) {
        const auto& p = parts[i];
        int pdg = std::abs(p.PDG);

        if (pdg == 23 || pdg == 22) {
            // access parents via index range
            int nParents = p.parents_end - p.parents_begin;

            // usually we should have two parents (e+ and e-)
            if (nParents == 2) {
                int parent1 = std::abs(parts[p.parents_begin].PDG);
                int parent2 = std::abs(parts[p.parents_begin + 1].PDG);
                if (parent1 == 11 && parent2 == 11) {
                    primaryBosonIndex = i;
                    break;
                }
            }
            // but sometimes apparently only one is stored (?)
            else if (nParents == 1) {
                if (std::abs(parts[p.parents_begin].PDG) == 11){
                    primaryBosonIndex = i;
                    break;
                }
            }

            // but in our files apparently no parent is stored (?)
            // (not sure if there is any guarantee that this is the primary boson)
            else if (nParents == 0) {
                primaryBosonIndex = i;
                break;
            }
        }
    }

    // return dummy value of no primary boson found
    if( primaryBosonIndex < 0){ return -1; }

    // Step 2: get quark daughters
    const auto& boson = parts[primaryBosonIndex];
    int nDaughters = boson.daughters_end - boson.daughters_begin;

    // return dummy value if no daughters found
    if( nDaughters <= 0 ){ return -2; }
    for (int j = boson.daughters_begin; j < boson.daughters_end; ++j) {
        int qid = std::abs(parts[j].PDG);
        if (qid >= 1 && qid <= 5) {
            return qid;
        }
    }

    // return dummy value if something else went wrong
    return -3;
}
