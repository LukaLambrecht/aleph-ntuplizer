#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticle.h"
#include <vector>
#include <cmath>


#ifndef GenEventType_H
#define GenEventType_H

namespace GenEventType{

int get_genEventType(const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles) {
    // note: does not work on our Aleph input files,
    //       since the parents and daughters of MC particles are apparently not stored...
    //       instead, a simpler version is used; see below.

    // initialization
    int primaryBosonIndex = -1;

    // Step 1: find the primary Z or gamma (parent e+e-)
    for (size_t i = 0; i < genParticles.size(); ++i) {
        const auto& p = genParticles[i];
        int pdg = std::abs(p.PDG);

        if (pdg == 23 || pdg == 22) {
            // access parents via index range
            int nParents = p.parents_end - p.parents_begin;

            // usually we should have two parents (e+ and e-)
            if (nParents == 2) {
                int parent1 = std::abs(genParticles[p.parents_begin].PDG);
                int parent2 = std::abs(genParticles[p.parents_begin + 1].PDG);
                if (parent1 == 11 && parent2 == 11) {
                    primaryBosonIndex = i;
                    break;
                }
            }
            // but sometimes apparently only one is stored (?)
            else if (nParents == 1) {
                if (std::abs(genParticles[p.parents_begin].PDG) == 11){
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
    const auto& boson = genParticles[primaryBosonIndex];
    int nDaughters = boson.daughters_end - boson.daughters_begin;

    // return dummy value if no daughters found
    if( nDaughters <= 0 ){ return -2; }
    for (int j = boson.daughters_begin; j < boson.daughters_end; ++j) {
        int qid = std::abs(genParticles[j].PDG);
        if (qid >= 1 && qid <= 5) {
            return qid;
        }
    }

    // return dummy value if something else went wrong
    return -3;
}


int get_genEventTypeFromFirstQuark(const ROOT::VecOps::RVec<edm4hep::MCParticleData>& genParticles) {
    // helper function for deriving the gen-level event type.
    // note: for now, only valid with qqbar simulations,
    //       where the event type is between 1 (d dbar) and 5 (b bbar) (see PDG numbering scheme).
    // note: the event type is derived simply from the first quark in the list of MCParticles;
    //       there is in principle no guarantee for any kind of ordering;
    //       we just assume the first quark PDG ID in the MCParticle collection
    //       is the one corresponding to the type of quarks produced in the hard scattering.
    //       to be checked and refined later.
    // note: in the original analyzer that served as source for this one,
    //       the event type needed not to be derived, as the simulation was split per quark flavour,
    //       so instead the event type was just derived from the file name.
    for (const auto& genParticle : genParticles) {
        int pdgid = std::abs(genParticle.PDG);
        if( (pdgid >= 1) && (pdgid <= 6) ){ return pdgid; }
    }
    return -1;
}

}

#endif
