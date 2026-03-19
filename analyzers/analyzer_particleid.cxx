#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/ParticleIDData.h"
#include <vector>

// Extract ParticleID types for jet constituents
// Returns a nested vector structure matching the jet constituents structure
ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> getParticleIDTypes(
    const ROOT::VecOps::RVec<edm4hep::ParticleIDData>& particle_ids,
    const std::vector<std::vector<int>>& constituents) {

    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> types;

    for (const auto& jet_constituents : constituents) {
        ROOT::VecOps::RVec<int> constituent_types;
        for (auto idx : jet_constituents) {
            // Get the type from ParticleID at this index
            // Type values: 0:Track, 1:Electron, 2:Muon, 3:Track from V0,
            //              4:EM, 5:Ecal hadron/residual, 6:Hcal element, 7:Lcal element
            int type = particle_ids[idx].type;
            constituent_types.push_back(type);
        }
        types.push_back(constituent_types);
    }

    return types;
}

// Extract dummy ParticleID types for jet constituents.
// Returns same structure as getParticleIDTypes,
// but using dummy particle IDs (in case the real ones are not available).
ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> getDummyParticleIDTypes(
    const std::vector<std::vector<int>>& constituents) {
    
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> types;

    for (const auto& jet_constituents : constituents) {
        ROOT::VecOps::RVec<int> constituent_types;
        for (auto idx : jet_constituents) {
            int type = 0;
            constituent_types.push_back(type);
        }
        types.push_back(constituent_types);
    }

    return types;
}

// Check if particle is a muon (type 2)
ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> get_isMu_from_type(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>& types) {
    
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> result;
    
    for (const auto& jet_types : types) {
        ROOT::VecOps::RVec<int> jet_flags;
        for (auto type : jet_types) {
            jet_flags.push_back(type == 2 ? 1 : 0);
        }
        result.push_back(jet_flags);
    }
    
    return result;
}

// Check if particle is an electron (type 1)
ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> get_isEl_from_type(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>& types) {
    
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> result;
    
    for (const auto& jet_types : types) {
        ROOT::VecOps::RVec<int> jet_flags;
        for (auto type : jet_types) {
            jet_flags.push_back(type == 1 ? 1 : 0);
        }
        result.push_back(jet_flags);
    }
    
    return result;
}

// Check if particle is a photon/gamma (type 4: EM)
ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> get_isGamma_from_type(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>& types) {
    
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> result;
    
    for (const auto& jet_types : types) {
        ROOT::VecOps::RVec<int> jet_flags;
        for (auto type : jet_types) {
            jet_flags.push_back(type == 4 ? 1 : 0);
        }
        result.push_back(jet_flags);
    }
    
    return result;
}

// Check if particle is a charged hadron (types 0: Track, 3: Track from V0)
ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> get_isChargedHad_from_type(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>& types) {
    
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> result;
    
    for (const auto& jet_types : types) {
        ROOT::VecOps::RVec<int> jet_flags;
        for (auto type : jet_types) {
            jet_flags.push_back((type == 0) ? 1 : 0);
        }
        result.push_back(jet_flags);
    }
    
    return result;
}

// Check if particle is a neutral hadron (types 5: Ecal hadron/residual, 6: Hcal element, 7: Lcal element)
ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> get_isNeutralHad_from_type(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>& types) {
    
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> result;
    
    for (const auto& jet_types : types) {
        ROOT::VecOps::RVec<int> jet_flags;
        for (auto type : jet_types) {
            jet_flags.push_back((type == 5) ? 1 : 0);
        }
        result.push_back(jet_flags);
    }
    
    return result;
}

