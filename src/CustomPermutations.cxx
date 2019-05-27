//
// Created by lennart on 24/04/19.
//
#include <fstream>
#include <iostream>
#include "KLFitter/CustomPermutations.h"

namespace KLFitter {

    CustomPermutations::CustomPermutations(KLFitter::Particles **p, KLFitter::Particles **pp) : fParticles(p)
            , fParticlesPermuted(pp)
            , fPermutationIndex(-1) {

        std::ifstream i("/home/lennart/cernbox/PhD_Sync/Packages/KLFitter/test_perm.json");
        i >> fJsonPermutations;

    }

    std::vector<std::vector<int> > *CustomPermutations::PermutationTable() {
        return nullptr;
    }

    int CustomPermutations::SetParticles(KLFitter::Particles *particles) {

        return 1;
    }

    int CustomPermutations::SetPermutation(int index) {
        // check index
        if (index < 0 || index >= NPermutations()) {
            std::cout << "KLFitter::Permutations::SetPermutation(). Index out of range." << std::endl;
            return 0;
        }

        if (!fParticlesPermuted) {
            std::cout << "KLFitter::Permutations::SetPermutation(). Pointer to permuted particles not available." << std::endl;
            return 0;
        }

        // set permutation
        (*fParticlesPermuted) = &fParticlesTable[index];

        // set permutation index
        fPermutationIndex = index;

        // no error
        return 1;
    }

    int CustomPermutations::CreatePermutations(int nPartonsInPermutations) {
        Reset();
        // create new table of particles
        fParticlesTable = std::vector<KLFitter::Particles>{};

        // create new table of permutations
        fPermutationTable = std::vector<std::vector<int> >{};

        fTablePartons = std::vector<std::vector<int> >{};
        fTableElectrons = std::vector<std::vector<int> >{};
        fTableMuons = std::vector<std::vector<int> >{};
        fTablePhotons = std::vector<std::vector<int> >{};
        fTableTracks = std::vector<std::vector<int> >{};

        // loop over all permutations stored in json
        for (auto& permutation : fJsonPermutations) {
            KLFitter::Particles particles{};

            std::vector<int> bJetParticleIndices = {};
            std::vector<int> lJetParticleIndices = {};

            std::vector<int> bJetPermutedParticleIndices = {};
            std::vector<int> lJetPermutedParticleIndices = {};

            for (int i = 0; i < (*fParticles)->NPartons(); i++) {
                if ((*fParticles)->IsBTagged(i)) bJetParticleIndices.emplace_back(i);
                else lJetParticleIndices.emplace_back(i);
            }

            const auto bJetPermutedIndices = permutation["b_indices"].get<std::vector<int>>();
            const auto lJetPermutedIndices = permutation["light_indices"].get<std::vector<int>>();

            if (bJetParticleIndices.size() !=4) {
                std::cout << "count other than 4 bjets is not yet implemented." << std::endl;
                continue;
            }
//            for_each(bJetPermutedIndices.begin(),bJetPermutedIndices.end(),[](auto x) {std::cout << x;});
//            std::cout << " : ";
//            for_each(lJetPermutedIndices.begin(),lJetPermutedIndices.end(),[](auto x) {std::cout << x;});
//            std::cout << std::endl;
//            std::cout << bJetPermutedIndices << " : " << lJetPerm
            
            for (auto j : bJetPermutedIndices) {
                int permutedB = bJetParticleIndices[j];
                bJetPermutedParticleIndices.emplace_back(permutedB);
                
                particles.AddParticle(
                        (*fParticles)->Parton(permutedB),
                        (*fParticles)->DetEta(permutedB, KLFitter::Particles::kParton),
                        KLFitter::Particles::kParton,
                        (*fParticles)->NameParticle(permutedB, KLFitter::Particles::kParton),
                        (*fParticles)->JetIndex(permutedB),
                        (*fParticles)->IsBTagged(permutedB),
                        (*fParticles)->BTaggingEfficiency(permutedB),
                        (*fParticles)->BTaggingRejection(permutedB),
                        (*fParticles)->TrueFlavor(permutedB),
                        (*fParticles)->BTagWeight(permutedB));
            }
            for (auto j : lJetPermutedIndices) {
                int permutedL = lJetParticleIndices[j];
                lJetPermutedParticleIndices.emplace_back(permutedL);


//                auto names = particles.ParticleNameContainer(KLFitter::Particles::kParton);
//                if ((*fParticles)->NameParticle(permutedL, KLFitter::Particles::kParton) == "Jet_4") {
//                    if (std::any_of(names->begin(), names->end(), [](std::string x) { return x == "Jet_4"; })) {
//                        std::cout << "THIS IS NOT GOOD." << std::endl;
//                    } else {
//                        std::cout << "ok." << std::endl;
//                    }
//                }


                particles.AddParticle(
                        (*fParticles)->Parton(permutedL),
                        (*fParticles)->DetEta(permutedL, KLFitter::Particles::kParton),
                        KLFitter::Particles::kParton,
                        (*fParticles)->NameParticle(permutedL, KLFitter::Particles::kParton),
                        (*fParticles)->JetIndex(permutedL),
                        (*fParticles)->IsBTagged(permutedL),
                        (*fParticles)->BTaggingEfficiency(permutedL),
                        (*fParticles)->BTaggingRejection(permutedL),
                        (*fParticles)->TrueFlavor(permutedL),
                        (*fParticles)->BTagWeight(permutedL));
            }

            for (int i = 0; i < (*fParticles)->NElectrons(); ++i) { //TODO: implement possibility for more leptons.
                if (i > 0 ) break;

                int permutedE = 0;
                particles.AddParticle((*fParticles)->Electron(permutedE),
                                      (*fParticles)->DetEta(permutedE, KLFitter::Particles::kElectron),
                                      KLFitter::Particles::kElectron,
                                      (*fParticles)->NameParticle(permutedE, KLFitter::Particles::kElectron),
                                      (*fParticles)->ElectronIndex(permutedE));
            }
            for (int i = 0; i < (*fParticles)->NMuons(); ++i) {
                if (i > 0)break;
                int permutedM = 0;
                particles.AddParticle((*fParticles)->Muon(permutedM),
                                      (*fParticles)->DetEta(permutedM, KLFitter::Particles::kMuon),
                                      KLFitter::Particles::kMuon,
                                      (*fParticles)->NameParticle(permutedM, KLFitter::Particles::kMuon),
                                      (*fParticles)->MuonIndex(permutedM));
            }

            fParticlesTable.emplace_back(particles);
//            std::cout << permutation << std::endl;
        }



        return 1;
    }

    int KLFitter::CustomPermutations::Reset() {
        // Clear particle and permutation tables.
        fParticlesTable.clear();
        fPermutationTable.clear();

        // no error
        return 1;
    }

    int CustomPermutations::CheckParticles() {
        return 0;
    }

    CustomPermutations::~CustomPermutations() = default;

    CustomPermutations::CustomPermutations(const CustomPermutations &o) = default;


}
