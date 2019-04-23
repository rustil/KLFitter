// c++ includes
#include <iostream>
#include <string>
#include <vector>

// KLFitter includes
#include "KLFitter/DetectorAtlas_8TeV.h"

#include "KLFitter/Fitter.h"
#include "KLFitter/LikelihoodFourTopLeptonJets.h"
#include "KLFitter/Permutations.h"
#include "KLFitter/PhysicsConstants.h"

// ROOT includes
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

int test() {
    KLFitter::Fitter fitter{};
    KLFitter::DetectorAtlas_8TeV detector{
            "/home/lennart/cernbox/PhD_Sync/Packages/KLFitter/data/transferfunctions/atlasmc12a/akt4_LCtopo_PP6"};
    if (!fitter.SetDetector(&detector)) {
        std::cerr << "ERROR: Failed to set detector! Aborting" << std::endl;
        return EXIT_FAILURE;
    }

    KLFitter::LikelihoodFourTopLeptonJets likelihood{};

    // Set the likelihood properties.
    likelihood.PhysicsConstants()->SetMassTop(172.5);  // mass in GeV
    // Other b-tagging mode options are e.g. 'kNotag' and 'kVeto'.
    // For all possible options, refer to the documentation.
    likelihood.SetBTagging(KLFitter::LikelihoodBase::BtaggingMethod::kWorkingPoint);
    likelihood.SetFlagTopMassFixed(true);

    if (!fitter.SetLikelihood(&likelihood)) {
        std::cerr << "ERROR: Failed to set likelihood. Aborting." << std::endl;
        return 1;
    }

    // auto inFile = TFile::Open("/home/lennart/Data/klftest/SM4t-212560/SSML/mc16e/inclusive/ttbar_nonAllHad.root",
    // "READ");
    auto inFile =
            TFile::Open("/home/lennart/Data/klftest/SM4t-212560/1LOS/mc16d/1los_4top/4topNLO.root", "READ");

    TTreeReader reader("nominal_Loose", inFile);

    // TTreeReaderValue<float> lep_0_pt(reader, "lep_0_pt");
    // TTreeReaderValue<float> lep_0_eta(reader, "lep_0_eta");
    // TTreeReaderValue<float> lep_0_phi(reader, "lep_0_phi");
    // TTreeReaderValue<float> lep_0_e(reader, "lep_0_e");
    // TTreeReaderValue<int> lep_0_pdgId(reader, "lep_0_pdgId");

    TTreeReaderValue<std::vector<float>> jet_pt(reader, "jet_pt");
    TTreeReaderValue<std::vector<float>> jet_eta(reader, "jet_eta");
    TTreeReaderValue<std::vector<float>> jet_phi(reader, "jet_phi");
    TTreeReaderValue<std::vector<float>> jet_e(reader, "jet_e");

    TTreeReaderValue<std::vector<float>> jet_mv2c10(reader, "jet_mv2c10");
//  TTreeReaderValue<std::vector<char>> jet_isbtagged_MV2c10_77(reader, "jet_isbtagged_MV2c10_77");
    TTreeReaderValue<int> nBTags_MV2c10_70(reader, "nBTags_MV2c10_70");

    TTreeReaderValue<std::vector<float>> el_pt(reader, "el_pt");
    TTreeReaderValue<std::vector<float>> el_eta(reader, "el_eta");
    TTreeReaderValue<std::vector<float>> el_phi(reader, "el_phi");
    TTreeReaderValue<std::vector<float>> el_e(reader, "el_e");

    TTreeReaderValue<std::vector<float>> mu_pt(reader, "mu_pt");
    TTreeReaderValue<std::vector<float>> mu_eta(reader, "mu_eta");
    TTreeReaderValue<std::vector<float>> mu_phi(reader, "mu_phi");
    TTreeReaderValue<std::vector<float>> mu_e(reader, "mu_e");

    TTreeReaderValue<float> met_met(reader, "met_met");
    TTreeReaderValue<float> met_phi(reader, "met_phi");

    unsigned int nEventsMax = 100;
    unsigned int eventInd = 0;

    /* Out section */
    TFile outFile("./output.root", "RECREATE");
    TTree outTree("nominal", "nominal");

//  std::vector<float> klf_bhad_pt;
//  std::vector<float> klf_bhad_eta;
//  std::vector<float> klf_bhad_phi;
//  std::vector<float> klf_bhad_e;
//  std::vector<unsigned int> klf_bhad_jet_index;
//  std::vector<float> klf_blep_pt;
//  std::vector<float> klf_blep_eta;
//  std::vector<float> klf_blep_phi;
//  std::vector<float> klf_blep_e;
//  std::vector<unsigned int> klf_blep_jet_index;
//  std::vector<float> klf_lquark1_pt;
//  std::vector<float> klf_lquark1_eta;
//  std::vector<float> klf_lquark1_phi;
//  std::vector<float> klf_lquark1_e;
//  std::vector<unsigned int> klf_lquark1_jet_index;
//  std::vector<float> klf_lquark2_pt;
//  std::vector<float> klf_lquark2_eta;
//  std::vector<float> klf_lquark2_phi;
//  std::vector<float> klf_lquark2_e;
//  std::vector<unsigned int> klf_lquark2_jet_index;
//  std::vector<float> klf_lepton_pt;
//  std::vector<float> klf_lepton_eta;
//  std::vector<float> klf_lepton_phi;
//  std::vector<float> klf_lepton_e;
//  std::vector<float> klf_neutrino_pt;
//  std::vector<float> klf_neutrino_eta;
//  std::vector<float> klf_neutrino_phi;
//  std::vector<float> klf_neutrino_e;
//  std::vector<double> klf_loglikelihood;
//  std::vector<double> klf_event_probability;
//  std::vector<char> klf_fit_minuit_did_not_converge;
//  std::vector<char> klf_fit_aborted_to_nan;
//  std::vector<char> klf_fit_parameter_at_limit;
//  std::vector<char> klf_fit_invalid_transfer_function;
//  unsigned int klf_highest_prob_index;
//
//  // Prepare the output variables as branches of the tree.
//  outTree.Branch("klf_bhad_pt", &klf_bhad_pt);
//  outTree.Branch("klf_bhad_eta", &klf_bhad_eta);
//  outTree.Branch("klf_bhad_phi", &klf_bhad_phi);
//  outTree.Branch("klf_bhad_e", &klf_bhad_e);
//  outTree.Branch("klf_bhad_jet_index", &klf_bhad_jet_index);
//  outTree.Branch("klf_blep_pt", &klf_blep_pt);
//  outTree.Branch("klf_blep_eta", &klf_blep_eta);
//  outTree.Branch("klf_blep_phi", &klf_blep_phi);
//  outTree.Branch("klf_blep_e", &klf_blep_e);
//  outTree.Branch("klf_blep_jet_index", &klf_blep_jet_index);
//  outTree.Branch("klf_lquark1_pt", &klf_lquark1_pt);
//  outTree.Branch("klf_lquark1_eta", &klf_lquark1_eta);
//  outTree.Branch("klf_lquark1_phi", &klf_lquark1_phi);
//  outTree.Branch("klf_lquark1_e", &klf_lquark1_e);
//  outTree.Branch("klf_lquark1_jet_index", &klf_lquark1_jet_index);
//  outTree.Branch("klf_lquark2_pt", &klf_lquark2_pt);
//  outTree.Branch("klf_lquark2_eta", &klf_lquark2_eta);
//  outTree.Branch("klf_lquark2_phi", &klf_lquark2_phi);
//  outTree.Branch("klf_lquark2_e", &klf_lquark2_e);
//  outTree.Branch("klf_lquark2_jet_index", &klf_lquark2_jet_index);
//  outTree.Branch("klf_lepton_pt", &klf_lepton_pt);
//  outTree.Branch("klf_lepton_eta", &klf_lepton_eta);
//  outTree.Branch("klf_lepton_phi", &klf_lepton_phi);
//  outTree.Branch("klf_lepton_e", &klf_lepton_e);
//  outTree.Branch("klf_neutrino_pt", &klf_neutrino_pt);
//  outTree.Branch("klf_neutrino_eta", &klf_neutrino_eta);
//  outTree.Branch("klf_neutrino_phi", &klf_neutrino_phi);
//  outTree.Branch("klf_neutrino_e", &klf_neutrino_e);
//  outTree.Branch("klf_loglikelihood", &klf_loglikelihood);
//  outTree.Branch("klf_event_probability", &klf_event_probability);
//  outTree.Branch("klf_fit_minuit_did_not_converge", &klf_fit_minuit_did_not_converge);
//  outTree.Branch("klf_fit_aborted_to_nan", &klf_fit_aborted_to_nan);
//  outTree.Branch("klf_fit_parameter_at_limit", &klf_fit_parameter_at_limit);
//  outTree.Branch("klf_fit_invalid_transfer_function", &klf_fit_invalid_transfer_function);
//
//  outTree.Branch("klf_highest_prob_index", &klf_highest_prob_index);



    while (reader.Next()) {
        unsigned int nLep = mu_pt->size() + el_pt->size();
        if (eventInd >= nEventsMax) break;
        else if (nLep == 0 || nLep > 1 || jet_pt->size() != 10 || *nBTags_MV2c10_70 != 4) continue;

        if ((eventInd % 1) == 0) printf("[%10i | %10i]\n", eventInd, nEventsMax);

        KLFitter::Particles particles{};

        TLorentzVector lep;
        unsigned int pdgId;

        if (mu_pt->empty() && !el_pt->empty()) {
            pdgId = 11;
        } else if (!mu_pt->empty() && el_pt->empty()) {
            pdgId = 13;
        } else {
            pdgId = mu_pt->at(0) > el_pt->at(0) ? 13 : 11;
        }

        switch (pdgId) {
            case 11:
                lep.SetPtEtaPhiE(el_pt->at(0) / 1e3, el_eta->at(0), el_phi->at(0), el_e->at(0) / 1e3);
                if (abs(lep.Eta()) > 1.37 && abs(lep.Eta()) < 1.52) continue;  // crack region
                likelihood.SetLeptonType(KLFitter::LikelihoodFourTopLeptonJets::kElectron);
                particles.AddParticle(lep, lep.Eta(), KLFitter::Particles::kElectron);
                break;
            case 13:
                lep.SetPtEtaPhiE(mu_pt->at(0) / 1e3, mu_eta->at(0), mu_phi->at(0), mu_e->at(0) / 1e3);
                likelihood.SetLeptonType(KLFitter::LikelihoodFourTopLeptonJets::kMuon);
                particles.AddParticle(lep, lep.Eta(), KLFitter::Particles::kMuon);
                break;

            default:
                continue;
        }

        for (unsigned int ijet = 0; ijet < 10; ijet++) {
            TLorentzVector jet;
            jet.SetPtEtaPhiE(jet_pt->at(ijet)/1e3, jet_eta->at(ijet), jet_phi->at(ijet), jet_e->at(ijet)/1e3);
            // Arguments are as follows:
            //  1) TLorentzVector of jet
            //  2) jet eta
            //  3) KLFitter particle type. kParton for jets
            //  4) Internal name
            //  5) Index of the jet
            //  6) Is the jet btagged?
            //  7) tagging effciency required for kWorkingPoint
            //  8) 1./tagging inefficiency required for kWorkingPoint
            //  9) true flavour type
            //  10) btag discriminant
            particles.AddParticle(jet, jet_eta->at(ijet), KLFitter::Particles::kParton, "", ijet,
                                  static_cast<int>(jet_mv2c10->at(ijet) > 0.64), 0.6, 145., KLFitter::Particles::kNone,
                                  jet_mv2c10->at(ijet));
        }

        // Add particles to the likelihood.
        if (!fitter.SetParticles(&particles)) {
            std::cerr << "ERROR: Failed to add particles to KLFitter. Aborting." << std::endl;
            return 1;
        }


        eventInd++;
    }

    // Go to the output file and write it.
    outFile.cd();
    std::cout << std::endl << "Writing into output root file: " << "output.root" << std::endl << std::endl;
    outTree.Write();

    // Close both input and output ROOT files.
    outFile.Close();
    inFile->Close();

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    return test();
    // return EXIT::SUCCES;
}
