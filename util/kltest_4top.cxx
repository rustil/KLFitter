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
#include <TString.h>
#include <TChain.h>

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
    likelihood.SetBTagging(KLFitter::LikelihoodBase::BtaggingMethod::kVetoNoFit);
    likelihood.SetFlagTopMassFixed(true);

    likelihood.setFPermMethod(KLFitter::LikelihoodBase::kCustomPermutations);

    if (!fitter.SetLikelihood(&likelihood)) {
        std::cerr << "ERROR: Failed to set likelihood. Aborting." << std::endl;
        return 1;
    }

    auto chain = new TChain("nominal_Loose");
    std::vector<std::string> fileNames = {
            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000002.output.root",
            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000001.output.root",
            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000004.output.root",
            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000006.output.root",
            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000005.output.root",
            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000003.output.root"
    };

    for (auto fn : fileNames) chain->AddFile(fn.c_str());

    // auto inFile = TFile::Open("/home/lennart/Data/klftest/SM4t-212560/SSML/mc16e/inclusive/ttbar_nonAllHad.root",
    // "READ");
//    auto inFile =
//            TFile::Open("/home/lennart/Data/klftest/SM4t-212560/1LOS/mc16d/1los_4top/4topNLO.root", "READ");


//auto inFile =
//            TFile::Open("/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000001.output.root", "READ");


    TTreeReader reader(chain);
//    TTreeReader reader("nominal_Loose", inFile);

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

    TTreeReaderValue<std::vector<int>> jet_parentghost_pdgId(reader, "jet_parentghost_pdgId");
    TTreeReaderValue<std::vector<char>> jet_isbtagged_MV2c10_77(reader, "jet_isbtagged_MV2c10_77");


    //TODO:!
    unsigned int nEventsMax = 100;
    unsigned int eventInd = 0;

    /* Out section */
    TFile outFile("./output1.root", "RECREATE");
//    TTree* outTree = (TTree*)inFile->Get("nominal_Loose")->Clone();
    TTree outTree("nominal", "nominal");

    std::vector<float> cp_el_pt;
    std::vector<float> cp_el_eta;
    std::vector<float> cp_el_phi;
    std::vector<float> cp_el_e;

    std::vector<float> cp_mu_pt;
    std::vector<float> cp_mu_eta;
    std::vector<float> cp_mu_phi;
    std::vector<float> cp_mu_e;

    std::vector<float> cp_jet_pt;
    std::vector<float> cp_jet_eta;
    std::vector<float> cp_jet_phi;
    std::vector<float> cp_jet_e;
    std::vector<float> cp_jet_mv2c10;

    float cp_met_met;
    float cp_met_phi;

    std::vector<float> klf_bhad1_pt;
    std::vector<float> klf_bhad1_eta;
    std::vector<float> klf_bhad1_phi;
    std::vector<float> klf_bhad1_e;
    std::vector<unsigned int> klf_bhad1_jet_index;

    std::vector<float> klf_bhad2_pt;
    std::vector<float> klf_bhad2_eta;
    std::vector<float> klf_bhad2_phi;
    std::vector<float> klf_bhad2_e;
    std::vector<unsigned int> klf_bhad2_jet_index;

    std::vector<float> klf_bhad3_pt;
    std::vector<float> klf_bhad3_eta;
    std::vector<float> klf_bhad3_phi;
    std::vector<float> klf_bhad3_e;
    std::vector<unsigned int> klf_bhad3_jet_index;
    
    std::vector<float> klf_blep_pt;
    std::vector<float> klf_blep_eta;
    std::vector<float> klf_blep_phi;
    std::vector<float> klf_blep_e;
    std::vector<unsigned int> klf_blep_jet_index;
    
    std::vector<float> klf_lquark1_pt;
    std::vector<float> klf_lquark1_eta;
    std::vector<float> klf_lquark1_phi;
    std::vector<float> klf_lquark1_e;
    std::vector<unsigned int> klf_lquark1_jet_index;
    
    std::vector<float> klf_lquark2_pt;
    std::vector<float> klf_lquark2_eta;
    std::vector<float> klf_lquark2_phi;
    std::vector<float> klf_lquark2_e;
    std::vector<unsigned int> klf_lquark2_jet_index;

    std::vector<float> klf_lquark3_pt;
    std::vector<float> klf_lquark3_eta;
    std::vector<float> klf_lquark3_phi;
    std::vector<float> klf_lquark3_e;
    std::vector<unsigned int> klf_lquark3_jet_index;

    std::vector<float> klf_lquark4_pt;
    std::vector<float> klf_lquark4_eta;
    std::vector<float> klf_lquark4_phi;
    std::vector<float> klf_lquark4_e;
    std::vector<unsigned int> klf_lquark4_jet_index;

    std::vector<float> klf_lquark5_pt;
    std::vector<float> klf_lquark5_eta;
    std::vector<float> klf_lquark5_phi;
    std::vector<float> klf_lquark5_e;
    std::vector<unsigned int> klf_lquark5_jet_index;

    std::vector<float> klf_lquark6_pt;
    std::vector<float> klf_lquark6_eta;
    std::vector<float> klf_lquark6_phi;
    std::vector<float> klf_lquark6_e;
    std::vector<unsigned int> klf_lquark6_jet_index;
    
    std::vector<float> klf_lepton_pt;
    std::vector<float> klf_lepton_eta;
    std::vector<float> klf_lepton_phi;
    std::vector<float> klf_lepton_e;
    std::vector<float> klf_neutrino_pt;
    std::vector<float> klf_neutrino_eta;
    std::vector<float> klf_neutrino_phi;
    std::vector<float> klf_neutrino_e;
    std::vector<double> klf_loglikelihood;
    std::vector<double> klf_event_probability;
    std::vector<char> klf_fit_minuit_did_not_converge;
    std::vector<char> klf_fit_aborted_to_nan;
    std::vector<char> klf_fit_parameter_at_limit;
    std::vector<char> klf_fit_invalid_transfer_function;
    unsigned int klf_highest_prob_index;

    // Prepare the output variables as branches of the tree.
    outTree.Branch("klf_bhad1_pt", &klf_bhad1_pt);
    outTree.Branch("klf_bhad1_eta", &klf_bhad1_eta);
    outTree.Branch("klf_bhad1_phi", &klf_bhad1_phi);
    outTree.Branch("klf_bhad1_e", &klf_bhad1_e);
    outTree.Branch("klf_bhad1_jet_index", &klf_bhad1_jet_index);

    outTree.Branch("klf_bhad2_pt", &klf_bhad2_pt);
    outTree.Branch("klf_bhad2_eta", &klf_bhad2_eta);
    outTree.Branch("klf_bhad2_phi", &klf_bhad2_phi);
    outTree.Branch("klf_bhad2_e", &klf_bhad2_e);
    outTree.Branch("klf_bhad2_jet_index", &klf_bhad2_jet_index);

    outTree.Branch("klf_bhad3_pt", &klf_bhad3_pt);
    outTree.Branch("klf_bhad3_eta", &klf_bhad3_eta);
    outTree.Branch("klf_bhad3_phi", &klf_bhad3_phi);
    outTree.Branch("klf_bhad3_e", &klf_bhad3_e);
    outTree.Branch("klf_bhad3_jet_index", &klf_bhad3_jet_index);
    
    
    outTree.Branch("klf_blep_pt", &klf_blep_pt);
    outTree.Branch("klf_blep_eta", &klf_blep_eta);
    outTree.Branch("klf_blep_phi", &klf_blep_phi);
    outTree.Branch("klf_blep_e", &klf_blep_e);
    outTree.Branch("klf_blep_jet_index", &klf_blep_jet_index);
    
    outTree.Branch("klf_lquark1_pt", &klf_lquark1_pt);
    outTree.Branch("klf_lquark1_eta", &klf_lquark1_eta);
    outTree.Branch("klf_lquark1_phi", &klf_lquark1_phi);
    outTree.Branch("klf_lquark1_e", &klf_lquark1_e);
    outTree.Branch("klf_lquark1_jet_index", &klf_lquark1_jet_index);
    outTree.Branch("klf_lquark2_pt", &klf_lquark2_pt);
    outTree.Branch("klf_lquark2_eta", &klf_lquark2_eta);
    outTree.Branch("klf_lquark2_phi", &klf_lquark2_phi);
    outTree.Branch("klf_lquark2_e", &klf_lquark2_e);
    outTree.Branch("klf_lquark2_jet_index", &klf_lquark2_jet_index);

    outTree.Branch("klf_lquark3_pt", &klf_lquark3_pt);
    outTree.Branch("klf_lquark3_eta", &klf_lquark3_eta);
    outTree.Branch("klf_lquark3_phi", &klf_lquark3_phi);
    outTree.Branch("klf_lquark3_e", &klf_lquark3_e);
    outTree.Branch("klf_lquark3_jet_index", &klf_lquark3_jet_index);

    outTree.Branch("klf_lquark4_pt", &klf_lquark4_pt);
    outTree.Branch("klf_lquark4_eta", &klf_lquark4_eta);
    outTree.Branch("klf_lquark4_phi", &klf_lquark4_phi);
    outTree.Branch("klf_lquark4_e", &klf_lquark4_e);
    outTree.Branch("klf_lquark4_jet_index", &klf_lquark4_jet_index);

    outTree.Branch("klf_lquark5_pt", &klf_lquark5_pt);
    outTree.Branch("klf_lquark5_eta", &klf_lquark5_eta);
    outTree.Branch("klf_lquark5_phi", &klf_lquark5_phi);
    outTree.Branch("klf_lquark5_e", &klf_lquark5_e);
    outTree.Branch("klf_lquark5_jet_index", &klf_lquark5_jet_index);

    outTree.Branch("klf_lquark6_pt", &klf_lquark6_pt);
    outTree.Branch("klf_lquark6_eta", &klf_lquark6_eta);
    outTree.Branch("klf_lquark6_phi", &klf_lquark6_phi);
    outTree.Branch("klf_lquark6_e", &klf_lquark6_e);
    outTree.Branch("klf_lquark6_jet_index", &klf_lquark6_jet_index);
    
    outTree.Branch("klf_lepton_pt", &klf_lepton_pt);
    outTree.Branch("klf_lepton_eta", &klf_lepton_eta);
    outTree.Branch("klf_lepton_phi", &klf_lepton_phi);
    outTree.Branch("klf_lepton_e", &klf_lepton_e);
    outTree.Branch("klf_neutrino_pt", &klf_neutrino_pt);
    outTree.Branch("klf_neutrino_eta", &klf_neutrino_eta);
    outTree.Branch("klf_neutrino_phi", &klf_neutrino_phi);
    outTree.Branch("klf_neutrino_e", &klf_neutrino_e);
    outTree.Branch("klf_loglikelihood", &klf_loglikelihood);
    outTree.Branch("klf_event_probability", &klf_event_probability);
    outTree.Branch("klf_fit_minuit_did_not_converge", &klf_fit_minuit_did_not_converge);
    outTree.Branch("klf_fit_aborted_to_nan", &klf_fit_aborted_to_nan);
    outTree.Branch("klf_fit_parameter_at_limit", &klf_fit_parameter_at_limit);
    outTree.Branch("klf_fit_invalid_transfer_function", &klf_fit_invalid_transfer_function);

    outTree.Branch("klf_highest_prob_index", &klf_highest_prob_index);

    outTree.Branch("el_pt", &cp_el_pt);
    outTree.Branch("el_eta", &cp_el_eta);
    outTree.Branch("el_phi", &cp_el_phi);
    outTree.Branch("el_e", &cp_el_e);

    outTree.Branch("mu_pt", &cp_mu_pt);
    outTree.Branch("mu_eta", &cp_mu_eta);
    outTree.Branch("mu_phi", &cp_mu_phi);
    outTree.Branch("mu_e", &cp_mu_e);

    outTree.Branch("jet_pt", &cp_jet_pt);
    outTree.Branch("jet_eta", &cp_jet_eta);
    outTree.Branch("jet_phi", &cp_jet_phi);
    outTree.Branch("jet_e", &cp_jet_e);
    outTree.Branch("jet_mv2c10", &cp_jet_mv2c10);

    outTree.Branch("met_met", &cp_met_met);
    outTree.Branch("met_phi", &cp_met_phi);

    
    while (reader.Next()) {
        KLFitter::Particles particles{};

        unsigned int nLep = mu_pt->size() + el_pt->size();
        if (eventInd >= nEventsMax) break;
        else if (nLep != 1 || jet_pt->size() != 10 || std::count_if(jet_isbtagged_MV2c10_77->begin(), jet_isbtagged_MV2c10_77->end(), [](char flag) {return flag;}) != 4) continue;

        int partonsMatcheable(0);
        for (auto pdgId : *jet_parentghost_pdgId) {
            if (abs(pdgId) == 6 || abs(pdgId) == 24) ++partonsMatcheable;
        }

        if (partonsMatcheable != 10) continue;

        if ((eventInd % 1) == 0) printf("[%10i | %10i]\n", eventInd, nEventsMax);

//        else if (nLep == 0 || nLep > 1 || jet_pt->size() != 10 || std::count_if(jet_mv2c10->begin(), jet_mv2c10->end(), [](float x){return x > 0.64;}) != 4) continue;



        cp_el_pt = *el_pt;
        cp_el_eta = *el_eta;
        cp_el_phi = *el_phi;
        cp_el_e = *el_e;

        cp_mu_pt = *mu_pt;
        cp_mu_eta = *mu_eta;
        cp_mu_phi = *mu_phi;
        cp_mu_e = *mu_e;

        cp_jet_pt = *jet_pt;
        cp_jet_eta = *jet_eta;
        cp_jet_phi = *jet_phi;
        cp_jet_e = *jet_e;
        cp_jet_mv2c10 = *jet_mv2c10;

        cp_met_met = *met_met;
        cp_met_phi = *met_phi;




        TLorentzVector lep;
        unsigned int pdgId;

        if (mu_pt->empty() && !el_pt->empty()) {
            pdgId = 11;
        } else if (!mu_pt->empty() && el_pt->empty()) {
            pdgId = 13;
        } else {
//            pdgId = mu_pt->at(0) > el_pt->at(0) ? 13 : 11;
            continue;
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
            particles.AddParticle(jet, jet_eta->at(ijet), KLFitter::Particles::kParton, Form("Jet_%i", ijet), ijet,
                                  jet_isbtagged_MV2c10_77->at(ijet), 0.6, 145., KLFitter::Particles::kNone,
                                  jet_mv2c10->at(ijet));
        }

        // Add particles to the likelihood.
        if (!fitter.SetParticles(&particles)) {
            std::cerr << "ERROR: Failed to add particles to KLFitter. Aborting." << std::endl;
            return 1;
        }

        // Add MET information
        const double met_ex = *met_met * cos(*met_phi)/1e3;
        const double met_ey = *met_met * sin(*met_phi)/1e3;
        if (!fitter.SetET_miss_XY_SumET(met_ex, met_ey, *met_met)) {
            std::cerr << "ERROR: Failed to add MET to fitter. Aborting." << std::endl;
            return 1;
        }

//        // Loop over all permutations.
        const int nperm = fitter.CustomPermutations()->NPermutations();
//
        bool isFirst = true;
//
        for (int iperm = 0; iperm < nperm; iperm++) {

            klf_bhad1_pt.clear();
            klf_bhad1_eta.clear();
            klf_bhad1_phi.clear();
            klf_bhad1_e.clear();
            klf_bhad1_jet_index.clear();

            klf_bhad2_pt.clear();
            klf_bhad2_eta.clear();
            klf_bhad2_phi.clear();
            klf_bhad2_e.clear();
            klf_bhad2_jet_index.clear();

            klf_bhad3_pt.clear();
            klf_bhad3_eta.clear();
            klf_bhad3_phi.clear();
            klf_bhad3_e.clear();
            klf_bhad3_jet_index.clear();
            
            klf_blep_pt.clear();
            klf_blep_eta.clear();
            klf_blep_phi.clear();
            klf_blep_e.clear();
            klf_blep_jet_index.clear();
            klf_lquark1_pt.clear();
            klf_lquark1_eta.clear();
            klf_lquark1_phi.clear();
            klf_lquark1_e.clear();
            klf_lquark1_jet_index.clear();
            
            klf_lquark2_pt.clear();
            klf_lquark2_eta.clear();
            klf_lquark2_phi.clear();
            klf_lquark2_e.clear();
            klf_lquark2_jet_index.clear();

            klf_lquark3_pt.clear();
            klf_lquark3_eta.clear();
            klf_lquark3_phi.clear();
            klf_lquark3_e.clear();
            klf_lquark3_jet_index.clear();

            klf_lquark4_pt.clear();
            klf_lquark4_eta.clear();
            klf_lquark4_phi.clear();
            klf_lquark4_e.clear();
            klf_lquark4_jet_index.clear();

            klf_lquark5_pt.clear();
            klf_lquark5_eta.clear();
            klf_lquark5_phi.clear();
            klf_lquark5_e.clear();
            klf_lquark5_jet_index.clear();

            klf_lquark6_pt.clear();
            klf_lquark6_eta.clear();
            klf_lquark6_phi.clear();
            klf_lquark6_e.clear();
            klf_lquark6_jet_index.clear();
            
            klf_lepton_pt.clear();
            klf_lepton_eta.clear();
            klf_lepton_phi.clear();
            klf_lepton_e.clear();
            klf_neutrino_pt.clear();
            klf_neutrino_eta.clear();
            klf_neutrino_phi.clear();
            klf_neutrino_e.clear();
            klf_loglikelihood.clear();
            klf_event_probability.clear();
            klf_fit_minuit_did_not_converge.clear();
            klf_fit_aborted_to_nan.clear();
            klf_fit_parameter_at_limit.clear();
            klf_fit_invalid_transfer_function.clear();


            fitter.Fit(iperm);

            // Read the output and convergence status of the fit.
            unsigned int ConvergenceStatusBitWord = fitter.ConvergenceStatus();
            bool MinuitDidNotConverge = (ConvergenceStatusBitWord & fitter.MinuitDidNotConvergeMask) != 0;
            bool FitAbortedDueToNaN = (ConvergenceStatusBitWord & fitter.FitAbortedDueToNaNMask) != 0;
            bool AtLeastOneFitParameterAtItsLimit =
                    (ConvergenceStatusBitWord & fitter.AtLeastOneFitParameterAtItsLimitMask) != 0;
            bool InvalidTransferFunctionAtConvergence =
                    (ConvergenceStatusBitWord & fitter.InvalidTransferFunctionAtConvergenceMask) != 0;

            // Get log likelihood and event probability values. Note
            // that the event probablity is _not_ normalized.
            double likelihood = fitter.Likelihood()->LogLikelihood(fitter.Likelihood()->GetBestFitParameters());
            double event_probability = std::exp(fitter.Likelihood()->LogEventProbability());

            // Get the values of all fitted variables.
            auto modelParticles = fitter.Likelihood()->ParticlesModel();
            auto permutedParticles = fitter.Likelihood()->PParticlesPermuted();

            // Hadronic b quark.
            float bhad1_pt = modelParticles->Parton(0)->Pt();
            float bhad1_eta = modelParticles->Parton(0)->Eta();
            float bhad1_phi = modelParticles->Parton(0)->Phi();
            float bhad1_e = modelParticles->Parton(0)->E();
            unsigned int bhad1_index = (*permutedParticles)->JetIndex(0);

            float bhad2_pt = modelParticles->Parton(1)->Pt();
            float bhad2_eta = modelParticles->Parton(1)->Eta();
            float bhad2_phi = modelParticles->Parton(1)->Phi();
            float bhad2_e = modelParticles->Parton(1)->E();
            unsigned int bhad2_index = (*permutedParticles)->JetIndex(1);

            float bhad3_pt = modelParticles->Parton(2)->Pt();
            float bhad3_eta = modelParticles->Parton(2)->Eta();
            float bhad3_phi = modelParticles->Parton(2)->Phi();
            float bhad3_e = modelParticles->Parton(2)->E();
            unsigned int bhad3_index = (*permutedParticles)->JetIndex(2);

            // Leptonic b quark.
            float blep_pt = modelParticles->Parton(3)->Pt();
            float blep_eta = modelParticles->Parton(3)->Eta();
            float blep_phi = modelParticles->Parton(3)->Phi();
            float blep_e = modelParticles->Parton(3)->E();
            unsigned int blep_index = (*permutedParticles)->JetIndex(3);

            // Light quark 1.
            float lquark1_pt = modelParticles->Parton(4)->Pt();
            float lquark1_eta = modelParticles->Parton(4)->Eta();
            float lquark1_phi = modelParticles->Parton(4)->Phi();
            float lquark1_e = modelParticles->Parton(4)->E();
            unsigned int lquark1_index = (*permutedParticles)->JetIndex(4);

            // Light quark 2.
            float lquark2_pt = modelParticles->Parton(5)->Pt();
            float lquark2_eta = modelParticles->Parton(5)->Eta();
            float lquark2_phi = modelParticles->Parton(5)->Phi();
            float lquark2_e = modelParticles->Parton(5)->E();
            unsigned int lquark2_index = (*permutedParticles)->JetIndex(5);

            // Light quark 3.
            float lquark3_pt = modelParticles->Parton(6)->Pt();
            float lquark3_eta = modelParticles->Parton(6)->Eta();
            float lquark3_phi = modelParticles->Parton(6)->Phi();
            float lquark3_e = modelParticles->Parton(6)->E();
            unsigned int lquark3_index = (*permutedParticles)->JetIndex(6);

            // Light quark 4.
            float lquark4_pt = modelParticles->Parton(7)->Pt();
            float lquark4_eta = modelParticles->Parton(7)->Eta();
            float lquark4_phi = modelParticles->Parton(7)->Phi();
            float lquark4_e = modelParticles->Parton(7)->E();
            unsigned int lquark4_index = (*permutedParticles)->JetIndex(7);

            // Light quark 5.
            float lquark5_pt = modelParticles->Parton(8)->Pt();
            float lquark5_eta = modelParticles->Parton(8)->Eta();
            float lquark5_phi = modelParticles->Parton(8)->Phi();
            float lquark5_e = modelParticles->Parton(8)->E();
            unsigned int lquark5_index = (*permutedParticles)->JetIndex(8);

            // Light quark 6.
            float lquark6_pt = modelParticles->Parton(9)->Pt();
            float lquark6_eta = modelParticles->Parton(9)->Eta();
            float lquark6_phi = modelParticles->Parton(9)->Phi();
            float lquark6_e = modelParticles->Parton(9)->E();
            unsigned int lquark6_index = (*permutedParticles)->JetIndex(9);

            float lepton_pt = -9999;
            float lepton_eta = -9999;
            float lepton_phi = -9999;
            float lepton_e = -9999;

            // Always check for lepton type or the code will crash.
            if (pdgId == 11) {
                lepton_pt = modelParticles->Electron(0)->Pt();
                lepton_eta = modelParticles->Electron(0)->Eta();
                lepton_phi = modelParticles->Electron(0)->Phi();
                lepton_e = modelParticles->Electron(0)->E();
            } else if (pdgId == 13) {
                lepton_pt = modelParticles->Muon(0)->Pt();
                lepton_eta = modelParticles->Muon(0)->Eta();
                lepton_phi = modelParticles->Muon(0)->Phi();
                lepton_e = modelParticles->Muon(0)->E();
            }

            // Neutrino parameters.
            float neutrino_pt = modelParticles->Neutrino(0)->Pt();
            float neutrino_eta = modelParticles->Neutrino(0)->Eta();
            float neutrino_phi = modelParticles->Neutrino(0)->Phi();
            float neutrino_e = modelParticles->Neutrino(0)->E();

//             if (isFirst) {
//               printf("----------------------------------------------------------------------------------------------\n");
//               printf("----------------------------------------Permutation %2i----------------------------------------\n",
//                      iperm);
//               printf("----------------------------------------------------------------------------------------------\n");
//               printf("                  | hadronic b quark | leptonic b quark  |  light quark 1   |  light quark 2  |\n");
//               printf("Jet index         | %16i | %17i | %16i | %15i |\n", bhad1_index, blep_index, lquark1_index,
//                      lquark2_index);
//               printf("----------------------------------------------------------------------------------------------\n");
//               printf("                  | hadronic b quark | leptonic b quark  |  light quark 1   |  light quark 2  |\n");
//               printf("Output Energies   | %16.2f | %17.2f | %16.2f | %15.2f |\n", bhad1_e, blep_e, lquark1_e, lquark2_e);
//               printf("----------------------------------------------------------------------------------------------\n");
//               printf("                  | lepton energy    | neutrino pz       | loglikelihood    |  probability    |\n");
//               printf("Other values      | %16.2f | %17.2f | %16.2f | %15.2e |\n", lepton_e, neutrino_pt, likelihood,
//                      event_probability);
//               printf("----------------------------------------------------------------------------------------------\n");
//               printf("                  | Minuit Not Conv. | Fit Aborted: NaN  | >=1 Par at Limit | Invalid TF@Conv.|\n");
//               printf("Status Code       | %16i | %17i | %16i | %15i |\n", MinuitDidNotConverge, FitAbortedDueToNaN,
//                      AtLeastOneFitParameterAtItsLimit, InvalidTransferFunctionAtConvergence);
//             }
            
            // Fill the vectors with the output variables. Note: it's
            // better/safer to store booleans as chars in ROOT files.
            klf_fit_minuit_did_not_converge.emplace_back(static_cast<char>(MinuitDidNotConverge));
            klf_fit_aborted_to_nan.emplace_back(static_cast<char>(FitAbortedDueToNaN));
            klf_fit_parameter_at_limit.emplace_back(static_cast<char>(AtLeastOneFitParameterAtItsLimit));
            klf_fit_invalid_transfer_function.emplace_back(static_cast<char>(InvalidTransferFunctionAtConvergence));

            klf_bhad1_pt.emplace_back(bhad1_pt);
            klf_bhad1_eta.emplace_back(bhad1_eta);
            klf_bhad1_phi.emplace_back(bhad1_phi);
            klf_bhad1_e.emplace_back(bhad1_e);
            klf_bhad1_jet_index.emplace_back(bhad1_index);

            klf_bhad2_pt.emplace_back(bhad2_pt);
            klf_bhad2_eta.emplace_back(bhad2_eta);
            klf_bhad2_phi.emplace_back(bhad2_phi);
            klf_bhad2_e.emplace_back(bhad2_e);
            klf_bhad2_jet_index.emplace_back(bhad2_index);

            klf_bhad3_pt.emplace_back(bhad3_pt);
            klf_bhad3_eta.emplace_back(bhad3_eta);
            klf_bhad3_phi.emplace_back(bhad3_phi);
            klf_bhad3_e.emplace_back(bhad3_e);
            klf_bhad3_jet_index.emplace_back(bhad3_index);
                        
            klf_blep_pt.emplace_back(blep_pt);
            klf_blep_eta.emplace_back(blep_eta);
            klf_blep_phi.emplace_back(blep_phi);
            klf_blep_e.emplace_back(blep_e);
            klf_blep_jet_index.emplace_back(blep_index);
            
            klf_lquark1_pt.emplace_back(lquark1_pt);
            klf_lquark1_eta.emplace_back(lquark1_eta);
            klf_lquark1_phi.emplace_back(lquark1_phi);
            klf_lquark1_e.emplace_back(lquark1_e);
            klf_lquark1_jet_index.emplace_back(lquark1_index);
            
            klf_lquark2_pt.emplace_back(lquark2_pt);
            klf_lquark2_eta.emplace_back(lquark2_eta);
            klf_lquark2_phi.emplace_back(lquark2_phi);
            klf_lquark2_e.emplace_back(lquark2_e);
            klf_lquark2_jet_index.emplace_back(lquark2_index);

            klf_lquark3_pt.emplace_back(lquark3_pt);
            klf_lquark3_eta.emplace_back(lquark3_eta);
            klf_lquark3_phi.emplace_back(lquark3_phi);
            klf_lquark3_e.emplace_back(lquark3_e);
            klf_lquark3_jet_index.emplace_back(lquark3_index);
            
            klf_lquark4_pt.emplace_back(lquark4_pt);
            klf_lquark4_eta.emplace_back(lquark4_eta);
            klf_lquark4_phi.emplace_back(lquark4_phi);
            klf_lquark4_e.emplace_back(lquark4_e);
            klf_lquark4_jet_index.emplace_back(lquark4_index);

            klf_lquark5_pt.emplace_back(lquark5_pt);
            klf_lquark5_eta.emplace_back(lquark5_eta);
            klf_lquark5_phi.emplace_back(lquark5_phi);
            klf_lquark5_e.emplace_back(lquark5_e);
            klf_lquark5_jet_index.emplace_back(lquark5_index);

            klf_lquark6_pt.emplace_back(lquark6_pt);
            klf_lquark6_eta.emplace_back(lquark6_eta);
            klf_lquark6_phi.emplace_back(lquark6_phi);
            klf_lquark6_e.emplace_back(lquark6_e);
            klf_lquark6_jet_index.emplace_back(lquark6_index);
            
            klf_lepton_pt.emplace_back(lepton_pt);
            klf_lepton_eta.emplace_back(lepton_eta);
            klf_lepton_phi.emplace_back(lepton_phi);
            klf_lepton_e.emplace_back(lepton_e);
            klf_neutrino_pt.emplace_back(neutrino_pt);
            klf_neutrino_eta.emplace_back(neutrino_eta);
            klf_neutrino_phi.emplace_back(neutrino_phi);
            klf_neutrino_e.emplace_back(neutrino_e);
            klf_loglikelihood.emplace_back(likelihood);
            klf_event_probability.emplace_back(event_probability);

            if (isFirst) isFirst = false;
        }


        klf_highest_prob_index = std::distance( klf_event_probability.begin(), max_element(klf_event_probability.begin(), klf_event_probability.end()));


        eventInd++;
        outTree.Fill();

    }

    // Go to the output file and write it.
    outFile.cd();
    std::cout << std::endl << "Writing into output root file: " << "output.root" << std::endl << std::endl;
    outTree.Write();

    // Close both input and output ROOT files.
    outFile.Close();
//    inFile->Close();
//    chain->Close();

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    return test();
    // return EXIT::SUCCES;
}
