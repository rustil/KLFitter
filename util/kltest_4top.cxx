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
//TODO: Switch to new and enhanced Lorentzvectors: https://root.cern.ch/doc/master/classROOT_1_1Math_1_1LorentzVector.html
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TRandom3.h"
#include "TMath.h"
#include "NLohmann/json.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <TH1F.h>
#include <TROOT.h>

template <typename T>
int ClearAndResize(std::vector<T>& vec, unsigned int const newSize, T defaultValue) {
    vec.clear();
    vec.resize(newSize, defaultValue);

    return 1;
}

std::tuple<int,float> deltaRMatch(TLorentzVector& child, std::vector<float>& jet_pt, std::vector<float>& jet_eta, std::vector<float>& jet_phi, std::vector<float>& jet_e) {
    float minDR = 999.;
    size_t minDR_Ind;
    for (size_t i = 0; i < jet_pt.size(); i++) {

        if(jet_pt.at(i) < 0) continue; // Needed for iterative matching

        TLorentzVector comp {
            jet_pt.at(i),
            jet_eta.at(i),
            jet_phi.at(i),
            jet_e.at(i)
        };

        float dR = child.DeltaR(comp);
        if (dR < minDR) {
            minDR = dR;
            minDR_Ind = i;
        }
    }
    return std::make_tuple(minDR_Ind, minDR);
}

int test(std::string config_path) {
    

    nlohmann::json conf;
    std::ifstream i(config_path);
    i >> conf;
    auto permutationsFilePath = conf["permutationsFilePath"].get<std::string>();
    KLFitter::Fitter fitter(permutationsFilePath);
    KLFitter::DetectorAtlas_8TeV detector{conf["detectorPath"].get<std::string>()};

    int nPartonsToBeMatched = -1;
    if (conf.find("nPartonsMatched") != conf.end())
        nPartonsToBeMatched = conf["nPartonsMatched"].get<int>();

    std::cout << "nPartons to be matched: " << nPartonsToBeMatched << std::endl;

    if (!fitter.SetDetector(&detector)) {
        std::cerr << "ERROR: Failed to set detector! Aborting" << std::endl;
        return EXIT_FAILURE;
    }

    KLFitter::LikelihoodFourTopLeptonJets likelihood{};

    // Set the likelihood properties.
    likelihood.PhysicsConstants()->SetMassTop(172.5);  // mass in GeV
    likelihood.SetBTagging(KLFitter::LikelihoodBase::BtaggingMethod::kVetoNoFit);
    likelihood.SetFlagTopMassFixed(true);
    likelihood.setFPermMethod(KLFitter::LikelihoodBase::kCustomPermutations);

    if (!fitter.SetLikelihood(&likelihood)) {
        std::cerr << "ERROR: Failed to set likelihood. Aborting." << std::endl;
        return 1;
    }

    std::vector<std::string> fileNames = conf["inputFiles"].get<std::vector<std::string>>();
    auto chain = new TChain("nominal_Loose");
    for (auto& fn : fileNames) chain->AddFile(fn.c_str());

    TTreeReader reader(chain);

    TTreeReaderValue<std::vector<int>> jet_parentghost_top_barcode(reader, "jet_parentghost_top_barcode");
    TTreeReaderValue<std::vector<int>> truth_pdgid(reader, "truth_pdgid");
    TTreeReaderValue<std::vector<int>> truth_barcode(reader, "truth_barcode");
    TTreeReaderValue<std::vector<float>> truth_pt(reader, "truth_pt");
    TTreeReaderValue<std::vector<float>> truth_eta(reader, "truth_eta");
    TTreeReaderValue<std::vector<float>> truth_phi(reader, "truth_phi");
    TTreeReaderValue<std::vector<float>> truth_m(reader, "truth_m");

    TTreeReaderValue<std::vector<float>> jet_pt(reader, "jet_pt");
    TTreeReaderValue<std::vector<float>> jet_eta(reader, "jet_eta");
    TTreeReaderValue<std::vector<float>> jet_phi(reader, "jet_phi");
    TTreeReaderValue<std::vector<float>> jet_e(reader, "jet_e");

    TTreeReaderValue<std::vector<float>> jet_mv2c10(reader, "jet_mv2c10");
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
    TTreeReaderValue<std::vector<char>> jet_isbtagged_MV2c10_77(reader, "jet_isbtagged_MV2c10_77");

    TTreeReaderValue<int> truth_nLepTop(reader, "truth_nLepTop");
    TTreeReaderValue<int> truth_nTop(reader, "truth_nTop");

    TTreeReaderValue<int> ejets_MV2c10(reader, "ejets_MV2c10");
    TTreeReaderValue<int> mujets_MV2c10(reader, "mujets_MV2c10");

    TTreeReaderValue<float> mu(reader, "mu");
    TTreeReaderValue<unsigned long long> eventNumber(reader, "eventNumber");
    TTreeReaderValue<unsigned int> runNumber(reader, "runNumber");

    TTreeReaderValue<float> weight_jvt(reader, "weight_jvt");
    TTreeReaderValue<float> weight_pileup(reader, "weight_pileup");
    TTreeReaderValue<float> weight_leptonSF(reader, "weight_leptonSF");
    TTreeReaderValue<float> weight_mc(reader, "weight_mc");
    TTreeReaderValue<float> weight_globalLeptonTriggerSF(reader, "weight_globalLeptonTriggerSF");
    TTreeReaderValue<float> weight_bTagSF_MV2c10_77(reader, "weight_bTagSF_MV2c10_77");

    // For 4t
    // TODO: In the current version eta and phi of children are mixed up. "Fix" this manually here.
    std::cout << "eta and phi still swapped!" << std::endl;
    TTreeReaderValue<float> truth_top1_pt(reader, "truth_top1_pt");
    TTreeReaderValue<float> truth_top1_eta(reader, "truth_top1_eta");
    TTreeReaderValue<float> truth_top1_phi(reader, "truth_top1_phi");
    TTreeReaderValue<int> truth_top1_isHad(reader, "truth_top1_isHad");

    TTreeReaderValue<float> truth_top1_child0_pt(reader, "truth_top1_child0_pt");
//    TTreeReaderValue<float> truth_top1_child0_eta(reader, "truth_top1_child0_eta");
//    TTreeReaderValue<float> truth_top1_child0_phi(reader, "truth_top1_child0_phi");
    TTreeReaderValue<float> truth_top1_child0_eta(reader, "truth_top1_child0_phi");
    TTreeReaderValue<float> truth_top1_child0_phi(reader, "truth_top1_child0_eta");
    TTreeReaderValue<float> truth_top1_child0_e(reader, "truth_top1_child0_e");
    TTreeReaderValue<int> truth_top1_child0_pdgid(reader, "truth_top1_child0_pdgid");
    TTreeReaderValue<float> truth_top1_child1_pt(reader, "truth_top1_child1_pt");
//    TTreeReaderValue<float> truth_top1_child1_eta(reader, "truth_top1_child1_eta");
//    TTreeReaderValue<float> truth_top1_child1_phi(reader, "truth_top1_child1_phi");
    TTreeReaderValue<float> truth_top1_child1_eta(reader, "truth_top1_child1_phi");
    TTreeReaderValue<float> truth_top1_child1_phi(reader, "truth_top1_child1_eta");
    TTreeReaderValue<float> truth_top1_child1_e(reader, "truth_top1_child1_e");
    TTreeReaderValue<int> truth_top1_child1_pdgid(reader, "truth_top1_child1_pdgid");
    TTreeReaderValue<float> truth_top1_child2_pt(reader, "truth_top1_child2_pt");
//    TTreeReaderValue<float> truth_top1_child2_eta(reader, "truth_top1_child2_eta");
//    TTreeReaderValue<float> truth_top1_child2_phi(reader, "truth_top1_child2_phi");
    TTreeReaderValue<float> truth_top1_child2_eta(reader, "truth_top1_child2_phi");
    TTreeReaderValue<float> truth_top1_child2_phi(reader, "truth_top1_child2_eta");
    TTreeReaderValue<float> truth_top1_child2_e(reader, "truth_top1_child2_e");
    TTreeReaderValue<int> truth_top1_child2_pdgid(reader, "truth_top1_child2_pdgid");


    TTreeReaderValue<float> truth_top2_pt(reader, "truth_top2_pt");
    TTreeReaderValue<float> truth_top2_eta(reader, "truth_top2_eta");
    TTreeReaderValue<float> truth_top2_phi(reader, "truth_top2_phi");
    TTreeReaderValue<int> truth_top2_isHad(reader, "truth_top2_isHad");

    TTreeReaderValue<float> truth_top2_child0_pt(reader, "truth_top2_child0_pt");
//    TTreeReaderValue<float> truth_top2_child0_eta(reader, "truth_top2_child0_eta");
//    TTreeReaderValue<float> truth_top2_child0_phi(reader, "truth_top2_child0_phi");
    TTreeReaderValue<float> truth_top2_child0_eta(reader, "truth_top2_child0_phi");
    TTreeReaderValue<float> truth_top2_child0_phi(reader, "truth_top2_child0_eta");
    TTreeReaderValue<float> truth_top2_child0_e(reader, "truth_top2_child0_e");
    TTreeReaderValue<int> truth_top2_child0_pdgid(reader, "truth_top2_child0_pdgid");
    TTreeReaderValue<float> truth_top2_child1_pt(reader, "truth_top2_child1_pt");
//    TTreeReaderValue<float> truth_top2_child1_eta(reader, "truth_top2_child1_eta");
//    TTreeReaderValue<float> truth_top2_child1_phi(reader, "truth_top2_child1_phi");
    TTreeReaderValue<float> truth_top2_child1_eta(reader, "truth_top2_child1_phi");
    TTreeReaderValue<float> truth_top2_child1_phi(reader, "truth_top2_child1_eta");
    TTreeReaderValue<float> truth_top2_child1_e(reader, "truth_top2_child1_e");
    TTreeReaderValue<int> truth_top2_child1_pdgid(reader, "truth_top2_child1_pdgid");
    TTreeReaderValue<float> truth_top2_child2_pt(reader, "truth_top2_child2_pt");
//    TTreeReaderValue<float> truth_top2_child2_eta(reader, "truth_top2_child2_eta");
//    TTreeReaderValue<float> truth_top2_child2_phi(reader, "truth_top2_child2_phi");
    TTreeReaderValue<float> truth_top2_child2_eta(reader, "truth_top2_child2_phi");
    TTreeReaderValue<float> truth_top2_child2_phi(reader, "truth_top2_child2_eta");
    TTreeReaderValue<float> truth_top2_child2_e(reader, "truth_top2_child2_e");
    TTreeReaderValue<int> truth_top2_child2_pdgid(reader, "truth_top2_child2_pdgid");


    TTreeReaderValue<float> truth_tbar1_pt(reader, "truth_tbar1_pt");
    TTreeReaderValue<float> truth_tbar1_eta(reader, "truth_tbar1_eta");
    TTreeReaderValue<float> truth_tbar1_phi(reader, "truth_tbar1_phi");
    TTreeReaderValue<int> truth_tbar1_isHad(reader, "truth_tbar1_isHad");

    TTreeReaderValue<float> truth_tbar1_child0_pt(reader, "truth_tbar1_child0_pt");
//    TTreeReaderValue<float> truth_tbar1_child0_eta(reader, "truth_tbar1_child0_eta");
//    TTreeReaderValue<float> truth_tbar1_child0_phi(reader, "truth_tbar1_child0_phi");
    TTreeReaderValue<float> truth_tbar1_child0_eta(reader, "truth_tbar1_child0_phi");
    TTreeReaderValue<float> truth_tbar1_child0_phi(reader, "truth_tbar1_child0_eta");
    TTreeReaderValue<float> truth_tbar1_child0_e(reader, "truth_tbar1_child0_e");
    TTreeReaderValue<int> truth_tbar1_child0_pdgid(reader, "truth_tbar1_child0_pdgid");
    TTreeReaderValue<float> truth_tbar1_child1_pt(reader, "truth_tbar1_child1_pt");
//    TTreeReaderValue<float> truth_tbar1_child1_eta(reader, "truth_tbar1_child1_eta");
//    TTreeReaderValue<float> truth_tbar1_child1_phi(reader, "truth_tbar1_child1_phi");
    TTreeReaderValue<float> truth_tbar1_child1_eta(reader, "truth_tbar1_child1_phi");
    TTreeReaderValue<float> truth_tbar1_child1_phi(reader, "truth_tbar1_child1_eta");
    TTreeReaderValue<float> truth_tbar1_child1_e(reader, "truth_tbar1_child1_e");
    TTreeReaderValue<int> truth_tbar1_child1_pdgid(reader, "truth_tbar1_child1_pdgid");
    TTreeReaderValue<float> truth_tbar1_child2_pt(reader, "truth_tbar1_child2_pt");
//    TTreeReaderValue<float> truth_tbar1_child2_eta(reader, "truth_tbar1_child2_eta");
//    TTreeReaderValue<float> truth_tbar1_child2_phi(reader, "truth_tbar1_child2_phi");
    TTreeReaderValue<float> truth_tbar1_child2_eta(reader, "truth_tbar1_child2_phi");
    TTreeReaderValue<float> truth_tbar1_child2_phi(reader, "truth_tbar1_child2_eta");
    TTreeReaderValue<float> truth_tbar1_child2_e(reader, "truth_tbar1_child2_e");
    TTreeReaderValue<int> truth_tbar1_child2_pdgid(reader, "truth_tbar1_child2_pdgid");


    TTreeReaderValue<float> truth_tbar2_pt(reader, "truth_tbar2_pt");
    TTreeReaderValue<float> truth_tbar2_eta(reader, "truth_tbar2_eta");
    TTreeReaderValue<float> truth_tbar2_phi(reader, "truth_tbar2_phi");
    TTreeReaderValue<int> truth_tbar2_isHad(reader, "truth_tbar2_isHad");
    TTreeReaderValue<float> truth_tbar2_child0_pt(reader, "truth_tbar2_child0_pt");
//    TTreeReaderValue<float> truth_tbar2_child0_eta(reader, "truth_tbar2_child0_eta");
//    TTreeReaderValue<float> truth_tbar2_child0_phi(reader, "truth_tbar2_child0_phi");
    TTreeReaderValue<float> truth_tbar2_child0_eta(reader, "truth_tbar2_child0_phi");
    TTreeReaderValue<float> truth_tbar2_child0_phi(reader, "truth_tbar2_child0_eta");
    TTreeReaderValue<float> truth_tbar2_child0_e(reader, "truth_tbar2_child0_e");
    TTreeReaderValue<int> truth_tbar2_child0_pdgid(reader, "truth_tbar2_child0_pdgid");
    TTreeReaderValue<float> truth_tbar2_child1_pt(reader, "truth_tbar2_child1_pt");
//    TTreeReaderValue<float> truth_tbar2_child1_eta(reader, "truth_tbar2_child1_eta");
//    TTreeReaderValue<float> truth_tbar2_child1_phi(reader, "truth_tbar2_child1_phi");
    TTreeReaderValue<float> truth_tbar2_child1_eta(reader, "truth_tbar2_child1_phi");
    TTreeReaderValue<float> truth_tbar2_child1_phi(reader, "truth_tbar2_child1_eta");
    TTreeReaderValue<float> truth_tbar2_child1_e(reader, "truth_tbar2_child1_e");
    TTreeReaderValue<int> truth_tbar2_child1_pdgid(reader, "truth_tbar2_child1_pdgid");
    TTreeReaderValue<float> truth_tbar2_child2_pt(reader, "truth_tbar2_child2_pt");
//    TTreeReaderValue<float> truth_tbar2_child2_eta(reader, "truth_tbar2_child2_eta");
//    TTreeReaderValue<float> truth_tbar2_child2_phi(reader, "truth_tbar2_child2_phi");
    TTreeReaderValue<float> truth_tbar2_child2_eta(reader, "truth_tbar2_child2_phi");
    TTreeReaderValue<float> truth_tbar2_child2_phi(reader, "truth_tbar2_child2_eta");
    TTreeReaderValue<float> truth_tbar2_child2_e(reader, "truth_tbar2_child2_e");
    TTreeReaderValue<int> truth_tbar2_child2_pdgid(reader, "truth_tbar2_child2_pdgid");
//-----------------------
    TTreeReaderValue<std::vector<float>> jet_firstghost_e(reader, "jet_firstghost_e");
    TTreeReaderValue<std::vector<float>> jet_firstghost_eta(reader, "jet_firstghost_eta");
    TTreeReaderValue<std::vector<float>> jet_firstghost_phi(reader, "jet_firstghost_phi");
    TTreeReaderValue<std::vector<float>> jet_firstghost_pt(reader, "jet_firstghost_pt");
    TTreeReaderValue<std::vector<int>> jet_firstghost_pdgId(reader, "jet_firstghost_pdgId");
    TTreeReaderValue<std::vector<float>> jet_parentghost_e(reader, "jet_parentghost_e");
    TTreeReaderValue<std::vector<float>> jet_parentghost_eta(reader, "jet_parentghost_eta");
    TTreeReaderValue<std::vector<float>> jet_parentghost_phi(reader, "jet_parentghost_phi");
    TTreeReaderValue<std::vector<float>> jet_parentghost_pt(reader, "jet_parentghost_pt");
    TTreeReaderValue<std::vector<int>> jet_parentghost_pdgId(reader, "jet_parentghost_pdgId");
    TTreeReaderValue<std::vector<int>> el_true_pdg(reader, "el_true_pdg");
    TTreeReaderValue<std::vector<float>> el_true_eta(reader, "el_true_eta");
    TTreeReaderValue<std::vector<float>> el_true_pt(reader, "el_true_pt");
    TTreeReaderValue<std::vector<int>> mu_true_pdg(reader, "mu_true_pdg");
    TTreeReaderValue<std::vector<float>> mu_true_eta(reader, "mu_true_eta");
    TTreeReaderValue<std::vector<float>> mu_true_pt(reader, "mu_true_pt");




    // For ttbar
//    TTreeReaderValue<float> truth_tbar_pt(reader, "truth_tbar_pt");
//    TTreeReaderValue<float> truth_top_pt(reader, "truth_top_pt");
//    TTreeReaderValue<float> truth_tbar_eta(reader, "truth_tbar_eta");
//    TTreeReaderValue<float> truth_top_eta(reader, "truth_top_eta");
//    TTreeReaderValue<float> truth_tbar_phi(reader, "truth_tbar_phi");
//    TTreeReaderValue<float> truth_top_phi(reader, "truth_top_phi");

//    bool isTtBar = false;
//    if (conf.find("isTtBar") != conf.end()){
//        isTtBar = conf["isTtBar"].get<bool>();
//    }





    //TODO:!
    unsigned int nEventsMax = conf["nEvents"].get<int>();
    unsigned int eventInd = 0;

    std::string outFilePath;
    if (conf.find("outputFile") != conf.end())
        outFilePath = conf["outputFile"].get<std::string>();

    /* Out section */
    TFile outFile(outFilePath.c_str(), "RECREATE");
//    TTree* outTree = (TTree*)inFile->Get("nominal_Loose")->Clone();
    TTree outTree("nominal", "nominal");
    TH1F cutFlow("CutFlow", "CutFlow", 8,0,8);
    cutFlow.Sumw2();
    TH1F cutFlowWeighted("CutFlowWeighted", "CutFlowWeighted", 8,0,8);
    cutFlowWeighted.Sumw2();

    cutFlow.GetXaxis()->SetBinLabel(1, "PreSel");
    cutFlow.GetXaxis()->SetBinLabel(2, "nLep == 1");
    cutFlow.GetXaxis()->SetBinLabel(3, "ljets Selection");
    cutFlow.GetXaxis()->SetBinLabel(4, "nJets == 10");
    cutFlow.GetXaxis()->SetBinLabel(5, "nBJets == 4");
    cutFlow.GetXaxis()->SetBinLabel(6, "4 truth tops matched");
    cutFlow.GetXaxis()->SetBinLabel(7, "10 reco jets matched");
    cutFlow.GetXaxis()->SetBinLabel(8, "triplet sizes == 301");

    cutFlowWeighted.GetXaxis()->SetBinLabel(1, "PreSel");
    cutFlowWeighted.GetXaxis()->SetBinLabel(2, "nLep == 1");
    cutFlowWeighted.GetXaxis()->SetBinLabel(3, "ljets Selection");
    cutFlowWeighted.GetXaxis()->SetBinLabel(4, "nJets == 10");
    cutFlowWeighted.GetXaxis()->SetBinLabel(5, "nBJets == 4");
    cutFlowWeighted.GetXaxis()->SetBinLabel(6, "4 truth tops matched");
    cutFlowWeighted.GetXaxis()->SetBinLabel(7, "10 reco jets matched");
    cutFlowWeighted.GetXaxis()->SetBinLabel(8, "triplet sizes == 301");

    int evt_NTruthTopMatches;
    int evt_NTruthMatchedJets;

    int klf_NCorrectlyMatchedTops;
    char klf_leptonicTopMatched;

    std::map<int, std::vector<int>> truth_tripletMap;
    std::map<int, TLorentzVector> truth_recoTopTruePermMap;
    std::vector<std::vector<int>> klf_KLFtriplets;
    std::map<int, std::vector<int>> klf_truthKLFMap;


    std::vector<float> klf_tops_pt; // ordering according to model: had123 lep
    std::vector<float> klf_tops_eta;
    std::vector<float> klf_tops_phi;
    std::vector<float> klf_tops_e;
    std::vector<char> klf_tops_matchesTruthTop;
    std::vector<int> klf_tops_matchedTruthTopOrigInd;
    std::vector<int> klf_tops_recoJetInds_0;
    std::vector<int> klf_tops_recoJetInds_1;
    std::vector<int> klf_tops_recoJetInds_2;
    
    std::vector<float> klf_matchedTruth_tops_pt;
    std::vector<float> klf_matchedTruth_tops_eta;
    std::vector<float> klf_matchedTruth_tops_phi;
    std::vector<float> klf_matchedTruth_tops_e;

    std::vector<float> klf_reco_tops_pt;
    std::vector<float> klf_reco_tops_eta;
    std::vector<float> klf_reco_tops_phi;
    std::vector<float> klf_reco_tops_e;

    std::vector<float> klf_firstghost_tops_pt;
    std::vector<float> klf_firstghost_tops_eta;
    std::vector<float> klf_firstghost_tops_phi;
    std::vector<float> klf_firstghost_tops_e;

    std::vector<float> truth_firstghost_tops_pt;
    std::vector<float> truth_firstghost_tops_eta;
    std::vector<float> truth_firstghost_tops_phi;
    std::vector<float> truth_firstghost_tops_e;

    std::vector<float> truth_hadronic_tops_pt;
    std::vector<float> truth_hadronic_tops_eta;
    std::vector<float> truth_hadronic_tops_phi;
    std::vector<float> truth_hadronic_tops_e;
    

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

    int cp_ejets_MV2c10;
    int cp_mujets_MV2c10;

    float cp_mu;
    float cp_weight_jvt;
    float cp_weight_pileup;
    float cp_weight_leptonSF;
    float cp_weight_mc;
    float cp_weight_globalLeptonTriggerSF;
    float cp_weight_bTagSF_MV2c10_77;

    unsigned long long cp_eventNumber;
    unsigned int cp_runNumber;

/*! Automated stuff */

    float cp_truth_top1_pt;
    float cp_truth_top1_eta;
    float cp_truth_top1_phi;
    float cp_truth_top1_child0_pt;
    float cp_truth_top1_child0_eta;
    float cp_truth_top1_child0_phi;
    float cp_truth_top1_child0_e;
    int cp_truth_top1_child0_pdgid;
    float cp_truth_top1_child1_pt;
    float cp_truth_top1_child1_eta;
    float cp_truth_top1_child1_phi;
    float cp_truth_top1_child1_e;
    int cp_truth_top1_child1_pdgid;
    float cp_truth_top1_child2_pt;
    float cp_truth_top1_child2_eta;
    float cp_truth_top1_child2_phi;
    float cp_truth_top1_child2_e;
    int cp_truth_top1_child2_pdgid;
    float cp_truth_top2_pt;
    float cp_truth_top2_eta;
    float cp_truth_top2_phi;
    float cp_truth_top2_child0_pt;
    float cp_truth_top2_child0_eta;
    float cp_truth_top2_child0_phi;
    float cp_truth_top2_child0_e;
    int cp_truth_top2_child0_pdgid;
    float cp_truth_top2_child1_pt;
    float cp_truth_top2_child1_eta;
    float cp_truth_top2_child1_phi;
    float cp_truth_top2_child1_e;
    int cp_truth_top2_child1_pdgid;
    float cp_truth_top2_child2_pt;
    float cp_truth_top2_child2_eta;
    float cp_truth_top2_child2_phi;
    float cp_truth_top2_child2_e;
    int cp_truth_top2_child2_pdgid;
    float cp_truth_tbar1_pt;
    float cp_truth_tbar1_eta;
    float cp_truth_tbar1_phi;
    float cp_truth_tbar1_child0_pt;
    float cp_truth_tbar1_child0_eta;
    float cp_truth_tbar1_child0_phi;
    float cp_truth_tbar1_child0_e;
    int cp_truth_tbar1_child0_pdgid;
    float cp_truth_tbar1_child1_pt;
    float cp_truth_tbar1_child1_eta;
    float cp_truth_tbar1_child1_phi;
    float cp_truth_tbar1_child1_e;
    int cp_truth_tbar1_child1_pdgid;
    float cp_truth_tbar1_child2_pt;
    float cp_truth_tbar1_child2_eta;
    float cp_truth_tbar1_child2_phi;
    float cp_truth_tbar1_child2_e;
    int cp_truth_tbar1_child2_pdgid;
    float cp_truth_tbar2_pt;
    float cp_truth_tbar2_eta;
    float cp_truth_tbar2_phi;
    float cp_truth_tbar2_child0_pt;
    float cp_truth_tbar2_child0_eta;
    float cp_truth_tbar2_child0_phi;
    float cp_truth_tbar2_child0_e;
    int cp_truth_tbar2_child0_pdgid;
    float cp_truth_tbar2_child1_pt;
    float cp_truth_tbar2_child1_eta;
    float cp_truth_tbar2_child1_phi;
    float cp_truth_tbar2_child1_e;
    int cp_truth_tbar2_child1_pdgid;
    float cp_truth_tbar2_child2_pt;
    float cp_truth_tbar2_child2_eta;
    float cp_truth_tbar2_child2_phi;
    float cp_truth_tbar2_child2_e;
    int cp_truth_tbar2_child2_pdgid;
//-----------------------
    std::vector<float> cp_jet_firstghost_e;
    std::vector<float> cp_jet_firstghost_eta;
    std::vector<float> cp_jet_firstghost_phi;
    std::vector<float> cp_jet_firstghost_pt;
    std::vector<int> cp_jet_firstghost_pdgId;
    std::vector<float> cp_jet_parentghost_e;
    std::vector<float> cp_jet_parentghost_eta;
    std::vector<float> cp_jet_parentghost_phi;
    std::vector<float> cp_jet_parentghost_pt;
    std::vector<int> cp_jet_parentghost_pdgId;
    std::vector<int> cp_el_true_pdg;
    std::vector<float> cp_el_true_eta;
    std::vector<float> cp_el_true_pt;
    std::vector<int> cp_mu_true_pdg;
    std::vector<float> cp_mu_true_eta;
    std::vector<float> cp_mu_true_pt;

    //----------------------------

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

    outTree.Branch("evt_NTruthTopMatches", &evt_NTruthTopMatches);
    outTree.Branch("evt_NTruthMatchedJets", &evt_NTruthMatchedJets);

    outTree.Branch("klf_NCorrectlyMatchedTops", &klf_NCorrectlyMatchedTops);
    outTree.Branch("klf_leptonicTopMatched", &klf_leptonicTopMatched);

    outTree.Branch("klf_tops_pt", &klf_tops_pt);
    outTree.Branch("klf_tops_eta", &klf_tops_eta);
    outTree.Branch("klf_tops_phi", &klf_tops_phi);
    outTree.Branch("klf_tops_e", &klf_tops_e);
    outTree.Branch("klf_tops_matchesTruthTop", &klf_tops_matchesTruthTop);
    outTree.Branch("klf_tops_matchedTruthTopOrigInd", &klf_tops_matchedTruthTopOrigInd);
    outTree.Branch("klf_tops_recoJetInds_0", &klf_tops_recoJetInds_0);
    outTree.Branch("klf_tops_recoJetInds_1", &klf_tops_recoJetInds_1);
    outTree.Branch("klf_tops_recoJetInds_2", &klf_tops_recoJetInds_2);

    outTree.Branch("klf_matchedTruth_tops_pt", &klf_matchedTruth_tops_pt);
    outTree.Branch("klf_matchedTruth_tops_eta", &klf_matchedTruth_tops_eta);
    outTree.Branch("klf_matchedTruth_tops_phi", &klf_matchedTruth_tops_phi);
    outTree.Branch("klf_matchedTruth_tops_e", &klf_matchedTruth_tops_e);

    outTree.Branch("klf_reco_tops_pt", &klf_reco_tops_pt);
    outTree.Branch("klf_reco_tops_eta", &klf_reco_tops_eta);
    outTree.Branch("klf_reco_tops_phi", &klf_reco_tops_phi);
    outTree.Branch("klf_reco_tops_e", &klf_reco_tops_e);

    outTree.Branch("klf_firstghost_tops_pt", &klf_firstghost_tops_pt);
    outTree.Branch("klf_firstghost_tops_eta", &klf_firstghost_tops_eta);
    outTree.Branch("klf_firstghost_tops_phi", &klf_firstghost_tops_phi);
    outTree.Branch("klf_firstghost_tops_e", &klf_firstghost_tops_e);

    outTree.Branch("truth_firstghost_tops_pt", &truth_firstghost_tops_pt);
    outTree.Branch("truth_firstghost_tops_eta", &truth_firstghost_tops_eta);
    outTree.Branch("truth_firstghost_tops_phi", &truth_firstghost_tops_phi);
    outTree.Branch("truth_firstghost_tops_e", &truth_firstghost_tops_e);

    outTree.Branch("truth_hadronic_tops_pt", &truth_hadronic_tops_pt);
    outTree.Branch("truth_hadronic_tops_eta", &truth_hadronic_tops_eta);
    outTree.Branch("truth_hadronic_tops_phi", &truth_hadronic_tops_phi);
    outTree.Branch("truth_hadronic_tops_e", &truth_hadronic_tops_e);

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

    outTree.Branch("ejets_MV2c10", &cp_ejets_MV2c10);
    outTree.Branch("mujets_MV2c10", &cp_mujets_MV2c10);

    outTree.Branch("mu", &cp_mu);
    outTree.Branch("weight_jvt", &cp_weight_jvt);
    outTree.Branch("weight_pileup", &cp_weight_pileup);
    outTree.Branch("weight_leptonSF", &cp_weight_leptonSF);
    outTree.Branch("weight_mc", &cp_weight_mc);
    outTree.Branch("weight_globalLeptonTriggerSF", &cp_weight_globalLeptonTriggerSF);
    outTree.Branch("weight_bTagSF_MV2c10_77", &cp_weight_bTagSF_MV2c10_77);

    outTree.Branch("eventNumber", &cp_eventNumber);
    outTree.Branch("runNumber", &cp_runNumber);

    int nPartonsMatcheable;
    outTree.Branch("nPartonsMatcheable", &nPartonsMatcheable);

    outTree.Branch("truth_top1_pt", &cp_truth_top1_pt);
    outTree.Branch("truth_top1_eta", &cp_truth_top1_eta);
    outTree.Branch("truth_top1_phi", &cp_truth_top1_phi);
    outTree.Branch("truth_top1_child0_pt", &cp_truth_top1_child0_pt);
    outTree.Branch("truth_top1_child0_eta", &cp_truth_top1_child0_eta);
    outTree.Branch("truth_top1_child0_phi", &cp_truth_top1_child0_phi);
    outTree.Branch("truth_top1_child0_e", &cp_truth_top1_child0_e);
    outTree.Branch("truth_top1_child0_pdgid", &cp_truth_top1_child0_pdgid);
    outTree.Branch("truth_top1_child1_pt", &cp_truth_top1_child1_pt);
    outTree.Branch("truth_top1_child1_eta", &cp_truth_top1_child1_eta);
    outTree.Branch("truth_top1_child1_phi", &cp_truth_top1_child1_phi);
    outTree.Branch("truth_top1_child1_e", &cp_truth_top1_child1_e);
    outTree.Branch("truth_top1_child1_pdgid", &cp_truth_top1_child1_pdgid);
    outTree.Branch("truth_top1_child2_pt", &cp_truth_top1_child2_pt);
    outTree.Branch("truth_top1_child2_eta", &cp_truth_top1_child2_eta);
    outTree.Branch("truth_top1_child2_phi", &cp_truth_top1_child2_phi);
    outTree.Branch("truth_top1_child2_e", &cp_truth_top1_child2_e);
    outTree.Branch("truth_top1_child2_pdgid", &cp_truth_top1_child2_pdgid);
    outTree.Branch("truth_top2_pt", &cp_truth_top2_pt);
    outTree.Branch("truth_top2_eta", &cp_truth_top2_eta);
    outTree.Branch("truth_top2_phi", &cp_truth_top2_phi);
    outTree.Branch("truth_top2_child0_pt", &cp_truth_top2_child0_pt);
    outTree.Branch("truth_top2_child0_eta", &cp_truth_top2_child0_eta);
    outTree.Branch("truth_top2_child0_phi", &cp_truth_top2_child0_phi);
    outTree.Branch("truth_top2_child0_e", &cp_truth_top2_child0_e);
    outTree.Branch("truth_top2_child0_pdgid", &cp_truth_top2_child0_pdgid);
    outTree.Branch("truth_top2_child1_pt", &cp_truth_top2_child1_pt);
    outTree.Branch("truth_top2_child1_eta", &cp_truth_top2_child1_eta);
    outTree.Branch("truth_top2_child1_phi", &cp_truth_top2_child1_phi);
    outTree.Branch("truth_top2_child1_e", &cp_truth_top2_child1_e);
    outTree.Branch("truth_top2_child1_pdgid", &cp_truth_top2_child1_pdgid);
    outTree.Branch("truth_top2_child2_pt", &cp_truth_top2_child2_pt);
    outTree.Branch("truth_top2_child2_eta", &cp_truth_top2_child2_eta);
    outTree.Branch("truth_top2_child2_phi", &cp_truth_top2_child2_phi);
    outTree.Branch("truth_top2_child2_e", &cp_truth_top2_child2_e);
    outTree.Branch("truth_top2_child2_pdgid", &cp_truth_top2_child2_pdgid);
    outTree.Branch("truth_tbar1_pt", &cp_truth_tbar1_pt);
    outTree.Branch("truth_tbar1_eta", &cp_truth_tbar1_eta);
    outTree.Branch("truth_tbar1_phi", &cp_truth_tbar1_phi);
    outTree.Branch("truth_tbar1_child0_pt", &cp_truth_tbar1_child0_pt);
    outTree.Branch("truth_tbar1_child0_eta", &cp_truth_tbar1_child0_eta);
    outTree.Branch("truth_tbar1_child0_phi", &cp_truth_tbar1_child0_phi);
    outTree.Branch("truth_tbar1_child0_e", &cp_truth_tbar1_child0_e);
    outTree.Branch("truth_tbar1_child0_pdgid", &cp_truth_tbar1_child0_pdgid);
    outTree.Branch("truth_tbar1_child1_pt", &cp_truth_tbar1_child1_pt);
    outTree.Branch("truth_tbar1_child1_eta", &cp_truth_tbar1_child1_eta);
    outTree.Branch("truth_tbar1_child1_phi", &cp_truth_tbar1_child1_phi);
    outTree.Branch("truth_tbar1_child1_e", &cp_truth_tbar1_child1_e);
    outTree.Branch("truth_tbar1_child1_pdgid", &cp_truth_tbar1_child1_pdgid);
    outTree.Branch("truth_tbar1_child2_pt", &cp_truth_tbar1_child2_pt);
    outTree.Branch("truth_tbar1_child2_eta", &cp_truth_tbar1_child2_eta);
    outTree.Branch("truth_tbar1_child2_phi", &cp_truth_tbar1_child2_phi);
    outTree.Branch("truth_tbar1_child2_e", &cp_truth_tbar1_child2_e);
    outTree.Branch("truth_tbar1_child2_pdgid", &cp_truth_tbar1_child2_pdgid);
    outTree.Branch("truth_tbar2_pt", &cp_truth_tbar2_pt);
    outTree.Branch("truth_tbar2_eta", &cp_truth_tbar2_eta);
    outTree.Branch("truth_tbar2_phi", &cp_truth_tbar2_phi);
    outTree.Branch("truth_tbar2_child0_pt", &cp_truth_tbar2_child0_pt);
    outTree.Branch("truth_tbar2_child0_eta", &cp_truth_tbar2_child0_eta);
    outTree.Branch("truth_tbar2_child0_phi", &cp_truth_tbar2_child0_phi);
    outTree.Branch("truth_tbar2_child0_e", &cp_truth_tbar2_child0_e);
    outTree.Branch("truth_tbar2_child0_pdgid", &cp_truth_tbar2_child0_pdgid);
    outTree.Branch("truth_tbar2_child1_pt", &cp_truth_tbar2_child1_pt);
    outTree.Branch("truth_tbar2_child1_eta", &cp_truth_tbar2_child1_eta);
    outTree.Branch("truth_tbar2_child1_phi", &cp_truth_tbar2_child1_phi);
    outTree.Branch("truth_tbar2_child1_e", &cp_truth_tbar2_child1_e);
    outTree.Branch("truth_tbar2_child1_pdgid", &cp_truth_tbar2_child1_pdgid);
    outTree.Branch("truth_tbar2_child2_pt", &cp_truth_tbar2_child2_pt);
    outTree.Branch("truth_tbar2_child2_eta", &cp_truth_tbar2_child2_eta);
    outTree.Branch("truth_tbar2_child2_phi", &cp_truth_tbar2_child2_phi);
    outTree.Branch("truth_tbar2_child2_e", &cp_truth_tbar2_child2_e);
    outTree.Branch("truth_tbar2_child2_pdgid", &cp_truth_tbar2_child2_pdgid);
//-----------------------
    outTree.Branch("jet_firstghost_e", &cp_jet_firstghost_e);
    outTree.Branch("jet_firstghost_eta", &cp_jet_firstghost_eta);
    outTree.Branch("jet_firstghost_phi", &cp_jet_firstghost_phi);
    outTree.Branch("jet_firstghost_pt", &cp_jet_firstghost_pt);
    outTree.Branch("jet_firstghost_pdgId", &cp_jet_firstghost_pdgId);
    outTree.Branch("jet_parentghost_e", &cp_jet_parentghost_e);
    outTree.Branch("jet_parentghost_eta", &cp_jet_parentghost_eta);
    outTree.Branch("jet_parentghost_phi", &cp_jet_parentghost_phi);
    outTree.Branch("jet_parentghost_pt", &cp_jet_parentghost_pt);
    outTree.Branch("jet_parentghost_pdgId", &cp_jet_parentghost_pdgId);
    outTree.Branch("el_true_pdg", &cp_el_true_pdg);
    outTree.Branch("el_true_eta", &cp_el_true_eta);
    outTree.Branch("el_true_pt", &cp_el_true_pt);
    outTree.Branch("mu_true_pdg", &cp_mu_true_pdg);
    outTree.Branch("mu_true_eta", &cp_mu_true_eta);
    outTree.Branch("mu_true_pt", &cp_mu_true_pt);

//    std::vector<float> globalChildRecoDR;
//    std::vector<float> iterativeChildRecoDR;
//    std::vector<int> globalChildRecoIndMult;
//    bool iterativeChildRecoDR_lt_01;
//
//    std::vector<float> globalChildFirstGhostDR;
//    std::vector<float> globalPDGChildFirstGhostDR;
//    std::vector<int> globalPDGChildFirstGhostIndMult;
//
//    std::vector<float> iterativePDGChildFirstGhostDRs;
//
//
//    outTree.Branch("evt_globalChildRecoDR", &globalChildRecoDR);
//    outTree.Branch("evt_iterativeChildRecoDR", &iterativeChildRecoDR);
//    outTree.Branch("evt_globalChildRecoIndMult", &globalChildRecoIndMult);
//    outTree.Branch("evt_iterativeChildRecoDR_lt_01", &iterativeChildRecoDR_lt_01);
//
//
//    outTree.Branch("evt_globalChildFirstGhostDR", &globalChildFirstGhostDR);
//    outTree.Branch("evt_globalPDGChildFirstGhostDR", &globalPDGChildFirstGhostDR);
//    outTree.Branch("evt_globalPDGChildFirstGhostIndMult", &globalPDGChildFirstGhostIndMult);
//
//    outTree.Branch("evt_iterativePDGChildFirstGhostDR", &iterativePDGChildFirstGhostDRs);

    time_t rawtime = time(NULL);
    struct tm * timeinfo;
    char buffer [80];



    while (reader.Next()) {

        if (nEventsMax >= 0 &&  eventInd >= nEventsMax) continue;

        evt_NTruthTopMatches = 0;
        evt_NTruthMatchedJets = 0;

        klf_NCorrectlyMatchedTops = 0;
        klf_leptonicTopMatched = 0;

        truth_tripletMap.clear();
        truth_recoTopTruePermMap.clear();
        klf_KLFtriplets.clear();
        klf_truthKLFMap.clear();

        int nLep = mu_pt->size() + el_pt->size();
        int nBJets = std::count_if(jet_isbtagged_MV2c10_77->begin(),
                                   jet_isbtagged_MV2c10_77->end(),
                                   [](char flag) {return flag;});

        // CUTS
        cutFlow.Fill(0);
        float weight = *weight_jvt * *weight_pileup * *weight_leptonSF * *weight_mc * *weight_bTagSF_MV2c10_77;
        cutFlowWeighted.Fill(0., weight);

        if(nLep != 1) continue;
        else {
            cutFlow.Fill(1);
            cutFlowWeighted.Fill(1, weight);
        }

        if (!(*ejets_MV2c10 || *mujets_MV2c10)) continue;
        else {
            cutFlow.Fill(2);
            cutFlowWeighted.Fill(2, weight);
        }

        if (jet_pt->size() != 10) continue;
        else {
            cutFlow.Fill(3);
            cutFlowWeighted.Fill(3, weight);
        }

        if (nBJets != 4) continue;
        else {
            cutFlow.Fill(4);
            cutFlowWeighted.Fill(4, weight);
        }


        for (int k =0 ; k < jet_parentghost_top_barcode->size(); k++) {
            auto key = jet_parentghost_top_barcode->at(k);
            if (key == 0) continue;
            auto search = truth_tripletMap.find(key);
            if (search == truth_tripletMap.end()) {
                truth_tripletMap[key] = {k};
            } else {
                truth_tripletMap[key].emplace_back(k);
            }
        }

        evt_NTruthTopMatches = truth_tripletMap.size();
        std::for_each(truth_tripletMap.begin(), truth_tripletMap.end(),
                      [&evt_NTruthMatchedJets] (auto& vec) { evt_NTruthMatchedJets += vec.second.size();}
        );

        // CUT
        if (truth_tripletMap.size() !=4) {
            continue;
        } else {
            cutFlow.Fill(5);
            cutFlowWeighted.Fill(5, weight);
        }

        if (evt_NTruthMatchedJets != nPartonsToBeMatched) {
            continue;
        } else {
            cutFlow.Fill(6);
            cutFlowWeighted.Fill(6, weight);
        }


        int tripletSizes = 0; // nTriplet *  100 + nDuplet * 10 +  nSinglet = 301 for 1LOS tree level FS. (3 tops with three jets and 1 w/ one)
        std::for_each(truth_tripletMap.begin(), truth_tripletMap.end(), [&tripletSizes] (auto& pair) {
            int size = pair.second.size();
            if (size == 3) tripletSizes += 100;
            else if (size == 2) tripletSizes += 10;
            else if (size == 1) tripletSizes += 1;
        });

        // CUT
        if (tripletSizes != 301) {
            continue;
        } else {
            cutFlow.Fill(7);
            cutFlowWeighted.Fill(7, weight);
        }


        KLFitter::Particles particles{};

        if ((eventInd % 1) == 0) {
            time(&rawtime);
            timeinfo = localtime(&rawtime);
            strftime(buffer, 80, "%FT%T", timeinfo);
            printf("%10s: [%10i | %10i]\n", buffer, eventInd, nEventsMax);
        }

        std::map<int, TLorentzVector> truthTopMap {};
        for (auto& pair: truth_tripletMap) {
            int truthTopInd = -1;
            for (int k = 0; k < truth_barcode->size(); k++) {
                if (pair.first == truth_barcode->at(k)) {
                    truthTopInd = k;
                    break;
                }
            }

            if (truthTopInd >= 0 && truthTopInd < truth_pt->size()) {
                TLorentzVector vec = {};
                vec.SetPtEtaPhiM(truth_pt->at(truthTopInd),
                                 truth_eta->at(truthTopInd),
                                 truth_phi->at(truthTopInd),
                                 truth_m->at(truthTopInd)
                );
                truthTopMap[pair.first] = vec;
            } else std::cerr << "Error for truthTopMap for barcode: " << pair.first << " in Event: " << *eventNumber << std::endl;
        }

        // Reco to Truth Map

        for (auto& pair: truth_tripletMap) {
            if (pair.second.size() == 3) { //hadronic
                std::vector<TLorentzVector> vecs = {};
                for (int ind : pair.second) {

                    if (ind >= 0 && ind < jet_pt->size()) {
                        TLorentzVector vec = {};
                        vec.SetPtEtaPhiE(
                                jet_pt->at(ind),
                                jet_eta->at(ind),
                                jet_phi->at(ind),
                                jet_e->at(ind)
                        );
                        vecs.emplace_back(vec);
                    } else std::cerr << "Error for Reco Truth Map for barcode: " << pair.first << " in Event: " << *eventNumber << std::endl;
                }
                truth_recoTopTruePermMap[pair.first] = (vecs.at(0) + vecs.at(1) + vecs.at(2));
            } else if (pair.second.size() == 1) { // leptonic - do nothing since MET is well.. missing.
//                TLorentzVector b = {};
//                b.SetPtEtaPhiE(jet_pt->at(pair.second.at(0)),
//                               jet_eta->at(pair.second.at(0)),
//                               jet_phi->at(pair.second.at(0)),
//                               jet_e->at(pair.second.at(0))
//                        );
                                } else std::cerr << "Unusual amount of matching children: " << pair.second.size() << std::endl;
            }


        ClearAndResize<float>(klf_tops_pt, 4, -99);
        ClearAndResize<float>(klf_tops_eta,4, -99);
        ClearAndResize<float>(klf_tops_phi ,4, -99);
        ClearAndResize<float>(klf_tops_e ,4, -99);
        ClearAndResize<char>(klf_tops_matchesTruthTop ,4, static_cast<char>(false));
        ClearAndResize<int>(klf_tops_matchedTruthTopOrigInd , 4, -99);
        ClearAndResize<int>(klf_tops_recoJetInds_0 ,4, -99);
        ClearAndResize<int>(klf_tops_recoJetInds_1 ,4, -99);
        ClearAndResize<int>(klf_tops_recoJetInds_2 ,4, -99);

        ClearAndResize<float>(klf_matchedTruth_tops_pt , 4, -99);
        ClearAndResize<float>(klf_matchedTruth_tops_eta, 4, -99);
        ClearAndResize<float>(klf_matchedTruth_tops_phi , 4, -99);
        ClearAndResize<float>(klf_matchedTruth_tops_e , 4, -99);

        ClearAndResize<float>(klf_reco_tops_pt ,4, -99);
        ClearAndResize<float>(klf_reco_tops_eta,4, -99);
        ClearAndResize<float>(klf_reco_tops_phi ,4, -99);
        ClearAndResize<float>(klf_reco_tops_e ,4, -99);

        ClearAndResize<float>(klf_firstghost_tops_pt ,4, -99);
        ClearAndResize<float>(klf_firstghost_tops_eta,4, -99);
        ClearAndResize<float>(klf_firstghost_tops_phi ,4, -99);
        ClearAndResize<float>(klf_firstghost_tops_e ,4, -99);

        ClearAndResize<float>(truth_firstghost_tops_pt ,3, -99);
        ClearAndResize<float>(truth_firstghost_tops_eta,3, -99);
        ClearAndResize<float>(truth_firstghost_tops_phi ,3, -99);
        ClearAndResize<float>(truth_firstghost_tops_e ,3, -99);

        ClearAndResize<float>(truth_hadronic_tops_pt ,3, -99);
        ClearAndResize<float>(truth_hadronic_tops_eta,3, -99);
        ClearAndResize<float>(truth_hadronic_tops_phi ,3, -99);
        ClearAndResize<float>(truth_hadronic_tops_e ,3, -99);
        
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

        // Automated stuff
        cp_truth_top1_pt = *truth_top1_pt;
        cp_truth_top1_eta = *truth_top1_eta;
        cp_truth_top1_phi = *truth_top1_phi;
        cp_truth_top1_child0_pt = *truth_top1_child0_pt;
        cp_truth_top1_child0_eta = *truth_top1_child0_eta;
        cp_truth_top1_child0_phi = *truth_top1_child0_phi;
        cp_truth_top1_child0_e = *truth_top1_child0_e;
        cp_truth_top1_child0_pdgid = *truth_top1_child0_pdgid;
        cp_truth_top1_child1_pt = *truth_top1_child1_pt;
        cp_truth_top1_child1_eta = *truth_top1_child1_eta;
        cp_truth_top1_child1_phi = *truth_top1_child1_phi;
        cp_truth_top1_child1_e = *truth_top1_child1_e;
        cp_truth_top1_child1_pdgid = *truth_top1_child1_pdgid;
        cp_truth_top1_child2_pt = *truth_top1_child2_pt;
        cp_truth_top1_child2_eta = *truth_top1_child2_eta;
        cp_truth_top1_child2_phi = *truth_top1_child2_phi;
        cp_truth_top1_child2_e = *truth_top1_child2_e;
        cp_truth_top1_child2_pdgid = *truth_top1_child2_pdgid;
        cp_truth_top2_pt = *truth_top2_pt;
        cp_truth_top2_eta = *truth_top2_eta;
        cp_truth_top2_phi = *truth_top2_phi;
        cp_truth_top2_child0_pt = *truth_top2_child0_pt;
        cp_truth_top2_child0_eta = *truth_top2_child0_eta;
        cp_truth_top2_child0_phi = *truth_top2_child0_phi;
        cp_truth_top2_child0_e = *truth_top2_child0_e;
        cp_truth_top2_child0_pdgid = *truth_top2_child0_pdgid;
        cp_truth_top2_child1_pt = *truth_top2_child1_pt;
        cp_truth_top2_child1_eta = *truth_top2_child1_eta;
        cp_truth_top2_child1_phi = *truth_top2_child1_phi;
        cp_truth_top2_child1_e = *truth_top2_child1_e;
        cp_truth_top2_child1_pdgid = *truth_top2_child1_pdgid;
        cp_truth_top2_child2_pt = *truth_top2_child2_pt;
        cp_truth_top2_child2_eta = *truth_top2_child2_eta;
        cp_truth_top2_child2_phi = *truth_top2_child2_phi;
        cp_truth_top2_child2_e = *truth_top2_child2_e;
        cp_truth_top2_child2_pdgid = *truth_top2_child2_pdgid;
        cp_truth_tbar1_pt = *truth_tbar1_pt;
        cp_truth_tbar1_eta = *truth_tbar1_eta;
        cp_truth_tbar1_phi = *truth_tbar1_phi;
        cp_truth_tbar1_child0_pt = *truth_tbar1_child0_pt;
        cp_truth_tbar1_child0_eta = *truth_tbar1_child0_eta;
        cp_truth_tbar1_child0_phi = *truth_tbar1_child0_phi;
        cp_truth_tbar1_child0_e = *truth_tbar1_child0_e;
        cp_truth_tbar1_child0_pdgid = *truth_tbar1_child0_pdgid;
        cp_truth_tbar1_child1_pt = *truth_tbar1_child1_pt;
        cp_truth_tbar1_child1_eta = *truth_tbar1_child1_eta;
        cp_truth_tbar1_child1_phi = *truth_tbar1_child1_phi;
        cp_truth_tbar1_child1_e = *truth_tbar1_child1_e;
        cp_truth_tbar1_child1_pdgid = *truth_tbar1_child1_pdgid;
        cp_truth_tbar1_child2_pt = *truth_tbar1_child2_pt;
        cp_truth_tbar1_child2_eta = *truth_tbar1_child2_eta;
        cp_truth_tbar1_child2_phi = *truth_tbar1_child2_phi;
        cp_truth_tbar1_child2_e = *truth_tbar1_child2_e;
        cp_truth_tbar1_child2_pdgid = *truth_tbar1_child2_pdgid;
        cp_truth_tbar2_pt = *truth_tbar2_pt;
        cp_truth_tbar2_eta = *truth_tbar2_eta;
        cp_truth_tbar2_phi = *truth_tbar2_phi;
        cp_truth_tbar2_child0_pt = *truth_tbar2_child0_pt;
        cp_truth_tbar2_child0_eta = *truth_tbar2_child0_eta;
        cp_truth_tbar2_child0_phi = *truth_tbar2_child0_phi;
        cp_truth_tbar2_child0_e = *truth_tbar2_child0_e;
        cp_truth_tbar2_child0_pdgid = *truth_tbar2_child0_pdgid;
        cp_truth_tbar2_child1_pt = *truth_tbar2_child1_pt;
        cp_truth_tbar2_child1_eta = *truth_tbar2_child1_eta;
        cp_truth_tbar2_child1_phi = *truth_tbar2_child1_phi;
        cp_truth_tbar2_child1_e = *truth_tbar2_child1_e;
        cp_truth_tbar2_child1_pdgid = *truth_tbar2_child1_pdgid;
        cp_truth_tbar2_child2_pt = *truth_tbar2_child2_pt;
        cp_truth_tbar2_child2_eta = *truth_tbar2_child2_eta;
        cp_truth_tbar2_child2_phi = *truth_tbar2_child2_phi;
        cp_truth_tbar2_child2_e = *truth_tbar2_child2_e;
        cp_truth_tbar2_child2_pdgid = *truth_tbar2_child2_pdgid;
//-----------------------
        cp_jet_firstghost_e = *jet_firstghost_e;
        cp_jet_firstghost_eta = *jet_firstghost_eta;
        cp_jet_firstghost_phi = *jet_firstghost_phi;
        cp_jet_firstghost_pt = *jet_firstghost_pt;
        cp_jet_firstghost_pdgId = *jet_firstghost_pdgId;
        cp_jet_parentghost_e = *jet_parentghost_e;
        cp_jet_parentghost_eta = *jet_parentghost_eta;
        cp_jet_parentghost_phi = *jet_parentghost_phi;
        cp_jet_parentghost_pt = *jet_parentghost_pt;
        cp_jet_parentghost_pdgId = *jet_parentghost_pdgId;
        cp_el_true_pdg = *el_true_pdg;
        cp_el_true_eta = *el_true_eta;
        cp_el_true_pt = *el_true_pt;
        cp_mu_true_pdg = *mu_true_pdg;
        cp_mu_true_eta = *mu_true_eta;
        cp_mu_true_pt = *mu_true_pt;

        //---

//        if (!isTtBar) {
//            cp_truth_top_pt = {*truth_top1_pt, *truth_top2_pt, *truth_tbar1_pt, *truth_tbar2_pt};
//            cp_truth_top_eta = {*truth_top1_eta, *truth_top2_eta, *truth_tbar1_eta, *truth_tbar2_eta};
//            cp_truth_top_phi = {*truth_top1_phi, *truth_top2_phi, *truth_tbar1_phi, *truth_tbar2_phi};
//            cp_truth_top_type = {1, 1, -1, -1};
//        } else {
//            cp_truth_top_pt = {*truth_top_pt, *truth_tbar_pt};
//            cp_truth_top_eta = {*truth_top_eta, *truth_tbar_eta};
//            cp_truth_top_phi = {*truth_top1_phi, *truth_tbar_phi};
//            cp_truth_top_type = {1, -1};
//        }

        cp_ejets_MV2c10 = *ejets_MV2c10;
        cp_mujets_MV2c10 = *mujets_MV2c10;

        cp_mu = *mu;
        cp_weight_jvt = *weight_jvt;
        cp_weight_pileup = *weight_pileup;
        cp_weight_leptonSF = *weight_leptonSF;
        cp_weight_mc = *weight_mc;
        cp_weight_globalLeptonTriggerSF = *weight_globalLeptonTriggerSF;
        cp_weight_bTagSF_MV2c10_77 = *weight_bTagSF_MV2c10_77;

        cp_eventNumber = *eventNumber;
        cp_runNumber = *runNumber;

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
        if (!fitter.SetET_miss_XY_SumET(met_ex, met_ey, *met_met)) { // TODO: This is wrong and should be sumET instead of metmet. Not in ntuple currently
            std::cerr << "ERROR: Failed to add MET to fitter. Aborting." << std::endl;
            return 1;
        }

//        // Loop over all permutations.
        const int nperm = fitter.CustomPermutations()->NPermutations();
        bool isFirst = true;

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

        for (int iperm = 0; iperm < nperm; iperm++) {

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



        auto highestPerm = fitter.CustomPermutations()->getFJsonPermutations().at(klf_highest_prob_index);
        // Build KLF best permutation map
        auto KLFJetIndices = highestPerm["light_indices"].get<std::vector<int>>();
        auto KLFBIndices = highestPerm["b_indices"].get<std::vector<int>>();

        // The indices in the permutation for KLF are not so much indices as they are the order for light and b jets. Since light and b-jets are together in the jet collection we need to map that.
        std::vector<int> recoBIndices {};
        std::vector<int> recoJetIndices {};
        for (int j = 0; j < jet_isbtagged_MV2c10_77->size(); j++) {
            if (jet_isbtagged_MV2c10_77->at(j) == 1) recoBIndices.emplace_back(j);
            else recoJetIndices.emplace_back(j);
        }


        for (int j = 0; j < recoJetIndices.size(); j += 2) {
            int permBtoRecoB = KLFBIndices.at(j / 2);
            std::vector<int> vec = {
                    recoBIndices.at(permBtoRecoB),
                    recoJetIndices.at(KLFJetIndices.at(j)),
                    recoJetIndices.at(KLFJetIndices.at(j+1)),
            };

            klf_KLFtriplets.emplace_back(vec);
        }

        // Compare leptonic b index
        int KLFLepBIndex = recoBIndices.at(KLFBIndices.at(KLFBIndices.size()-1));

        // Add Leptonic bjet to 'triplets'
        std::vector<int> lepTriplet = {KLFLepBIndex};
        klf_KLFtriplets.emplace_back(lepTriplet);

        int TrueLepBIndex = -1;
//        true leptonic b
        for (auto& pair : truth_tripletMap) {
            if (pair.second.size() == 1) {
                TrueLepBIndex = pair.second.at(0);
                break;
            }
        }

        // TODO: Match truth top to KLF top. Get Number of matched tops etc.

// Sort perms. Then look for match.

        for (auto& triplet : klf_KLFtriplets) {
            std::sort(triplet.begin(), triplet.end());
        }

        for (auto& pair : truth_tripletMap) {
            std::sort(pair.second.begin(), pair.second.end());

            if (pair.second.size() == 1) { //leptonic
                if (KLFLepBIndex == TrueLepBIndex) klf_truthKLFMap[pair.first] = {KLFLepBIndex};
            } else { //hadronic
                for (auto& triplet : klf_KLFtriplets) {
                    if ( triplet == pair.second) klf_truthKLFMap[pair.first] = triplet;
                }
            }
        }


        klf_NCorrectlyMatchedTops = klf_truthKLFMap.size();
        klf_leptonicTopMatched = static_cast<char>(KLFLepBIndex == TrueLepBIndex);

        /*****************
         * Add KLF tops in KLF order.
         */
        TLorentzVector t1, t1_1, t1_2, t1_3 = {};
        t1_1.SetPtEtaPhiE(
                klf_bhad1_pt.at(klf_highest_prob_index),
                klf_bhad1_eta.at(klf_highest_prob_index),
                klf_bhad1_phi.at(klf_highest_prob_index),
                klf_bhad1_e.at(klf_highest_prob_index)
        );

        t1_2.SetPtEtaPhiE(
                klf_lquark1_pt.at(klf_highest_prob_index),
                klf_lquark1_eta.at(klf_highest_prob_index),
                klf_lquark1_phi.at(klf_highest_prob_index),
                klf_lquark1_e.at(klf_highest_prob_index)
        );

        t1_3.SetPtEtaPhiE(
                klf_lquark2_pt.at(klf_highest_prob_index),
                klf_lquark2_eta.at(klf_highest_prob_index),
                klf_lquark2_phi.at(klf_highest_prob_index),
                klf_lquark2_e.at(klf_highest_prob_index)
        );

        t1 = t1_1 + t1_2 + t1_3;

        TLorentzVector t2, t2_1, t2_2, t2_3 = {};
        t2_1.SetPtEtaPhiE(
                klf_bhad2_pt.at(klf_highest_prob_index),
                klf_bhad2_eta.at(klf_highest_prob_index),
                klf_bhad2_phi.at(klf_highest_prob_index),
                klf_bhad2_e.at(klf_highest_prob_index)
        );

        t2_2.SetPtEtaPhiE(
                klf_lquark3_pt.at(klf_highest_prob_index),
                klf_lquark3_eta.at(klf_highest_prob_index),
                klf_lquark3_phi.at(klf_highest_prob_index),
                klf_lquark3_e.at(klf_highest_prob_index)
        );

        t2_3.SetPtEtaPhiE(
                klf_lquark4_pt.at(klf_highest_prob_index),
                klf_lquark4_eta.at(klf_highest_prob_index),
                klf_lquark4_phi.at(klf_highest_prob_index),
                klf_lquark4_e.at(klf_highest_prob_index)
        );

        t2 = t2_1 + t2_2 + t2_3;

        TLorentzVector t3, t3_1, t3_2, t3_3 = {};
        t3_1.SetPtEtaPhiE(
                klf_bhad3_pt.at(klf_highest_prob_index),
                klf_bhad3_eta.at(klf_highest_prob_index),
                klf_bhad3_phi.at(klf_highest_prob_index),
                klf_bhad3_e.at(klf_highest_prob_index)
        );

        t3_2.SetPtEtaPhiE(
                klf_lquark5_pt.at(klf_highest_prob_index),
                klf_lquark5_eta.at(klf_highest_prob_index),
                klf_lquark5_phi.at(klf_highest_prob_index),
                klf_lquark5_e.at(klf_highest_prob_index)
        );

        t3_3.SetPtEtaPhiE(
                klf_lquark6_pt.at(klf_highest_prob_index),
                klf_lquark6_eta.at(klf_highest_prob_index),
                klf_lquark6_phi.at(klf_highest_prob_index),
                klf_lquark6_e.at(klf_highest_prob_index)
        );

        t3 = t3_1 + t3_2 + t3_3;

        TLorentzVector t4, t4_1, t4_2, t4_3 = {};
        t4_1.SetPtEtaPhiE(
               klf_blep_pt.at(klf_highest_prob_index),
               klf_blep_eta.at(klf_highest_prob_index),
               klf_blep_phi.at(klf_highest_prob_index),
               klf_blep_e.at(klf_highest_prob_index)
        );

        t4_2.SetPtEtaPhiE(
                klf_lepton_pt.at(klf_highest_prob_index),
                klf_lepton_eta.at(klf_highest_prob_index),
                klf_lepton_phi.at(klf_highest_prob_index),
                klf_lepton_e.at(klf_highest_prob_index)
        );

        t4_3.SetPtEtaPhiE(
                klf_neutrino_pt.at(klf_highest_prob_index),
                klf_neutrino_eta.at(klf_highest_prob_index),
                klf_neutrino_phi.at(klf_highest_prob_index),
                klf_neutrino_e.at(klf_highest_prob_index)
        );

        t4 = t4_1 + t4_2 + t4_3;

        // Leptonic top not available at reco due to MET.
        std::vector<TLorentzVector> klfTops = {t1, t2, t3, t4};
        for (size_t ind =0 ; ind<4; ind++) {
            klf_tops_pt.at(ind) = klfTops.at(ind).Pt();
            klf_tops_eta.at(ind) = klfTops.at(ind).Eta();
            klf_tops_phi.at(ind) = klfTops.at(ind).Phi();
            klf_tops_e.at(ind) = klfTops.at(ind).E();
        }


        // What do I have:
        // - jet reco
        // - first and parentghosts to jet reco
        // - truth tops (and jets)
        // - klf correct permutation and tops

        // I want to match:
        // - Reco to KLF (and thus klf to first and parentghost)
        //      - This is given by the model particles and the jet_index variable.
        // - KLF to truth tops
        //      -> make a truth four top collection and add references between the klf and the truth top collection just as for the reco case
        // - firstghost to truth tops to show how bad we are.




        for (size_t klfInd = 0; klfInd < klf_KLFtriplets.size(); klfInd++) {
            auto &triplet = klf_KLFtriplets.at(klfInd);
            if (triplet.size() == 3) {
                klf_tops_recoJetInds_0.at(klfInd) = triplet.at(0);
                klf_tops_recoJetInds_1.at(klfInd) = triplet.at(1);
                klf_tops_recoJetInds_2.at(klfInd) = triplet.at(2);
            
                TLorentzVector tReco, tReco1, tReco2, tReco3 = {};
                tReco1.SetPtEtaPhiE(
                        jet_pt->at(triplet.at(0)),
                        jet_eta->at(triplet.at(0)),
                        jet_phi->at(triplet.at(0)),
                        jet_e->at(triplet.at(0))               
                        );
                tReco2.SetPtEtaPhiE(
                        jet_pt->at(triplet.at(1)),
                        jet_eta->at(triplet.at(1)),
                        jet_phi->at(triplet.at(1)),
                        jet_e->at(triplet.at(1))
                );
                tReco3.SetPtEtaPhiE(
                        jet_pt->at(triplet.at(2)),
                        jet_eta->at(triplet.at(2)),
                        jet_phi->at(triplet.at(2)),
                        jet_e->at(triplet.at(2))
                );

                tReco = tReco1 + tReco2 + tReco3;
                
                klf_reco_tops_pt.at(klfInd) = tReco.Pt();
                klf_reco_tops_eta.at(klfInd) = tReco.Eta();
                klf_reco_tops_phi.at(klfInd) = tReco.Phi();
                klf_reco_tops_e.at(klfInd) = tReco.E();

                TLorentzVector tFirstghost, tFirstghost1, tFirstghost2, tFirstghost3 = {};
                tFirstghost1.SetPtEtaPhiE(
                        jet_firstghost_pt->at(triplet.at(0)),
                        jet_firstghost_eta->at(triplet.at(0)),
                        jet_firstghost_phi->at(triplet.at(0)),
                        jet_firstghost_e->at(triplet.at(0))
                );
                tFirstghost2.SetPtEtaPhiE(
                        jet_firstghost_pt->at(triplet.at(1)),
                        jet_firstghost_eta->at(triplet.at(1)),
                        jet_firstghost_phi->at(triplet.at(1)),
                        jet_firstghost_e->at(triplet.at(1))
                );
                tFirstghost3.SetPtEtaPhiE(
                        jet_firstghost_pt->at(triplet.at(2)),
                        jet_firstghost_eta->at(triplet.at(2)),
                        jet_firstghost_phi->at(triplet.at(2)),
                        jet_firstghost_e->at(triplet.at(2))
                );

                tFirstghost = tFirstghost1 + tFirstghost2 + tFirstghost3;

                klf_firstghost_tops_pt.at(klfInd) = tFirstghost.Pt();
                klf_firstghost_tops_eta.at(klfInd) = tFirstghost.Eta();
                klf_firstghost_tops_phi.at(klfInd) = tFirstghost.Phi();
                klf_firstghost_tops_e.at(klfInd) = tFirstghost.E();

                // no leptonic reco top because of met. Also true for true.. don't have true met for whatever reason :(.

            } else if (triplet.size() == 1) {
                klf_tops_recoJetInds_0.at(klfInd) =  triplet.at(0);

            }

            // go through KLFtruth map again and find out which barcode belongs to klfInd by comparing the triplets again. Then used the barcode to find the truth top index in truth_barcode which can then be used to build the truth top.
//            bool foundMatch = false;
            for (auto& pair : klf_truthKLFMap) {
                if (triplet == pair.second) { // we have a match
//                    foundMatch = true;
                    auto searchResult = std::find_if(truth_barcode->begin(), truth_barcode->end(), [&pair] (int bar){
                        return pair.first == bar;
                    });
                    int truthInd = -1;
                    if (searchResult != truth_barcode->end()) {
                        truthInd = std::distance(truth_barcode->begin(), searchResult);
                    }

                    if (truthInd == -1 || truthInd > truth_pt->size()) {
                        std::cerr << "Did not find barcode in barcodes. This should never happen. ind: " << truthInd << " for barcode: " << pair.first << std::endl;
                        continue;
                    }

                    TLorentzVector tTruth = {};
                    tTruth.SetPtEtaPhiM(
                            truth_pt->at(truthInd),
                            truth_eta->at(truthInd),
                            truth_phi->at(truthInd),
                            truth_m->at(truthInd)
                            );

                    klf_matchedTruth_tops_pt.at(klfInd) = tTruth.Pt();
                    klf_matchedTruth_tops_eta.at(klfInd) = tTruth.Eta();
                    klf_matchedTruth_tops_phi.at(klfInd) = tTruth.Phi();
                    klf_matchedTruth_tops_e.at(klfInd) = tTruth.E();

                    klf_tops_matchesTruthTop.at(klfInd) = static_cast<char>(true);
                    klf_tops_matchedTruthTopOrigInd.at(klfInd) = truthInd; // I'm an idiot and should have base everything around that index instead of barcode.
                }
            }

        }

        // Go through all truth tops and match them to firstghosts
        int l = 0;
        for (auto& pair : truth_tripletMap) {
            if (pair.second.size() != 3 || l < 0 || l > 2) continue;

            TLorentzVector ghost;
            int jetc = 0;
            for (int jeti : pair.second) {
                TLorentzVector ghostFS = {};
                ghostFS.SetPtEtaPhiE(
                        jet_firstghost_pt->at(jeti),
                        jet_firstghost_eta->at(jeti),
                        jet_firstghost_phi->at(jeti),
                        jet_firstghost_e->at(jeti)
                        );
                if (jetc == 0) ghost = ghostFS;
                else ghost += ghostFS;

                jetc++;
            }

            truth_firstghost_tops_pt.at(l) = ghost.Pt();
            truth_firstghost_tops_eta.at(l) = ghost.Eta();
            truth_firstghost_tops_phi.at(l) = ghost.Phi();
            truth_firstghost_tops_e.at(l) = ghost.E();

            if ( truthTopMap.find(pair.first) != truthTopMap.end() ) {

                truth_hadronic_tops_pt.at(l) = truthTopMap[pair.first].Pt();
                truth_hadronic_tops_eta.at(l) = truthTopMap[pair.first].Eta();
                truth_hadronic_tops_phi.at(l) = truthTopMap[pair.first].Phi();
                truth_hadronic_tops_e.at(l) = truthTopMap[pair.first].E();
            }
            
            l++;



        }




        eventInd++;
        outTree.Fill();

    }

    // Go to the output file and write it.
    outFile.cd();
    std::cout << std::endl << "Writing into output root file: " << outFilePath << std::endl << std::endl;
    cutFlow.Write();
    cutFlowWeighted.Write();
    outTree.Write();

    // Close both input and output ROOT files.
    outFile.Close();
//    inFile->Close();
//    chain->Close();

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    std::string config_path = "./config.json";
    if (argc == 2)
        config_path = argv[1];

    return test(config_path);
    // return EXIT::SUCCES;
}



