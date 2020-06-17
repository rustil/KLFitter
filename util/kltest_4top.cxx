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

struct topChild {
    int recoJetIndex_iterativeChildReco;
    int recoJetIndex_globalChildReco;
    int recoJetIndex_globalChildFirstGhost;
    int recoJetIndex_PDGglobalChildFirstGhost;
    int recoJetIndex_PDGiterativeChildFirstGhost;
    int childPDGId;
    int ghostPDGId;
    int PDGGhostPDGId;
    float dR_iterativeChildReco;
    float dR_globalChildReco;
    float dR_globalChildFirstGhost;
    float dR_PDGglobalChildFirstGhost;
    float dR_PDGiterativeChildFirstGhost;
    TLorentzVector childVector;
//    TLorentzVector truthTop;
    bool isHadronic;
};

struct top {
    float truthBPt;
    int truthBChildInd;
    float recoBPt;
    int recoBInd;
    bool isHadronic;
    TLorentzVector truthTop;
    std::vector<topChild> children;
    std::vector<int> recoChildrenIndices;
    TLorentzVector SumChildren;
    TLorentzVector SumIterativePDGMatchedFirstGhosts;
    float dR_top_SumChildren;
    float dR_top_SumIterative;
    float dR_SumChildren_SumIterative;

};

void WMatchs(std::vector<std::tuple<std::tuple<int,int>, TLorentzVector>>& Ws, std::vector<std::tuple<int, TLorentzVector>> Wvecs) {
    while (Wvecs.size() > 1) {
        auto back = Wvecs.back();
        Wvecs.pop_back();

        for (auto& W : Wvecs) {
            if (std::get<1>(W).M() == std::get<1>(back).M() && !(std::get<0>(W) == std::get<0>(back))) {
                Ws.emplace_back(
                        std::make_tuple(std::make_tuple(std::get<0>(W), std::get<0>(back)), std::get<1>(W))
                );
                std::remove(Wvecs.begin(), Wvecs.end(), W);
                break;
            }
        }
    }
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

std::tuple<int,float> deltaRMatchWPDGId(TLorentzVector& child, int childPDGId, std::vector<int>& jet_PDGIDs,  std::vector<float>& jet_pt, std::vector<float>& jet_eta, std::vector<float>& jet_phi, std::vector<float>& jet_e) {
    float minDR = 999.;
    size_t minDR_Ind;
    for (size_t i = 0; i < jet_pt.size(); i++) {

        if (childPDGId != jet_PDGIDs.at(i)) continue;

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

    int nPartonsMatched = -1;
    if (conf.find("nPartonsMatched") != conf.end())
     nPartonsMatched = conf["nPartonsMatched"].get<int>();

    std::cout << "nPartonsMatched: " << nPartonsMatched << std::endl;

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


    std::vector<std::string> fileNames = conf["inputFiles"].get<std::vector<std::string>>();

    auto chain = new TChain("nominal_Loose");
//    std::vector<std::string> fileNames = {
//            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000002.output.root",
//            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000001.output.root",
//            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000004.output.root",
//            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000006.output.root",
//            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000005.output.root",
//            "/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000003.output.root"
//    };

    for (auto fn : fileNames) chain->AddFile(fn.c_str());

    // auto inFile = TFile::Open("/home/lennart/Data/klftest/SM4t-212560/SSML/mc16e/inclusive/ttbar_nonAllHad.root",
    // "READ");
//    auto inFile =
//            TFile::Open("/home/lennart/Data/klftest/SM4t-212560/1LOS/mc16d/1los_4top/4topNLO.root", "READ");


//auto inFile =
//            TFile::Open("/home/lennart/Data/klftest/wTruth/user.lrustige.412043.aMcAtNloPythia8EvtGen.DAOD_TOPQ1.e7101_a875_r10724_p3629.SM4t-212560_mc16e_KLFitter_v1_output_root/user.lrustige.18026714._000001.output.root", "READ");


    TTreeReader reader(chain);
//    TTreeReader reader("nominal_Loose", inFile);

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

    bool isTtBar = false;
    if (conf.find("isTtBar") != conf.end()){
        isTtBar = conf["isTtBar"].get<bool>();
    }





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

//    std::vector<float> cp_truth_top_pt;
//    std::vector<float> cp_truth_top_eta;
//    std::vector<float> cp_truth_top_phi;
//    std::vector<int> cp_truth_top_type;

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

    std::vector<float> globalChildRecoDR;
    std::vector<float> iterativeChildRecoDR;
    std::vector<int> globalChildRecoIndMult;
    bool iterativeChildRecoDR_lt_01;

    std::vector<float> globalChildFirstGhostDR;
    std::vector<float> globalPDGChildFirstGhostDR;
    std::vector<int> globalPDGChildFirstGhostIndMult;

    std::vector<float> iterativePDGChildFirstGhostDRs;


    outTree.Branch("evt_globalChildRecoDR", &globalChildRecoDR);
    outTree.Branch("evt_iterativeChildRecoDR", &iterativeChildRecoDR);
    outTree.Branch("evt_globalChildRecoIndMult", &globalChildRecoIndMult);
    outTree.Branch("evt_iterativeChildRecoDR_lt_01", &iterativeChildRecoDR_lt_01);


    outTree.Branch("evt_globalChildFirstGhostDR", &globalChildFirstGhostDR);
    outTree.Branch("evt_globalPDGChildFirstGhostDR", &globalPDGChildFirstGhostDR);
    outTree.Branch("evt_globalPDGChildFirstGhostIndMult", &globalPDGChildFirstGhostIndMult);

    outTree.Branch("evt_iterativePDGChildFirstGhostDR", &iterativePDGChildFirstGhostDRs);

    time_t rawtime = time(NULL);
    struct tm * timeinfo;
    char buffer [80];



    while (reader.Next()) {

        globalChildRecoDR.clear();
        globalChildRecoIndMult.clear();
        iterativeChildRecoDR.clear();

        globalChildFirstGhostDR.clear();
        globalPDGChildFirstGhostDR.clear();
        globalPDGChildFirstGhostIndMult.clear();

        iterativePDGChildFirstGhostDRs.clear();

        KLFitter::Particles particles{};

        unsigned int nLep = mu_pt->size() + el_pt->size();
        if (eventInd >= nEventsMax) break;
        else if (nLep != 1 || jet_pt->size() != 10 || std::count_if(jet_isbtagged_MV2c10_77->begin(), jet_isbtagged_MV2c10_77->end(), [](char flag) {return flag;}) != 4) continue;


            nPartonsMatcheable = 0;

            for (auto pdgId : *jet_parentghost_pdgId) {
                if (abs(pdgId) == 6 || abs(pdgId) == 24) ++nPartonsMatcheable;
            }

        if (nPartonsMatched != -1 && nPartonsMatcheable != nPartonsMatched) continue;

        if ((eventInd % 1) == 0) {
            time(&rawtime);
            timeinfo = localtime(&rawtime);
            strftime(buffer, 80, "%FT%T", timeinfo);
            printf("%10s: [%10i | %10i]\n", buffer, eventInd, nEventsMax);
        }

//        else if (nLep == 0 || nLep > 1 || jet_pt->size() != 10 || std::count_if(jet_mv2c10->begin(), jet_mv2c10->end(), [](float x){return x > 0.64;}) != 4) continue;

        //TODO: not exactly optimal here. Match truth children from top with truth matched to reco to find correct permutation.
        // Compare with permutation given by klf. In principle I'm only interested in the correct association of qqb and
        // lvb to "a" top (maybe of the correct charge) and I don't care which top in a first step. After that I could ask
        // for the right top as well (which should be always the case anyhow?). KLFitter doesn't know the charge of the tops (I think).

//        std::vector<TLorentzVector> jets(10);
//        for(std::size_t k = 0; k < jets.size(); ++k) {
//            jets.at(k).SetPtEtaPhiE(jet_pt->at(k), jet_eta->at(k), jet_phi->at(k), jet_e->at(k));
//        }

// TODO: use true bjets to define order of tops and then permutation accordingly

//        std::vector<TLorentzVector> truthBJets(4);
//        std::vector<int> truthBJetIndices(4);
//        std::size_t j = 0;
//        for(std::size_t k = 0; k < jet_firstghost_pt->size(); ++k) {
////            if (jet_firstghost_pt->size() < 10) break;
//
//            if (TMath::Abs(jet_firstghost_pdgId->at(k)) == 5) {
//                truthBJets.at(j).SetPtEtaPhiE(jet_firstghost_pt->at(k), jet_firstghost_eta->at(k), jet_firstghost_phi->at(k), jet_firstghost_e->at(k));
//                truthBJetIndices.at(j) = k;
//                j++;
//            }
//        }

//        std::sort(truthBJets.begin(), truthBJets.end(), [](TLorentzVector v1, TLorentzVector v2) {v1.Pt() > v2.Pt();}); // They should already be pt ordered as the jet collection is
//        std::sort(truthBJets.begin(), truthBJets.end(), [](TLorentzVector v1, TLorentzVector v2) {v1.Pt() > v2.Pt();});
//        std::for_each(truthBJets.begin(), truthBJets.end(), [](TLorentzVector v1) {std::cout << v1.Pt() << ", ";}); // Indeed they are.

        // Build fourtop vector in top1 tbar 1 top 2 tbar 2

        std::vector<topChild> children;


        std::vector<TLorentzVector> truthTops (4);
        truthTops.at(0).SetPtEtaPhiM(*truth_top1_pt, *truth_top1_eta, *truth_top1_phi, 172.5*1e3);
        truthTops.at(1).SetPtEtaPhiM(*truth_top2_pt, *truth_top2_eta, *truth_top2_phi, 172.5*1e3);
        truthTops.at(2).SetPtEtaPhiM(*truth_tbar1_pt, *truth_tbar1_eta, *truth_tbar1_phi, 172.5*1e3);
        truthTops.at(3).SetPtEtaPhiM(*truth_tbar2_pt, *truth_tbar2_eta, *truth_tbar2_phi, 172.5*1e3);

        std::vector<bool> truthTops_areHad {
            *truth_top1_isHad == 1,
            *truth_top2_isHad == 1,
            *truth_tbar1_isHad == 1,
            *truth_tbar2_isHad == 1
        };

        TLorentzVector child0_0_Vector = {};
        child0_0_Vector.SetPtEtaPhiE(*truth_top1_child0_pt, *truth_top1_child0_eta, *truth_top1_child0_phi, *truth_top1_child0_e);
        TLorentzVector child0_1_Vector = {};
        child0_1_Vector.SetPtEtaPhiE(*truth_top1_child1_pt, *truth_top1_child1_eta, *truth_top1_child1_phi, *truth_top1_child1_e);
        TLorentzVector child0_2_Vector = {};
        child0_2_Vector.SetPtEtaPhiE(*truth_top1_child2_pt, *truth_top1_child2_eta, *truth_top1_child2_phi, *truth_top1_child2_e);

        TLorentzVector child1_0_Vector = {};
        child1_0_Vector.SetPtEtaPhiE(*truth_top2_child0_pt, *truth_top2_child0_eta, *truth_top2_child0_phi, *truth_top2_child0_e);
        TLorentzVector child1_1_Vector = {};
        child1_1_Vector.SetPtEtaPhiE(*truth_top2_child1_pt, *truth_top2_child1_eta, *truth_top2_child1_phi, *truth_top2_child1_e);
        TLorentzVector child1_2_Vector = {};
        child1_2_Vector.SetPtEtaPhiE(*truth_top2_child2_pt, *truth_top2_child2_eta, *truth_top2_child2_phi, *truth_top2_child2_e);

        TLorentzVector child2_0_Vector = {};
        child2_0_Vector.SetPtEtaPhiE(*truth_tbar1_child0_pt, *truth_tbar1_child0_eta, *truth_tbar1_child0_phi, *truth_tbar1_child0_e);
        TLorentzVector child2_1_Vector = {};
        child2_1_Vector.SetPtEtaPhiE(*truth_tbar1_child1_pt, *truth_tbar1_child1_eta, *truth_tbar1_child1_phi, *truth_tbar1_child1_e);
        TLorentzVector child2_2_Vector = {};
        child2_2_Vector.SetPtEtaPhiE(*truth_tbar1_child2_pt, *truth_tbar1_child2_eta, *truth_tbar1_child2_phi, *truth_tbar1_child2_e);

        TLorentzVector child3_0_Vector = {};
        child3_0_Vector.SetPtEtaPhiE(*truth_tbar2_child0_pt, *truth_tbar2_child0_eta, *truth_tbar2_child0_phi, *truth_tbar2_child0_e);
        TLorentzVector child3_1_Vector = {};
        child3_1_Vector.SetPtEtaPhiE(*truth_tbar2_child1_pt, *truth_tbar2_child1_eta, *truth_tbar2_child1_phi, *truth_tbar2_child1_e);
        TLorentzVector child3_2_Vector = {};
        child3_2_Vector.SetPtEtaPhiE(*truth_tbar2_child2_pt, *truth_tbar2_child2_eta, *truth_tbar2_child2_phi, *truth_tbar2_child2_e);

        std::vector<TLorentzVector> childrenVectors = {
                child0_0_Vector, child0_1_Vector, child0_2_Vector,
                child1_0_Vector, child1_1_Vector, child1_2_Vector,
                child2_0_Vector, child2_1_Vector, child2_2_Vector,
                child3_0_Vector, child3_1_Vector, child3_2_Vector,
        };

        std::vector<int> childrenPdgs = {
            *truth_top1_child0_pdgid, *truth_top1_child1_pdgid, *truth_top1_child2_pdgid,
            *truth_top2_child0_pdgid, *truth_top2_child1_pdgid, *truth_top2_child2_pdgid,
            *truth_tbar1_child0_pdgid, *truth_tbar1_child1_pdgid, *truth_tbar1_child2_pdgid,
            *truth_tbar2_child0_pdgid, *truth_tbar2_child1_pdgid, *truth_tbar2_child2_pdgid,
        };

        std::vector<top> enhancedTops(4);


        std::vector<int> globalChildRecoInds = {};
        std::vector<int> globalPDGChildFirstGhostInds = {};
        std::vector<float> mutable_jet_pt = *jet_pt;
        std::vector<float> mutable_jet_eta = *jet_eta;
        std::vector<float> mutable_jet_phi = *jet_phi;
        std::vector<float> mutable_jet_e = *jet_e;

        std::vector<float> mutable_jet_firstghost_pt = *jet_firstghost_pt;
        std::vector<float> mutable_jet_firstghost_eta = *jet_firstghost_eta;
        std::vector<float> mutable_jet_firstghost_phi = *jet_firstghost_phi;
        std::vector<float> mutable_jet_firstghost_e = *jet_firstghost_e;
        std::vector<int> mutable_jet_firstghost_pdgId = *jet_firstghost_pdgId;

        std::vector<std::tuple<int, float, float, TLorentzVector>> parentGhostVectors = {};
        for (size_t k =0 ; k< jet_parentghost_pt->size(); k++) {
            TLorentzVector v = {jet_parentghost_pt->at(k),
                                jet_parentghost_eta->at(k),
                                jet_parentghost_phi->at(k),
                                jet_parentghost_e->at(k),
                                };
            auto tup = std::make_tuple(jet_parentghost_pdgId->at(k), v.Pt(), v.M(), v);
            parentGhostVectors.emplace_back(tup);
        }


        std::vector<std::tuple<int, TLorentzVector>> Wvecs = {};
        for (size_t k =0 ; k< jet_parentghost_pt->size(); k++) {
            if (TMath::Abs(jet_parentghost_pdgId->at(k)) == 24) {
            TLorentzVector Wv = {jet_parentghost_pt->at(k),
                                jet_parentghost_eta->at(k),
                                jet_parentghost_phi->at(k),
                                jet_parentghost_e->at(k)
            };
            auto tup = std::make_tuple(k, Wv);
            Wvecs.emplace_back(tup);
            }
        }

//        std::vector<std::tuple<int, TLorentzVector>> mutable_Wvecs = Wvecs;
        std::vector<std::tuple<std::tuple<int,int>, TLorentzVector>> res = {};
        WMatchs(res, Wvecs);



//        std::vector<std::tuple<int, float, float, TLorentzVector>> condensed_Wvecs = {};
//        for (auto W : Wvecs) {
//            if (mutable_Wvecs.empty()) break;
//
//            for (auto mutW : mutable_Wvecs) {
//                if (std::get<0>(W) == std::get<0>(mutW)) continue;
//                else if (std::get<1>(W) == std::get<1>(mutW))
//                    std::tuple<int, int> indices = std::make_tuple(std::get<0>(W), std::get<1>(mutW));
//                    condensed_Wvecs.emplace_back()
//
//
//            }
//
//        }


// lets test something..

// AWESOME AS FUCK this seems to work. It seems I can use the parentghost top barcode to match reoc jets to truth objects where the truth objects are then in the truth container.
// This way I would need to check which jets have the same top barcode and should get my triplets. Let's do this.
    std::map<int, std::vector<int>> tripletMap = {};
        for (int k =0 ; k < jet_parentghost_top_barcode->size(); k++) {
            auto key = jet_parentghost_top_barcode->at(k);
            auto search = tripletMap.find(key);
            if (search == tripletMap.end()) {
                tripletMap[key] = {k};
            } else {
                tripletMap[key].emplace_back(k);
            }
        }

    std::cout << "help";
//        for (size_t k =0 ; k< jet_parentghost_top_barcode->size(); k++) {
//            //match parentghost to truth
//            int truthInd = -1;
//            for( size_t j=0; j< truth_barcode->size() ; j++) {
//                if (jet_parentghost_top_barcode->at(k) == truth_barcode->at(j)) {
//                    truthInd = j;
//                    break;
//                }
//            }
//
//               if (truthInd < 0 ) {
//                   std::cerr << "not good." << std::endl;
//                   return -1;
//               }
//
//               if (TMath::Abs(truth_pdgid->at(truthInd)) !=6) {
//                   std::cerr << "not good." << std::endl;
//                   continue;
//               }
//
//               TLorentzVector truthTop = {};
//               truthTop.SetPtEtaPhiM(truth_pt->at(truthInd),
//                                     truth_eta->at(truthInd),
//                                     truth_phi->at(truthInd),
//                                     truth_m->at(truthInd)
//                       );
//
//               std::cout<< "Real Truth: " << truthTop.M() << " compared to: ";
//               for (auto& top: truthTops) {
//                   std::cout << truthTop.DeltaR(top) << ", ";
//               }
//               std::cout << std::endl;
//
//
//            }



//        std::vector<std::tuple<int, float, float, TLorentzVector>> bvecs = {};
//        std::vector<std::tuple<int, int, int, float, TLorentzVector, TLorentzVector>> pottops = {};
//        for (size_t k =0 ; k< jet_parentghost_pt->size(); k++) {
//            if (TMath::Abs(jet_parentghost_pdgId->at(k)) == 6) {
//                TLorentzVector bv = {jet_parentghost_pt->at(k),
//                                     jet_parentghost_eta->at(k),
//                                     jet_parentghost_phi->at(k),
//                                     jet_parentghost_e->at(k)
//                };
//                auto tup = std::make_tuple(k, bv.Pt(), bv.M(), bv);
//                bvecs.emplace_back(tup);
//
//                TLorentzVector firstB = {jet_firstghost_pt->at(k),
//                                         jet_firstghost_eta->at(k),
//                                         jet_firstghost_phi->at(k),
//                                         jet_firstghost_e->at(k)
//                };
//
////                std::cout << bv.E() << ", " << bv.Pt() << ", " << bv.M() << std::endl;
//
//                std::cout << "firstB: " << firstB.E() << " childrenB: " <<
//                ", " << child0_0_Vector.E() <<
//                ", " << child1_0_Vector.E() <<
//                ", " << child2_0_Vector.E() <<
//                ", " << child3_0_Vector.E() << std::endl;
//
////                std::cout << "parent: " << bv.M() << ", Children sum M: " << childrenVectors
//
//                for (size_t j = 0; j < res.size(); j++) {
//                    TLorentzVector potTop = (std::get<1>(res.at(j)) + firstB);
//
////                    std::cout << potTop.E() - bv.E() << ", ";
//                    std::cout << potTop.M() << ", ";
//                }
//                std::cout << std::endl;
//
////                for (size_t z = 0 ; z <truthTops.size();z++) {
////                    std::cout << truthTops.at(z).E() - bv.E() << ", ";
////                }
////                std::cout << std::endl;
//            }
//        }
//
////                    auto potTop = bv;
////
////                    int minInd = -1;
////                    float mindR = 1000;
////                    for (size_t z = 0 ; z <truthTops.size();z++) {
////                        float dR = potTop.DeltaR(truthTops.at(z));
////                        if (dR < mindR) {
////                            mindR = dR;
////                            minInd = z;
////                        }
////                    }
////
//////                    auto bWtInd = std::make_tuple(k, j, minInd);
//////                    auto temp = std::make_tuple(k,j,minInd, mindR, potTop, truthTops.at(minInd));
////                    auto temp = std::make_tuple(k,-1,minInd, mindR, potTop, truthTops.at(minInd));
////                    pottops.emplace_back(temp);
////
//////                    pottops.emplace_back(std::make_tuple(k, j, ));
////                }
////            }
//
//        std::cout << "Well.." << std::endl;
//        for (auto& top: truthTops) {
////            std::cout << top.E() << ", " << top.Pt() << ", " << top.M() << std::endl;
//            std::cout << top.M() << ", ";
//        }
//        std::cout << std::endl << std::endl;
//
//
//
//
//
//
//
        std::cout << "Halt!" << std::endl;



        for (size_t k = 0; k < childrenVectors.size(); k++) {
            int jetInd = -2;
            int globalJetInd = -2;
            float globaldeltaR = 1000;
            float deltaR = 1000;
            float childFirstGhostDR = 1000;
            int childFirstGhostInd = -2;
            float PDGChildFirstGhostDR = 1000;
            int PDGChildFirstGhostInd = -2;
            int iterativePDGChildFirstGhostInd = -2;
            float iterativePDGChildFirstGhostDR = 1000;

            if (TMath::Abs(childrenPdgs.at(k)) < 6) {
                // Rematch truth to reco.. does not look good on purely dR!

                // This is the real minimum..and values are overwritten. Next is the iterative minimum
                auto results = deltaRMatch(childrenVectors.at(k), *jet_pt, *jet_eta, *jet_phi, *jet_e);
                globalJetInd = std::get<0> (results);
                globaldeltaR = std::get<1> (results);
                if (jetInd != -2) globalChildRecoInds.emplace_back(jetInd);
                globalChildRecoDR.emplace_back(deltaR);

                // Iterative Minimum. First come, first serve
                results = deltaRMatch(childrenVectors.at(k), mutable_jet_pt,mutable_jet_eta,mutable_jet_phi,mutable_jet_e);
                jetInd = std::get<0> (results);
                deltaR = std::get<1> (results);
                mutable_jet_pt.at(jetInd) = -1;

                iterativeChildRecoDR.emplace_back(deltaR);


                // Now match child to firstghost. Should be 0. -> but isn't
                results = deltaRMatch(childrenVectors.at(k), *jet_firstghost_pt,*jet_firstghost_eta,*jet_firstghost_phi,*jet_firstghost_e);
                childFirstGhostInd = std::get<0> (results);
                childFirstGhostDR = std::get<1> (results);

                globalChildFirstGhostDR.emplace_back(childFirstGhostDR);

                // Now match child to firstghost with PDG
                results = deltaRMatchWPDGId(childrenVectors.at(k),childrenPdgs.at(k), *jet_firstghost_pdgId,  *jet_firstghost_pt,*jet_firstghost_eta,*jet_firstghost_phi,*jet_firstghost_e);
                PDGChildFirstGhostInd = std::get<0> (results);
                PDGChildFirstGhostDR = std::get<1> (results);

                globalPDGChildFirstGhostDR.emplace_back(childFirstGhostDR);
                globalPDGChildFirstGhostInds.emplace_back(PDGChildFirstGhostInd);


                // Now match child to firstghost with PDG and iteratively
                results = deltaRMatchWPDGId(childrenVectors.at(k),childrenPdgs.at(k), mutable_jet_firstghost_pdgId,  mutable_jet_firstghost_pt,mutable_jet_firstghost_eta,mutable_jet_firstghost_phi,mutable_jet_firstghost_e);

                iterativePDGChildFirstGhostInd = std::get<0> (results);
                iterativePDGChildFirstGhostDR = std::get<1> (results);

                mutable_jet_firstghost_pt.at(iterativePDGChildFirstGhostInd) = -1;
                iterativePDGChildFirstGhostDRs.emplace_back(iterativePDGChildFirstGhostDR);


            }

            topChild child = {
                        jetInd,
                        globalJetInd,
                        childFirstGhostInd,
                        PDGChildFirstGhostInd,
                        iterativePDGChildFirstGhostInd,
                        childrenPdgs.at(k),
                        childFirstGhostInd >= 0 ? jet_firstghost_pdgId->at(childFirstGhostInd) : -99,
                        PDGChildFirstGhostInd >= 0 ? jet_firstghost_pdgId->at(PDGChildFirstGhostInd) : -99,
                        deltaR,
                        globaldeltaR,
                        childFirstGhostDR,
                        PDGChildFirstGhostDR,
                        iterativePDGChildFirstGhostDR,
                        childrenVectors.at(k),
                        truthTops_areHad.at( k / 3)
            };
        children.emplace_back( child);

        }

        // Measure to show just how bad the matching is!
        std::map<int,int> multiplicities;
        std::for_each(globalChildRecoInds.begin(), globalChildRecoInds.end(), [&multiplicities](int val) {multiplicities[val]++;});
        std::for_each(multiplicities.begin(), multiplicities.end(), [&globalChildRecoIndMult](std::pair<int,int> val) {globalChildRecoIndMult.emplace_back(val.second);});

        int count = std::count_if(iterativeChildRecoDR.begin(), iterativeChildRecoDR.end(), [] (float val) {return val >= 0.01;});
        iterativeChildRecoDR_lt_01 = count == 0;

        // Same for global PDG ChildFirstghost
        multiplicities.clear();
        std::for_each(globalPDGChildFirstGhostInds.begin(), globalPDGChildFirstGhostInds.end(), [&multiplicities](int val) {multiplicities[val]++;});
        std::for_each(multiplicities.begin(), multiplicities.end(), [&globalPDGChildFirstGhostIndMult](std::pair<int,int> val) {globalPDGChildFirstGhostIndMult.emplace_back(val.second);});

        for (size_t k = 0; k < enhancedTops.size(); k++) {
            TLorentzVector bChild;
            int bJetRecoInd = -2;
            int childBInd = -2;

            for (size_t j = 0; j < 3; j++) {
                if(TMath::Abs(children.at(j + 3*k).childPDGId) == 5) {
                    bChild = children.at(j + 3*k).childVector;
                    bJetRecoInd = children.at(j + 3*k).recoJetIndex_globalChildFirstGhost;
                    childBInd = j + 3*k;
                    break;
                }
            }

            std::vector<topChild> ownChildren {
                children.at(3*k + 0),
                children.at(3*k + 1),
                children.at(3*k + 2)
            };

            std::vector<int> recoChildIndices {
                    children.at(3*k + 0).recoJetIndex_globalChildFirstGhost,
                    children.at(3*k + 1).recoJetIndex_globalChildFirstGhost,
                    children.at(3*k + 2).recoJetIndex_globalChildFirstGhost
            };

            if (bJetRecoInd < 0) {
                std::cerr << "No Bjet found for top" << std::endl;
                break;
            }

            TLorentzVector SumChildren = children.at(3*k + 0).childVector
                    + children.at(3*k + 1).childVector
                    + children.at(3*k + 2).childVector;

            std::vector<TLorentzVector> iterativePDGMatchedFirstGhosts {};
            for (size_t j = 0; j< 3; j++ ) {
                topChild child = children.at(3 * k + j);
                int ind = child.recoJetIndex_PDGiterativeChildFirstGhost;
                if (ind < 0) {
                    continue;
                } else {

                    TLorentzVector vecci = {
                            jet_firstghost_pt->at(ind),
                            jet_firstghost_eta->at(ind),
                            jet_firstghost_phi->at(ind),
                            jet_firstghost_e->at(ind),
                    };
                    iterativePDGMatchedFirstGhosts.emplace_back(vecci);
                }
            }

            TLorentzVector sumPDGetc = {};
            if (iterativePDGMatchedFirstGhosts.size() == 3) {
                sumPDGetc = iterativePDGMatchedFirstGhosts.at(0)
                            + iterativePDGMatchedFirstGhosts.at(1)
                            + iterativePDGMatchedFirstGhosts.at(2);
            }

            top enhancedTop = {
                    static_cast<float>(bChild.Pt()),
                    childBInd,
                    jet_pt->at(bJetRecoInd),
                    bJetRecoInd,
                    truthTops_areHad.at(k),
                    truthTops.at(k),
                    ownChildren,
                    recoChildIndices,
                    SumChildren,
                    sumPDGetc,
                    static_cast<float>(truthTops.at(k).DeltaR(SumChildren)),
                    static_cast<float>(truthTops.at(k).DeltaR(sumPDGetc)),
                    static_cast<float>(SumChildren.DeltaR(sumPDGetc))
            };

            enhancedTops.at(k) = enhancedTop;

        }



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
        if (!fitter.SetET_miss_XY_SumET(met_ex, met_ey, *met_met)) {
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
        auto jets = highestPerm["light_indices"].get<std::vector<int>>();;
//        std::for_each(highestPerm.begin(), highestPerm.end(), [](int x) {std::cout << x;});
//        std::cout << std::endl;


        eventInd++;
        outTree.Fill();

    }

    // Go to the output file and write it.
    outFile.cd();
    std::cout << std::endl << "Writing into output root file: " << outFilePath << std::endl << std::endl;
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
