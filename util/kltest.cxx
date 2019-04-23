// c++ includes
#include <iostream>
#include <string>
#include <vector>

// KLFitter includes
#include "KLFitter/DetectorAtlas_8TeV.h"

#include "KLFitter/Fitter.h"
#include "KLFitter/LikelihoodTopLeptonJets.h"
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

  KLFitter::LikelihoodTopLeptonJets likelihood{};

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
      TFile::Open("/home/lennart/Data/klftest/SM4t-212560/1LOS/mc16d/1los_ttbar/ttbar_nonAllHad.root", "READ");

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

  unsigned int nEventsMax = 1000;
  unsigned int eventInd = 0;

  /* Out section */
  TFile outFile("./output.root", "RECREATE");
  TTree outTree("nominal", "nominal");

  std::vector<float> klf_bhad_pt;
  std::vector<float> klf_bhad_eta;
  std::vector<float> klf_bhad_phi;
  std::vector<float> klf_bhad_e;
  std::vector<unsigned int> klf_bhad_jet_index;
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
  outTree.Branch("klf_bhad_pt", &klf_bhad_pt);
  outTree.Branch("klf_bhad_eta", &klf_bhad_eta);
  outTree.Branch("klf_bhad_phi", &klf_bhad_phi);
  outTree.Branch("klf_bhad_e", &klf_bhad_e);
  outTree.Branch("klf_bhad_jet_index", &klf_bhad_jet_index);
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



  while (reader.Next()) {
    if (eventInd >= nEventsMax) break;

    if ((eventInd % 500) == 0) printf("[%10i | %10i]\n", eventInd, nEventsMax);

    if (mu_pt->size() > 0 && el_pt->size() > 0) continue;

    // Clear all vectors for of the variables.
    klf_bhad_pt.clear();
    klf_bhad_eta.clear();
    klf_bhad_phi.clear();
    klf_bhad_e.clear();
    klf_bhad_jet_index.clear();
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

    KLFitter::Particles particles{};

    // if (jet_pt->size() < 4) continue;

    TLorentzVector lep;
    unsigned int pdgId;

    if (mu_pt->size() == 0 && el_pt->size() > 0) {
      pdgId = 11;
    } else if (mu_pt->size() > 0 && el_pt->size() == 0) {
      pdgId = 13;
    } else {
      pdgId = mu_pt->at(0) > el_pt->at(0) ? 13 : 11;
    }

    switch (pdgId) {
      case 11:
        lep.SetPtEtaPhiE(el_pt->at(0)/1e3, el_eta->at(0), el_phi->at(0), el_e->at(0)/1e3);
        if (abs(lep.Eta()) > 1.37 && abs(lep.Eta()) < 1.52) continue;  // crack region
        likelihood.SetLeptonType(KLFitter::LikelihoodTopLeptonJets::kElectron);
        particles.AddParticle(lep, lep.Eta(), KLFitter::Particles::kElectron);
        break;

      case 13:
        lep.SetPtEtaPhiE(mu_pt->at(0)/1e3, mu_eta->at(0), mu_phi->at(0), mu_e->at(0)/1e3);
        likelihood.SetLeptonType(KLFitter::LikelihoodTopLeptonJets::kMuon);
        particles.AddParticle(lep, lep.Eta(), KLFitter::Particles::kMuon);
        break;

      default:
        continue;
    }

    for (unsigned int ijet = 0; ijet < 4; ijet++) {
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

    // Add MET information
    const double met_ex = *met_met * cos(*met_phi)/1e3;
    const double met_ey = *met_met * sin(*met_phi)/1e3;
    if (!fitter.SetET_miss_XY_SumET(met_ex, met_ey, *met_met)) {
      std::cerr << "ERROR: Failed to add MET to fitter. Aborting." << std::endl;
      return 1;
    }

    // Loop over all permutations.
    const int nperm = fitter.Permutations()->NPermutations();

    bool isFirst = true;

    for (int iperm = 0; iperm < nperm; iperm++) {
      // Do the fitting magic.
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
      float bhad_pt = modelParticles->Parton(0)->Pt();
      float bhad_eta = modelParticles->Parton(0)->Eta();
      float bhad_phi = modelParticles->Parton(0)->Phi();
      float bhad_e = modelParticles->Parton(0)->E();
      unsigned int bhad_index = (*permutedParticles)->JetIndex(0);

      // Leptonic b quark.
      float blep_pt = modelParticles->Parton(1)->Pt();
      float blep_eta = modelParticles->Parton(1)->Eta();
      float blep_phi = modelParticles->Parton(1)->Phi();
      float blep_e = modelParticles->Parton(1)->E();
      unsigned int blep_index = (*permutedParticles)->JetIndex(1);

      // Light quark 1.
      float lquark1_pt = modelParticles->Parton(2)->Pt();
      float lquark1_eta = modelParticles->Parton(2)->Eta();
      float lquark1_phi = modelParticles->Parton(2)->Phi();
      float lquark1_e = modelParticles->Parton(2)->E();
      unsigned int lquark1_index = (*permutedParticles)->JetIndex(2);

      // Light quark 2.
      float lquark2_pt = modelParticles->Parton(3)->Pt();
      float lquark2_eta = modelParticles->Parton(3)->Eta();
      float lquark2_phi = modelParticles->Parton(3)->Phi();
      float lquark2_e = modelParticles->Parton(3)->E();
      unsigned int lquark2_index = (*permutedParticles)->JetIndex(3);

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

      // if (isFirst) {
      //   printf("----------------------------------------------------------------------------------------------\n");
      //   printf("----------------------------------------Permutation %2i----------------------------------------\n",
      //          iperm);
      //   printf("----------------------------------------------------------------------------------------------\n");
      //   printf("                  | hadronic b quark | leptonic b quark  |  light quark 1   |  light quark 2  |\n");
      //   printf("Jet index         | %16i | %17i | %16i | %15i |\n", bhad_index, blep_index, lquark1_index,
      //          lquark2_index);
      //   printf("----------------------------------------------------------------------------------------------\n");
      //   printf("                  | hadronic b quark | leptonic b quark  |  light quark 1   |  light quark 2  |\n");
      //   printf("Output Energies   | %16.2f | %17.2f | %16.2f | %15.2f |\n", bhad_e, blep_e, lquark1_e, lquark2_e);
      //   printf("----------------------------------------------------------------------------------------------\n");
      //   printf("                  | lepton energy    | neutrino pz       | loglikelihood    |  probability    |\n");
      //   printf("Other values      | %16.2f | %17.2f | %16.2f | %15.2e |\n", lepton_e, neutrino_pt, likelihood,
      //          event_probability);
      //   printf("----------------------------------------------------------------------------------------------\n");
      //   printf("                  | Minuit Not Conv. | Fit Aborted: NaN  | >=1 Par at Limit | Invalid TF@Conv.|\n");
      //   printf("Status Code       | %16i | %17i | %16i | %15i |\n", MinuitDidNotConverge, FitAbortedDueToNaN,
      //          AtLeastOneFitParameterAtItsLimit, InvalidTransferFunctionAtConvergence);
      // }

      // Fill the vectors with the output variables. Note: it's
      // better/safer to store booleans as chars in ROOT files.
      klf_fit_minuit_did_not_converge.emplace_back(static_cast<char>(MinuitDidNotConverge));
      klf_fit_aborted_to_nan.emplace_back(static_cast<char>(FitAbortedDueToNaN));
      klf_fit_parameter_at_limit.emplace_back(static_cast<char>(AtLeastOneFitParameterAtItsLimit));
      klf_fit_invalid_transfer_function.emplace_back(static_cast<char>(InvalidTransferFunctionAtConvergence));

      klf_bhad_pt.emplace_back(bhad_pt);
      klf_bhad_eta.emplace_back(bhad_eta);
      klf_bhad_phi.emplace_back(bhad_phi);
      klf_bhad_e.emplace_back(bhad_e);
      klf_bhad_jet_index.emplace_back(bhad_index);
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

    }
    if (isFirst) isFirst = false;

    klf_highest_prob_index = std::distance( klf_event_probability.begin(), max_element(klf_event_probability.begin(), klf_event_probability.end()));

    outTree.Fill();
    eventInd++;
  }

  // Go to the output file and write it.
  outFile.cd();
  std::cout << std::endl << "Writing into output root file: " << "top-ljets-output.root" << std::endl << std::endl;
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
