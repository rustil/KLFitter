#include <memory>

//
// Created by lennart on 16/04/19.
//

#include "KLFitter/LikelihoodFourTopLeptonJets.h"

#include <algorithm>
#include <iostream>

#include "BAT/BCMath.h"
#include "BAT/BCParameter.h"
#include "KLFitter/DetectorBase.h"
#include "KLFitter/Particles.h"
#include "KLFitter/Permutations.h"
#include "KLFitter/PhysicsConstants.h"
#include "KLFitter/ResolutionBase.h"
#include "TLorentzVector.h"

namespace KLFitter {

    LikelihoodFourTopLeptonJets::LikelihoodFourTopLeptonJets()
            : LikelihoodBase::LikelihoodBase(), m_flag_top_mass_fixed(false), m_flag_get_par_sigmas_from_TFs(false),
              m_et_miss_x(0.), m_et_miss_y(0.), m_et_miss_sum(0.), m_lepton_type(kElectron) {
        // define model particles
        this->DefineModelParticles();

        // define parameters
        this->DefineParameters();
    }

// ---------------------------------------------------------
    LikelihoodFourTopLeptonJets::~LikelihoodFourTopLeptonJets() = default;

// ---------------------------------------------------------
    int LikelihoodFourTopLeptonJets::SetET_miss_XY_SumET(double etx, double ety, double sumet) {
        // set missing ET x and y component and the m_et_miss_sum
        m_et_miss_x = etx;
        m_et_miss_y = ety;
        m_et_miss_sum = sumet;

        // no error
        return 1;
    }

// ---------------------------------------------------------
    void LikelihoodFourTopLeptonJets::RequestResolutionFunctions() {
        (*fDetector)->RequestResolutionType(ResolutionType::EnergyLightJet);
        (*fDetector)->RequestResolutionType(ResolutionType::EnergyBJet);
        (*fDetector)->RequestResolutionType(ResolutionType::EnergyElectron);
        (*fDetector)->RequestResolutionType(ResolutionType::EnergyMuon);
        (*fDetector)->RequestResolutionType(ResolutionType::MissingET);
    }

// ---------------------------------------------------------
    void LikelihoodFourTopLeptonJets::SetLeptonType(LeptonType leptontype) {
        if (leptontype != kElectron && leptontype != kMuon) {
            std::cout << "KLFitter::SetLeptonType(). Warning: lepton type not defined. Set electron as lepton type."
                      << std::endl;
            m_lepton_type = kElectron;
        } else {
            m_lepton_type = leptontype;
        }

        // define model particles
        DefineModelParticles();
    }

// ---------------------------------------------------------
    void LikelihoodFourTopLeptonJets::SetLeptonType(int leptontype) {
        if (leptontype != 1 && leptontype != 2) {
            std::cout << "KLFitter::SetLeptonType(). Warning: lepton type not defined. Set electron as lepton type."
                      << std::endl;
            leptontype = 1;
        }

        if (leptontype == 1) {
            SetLeptonType(kElectron);
        } else if (leptontype == 2) {
            SetLeptonType(kMuon);
        }
    }

// ---------------------------------------------------------
    int LikelihoodFourTopLeptonJets::DefineModelParticles() {
        // create the particles of the model
        fParticlesModel = std::make_unique<Particles>();

        // add model particles
        TLorentzVector dummy{0, 0, 0, 0};
        fParticlesModel->AddParticle(&dummy,
                                     Particles::kParton,  // type
                                     "hadronic b quark 1",  // name
                                     0,                   // index of corresponding particle
                                     Particles::kB);      // b jet (truth)

        fParticlesModel->AddParticle(&dummy,
                                     Particles::kParton,  // type
                                     "hadronic b quark 2",  // name
                                     1,                   // index of corresponding particle
                                     Particles::kB);      // b jet (truth)

        fParticlesModel->AddParticle(&dummy,
                                     Particles::kParton,  // type
                                     "hadronic b quark 3",  // name
                                     2,                   // index of corresponding particle
                                     Particles::kB);      // b jet (truth)

        fParticlesModel->AddParticle(&dummy,
                                     Particles::kParton,
                                     "leptonic b quark",
                                     3,                   // index of corresponding particle
                                     Particles::kB);      // b jet (truth)

        fParticlesModel->AddParticle(&dummy,
                                     Particles::kParton,
                                     "light quark 1",
                                     4,                   // index of corresponding particle
                                     Particles::kLight);  // light jet (truth)

        fParticlesModel->AddParticle(&dummy,
                                     Particles::kParton,
                                     "light quark 2",
                                     5,                   // index of corresponding particle
                                     Particles::kLight);  // light jet (truth)
        fParticlesModel->AddParticle(&dummy,
                                     Particles::kParton,
                                     "light quark 3",
                                     6,                   // index of corresponding particle
                                     Particles::kLight);  // light jet (truth)
        fParticlesModel->AddParticle(&dummy,
                                     Particles::kParton,
                                     "light quark 4",
                                     7,                   // index of corresponding particle
                                     Particles::kLight);  // light jet (truth)
        fParticlesModel->AddParticle(&dummy,
                                     Particles::kParton,
                                     "light quark 5",
                                     8,                   // index of corresponding particle
                                     Particles::kLight);  // light jet (truth)
        fParticlesModel->AddParticle(&dummy,
                                     Particles::kParton,
                                     "light quark 6",
                                     9,                   // index of corresponding particle
                                     Particles::kLight);  // light jet (truth)

        if (m_lepton_type == kElectron) {
            fParticlesModel->AddParticle(&dummy, Particles::kElectron, "electron");
        } else if (m_lepton_type == kMuon) {
            fParticlesModel->AddParticle(&dummy, Particles::kMuon, "muon");
        }

        fParticlesModel->AddParticle(&dummy, Particles::kNeutrino, "neutrino");

        fParticlesModel->AddParticle(&dummy, Particles::kBoson, "hadronic W 1");
        fParticlesModel->AddParticle(&dummy, Particles::kBoson, "hadronic W 2");
        fParticlesModel->AddParticle(&dummy, Particles::kBoson, "hadronic W 3");
        fParticlesModel->AddParticle(&dummy, Particles::kBoson, "leptonic W");

        fParticlesModel->AddParticle(&dummy, Particles::kParton, "hadronic top 1");
        fParticlesModel->AddParticle(&dummy, Particles::kParton, "hadronic top 2");
        fParticlesModel->AddParticle(&dummy, Particles::kParton, "hadronic top 3");
        fParticlesModel->AddParticle(&dummy, Particles::kParton, "leptonic top");

        // no error
        return 1;
    }

    // ---------------------------------------------------------
    void LikelihoodFourTopLeptonJets::DefineParameters() {
        // add parameters of model
        AddParameter("energy hadronic b 1", fPhysicsConstants.MassBottom(), 1000.0);  // parBhadE
        AddParameter("energy hadronic b 2", fPhysicsConstants.MassBottom(), 1000.0);  // parBhadE
        AddParameter("energy hadronic b 3", fPhysicsConstants.MassBottom(), 1000.0);  // parBhadE
        AddParameter("energy leptonic b", fPhysicsConstants.MassBottom(), 1000.0);  // parBlepE
        AddParameter("energy light quark 1", 0.0, 1000.0);                          // parLQ1E
        AddParameter("energy light quark 2", 0.0, 1000.0);                          // parLQ2E
        AddParameter("energy light quark 3", 0.0, 1000.0);                          // parLQ2E
        AddParameter("energy light quark 4", 0.0, 1000.0);                          // parLQ2E
        AddParameter("energy light quark 5", 0.0, 1000.0);                          // parLQ2E
        AddParameter("energy light quark 6", 0.0, 1000.0);                          // parLQ2E
        AddParameter("energy lepton", 0.0, 1000.0);                                 // parLepE
        AddParameter("p_x neutrino", -1000.0, 1000.0);                              // parNuPx
        AddParameter("p_y neutrino", -1000.0, 1000.0);                              // parNuPy
        AddParameter("p_z neutrino", -1000.0, 1000.0);                              // parNuPz
        AddParameter("top mass", 100.0, 1000.0);                                   // parTopM
    }

    // ---------------------------------------------------------
    int LikelihoodFourTopLeptonJets::CalculateLorentzVectors(std::vector <double> const& parameters) {
        double scale;
        
        double thad1_fit_e;
        double thad1_fit_px;
        double thad1_fit_py;
        double thad1_fit_pz;

        double thad2_fit_e;
        double thad2_fit_px;
        double thad2_fit_py;
        double thad2_fit_pz;

        double thad3_fit_e;
        double thad3_fit_px;
        double thad3_fit_py;
        double thad3_fit_pz;
        
        
        double tlep_fit_e;
        double tlep_fit_px;
        double tlep_fit_py;
        double tlep_fit_pz;

        // hadronic b quark
        m_bhad1_fit_e = parameters[parBhad1E];
        scale = sqrt(m_bhad1_fit_e * m_bhad1_fit_e - m_bhad1_meas_m * m_bhad1_meas_m) / m_bhad1_meas_p;
        m_bhad1_fit_px = scale * m_bhad1_meas_px;
        m_bhad1_fit_py = scale * m_bhad1_meas_py;
        m_bhad1_fit_pz = scale * m_bhad1_meas_pz;

        // hadronic b quark
        m_bhad2_fit_e = parameters[parBhad2E];
        scale = sqrt(m_bhad2_fit_e * m_bhad2_fit_e - m_bhad2_meas_m * m_bhad2_meas_m) / m_bhad2_meas_p;
        m_bhad2_fit_px = scale * m_bhad2_meas_px;
        m_bhad2_fit_py = scale * m_bhad2_meas_py;
        m_bhad2_fit_pz = scale * m_bhad2_meas_pz;

        // hadronic b quark
        m_bhad3_fit_e = parameters[parBhad3E];
        scale = sqrt(m_bhad3_fit_e * m_bhad3_fit_e - m_bhad3_meas_m * m_bhad3_meas_m) / m_bhad3_meas_p;
        m_bhad3_fit_px = scale * m_bhad3_meas_px;
        m_bhad3_fit_py = scale * m_bhad3_meas_py;
        m_bhad3_fit_pz = scale * m_bhad3_meas_pz;

        // leptonic b quark
        m_blep_fit_e = parameters[parBlepE];
        scale = sqrt(m_blep_fit_e * m_blep_fit_e - m_blep_meas_m * m_blep_meas_m) / m_blep_meas_p;
        m_blep_fit_px = scale * m_blep_meas_px;
        m_blep_fit_py = scale * m_blep_meas_py;
        m_blep_fit_pz = scale * m_blep_meas_pz;

        // light quark 1
        m_lq1_fit_e = parameters[parLQ1E];
        scale = sqrt(m_lq1_fit_e * m_lq1_fit_e - m_lq1_meas_m * m_lq1_meas_m) / m_lq1_meas_p;
        m_lq1_fit_px = scale * m_lq1_meas_px;
        m_lq1_fit_py = scale * m_lq1_meas_py;
        m_lq1_fit_pz = scale * m_lq1_meas_pz;

        // light quark 2
        m_lq2_fit_e = parameters[parLQ2E];
        scale = sqrt(m_lq2_fit_e * m_lq2_fit_e - m_lq2_meas_m * m_lq2_meas_m) / m_lq2_meas_p;
        m_lq2_fit_px  = scale * m_lq2_meas_px;
        m_lq2_fit_py  = scale * m_lq2_meas_py;
        m_lq2_fit_pz  = scale * m_lq2_meas_pz;

        // light quark 3
        m_lq3_fit_e = parameters[parLQ3E];
        scale = sqrt(m_lq3_fit_e * m_lq3_fit_e - m_lq3_meas_m * m_lq3_meas_m) / m_lq3_meas_p;
        m_lq3_fit_px  = scale * m_lq3_meas_px;
        m_lq3_fit_py  = scale * m_lq3_meas_py;
        m_lq3_fit_pz  = scale * m_lq3_meas_pz;

        // light quark 4
        m_lq4_fit_e = parameters[parLQ4E];
        scale = sqrt(m_lq4_fit_e * m_lq4_fit_e - m_lq4_meas_m * m_lq4_meas_m) / m_lq4_meas_p;
        m_lq4_fit_px  = scale * m_lq4_meas_px;
        m_lq4_fit_py  = scale * m_lq4_meas_py;
        m_lq4_fit_pz  = scale * m_lq4_meas_pz;

        // light quark 5
        m_lq5_fit_e = parameters[parLQ5E];
        scale = sqrt(m_lq5_fit_e * m_lq5_fit_e - m_lq5_meas_m * m_lq5_meas_m) / m_lq5_meas_p;
        m_lq5_fit_px  = scale * m_lq5_meas_px;
        m_lq5_fit_py  = scale * m_lq5_meas_py;
        m_lq5_fit_pz  = scale * m_lq5_meas_pz;

        // light quark 6
        m_lq6_fit_e = parameters[parLQ6E];
        scale = sqrt(m_lq6_fit_e * m_lq6_fit_e - m_lq6_meas_m * m_lq6_meas_m) / m_lq6_meas_p;
        m_lq6_fit_px  = scale * m_lq6_meas_px;
        m_lq6_fit_py  = scale * m_lq6_meas_py;
        m_lq6_fit_pz  = scale * m_lq6_meas_pz;

        // lepton
        m_lep_fit_e = parameters[parLepE];
        scale = m_lep_fit_e / m_lep_meas_e;
        m_lep_fit_px = scale * m_lep_meas_px;
        m_lep_fit_py = scale * m_lep_meas_py;
        m_lep_fit_pz = scale * m_lep_meas_pz;

        // neutrino
        m_nu_fit_px = parameters[parNuPx];
        m_nu_fit_py = parameters[parNuPy];
        m_nu_fit_pz = parameters[parNuPz];
        m_nu_fit_e  = sqrt(m_nu_fit_px * m_nu_fit_px + m_nu_fit_py * m_nu_fit_py + m_nu_fit_pz * m_nu_fit_pz);

        // hadronic W1
        m_whad1_fit_e  = m_lq1_fit_e + m_lq2_fit_e;
        m_whad1_fit_px = m_lq1_fit_px + m_lq2_fit_px;
        m_whad1_fit_py = m_lq1_fit_py + m_lq2_fit_py;
        m_whad1_fit_pz = m_lq1_fit_pz + m_lq2_fit_pz;
        m_whad1_fit_m = sqrt(m_whad1_fit_e * m_whad1_fit_e - (m_whad1_fit_px * m_whad1_fit_px + m_whad1_fit_py * m_whad1_fit_py + m_whad1_fit_pz * m_whad1_fit_pz));

        // hadronic W2
        m_whad2_fit_e  = m_lq3_fit_e + m_lq4_fit_e;
        m_whad2_fit_px = m_lq3_fit_px + m_lq4_fit_px;
        m_whad2_fit_py = m_lq3_fit_py + m_lq4_fit_py;
        m_whad2_fit_pz = m_lq3_fit_pz + m_lq4_fit_pz;
        m_whad2_fit_m = sqrt(m_whad2_fit_e * m_whad2_fit_e - (m_whad2_fit_px * m_whad2_fit_px + m_whad2_fit_py * m_whad2_fit_py + m_whad2_fit_pz * m_whad2_fit_pz));

        // hadronic W3
        m_whad3_fit_e  = m_lq5_fit_e + m_lq6_fit_e;
        m_whad3_fit_px = m_lq5_fit_px + m_lq6_fit_px;
        m_whad3_fit_py = m_lq5_fit_py + m_lq6_fit_py;
        m_whad3_fit_pz = m_lq5_fit_pz + m_lq6_fit_pz;
        m_whad3_fit_m = sqrt(m_whad3_fit_e * m_whad3_fit_e - (m_whad3_fit_px * m_whad3_fit_px + m_whad3_fit_py * m_whad3_fit_py + m_whad3_fit_pz * m_whad3_fit_pz));

        // leptonic W
        m_wlep_fit_e  = m_lep_fit_e + m_nu_fit_e;
        m_wlep_fit_px = m_lep_fit_px + m_nu_fit_px;
        m_wlep_fit_py = m_lep_fit_py + m_nu_fit_py;
        m_wlep_fit_pz = m_lep_fit_pz + m_nu_fit_pz;
        m_wlep_fit_m = sqrt(m_wlep_fit_e * m_wlep_fit_e - (m_wlep_fit_px * m_wlep_fit_px + m_wlep_fit_py * m_wlep_fit_py + m_wlep_fit_pz * m_wlep_fit_pz));

        // hadronic top
        thad1_fit_e = m_whad1_fit_e + m_bhad1_fit_e;
        thad1_fit_px = m_whad1_fit_px + m_bhad1_fit_px;
        thad1_fit_py = m_whad1_fit_py + m_bhad1_fit_py;
        thad1_fit_pz = m_whad1_fit_pz + m_bhad1_fit_pz;
        m_thad1_fit_m = sqrt(thad1_fit_e * thad1_fit_e - (thad1_fit_px * thad1_fit_px + thad1_fit_py * thad1_fit_py + thad1_fit_pz * thad1_fit_pz));

        // hadronic top
        thad2_fit_e = m_whad2_fit_e + m_bhad2_fit_e;
        thad2_fit_px = m_whad2_fit_px + m_bhad2_fit_px;
        thad2_fit_py = m_whad2_fit_py + m_bhad2_fit_py;
        thad2_fit_pz = m_whad2_fit_pz + m_bhad2_fit_pz;
        m_thad2_fit_m = sqrt(thad2_fit_e * thad2_fit_e - (thad2_fit_px * thad2_fit_px + thad2_fit_py * thad2_fit_py + thad2_fit_pz * thad2_fit_pz));

        // hadronic top
        thad3_fit_e = m_whad3_fit_e + m_bhad3_fit_e;
        thad3_fit_px = m_whad3_fit_px + m_bhad3_fit_px;
        thad3_fit_py = m_whad3_fit_py + m_bhad3_fit_py;
        thad3_fit_pz = m_whad3_fit_pz + m_bhad3_fit_pz;
        m_thad3_fit_m = sqrt(thad3_fit_e * thad3_fit_e - (thad3_fit_px * thad3_fit_px + thad3_fit_py * thad3_fit_py + thad3_fit_pz * thad3_fit_pz));

        // leptonic top
        tlep_fit_e = m_wlep_fit_e + m_blep_fit_e;
        tlep_fit_px = m_wlep_fit_px + m_blep_fit_px;
        tlep_fit_py = m_wlep_fit_py + m_blep_fit_py;
        tlep_fit_pz = m_wlep_fit_pz + m_blep_fit_pz;
        m_tlep_fit_m = sqrt(tlep_fit_e * tlep_fit_e - (tlep_fit_px * tlep_fit_px + tlep_fit_py * tlep_fit_py + tlep_fit_pz * tlep_fit_pz));

        // no error
        return 1;
    }

// ---------------------------------------------------------
    int LikelihoodFourTopLeptonJets::RemoveInvariantParticlePermutations() {
        // error code
        int err = 1;
        // return error code
        return err;
    }

    // ---------------------------------------------------------
    int LikelihoodFourTopLeptonJets::AdjustParameterRanges() {
        // adjust limits
        double nsigmas_jet    = m_flag_get_par_sigmas_from_TFs ? 10 : 7;
        double nsigmas_lepton = m_flag_get_par_sigmas_from_TFs ? 10 : 2;
        double nsigmas_met    = m_flag_get_par_sigmas_from_TFs ? 10 : 1;

        double E = (*fParticlesPermuted)->Parton(0)->E();
        double m = fPhysicsConstants.MassBottom();
        if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(0)->M());
        double sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_bhad1->GetSigma(E) : sqrt(E);
        double Emin = std::max(m, E - nsigmas_jet * sigma);
        double Emax = E + nsigmas_jet * sigma;

        SetParameterRange(parBhad1E, Emin, Emax);

        E = (*fParticlesPermuted)->Parton(1)->E();
        m = fPhysicsConstants.MassBottom();
        if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(1)->M());
        sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_bhad2->GetSigma(E) : sqrt(E);
        Emin = std::max(m, E - nsigmas_jet * sigma);
        Emax = E + nsigmas_jet * sigma;

        SetParameterRange(parBhad2E, Emin, Emax);

        E = (*fParticlesPermuted)->Parton(2)->E();
        m = fPhysicsConstants.MassBottom();
        if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(2)->M());
        sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_bhad3->GetSigma(E) : sqrt(E);
        Emin = std::max(m, E - nsigmas_jet * sigma);
        Emax = E + nsigmas_jet * sigma;

        SetParameterRange(parBhad3E, Emin, Emax);

        E = (*fParticlesPermuted)->Parton(3)->E();
        m = fPhysicsConstants.MassBottom();
        if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(3)->M());
        sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_blep->GetSigma(E) : sqrt(E);
        Emin = std::max(m, E - nsigmas_jet * sigma);
        Emax = E + nsigmas_jet * sigma;

        SetParameterRange(parBlepE, Emin, Emax);

        E = (*fParticlesPermuted)->Parton(4)->E();
        m = 0.001;
        if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(4)->M());
        sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq1->GetSigma(E) : sqrt(E);
        Emin = std::max(m, E - nsigmas_jet * sigma);
        Emax = E + nsigmas_jet * sigma;
        SetParameterRange(parLQ1E, Emin, Emax);

        E = (*fParticlesPermuted)->Parton(5)->E();
        m = 0.001;
        if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(5)->M());
        sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq2->GetSigma(E) : sqrt(E);
        Emin = std::max(m, E - nsigmas_jet * sigma);
        Emax = E + nsigmas_jet * sigma;
        SetParameterRange(parLQ2E, Emin, Emax);

        E = (*fParticlesPermuted)->Parton(6)->E();
        m = 0.001;
        if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(6)->M());
        sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq3->GetSigma(E) : sqrt(E);
        Emin = std::max(m, E - nsigmas_jet * sigma);
        Emax = E + nsigmas_jet * sigma;
        SetParameterRange(parLQ3E, Emin, Emax);

        E = (*fParticlesPermuted)->Parton(7)->E();
        m = 0.001;
        if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(7)->M());
        sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq4->GetSigma(E) : sqrt(E);
        Emin = std::max(m, E - nsigmas_jet * sigma);
        Emax = E + nsigmas_jet * sigma;
        SetParameterRange(parLQ4E, Emin, Emax);

        E = (*fParticlesPermuted)->Parton(8)->E();
        m = 0.001;
        if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(8)->M());
        sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq5->GetSigma(E) : sqrt(E);
        Emin = std::max(m, E - nsigmas_jet * sigma);
        Emax = E + nsigmas_jet * sigma;
        SetParameterRange(parLQ5E, Emin, Emax);

        E = (*fParticlesPermuted)->Parton(9)->E();
        m = 0.001;
        if (fFlagUseJetMass) m = std::max(0.0, (*fParticlesPermuted)->Parton(9)->M());
        sigma = m_flag_get_par_sigmas_from_TFs ? m_res_energy_lq6->GetSigma(E) : sqrt(E);
        Emin = std::max(m, E - nsigmas_jet * sigma);
        Emax = E + nsigmas_jet * sigma;
        SetParameterRange(parLQ6E, Emin, Emax);

        if (m_lepton_type == kElectron) {
            E = (*fParticlesPermuted)->Electron(0)->E();
            sigma = m_flag_get_par_sigmas_from_TFs ? m_res_lepton->GetSigma(E) : sqrt(E);
            Emin = std::max(0.001, E - nsigmas_lepton * sigma);
            Emax = E + nsigmas_lepton * sigma;
        } else if (m_lepton_type == kMuon) {
            E = (*fParticlesPermuted)->Muon(0)->E();
            double sintheta = sin((*fParticlesPermuted)->Muon(0)->Theta());
            sigma = m_flag_get_par_sigmas_from_TFs ? m_res_lepton->GetSigma(E * sintheta) / sintheta : E * E * sintheta;
            double sigrange = nsigmas_lepton * sigma;
            Emin = std::max(0.001, E - sigrange);
            Emax = E + sigrange;
        }
        SetParameterRange(parLepE, Emin, Emax);

        sigma = m_flag_get_par_sigmas_from_TFs ? m_res_met->GetSigma(m_et_miss_sum) : 100;
        double sigrange = nsigmas_met * sigma;
        SetParameterRange(parNuPx, m_et_miss_x - sigrange, m_et_miss_x + sigrange);
        SetParameterRange(parNuPy, m_et_miss_y - sigrange, m_et_miss_y + sigrange);

        if (m_flag_top_mass_fixed)
            SetParameterRange(parTopM, fPhysicsConstants.MassTop(), fPhysicsConstants.MassTop());

        // no error
        return 1;
    }

    // ---------------------------------------------------------
    double LikelihoodFourTopLeptonJets::LogLikelihood(const std::vector<double> & parameters) {
        // calculate 4-vectors
        CalculateLorentzVectors(parameters);

        // define log of likelihood
        double logprob(0.);

        // temporary flag for a safe use of the transfer functions
        bool TFgoodTmp(true);

        // jet energy resolution terms
        logprob += m_res_energy_bhad1->logp(m_bhad1_fit_e, m_bhad1_meas_e, &TFgoodTmp);
        if (!TFgoodTmp) fTFgood = false;
        
        logprob += m_res_energy_bhad2->logp(m_bhad2_fit_e, m_bhad2_meas_e, &TFgoodTmp);
        if (!TFgoodTmp) fTFgood = false;
        
        logprob += m_res_energy_bhad3->logp(m_bhad3_fit_e, m_bhad3_meas_e, &TFgoodTmp);
        if (!TFgoodTmp) fTFgood = false;

        logprob += m_res_energy_blep->logp(m_blep_fit_e, m_blep_meas_e, &TFgoodTmp);
        if (!TFgoodTmp) fTFgood = false;

        logprob += m_res_energy_lq1->logp(m_lq1_fit_e, m_lq1_meas_e, &TFgoodTmp);
        if (!TFgoodTmp) fTFgood = false;

        logprob += m_res_energy_lq2->logp(m_lq2_fit_e, m_lq2_meas_e, &TFgoodTmp);
        if (!TFgoodTmp) fTFgood = false;

        logprob += m_res_energy_lq3->logp(m_lq3_fit_e, m_lq3_meas_e, &TFgoodTmp);
        if (!TFgoodTmp) fTFgood = false;

        logprob += m_res_energy_lq4->logp(m_lq4_fit_e, m_lq4_meas_e, &TFgoodTmp);
        if (!TFgoodTmp) fTFgood = false;

        logprob += m_res_energy_lq5->logp(m_lq5_fit_e, m_lq5_meas_e, &TFgoodTmp);
        if (!TFgoodTmp) fTFgood = false;

        logprob += m_res_energy_lq6->logp(m_lq6_fit_e, m_lq6_meas_e, &TFgoodTmp);
        if (!TFgoodTmp) fTFgood = false;

        // lepton energy resolution terms
        if (m_lepton_type == kElectron) {
            logprob += m_res_lepton->logp(m_lep_fit_e, m_lep_meas_e, &TFgoodTmp);
        } else if (m_lepton_type == kMuon) {
            logprob += m_res_lepton->logp(m_lep_fit_e * m_lep_meas_sintheta, m_lep_meas_pt, &TFgoodTmp);
        }
        if (!TFgoodTmp) fTFgood = false;

        // neutrino px and py
        logprob += m_res_met->logp(m_nu_fit_px, m_et_miss_x, &TFgoodTmp, m_et_miss_sum);
        if (!TFgoodTmp) fTFgood = false;

        logprob += m_res_met->logp(m_nu_fit_py, m_et_miss_y, &TFgoodTmp, m_et_miss_sum);
        if (!TFgoodTmp) fTFgood = false;

        // physics constants
        double massW = fPhysicsConstants.MassW();
        double gammaW = fPhysicsConstants.GammaW();
        double gammaTop = fPhysicsConstants.GammaTop();

        // Breit-Wigner of hadronically decaying W-boson
        logprob += BCMath::LogBreitWignerRel(m_whad1_fit_m, massW, gammaW);
        logprob += BCMath::LogBreitWignerRel(m_whad2_fit_m, massW, gammaW);
        logprob += BCMath::LogBreitWignerRel(m_whad3_fit_m, massW, gammaW);

        // Breit-Wigner of leptonically decaying W-boson
        logprob += BCMath::LogBreitWignerRel(m_wlep_fit_m, massW, gammaW);

        // Breit-Wigner of hadronically decaying top quark
        logprob += BCMath::LogBreitWignerRel(m_thad1_fit_m, parameters[parTopM], gammaTop);
        logprob += BCMath::LogBreitWignerRel(m_thad2_fit_m, parameters[parTopM], gammaTop);
        logprob += BCMath::LogBreitWignerRel(m_thad3_fit_m, parameters[parTopM], gammaTop);

        // Breit-Wigner of leptonically decaying top quark
        logprob += BCMath::LogBreitWignerRel(m_tlep_fit_m, parameters[parTopM], gammaTop);

        // return log of likelihood
        return logprob;
    }

    // ---------------------------------------------------------
    std::vector<double> LikelihoodFourTopLeptonJets::GetInitialParameters() {
        std::vector<double> values(GetInitialParametersWoNeutrinoPz());

        // check second neutrino solution
        std::vector<double> neutrino_pz_solutions = GetNeutrinoPzSolutions();
        if (neutrino_pz_solutions.size() == 1) {
            values[parNuPz] = neutrino_pz_solutions[0];
        } else if (neutrino_pz_solutions.size() == 2) {
            double sol1, sol2;
            values[parNuPz] = neutrino_pz_solutions[0];
            sol1 = LogLikelihood(values);
            values[parNuPz] = neutrino_pz_solutions[1];
            sol2 = LogLikelihood(values);

            if (sol1 > sol2) values[parNuPz] = neutrino_pz_solutions[0];
        }

        return values;
    }

    // ---------------------------------------------------------
    std::vector<double> LikelihoodFourTopLeptonJets::GetInitialParametersWoNeutrinoPz() {
        std::vector<double> values(GetNParameters());

        // energies of the quarks
        values[parBhad1E] = m_bhad1_meas_e;
        values[parBhad2E] = m_bhad2_meas_e;
        values[parBhad3E] = m_bhad3_meas_e;
        values[parBlepE] = m_blep_meas_e;
        values[parLQ1E]  = m_lq1_meas_e;
        values[parLQ2E]  = m_lq2_meas_e;
        values[parLQ3E]  = m_lq3_meas_e;
        values[parLQ4E]  = m_lq4_meas_e;
        values[parLQ5E]  = m_lq5_meas_e;
        values[parLQ6E]  = m_lq6_meas_e;

        // energy of the lepton
        if (m_lepton_type == kElectron) {
            values[parLepE] = (*fParticlesPermuted)->Electron(0)->E();
        } else if (m_lepton_type == kMuon) {
            values[parLepE] = (*fParticlesPermuted)->Muon(0)->E();
        }

        // missing px and py
        values[parNuPx] = m_et_miss_x;
        values[parNuPy] = m_et_miss_y;

        // pz of the neutrino
        values[parNuPz] = 0.;

        // top mass
        double mtop = (*(*fParticlesPermuted)->Parton(0) + *(*fParticlesPermuted)->Parton(4) + *(*fParticlesPermuted)->Parton(5)).M(); // 1 -3 are other bquarks
        if (mtop < GetParameter(parTopM)->GetLowerLimit()) {
            mtop = GetParameter(parTopM)->GetLowerLimit();
        } else if (mtop > GetParameter(parTopM)->GetUpperLimit()) {
            mtop = GetParameter(parTopM)->GetUpperLimit();
        }
        values[parTopM] = mtop;

        // return the vector
        return values;
    }
// ---------------------------------------------------------
    std::vector<double> LikelihoodFourTopLeptonJets::GetNeutrinoPzSolutions() {
        return CalculateNeutrinoPzSolutions();
    }

// ---------------------------------------------------------
    std::vector<double> LikelihoodFourTopLeptonJets::CalculateNeutrinoPzSolutions(TLorentzVector* additionalParticle) {
        std::vector<double> pz;

        class PhysicsConstants constants;
        // electron mass
        double mE = 0.;

        double px_c = 0.0;
        double py_c = 0.0;
        double pz_c = 0.0;
        double Ec = 0.0;

        if (m_lepton_type == kElectron) {
            px_c = (*fParticlesPermuted)->Electron(0)->Px();
            py_c = (*fParticlesPermuted)->Electron(0)->Py();
            pz_c = (*fParticlesPermuted)->Electron(0)->Pz();
            Ec = (*fParticlesPermuted)->Electron(0)->E();
        } else if (m_lepton_type == kMuon) {
            px_c = (*fParticlesPermuted)->Muon(0)->Px();
            py_c = (*fParticlesPermuted)->Muon(0)->Py();
            pz_c = (*fParticlesPermuted)->Muon(0)->Pz();
            Ec = (*fParticlesPermuted)->Muon(0)->E();
        }

        // add additional particle to "charged lepton" 4-vector
        if (additionalParticle) {
            px_c += additionalParticle->Px();
            py_c += additionalParticle->Py();
            pz_c += additionalParticle->Pz();
            Ec += additionalParticle->E();
        }

        double px_nu = m_et_miss_x;
        double py_nu = m_et_miss_y;
        double alpha = constants.MassW() * constants.MassW() - mE * mE + 2 * (px_c * px_nu + py_c * py_nu);

        double a = pz_c * pz_c - Ec * Ec;
        double b = alpha * pz_c;
        double c = - Ec * Ec * (px_nu * px_nu + py_nu * py_nu) + alpha * alpha / 4.;

        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0.) return pz;

        double pz_offset = - b / (2 * a);

        double squareRoot = sqrt(discriminant);
        if (squareRoot < 1.e-6) {
            pz.push_back(pz_offset);
        } else {
            pz.push_back(pz_offset + squareRoot / (2 * a));
            pz.push_back(pz_offset - squareRoot / (2 * a));
        }

        return pz;
    }

    // ---------------------------------------------------------
    int LikelihoodFourTopLeptonJets::SavePermutedParticles() {
        m_bhad1_meas_e      = (*fParticlesPermuted)->Parton(0)->E();
        m_bhad1_meas_deteta = (*fParticlesPermuted)->DetEta(0, Particles::kParton);
        m_bhad1_meas_px     = (*fParticlesPermuted)->Parton(0)->Px();
        m_bhad1_meas_py     = (*fParticlesPermuted)->Parton(0)->Py();
        m_bhad1_meas_pz     = (*fParticlesPermuted)->Parton(0)->Pz();
        m_bhad1_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(0)->M(), fPhysicsConstants.MassBottom(), &m_bhad1_meas_px, &m_bhad1_meas_py, &m_bhad1_meas_pz, m_bhad1_meas_e);
        m_bhad1_meas_p      = sqrt(m_bhad1_meas_e*m_bhad1_meas_e - m_bhad1_meas_m*m_bhad1_meas_m);

        m_bhad2_meas_e      = (*fParticlesPermuted)->Parton(1)->E();
        m_bhad2_meas_deteta = (*fParticlesPermuted)->DetEta(1, Particles::kParton);
        m_bhad2_meas_px     = (*fParticlesPermuted)->Parton(1)->Px();
        m_bhad2_meas_py     = (*fParticlesPermuted)->Parton(1)->Py();
        m_bhad2_meas_pz     = (*fParticlesPermuted)->Parton(1)->Pz();
        m_bhad2_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(1)->M(), fPhysicsConstants.MassBottom(), &m_bhad2_meas_px, &m_bhad2_meas_py, &m_bhad2_meas_pz, m_bhad2_meas_e);
        m_bhad2_meas_p      = sqrt(m_bhad2_meas_e*m_bhad2_meas_e - m_bhad2_meas_m*m_bhad2_meas_m);

        m_bhad3_meas_e      = (*fParticlesPermuted)->Parton(2)->E();
        m_bhad3_meas_deteta = (*fParticlesPermuted)->DetEta(2, Particles::kParton);
        m_bhad3_meas_px     = (*fParticlesPermuted)->Parton(2)->Px();
        m_bhad3_meas_py     = (*fParticlesPermuted)->Parton(2)->Py();
        m_bhad3_meas_pz     = (*fParticlesPermuted)->Parton(2)->Pz();
        m_bhad3_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(2)->M(), fPhysicsConstants.MassBottom(), &m_bhad3_meas_px, &m_bhad3_meas_py, &m_bhad3_meas_pz, m_bhad3_meas_e);
        m_bhad3_meas_p      = sqrt(m_bhad3_meas_e*m_bhad3_meas_e - m_bhad3_meas_m*m_bhad3_meas_m);


        m_blep_meas_e      = (*fParticlesPermuted)->Parton(3)->E();
        m_blep_meas_deteta = (*fParticlesPermuted)->DetEta(3, Particles::kParton);
        m_blep_meas_px     = (*fParticlesPermuted)->Parton(3)->Px();
        m_blep_meas_py     = (*fParticlesPermuted)->Parton(3)->Py();
        m_blep_meas_pz     = (*fParticlesPermuted)->Parton(3)->Pz();
        m_blep_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(3)->M(), fPhysicsConstants.MassBottom(), &m_blep_meas_px, &m_blep_meas_py, &m_blep_meas_pz, m_blep_meas_e);
        m_blep_meas_p      = sqrt(m_blep_meas_e*m_blep_meas_e - m_blep_meas_m*m_blep_meas_m);

        m_lq1_meas_e      = (*fParticlesPermuted)->Parton(4)->E();
        m_lq1_meas_deteta = (*fParticlesPermuted)->DetEta(4, Particles::kParton);
        m_lq1_meas_px     = (*fParticlesPermuted)->Parton(4)->Px();
        m_lq1_meas_py     = (*fParticlesPermuted)->Parton(4)->Py();
        m_lq1_meas_pz     = (*fParticlesPermuted)->Parton(4)->Pz();
        m_lq1_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(4)->M(), 0., &m_lq1_meas_px, &m_lq1_meas_py, &m_lq1_meas_pz, m_lq1_meas_e);
        m_lq1_meas_p      = sqrt(m_lq1_meas_e*m_lq1_meas_e - m_lq1_meas_m*m_lq1_meas_m);

        m_lq2_meas_e      = (*fParticlesPermuted)->Parton(5)->E();
        m_lq2_meas_deteta = (*fParticlesPermuted)->DetEta(5, Particles::kParton);
        m_lq2_meas_px     = (*fParticlesPermuted)->Parton(5)->Px();
        m_lq2_meas_py     = (*fParticlesPermuted)->Parton(5)->Py();
        m_lq2_meas_pz     = (*fParticlesPermuted)->Parton(5)->Pz();
        m_lq2_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(5)->M(), 0., &m_lq2_meas_px, &m_lq2_meas_py, &m_lq2_meas_pz, m_lq2_meas_e);
        m_lq2_meas_p      = sqrt(m_lq2_meas_e*m_lq2_meas_e - m_lq2_meas_m*m_lq2_meas_m);

        m_lq3_meas_e      = (*fParticlesPermuted)->Parton(6)->E();
        m_lq3_meas_deteta = (*fParticlesPermuted)->DetEta(6, Particles::kParton);
        m_lq3_meas_px     = (*fParticlesPermuted)->Parton(6)->Px();
        m_lq3_meas_py     = (*fParticlesPermuted)->Parton(6)->Py();
        m_lq3_meas_pz     = (*fParticlesPermuted)->Parton(6)->Pz();
        m_lq3_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(6)->M(), 0., &m_lq3_meas_px, &m_lq3_meas_py, &m_lq3_meas_pz, m_lq3_meas_e);
        m_lq3_meas_p      = sqrt(m_lq3_meas_e*m_lq3_meas_e - m_lq3_meas_m*m_lq3_meas_m);

        m_lq4_meas_e      = (*fParticlesPermuted)->Parton(7)->E();
        m_lq4_meas_deteta = (*fParticlesPermuted)->DetEta(7, Particles::kParton);
        m_lq4_meas_px     = (*fParticlesPermuted)->Parton(7)->Px();
        m_lq4_meas_py     = (*fParticlesPermuted)->Parton(7)->Py();
        m_lq4_meas_pz     = (*fParticlesPermuted)->Parton(7)->Pz();
        m_lq4_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(7)->M(), 0., &m_lq4_meas_px, &m_lq4_meas_py, &m_lq4_meas_pz, m_lq4_meas_e);
        m_lq4_meas_p      = sqrt(m_lq4_meas_e*m_lq4_meas_e - m_lq4_meas_m*m_lq4_meas_m);

        m_lq5_meas_e      = (*fParticlesPermuted)->Parton(8)->E();
        m_lq5_meas_deteta = (*fParticlesPermuted)->DetEta(8, Particles::kParton);
        m_lq5_meas_px     = (*fParticlesPermuted)->Parton(8)->Px();
        m_lq5_meas_py     = (*fParticlesPermuted)->Parton(8)->Py();
        m_lq5_meas_pz     = (*fParticlesPermuted)->Parton(8)->Pz();
        m_lq5_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(8)->M(), 0., &m_lq5_meas_px, &m_lq5_meas_py, &m_lq5_meas_pz, m_lq5_meas_e);
        m_lq5_meas_p      = sqrt(m_lq5_meas_e*m_lq5_meas_e - m_lq5_meas_m*m_lq5_meas_m);

        m_lq6_meas_e      = (*fParticlesPermuted)->Parton(9)->E();
        m_lq6_meas_deteta = (*fParticlesPermuted)->DetEta(9, Particles::kParton);
        m_lq6_meas_px     = (*fParticlesPermuted)->Parton(9)->Px();
        m_lq6_meas_py     = (*fParticlesPermuted)->Parton(9)->Py();
        m_lq6_meas_pz     = (*fParticlesPermuted)->Parton(9)->Pz();
        m_lq6_meas_m      = SetPartonMass((*fParticlesPermuted)->Parton(9)->M(), 0., &m_lq6_meas_px, &m_lq6_meas_py, &m_lq6_meas_pz, m_lq6_meas_e);
        m_lq6_meas_p      = sqrt(m_lq6_meas_e*m_lq6_meas_e - m_lq6_meas_m*m_lq6_meas_m);

        TLorentzVector* lepton(0);
        if (m_lepton_type == kElectron) {
            lepton = (*fParticlesPermuted)->Electron(0);
            m_lep_meas_deteta = (*fParticlesPermuted)->DetEta(0, Particles::kElectron);
        } else {
            lepton = (*fParticlesPermuted)->Muon(0);
            m_lep_meas_deteta = (*fParticlesPermuted)->DetEta(0, Particles::kMuon);
        }
        m_lep_meas_e        = lepton->E();
        m_lep_meas_sintheta = sin(lepton->Theta());
        m_lep_meas_pt       = lepton->Pt();
        m_lep_meas_px       = lepton->Px();
        m_lep_meas_py       = lepton->Py();
        m_lep_meas_pz       = lepton->Pz();

        // no error
        return 1;
    }

    // ---------------------------------------------------------
    int LikelihoodFourTopLeptonJets::SaveResolutionFunctions() {
        m_res_energy_bhad1 = (*fDetector)->ResEnergyBJet(m_bhad1_meas_deteta);
        m_res_energy_bhad2 = (*fDetector)->ResEnergyBJet(m_bhad2_meas_deteta);
        m_res_energy_bhad3 = (*fDetector)->ResEnergyBJet(m_bhad3_meas_deteta);
        m_res_energy_blep = (*fDetector)->ResEnergyBJet(m_blep_meas_deteta);
        m_res_energy_lq1  = (*fDetector)->ResEnergyLightJet(m_lq1_meas_deteta);
        m_res_energy_lq2  = (*fDetector)->ResEnergyLightJet(m_lq2_meas_deteta);
        m_res_energy_lq3  = (*fDetector)->ResEnergyLightJet(m_lq3_meas_deteta);
        m_res_energy_lq4  = (*fDetector)->ResEnergyLightJet(m_lq4_meas_deteta);
        m_res_energy_lq5  = (*fDetector)->ResEnergyLightJet(m_lq5_meas_deteta);
        m_res_energy_lq6  = (*fDetector)->ResEnergyLightJet(m_lq6_meas_deteta);
        if (m_lepton_type == kElectron) {
            m_res_lepton = (*fDetector)->ResEnergyElectron(m_lep_meas_deteta);
        } else if (m_lepton_type == kMuon) {
            m_res_lepton = (*fDetector)->ResEnergyMuon(m_lep_meas_deteta);
        }
        m_res_met = (*fDetector)->ResMissingET();

        // no error
        return 1;
    }

    // ---------------------------------------------------------
    int LikelihoodFourTopLeptonJets::BuildModelParticles() {
        if (!GetBestFitParameters().empty()) CalculateLorentzVectors(GetBestFitParameters());

        TLorentzVector* bhad1 = fParticlesModel->Parton(0);
        TLorentzVector* bhad2 = fParticlesModel->Parton(1);
        TLorentzVector* bhad3 = fParticlesModel->Parton(2);
        TLorentzVector* blep = fParticlesModel->Parton(3);
        TLorentzVector* lq1  = fParticlesModel->Parton(4);
        TLorentzVector* lq2  = fParticlesModel->Parton(5);
        TLorentzVector* lq3  = fParticlesModel->Parton(6);
        TLorentzVector* lq4  = fParticlesModel->Parton(7);
        TLorentzVector* lq5  = fParticlesModel->Parton(8);
        TLorentzVector* lq6  = fParticlesModel->Parton(9);
        TLorentzVector* lep(0);
        if (m_lepton_type == kElectron) {
            lep  = fParticlesModel->Electron(0);
        } else if (m_lepton_type == kMuon) {
            lep  = fParticlesModel->Muon(0);
        } else {
            std::cerr << "Lepton type not set in LikelihoodFourTopLeptonJets::BuildModelParticles" << std::endl;
            return 0;
        }
        TLorentzVector* nu   = fParticlesModel->Neutrino(0);
        TLorentzVector* whad1  = fParticlesModel->Boson(0);
        TLorentzVector* whad2  = fParticlesModel->Boson(1);
        TLorentzVector* whad3  = fParticlesModel->Boson(2);
        TLorentzVector* wlep  = fParticlesModel->Boson(3);
        TLorentzVector* thad1  = fParticlesModel->Parton(10);
        TLorentzVector* thad2  = fParticlesModel->Parton(11);
        TLorentzVector* thad3  = fParticlesModel->Parton(12);
        TLorentzVector* tlep  = fParticlesModel->Parton(13);

        bhad1->SetPxPyPzE(m_bhad1_fit_px, m_bhad1_fit_py, m_bhad1_fit_pz, m_bhad1_fit_e);
        bhad2->SetPxPyPzE(m_bhad2_fit_px, m_bhad2_fit_py, m_bhad2_fit_pz, m_bhad2_fit_e);
        bhad3->SetPxPyPzE(m_bhad3_fit_px, m_bhad3_fit_py, m_bhad3_fit_pz, m_bhad3_fit_e);
        blep->SetPxPyPzE(m_blep_fit_px, m_blep_fit_py, m_blep_fit_pz, m_blep_fit_e);
        lq1 ->SetPxPyPzE(m_lq1_fit_px,  m_lq1_fit_py,  m_lq1_fit_pz,  m_lq1_fit_e);
        lq2 ->SetPxPyPzE(m_lq2_fit_px,  m_lq2_fit_py,  m_lq2_fit_pz,  m_lq2_fit_e);
        lq3 ->SetPxPyPzE(m_lq3_fit_px,  m_lq3_fit_py,  m_lq3_fit_pz,  m_lq3_fit_e);
        lq4 ->SetPxPyPzE(m_lq4_fit_px,  m_lq4_fit_py,  m_lq4_fit_pz,  m_lq4_fit_e);
        lq5 ->SetPxPyPzE(m_lq5_fit_px,  m_lq5_fit_py,  m_lq5_fit_pz,  m_lq5_fit_e);
        lq6 ->SetPxPyPzE(m_lq6_fit_px,  m_lq6_fit_py,  m_lq6_fit_pz,  m_lq6_fit_e);
        lep ->SetPxPyPzE(m_lep_fit_px,  m_lep_fit_py,  m_lep_fit_pz,  m_lep_fit_e);
        nu  ->SetPxPyPzE(m_nu_fit_px,   m_nu_fit_py,   m_nu_fit_pz,   m_nu_fit_e);

        (*whad1) = (*lq1)  + (*lq2);
        (*whad2) = (*lq3)  + (*lq4);
        (*whad3) = (*lq5)  + (*lq6);
        (*wlep) = (*lep)  + (*nu);
        (*thad1) = (*whad1) + (*bhad1);
        (*thad2) = (*whad2) + (*bhad2);
        (*thad3) = (*whad3) + (*bhad3);
        (*tlep) = (*wlep) + (*blep);

        // no error
        return 1;
    }


    // ---------------------------------------------------------
    std::vector<double> LikelihoodFourTopLeptonJets::LogLikelihoodComponents(std::vector<double> parameters) {
        std::vector<double> vecci;

        // calculate 4-vectors
        CalculateLorentzVectors(parameters);

        // temporary flag for a safe use of the transfer functions
        bool TFgoodTmp(true);

        // jet energy resolution terms
        vecci.push_back(m_res_energy_bhad1->logp(m_bhad1_fit_e, m_bhad1_meas_e, &TFgoodTmp));  // comp0
        if (!TFgoodTmp) fTFgood = false;

        vecci.push_back(m_res_energy_bhad2->logp(m_bhad2_fit_e, m_bhad2_meas_e, &TFgoodTmp));  // comp0
        if (!TFgoodTmp) fTFgood = false;

        vecci.push_back(m_res_energy_bhad3->logp(m_bhad3_fit_e, m_bhad3_meas_e, &TFgoodTmp));  // comp0
        if (!TFgoodTmp) fTFgood = false;

        vecci.push_back(m_res_energy_blep->logp(m_blep_fit_e, m_blep_meas_e, &TFgoodTmp));  // comp1
        if (!TFgoodTmp) fTFgood = false;

        vecci.push_back(m_res_energy_lq1->logp(m_lq1_fit_e, m_lq1_meas_e, &TFgoodTmp));  // comp2
        if (!TFgoodTmp) fTFgood = false;

        vecci.push_back(m_res_energy_lq2->logp(m_lq2_fit_e, m_lq2_meas_e, &TFgoodTmp));  // comp3
        if (!TFgoodTmp) fTFgood = false;

        vecci.push_back(m_res_energy_lq3->logp(m_lq3_fit_e, m_lq3_meas_e, &TFgoodTmp));  // comp3
        if (!TFgoodTmp) fTFgood = false;

        vecci.push_back(m_res_energy_lq4->logp(m_lq4_fit_e, m_lq4_meas_e, &TFgoodTmp));  // comp3
        if (!TFgoodTmp) fTFgood = false;

        vecci.push_back(m_res_energy_lq5->logp(m_lq5_fit_e, m_lq5_meas_e, &TFgoodTmp));  // comp3
        if (!TFgoodTmp) fTFgood = false;

        vecci.push_back(m_res_energy_lq6->logp(m_lq6_fit_e, m_lq6_meas_e, &TFgoodTmp));  // comp3
        if (!TFgoodTmp) fTFgood = false;
        
        // lepton energy resolution terms
        if (m_lepton_type == kElectron) {
            vecci.push_back(m_res_lepton->logp(m_lep_fit_e, m_lep_meas_e, &TFgoodTmp));  // comp4
        } else if (m_lepton_type == kMuon) {
            vecci.push_back(m_res_lepton->logp(m_lep_fit_e* m_lep_meas_sintheta, m_lep_meas_pt, &TFgoodTmp));  // comp4
        }
        if (!TFgoodTmp) fTFgood = false;

        // neutrino px and py
        vecci.push_back(m_res_met->logp(m_nu_fit_px, m_et_miss_x, &TFgoodTmp, m_et_miss_sum));  // comp5
        if (!TFgoodTmp) fTFgood = false;

        vecci.push_back(m_res_met->logp(m_nu_fit_py, m_et_miss_y, &TFgoodTmp, m_et_miss_sum));  // comp6
        if (!TFgoodTmp) fTFgood = false;

        // physics constants
        double massW = fPhysicsConstants.MassW();
        double gammaW = fPhysicsConstants.GammaW();
        double gammaTop = fPhysicsConstants.GammaTop();

        // Breit-Wigner of hadronically decaying W-boson
        vecci.push_back(BCMath::LogBreitWignerRel(m_whad1_fit_m, massW, gammaW));  // comp7
        vecci.push_back(BCMath::LogBreitWignerRel(m_whad2_fit_m, massW, gammaW));  // comp7
        vecci.push_back(BCMath::LogBreitWignerRel(m_whad3_fit_m, massW, gammaW));  // comp7

        // Breit-Wigner of leptonically decaying W-boson
        vecci.push_back(BCMath::LogBreitWignerRel(m_wlep_fit_m, massW, gammaW));  // comp8

        // Breit-Wigner of hadronically decaying top quark
        vecci.push_back(BCMath::LogBreitWignerRel(m_thad1_fit_m, parameters[parTopM], gammaTop));  // comp9
        vecci.push_back(BCMath::LogBreitWignerRel(m_thad2_fit_m, parameters[parTopM], gammaTop));  // comp9
        vecci.push_back(BCMath::LogBreitWignerRel(m_thad3_fit_m, parameters[parTopM], gammaTop));  // comp9

        // Breit-Wigner of leptonically decaying top quark
        vecci.push_back(BCMath::LogBreitWignerRel(m_tlep_fit_m, parameters[parTopM], gammaTop));  // comp10

        // return log of likelihood
        return vecci;
    }


}