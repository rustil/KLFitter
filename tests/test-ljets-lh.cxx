#include "KLFitter/LikelihoodTopLeptonJets.h"

#include "TLorentzVector.h"
#include <iostream>
#include <memory>

namespace {
  std::unique_ptr<KLFitter::Particles> getExampleParticles(float tag_eff, float tag_ineff) {
    TLorentzVector jet1{};
    jet1.SetPtEtaPhiE(133.56953, 0.2231264, 1.7798618, 137.56292);
    const float jet1_btag_weight{0.6868029};
    const bool jet1_has_btag{false};

    TLorentzVector jet2{};
    jet2.SetPtEtaPhiE(77.834281, 0.8158330, -1.533635, 105.72334);
    const float jet2_btag_weight{-0.869940};
    const bool jet2_has_btag{false};

    TLorentzVector jet3{};
    jet3.SetPtEtaPhiE(49.327293, 1.9828589, -1.878274, 182.64006);
    const float jet3_btag_weight{0.9999086};
    const bool jet3_has_btag{true};

    TLorentzVector jet4{};
    jet4.SetPtEtaPhiE(43.140816, 0.4029131, -0.472721, 47.186804);
    const float jet4_btag_weight{-0.223728};
    const bool jet4_has_btag{false};

    TLorentzVector lep{};
    lep.SetPtEtaPhiE(30.501886, 0.4483959, 2.9649317, 33.620113);

    auto particles = std::make_unique<KLFitter::Particles>();
    particles->AddParticle(&jet1, jet1.Eta(), KLFitter::Particles::kParton, "", 0, jet1_has_btag, tag_eff, tag_ineff, KLFitter::Particles::kNone, jet1_btag_weight);
    particles->AddParticle(&jet2, jet2.Eta(), KLFitter::Particles::kParton, "", 1, jet2_has_btag, tag_eff, tag_ineff, KLFitter::Particles::kNone, jet2_btag_weight);
    particles->AddParticle(&jet3, jet3.Eta(), KLFitter::Particles::kParton, "", 2, jet3_has_btag, tag_eff, tag_ineff, KLFitter::Particles::kNone, jet3_btag_weight);
    particles->AddParticle(&jet4, jet4.Eta(), KLFitter::Particles::kParton, "", 3, jet4_has_btag, tag_eff, tag_ineff, KLFitter::Particles::kNone, jet4_btag_weight);
    particles->AddParticle(&lep, lep.Eta(), KLFitter::Particles::kMuon, "", 0);
    return particles;
  }
}

int main() {
  // Get one set of example particles. Assume the following
  // efficiencies for the jet b-tagging algorithm:
  //   - 0.7 tag rate
  //   - 1/125 type-II error (false positives)
  const auto particles = getExampleParticles(0.7, 125);

  KLFitter::LikelihoodTopLeptonJets lh{};

  std::cout << "Hello World" << std::endl;
  return 0;
}
