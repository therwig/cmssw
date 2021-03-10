#include <vector>
#include <numeric>

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
// #include "DataFormats/L1TParticleFlow/interface/PFJet.h"
// #include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/L1TParticleFlow/interface/PFMet.h"

class L1MetPFHWProducer : public edm::global::EDProducer<> {
public:
  explicit L1MetPFHWProducer(const edm::ParameterSet&);
  ~L1MetPFHWProducer() override;

private:
  /// ///////////////// ///
  /// MANDATORY METHODS ///
  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;
  /// ///////////////// ///

  unsigned _maxCands;
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> _l1PFToken;

    //static l1t::PFJet makeJet(const std::vector<edm::Ptr<l1t::PFCandidate>>& parts);
};

// l1t::PFJet L1MetPFHWProducer::makeJet(const std::vector<edm::Ptr<l1t::PFCandidate>>& parts) {
//   l1t::PFCandidate seed = *parts.at(0);

//   auto sumpt = [](float a, const edm::Ptr<l1t::PFCandidate>& b) { return a + b->pt(); };

//   // Sum the pt
//   float pt = std::accumulate(parts.begin(), parts.end(), 0., sumpt);

//   // pt weighted d eta
//   std::vector<float> pt_deta;
//   pt_deta.resize(parts.size());
//   std::transform(parts.begin(), parts.end(), pt_deta.begin(), [&seed, &pt](const edm::Ptr<l1t::PFCandidate>& part) {
//     return (part->pt() / pt) * (part->eta() - seed.eta());
//   });
//   // Accumulate the pt weighted etas. Init to the seed eta, start accumulating at begin()+1 to skip seed
//   float eta = std::accumulate(pt_deta.begin() + 1, pt_deta.end(), seed.eta());

//   // pt weighted d phi
//   std::vector<float> pt_dphi;
//   pt_dphi.resize(parts.size());
//   std::transform(parts.begin(), parts.end(), pt_dphi.begin(), [&seed, &pt](const edm::Ptr<l1t::PFCandidate>& part) {
//     return (part->pt() / pt) * (part->phi() - seed.phi());
//   });
//   // Accumulate the pt weighted phis. Init to the seed phi, start accumulating at begin()+1 to skip seed
//   float phi = std::accumulate(pt_dphi.begin() + 1, pt_dphi.end(), seed.phi());

//   l1t::PFJet jet(pt, eta, phi);
//   for (auto it = parts.begin(); it != parts.end(); it++) {
//     jet.addConstituent(*it);
//   }

//   return jet;
// }

L1MetPFHWProducer::L1MetPFHWProducer(const edm::ParameterSet& cfg)
    : _maxCands(cfg.getParameter<unsigned>("maxCands")),
      _l1PFToken(consumes<std::vector<l1t::PFCandidate>>(cfg.getParameter<edm::InputTag>("src"))) {
    produces<l1t::PFCandidateCollection>();
    // produces<l1t::PFJetCollection>(); // FIXME
    // produces<reco::PFMETCollection>();
    // produces<reco::CandidateCollection>();
}

void L1MetPFHWProducer::produce(edm::StreamID /*unused*/,
                                      edm::Event& iEvent,
                                      const edm::EventSetup& iSetup) const {
  //std::unique_ptr<l1t::PFJetCollection> newPFJetCollection(new l1t::PFJetCollection);

  edm::Handle<l1t::PFCandidateCollection> l1PFCandidates;
  iEvent.getByToken(_l1PFToken, l1PFCandidates);

  float px=0;
  float py=0;
  for (unsigned i = 0; i < (*l1PFCandidates).size() && i<_maxCands; i++) {
      auto cand = (*l1PFCandidates)[i];
      // do stuff...
      px += cand.px();
      py += cand.py();
  }
  float pt = sqrt(px*px+py*py);
  float phi = atan2(py,px);

  std::cout << "MET = " << pt << " from " << (*l1PFCandidates).size() << " cands" << endl;

  // l1tk::LorentzVector p(pt*sin(phi), pt*cos(phi),0,0);
  //LorentzVector p; //(pt*sin(phi), pt*cos(phi),0,0);
  //LorentzVector p(pt*sin(phi), pt*cos(phi),0,0);

  l1t::PFMet met(pt, phi, int(pt), int(phi)); //0,0,p,0,pt,0,phi);
  auto metcoll = std::make_unique<l1t::PFMetCollection>();
  // auto metcoll = std::make_unique<l1t::PFCandidateCollection>();
  metcoll->push_back(met);
  iEvent.put(std::move(metcoll));

   //L1Candidate(p4, hwpt, hweta, hwphi, /*hwQuality=*/0), rawPt_(p4.Pt())
  
  //met.Set
  // reco::Candidate met;

  // reco::PFMET pfmet(specific, commonMETdata.sumet, p4, vtx); //...
  // auto pfmetcoll = std::make_unique<reco::PFMETCollection>();
  // pfmetcoll->push_back(pfmet);
  // event.put(std::move(pfmetcoll));


  // newPFJetCollection->swap(jets);
  // iEvent.put(std::move(newPFJetCollection));

  // std::vector<edm::Ptr<l1t::PFCandidate>> work;
  // for (unsigned i = 0; i < (*l1PFCandidates).size(); i++) {
  //   work.push_back(edm::Ptr<l1t::PFCandidate>(l1PFCandidates, i));
  // }

  // std::sort(work.begin(), work.end(), [](edm::Ptr<l1t::PFCandidate> i, edm::Ptr<l1t::PFCandidate> j) {
  //   return (i->pt() > j->pt());
  // });

  // float px=0;
  // float py=0;
  // for(int i=0; i< work.size() && i<_maxCands; i++){
  //     // met
  //     edm::Ptr<l1t::PFCandidate> seed = work.at(0);
  // }

  // std::vector<l1t::PFJet> jets;
  // jets.reserve(_nJets);

  // while (!work.empty() && jets.size() < _nJets) {
  //   // Take the first (highest pt) candidate as a seed
  //   edm::Ptr<l1t::PFCandidate> seed = work.at(0);
  //   // Get the particles within a _coneSize of the seed
  //   std::vector<edm::Ptr<l1t::PFCandidate>> particlesInCone;
  //   std::copy_if(
  //       work.begin(), work.end(), std::back_inserter(particlesInCone), [&](const edm::Ptr<l1t::PFCandidate>& part) {
  //         return reco::deltaR<l1t::PFCandidate, l1t::PFCandidate>(*seed, *part) <= _coneSize;
  //       });
  //   jets.push_back(makeJet(particlesInCone));
  //   // remove the clustered particles
  //   work.erase(std::remove_if(work.begin(),
  //                             work.end(),
  //                             [&](const edm::Ptr<l1t::PFCandidate>& part) {
  //                               return reco::deltaR<l1t::PFCandidate, l1t::PFCandidate>(*seed, *part) <= _coneSize;
  //                             }),
  //              work.end());
  // }
  // std::sort(jets.begin(), jets.end(), [](l1t::PFJet i, l1t::PFJet j) { return (i.pt() > j.pt()); });
  // newPFJetCollection->swap(jets);
  // iEvent.put(std::move(newPFJetCollection));
}

/////////////
// DESTRUCTOR
L1MetPFHWProducer::~L1MetPFHWProducer() {}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1MetPFHWProducer);
