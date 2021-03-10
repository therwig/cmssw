#ifndef DataFormats_L1TParticleFlow_PFMet_h
#define DataFormats_L1TParticleFlow_PFMet_h

#include <vector>
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Common/interface/Ptr.h"

namespace l1t {

  class PFMet : public L1Candidate {
  public:

    PFMet() {}
    PFMet(float pt, float phi, int hwpt = 0, int hwphi = 0)
        : L1Candidate(PolarLorentzVector(pt, 0, phi, 0), hwpt, 0, hwphi, /*hwQuality=*/0) {}

    PFMet(const LorentzVector& p4, int hwpt = 0, int hwphi = 0)
        : L1Candidate(p4, hwpt, hwphi, /*hwQuality=*/0) {}

  };

  typedef std::vector<l1t::PFMet> PFMetCollection;
  typedef edm::Ref<l1t::PFMetCollection> PFMetRef;
  typedef std::vector<l1t::PFMetRef> PFMetVectorRef;
}  // namespace l1t
#endif
