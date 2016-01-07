#ifndef VME_TDCTrigger_h
#define VME_TDCTrigger_h

#include "VME_TDCEvent.h"
#include <vector>
#include <map>

namespace VME
{
  class TDCTrigger
  {
    public:
      inline TDCTrigger() :
        fETTT(0), fNumEdges(0) {;}
      inline ~TDCTrigger() {;}

      typedef std::pair<unsigned int, unsigned int> Hit;
      inline void AddHit(unsigned int channel, unsigned int lead, unsigned int trail) { fEdges[channel].push_back(Hit(lead, trail)); }
      inline std::vector<Hit> GetHits(unsigned int channel) { return fEdges[channel]; }

      inline void AddError(const VME::TDCErrorFlag& err) { fErrors.push_back(err); }
      typedef std::vector<VME::TDCErrorFlag> Errors;
      inline Errors GetErrors() const { return fErrors; }

      inline void SetETTT(unsigned int ettt) { fETTT = ettt; }
      inline unsigned int GetETTT() const { return fETTT; }

    private:
      std::map<unsigned int, std::vector<Hit> > fEdges;
      Errors fErrors;
      unsigned int fETTT;
      unsigned int fNumEdges;
  };
}

#endif
