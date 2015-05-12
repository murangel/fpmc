#ifndef FpmcInterface_Fpmc_h
#define FpmcInterface_Fpmc_h

/** \class Fpmc
 *
 * Generates Fpmc HepMC events
 *
 ***************************************/

#include "HepMC/GenEvent.h"
#include <HepMC/IO_HERWIG.h>
//#include "CLHEP/Random/JamesRandom.h"
//#include "CLHEP/Random/RandFlat.h"

#include <vector>
#include <string>
#include <ostream>

namespace fpmc
{
  class Fpmc{
  public:
    Fpmc(double, long int, std::vector<std::string> const&);
    ~Fpmc();

    void begin();
    bool run();
    void end();
    void write(std::ostream&);

    const HepMC::GenEvent* event() const { return hepMCEvt_; } 

  private:

    HepMC::IO_HERWIG  conv_;
    HepMC::GenEvent  *hepMCEvt_;
 
    /// HERWIG verbosity
    unsigned int herwigVerbosity_;
    /// HepMC verbosity
    bool hepMCVerbosity_;
    /// Events to print if verbosity
    unsigned int maxEventsToPrint_;    

    unsigned int event_;
    double comEnergy_;
    // Not used temporarily (take from datacard)
    long int seed_;

    std::vector<std::string> params_;

    bool hadronize_;
  
    bool debug_;

    //CLHEP::HepRandomEngine* fRandomEngine;
    //CLHEP::RandFlat*        fRandomGenerator; 
  };
} 

#endif
