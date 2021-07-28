#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include <PlotUtils/MnvH1D.h>

typedef unsigned int uint;

class CCQEXSec : public XSec
{
public:
  CCQEXSec(const char* name)
    : XSec(name)
    {}

  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(ChainWrapper& chw, int entry)
  {
    if((int)chw.GetValue("mc_intType", entry)!=1) return false;
    if((int)chw.GetValue("mc_current", entry)!=1) return false;
    if((bool)chw.GetValue("mc_charm", entry)) return false;
    return true;
  }
};

void runXSecLooper()
{
  // Create the XSecLooper and tell it the input files
  XSecLooper loop("/minerva/data/users/rodriges/anatruth/*.root");
  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(14);
  
  //Setting the number of Universes in the GENIE error band, default 100 universes put 0 if you do not want universes to be included
  loop.setNumUniv(100); 
  
  // Add the total CCQE cross section as a function of energy
  CCQEXSec* ccqeEnu=new CCQEXSec("ccqeEnu");
  ccqeEnu->setUniformBinning(100, 0, 10);
  ccqeEnu->setVariable(XSec::kENu);
  ccqeEnu->setIsFluxIntegrated(false);
  ccqeEnu->setNormalizationType(XSec::kPerNeutron);
  ccqeEnu->setUniverses(100);//default value 100 universes, put 0 if you do not want universes to be included.

  loop.addXSec(ccqeEnu);

  // Add the CCQE differential cross section dsigma/dQ^2
  // flux-integrated over the range 1.5 to 12 GeV
  CCQEXSec* ccqeQ2=new CCQEXSec("ccqeQ2");
  ccqeQ2->setUniformBinning(50, 0, 2);
  ccqeQ2->setVariable(XSec::kQ2);
  ccqeQ2->setIsFluxIntegrated(true);
  ccqeQ2->setFluxIntLimits(1.5, 12);
  ccqeQ2->setNormalizationType(XSec::kPerNeutron);
  ccqeQ2->setUniverses(100);//default value, put 0 if you do not want universes to be included.
  
  loop.addXSec(ccqeQ2);

  // Once everything's set up, actually run the thing
  loop.runLoop();

  // Get the output histograms and save them to file
  TFile fout("newCodeTest.root", "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i){
    loop.getXSecs()[i]->getXSecHist()->Write();
  }
}

int main()
{
  runXSecLooper();
  return 0;
}

