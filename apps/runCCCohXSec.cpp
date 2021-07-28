#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include <PlotUtils/MnvH1D.h>

#include <cstdlib>
typedef unsigned int uint;

class CCCohXSec : public XSec
{
public:
  CCCohXSec(const char* name)
    : XSec(name)
    {}

  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(ChainWrapper& chw, int entry)
  {
    if((int)chw.GetValue("mc_intType", entry)!=4) return false;
    if((int)chw.GetValue("mc_current", entry)!=1) return false;
    return true;
  }
};

void runCCCohXSec()
{
  // Create the XSecLooper and tell it the input files
  XSecLooper loop("/minerva/data/users/higuera/mc_production/grid/central_value/minerva/ana/v10r7p1/00/00/00/04/*root");

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(14);

  //Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
   loop.setNumUniv(100); 

  // Add the differential cross section dsigma/theta
  // flux-integrated over the range 1.5 to 20 GeV
  CCCohXSec* ds_dtheta=new CCCohXSec("ds_dthetapi");
  double theta_edges[] = { 0.0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
                     60, 70, 80, 90  };
  int theta_nbins = 14; 
  ds_dtheta->setBinEdges(theta_nbins,theta_edges);
  ds_dtheta->setVariable(XSec::kThetaPiPlus);
  ds_dtheta->setIsFluxIntegrated(true);
  ds_dtheta->setFluxIntLimits(1.5, 20);
  ds_dtheta->setNormalizationType(XSec::kPerNucleon);
  
  loop.addXSec(ds_dtheta);

  // Once everything's set up, actually run the thing
  loop.runLoop();


  // Get the output histograms and save them to file
  string geniefilename =  "test.root";
  TFile fout(geniefilename.c_str(), "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i){
    loop.getXSecs()[i]->getXSecHist()->Write();
  }
}

int main()
{
  runCCCohXSec();
  return 0;
}

