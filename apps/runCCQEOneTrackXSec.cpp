#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include <PlotUtils/MnvH1D.h>

#include <cstdlib>
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

void runCCQEOneTrackXSec()
{
  // Create the XSecLooper and tell it the input files
  XSecLooper loop("/minerva/data/users/arturo/mc_production/grid/central_value/minerva/ana/v10r6/00/01/02/*/*.root");

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(14);

  //Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
   loop.setNumUniv(100); 

  // Add the total CCQE cross section as a function of energy
  CCQEXSec* ccqeEnu=new CCQEXSec("ccqe_enu_per_neutron");
  double enu_edges[] = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.5, 10.0 };
  int enu_nbins = 13; 
  ccqeEnu->setBinEdges(enu_nbins,enu_edges);
  ccqeEnu->setVariable(XSec::kENu);
  ccqeEnu->setIsFluxIntegrated(false);
  ccqeEnu->setNormalizationType(XSec::kPerNeutron);
  ccqeEnu->setUniverses(100);//default value, put 0 if you do not want universes to be included

  loop.addXSec(ccqeEnu);

  // Add the CCQE differential cross section dsigma/dQ^2
  // flux-integrated over the range 1.5 to 12 GeV
  CCQEXSec* ccqeQ2=new CCQEXSec("ccqe_q2_per_neutron");
  double q2_edges[] = { 0.0, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 2.0 };
  int q2_nbins = 8; 
  ccqeQ2->setBinEdges(q2_nbins,q2_edges);
  ccqeQ2->setVariable(XSec::kQ2);
  ccqeQ2->setIsFluxIntegrated(true);
  ccqeQ2->setFluxIntLimits(1.5, 12);
  ccqeQ2->setNormalizationType(XSec::kPerNeutron);
  ccqeQ2->setUniverses(100);

  loop.addXSec(ccqeQ2);

  // Add the CCQE differential cross section dsigma/dQ^2
  // flux-integrated over the range 1.5 to 12 GeV
  CCQEXSec* ccqeQ2QE=new CCQEXSec("ccqe_q2qe_per_neutron");
  ccqeQ2QE->setBinEdges(q2_nbins,q2_edges);
  ccqeQ2QE->setVariable(XSec::kQ2QE);
  ccqeQ2QE->setIsFluxIntegrated(true);
  ccqeQ2QE->setFluxIntLimits(1.5, 12);
  ccqeQ2QE->setNormalizationType(XSec::kPerNeutron);
  ccqeQ2QE->setUniverses(100);

  loop.addXSec(ccqeQ2QE);

  // Once everything's set up, actually run the thing
  loop.runLoop();

  // Get the output histograms and save them to file
  string geniefilename = Form( 
      "%s/ana/rootfiles/genie_minerva_hists_numu_v4.root", 
      getenv( "CCQUASIELASTICROOT" ) );
  TFile fout(geniefilename.c_str(), "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i){
    loop.getXSecs()[i]->getXSecHist()->Write();
  }
}

int main()
{
  runCCQEOneTrackXSec();
  return 0;
}

