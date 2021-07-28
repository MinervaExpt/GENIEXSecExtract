#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include <PlotUtils/MnvH1D.h>

#include <cstdlib>
typedef unsigned int uint;

class CCNuPionIncXSec : public XSec
{
public:
  CCNuPionIncXSec(const char* name, bool do1PiAna = true)
    : XSec(name),
    m_do1PiAna(do1PiAna)
    {};

  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(ChainWrapper& chw, int entry)
  {
    if((int)chw.GetValue("mc_incoming", entry)!=14) return false;
    if((int)chw.GetValue("mc_current", entry)!=1) return false;
    if((double)chw.GetValue("mc_w",entry) >= 1800.0) return false;
    
    int n_qpi = (int)chw.GetValue("truth_N_pip", entry) + (int)chw.GetValue("truth_N_pim", entry);
    if (n_qpi ==0) return false;
    
    if (m_do1PiAna) {
      if (n_qpi != 1) return false;
      if ((double)chw.GetValue("mc_w",entry) >= 1400.0) return false;
    }
    
    return true;
  }
  
  //input parameters
  bool m_do1PiAna;
};

void runNuECCQELikeXSec()
{
  // Create the XSecLooper and tell it the input files
  XSecLooper loop("/minerva/data/users/eberly/processing_nov26/mc_ntuples/*.root");

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(14);

  //Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
   loop.setNumUniv(0); 

  // Add the differential cross section dsigma/theta
  // flux-integrated over the range 1.5 to 10 GeV
  CCNuPionIncXSec* ds_dtheta = new CCNuPionIncXSec("ds_pi_theta",true);
  double theta_edges[] = { 0.0, 15.0, 22.0, 29.0, 36.0, 43.0, 50.0, 57.0, 72.0, 108.0, 130.0, 140.0, 150.0, 165.0 };
  int theta_nbins = 13; 
  ds_dtheta->setBinEdges(theta_nbins,theta_edges);
  ds_dtheta->setVariable(XSec::kThetaChargedPion);
  ds_dtheta->setIsFluxIntegrated(true);
  ds_dtheta->setFluxIntLimits(1.5, 10.0);
  ds_dtheta->setNormalizationType(XSec::kPerNucleon);  
  loop.addXSec(ds_dtheta);
  
  // Add the differential cross section dsigma/TPi
  // flux-integrated over the range 1.5 to 10 GeV
  CCNuPionIncXSec* ds_dKE = new CCNuPionIncXSec("ds_pi_ke",true);
  double KE_edges[] = { 60.0, 85.0, 105.0, 130.0, 160.0, 210.0, 300.0 };
  int KE_nbins = 6; 
  ds_dKE->setBinEdges(KE_nbins,KE_edges);
  ds_dKE->setVariable(XSec::kTChargedPion);
  ds_dKE->setIsFluxIntegrated(true);
  ds_dKE->setFluxIntLimits(1.5, 10.0);
  ds_dKE->setNormalizationType(XSec::kPerNucleon);
  loop.addXSec(ds_dKE);

  // Once everything's set up, actually run the thing
  loop.runLoop();


  // Get the output histograms and save them to file
  string geniefilename =  "/minerva/data/users/eberly/analysis/histos/genie_crossSections.root";
  TFile fout(geniefilename.c_str(), "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i){
    loop.getXSecs()[i]->getXSecHist()->Write();
  }
}

int main()
{
  runNuECCQELikeXSec();
  return 0;
}

