#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>
#include "TMath.h"
#include <cstdlib>
typedef unsigned int uint;

class CCQEAntiNu2DXSec : public XSec
{
public:
  CCQEAntiNu2DXSec(const char* name)
    : XSec(name)
    {};

bool isQELikeSignal( ChainWrapper& chw, int entry ) {

  bool qelike=(bool)chw.GetValue("truth_qelike",entry);
  double proton_ke =(double)chw.GetValue("truth_leading_proton_ke",entry);
  double muon_angle=(double)chw.GetValue("truth_muon_theta_tmk",entry)*180./TMath::Pi();
  bool fiducial=(bool)chw.GetValue("truth_is_fiducial",entry);
  bool charm=(bool)chw.GetValue("mc_charm",entry);
  // CCQE-like: 1 muon (from antineutrino) and no mesons/heavy baryons in final state. 
  // Any number of final state neutrons allowed.
  // Protons are only allowed if they have kinetic energy less than 120 MeV 
  // Photons from nuclear de-excitation are kept. These tend to be < 10 MeV. Events with photons from other sources are excluded. 
  // GENIE simulates nuclear de-excitations only for Oxygen atoms at present.
  // Muon angle cut at maximum 20 degrees
   if (!charm && fiducial && qelike && proton_ke <120.&& muon_angle <20.) return true;

   //if (qelike && proton_ke <=120.) return true;
  return false;

}
  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(ChainWrapper& chw, int entry)
  {
    if((int)chw.GetValue("mc_incoming", entry)!=-14) return false;
    if((int)chw.GetValue("mc_current", entry)!=1) return false;
    if(!isQELikeSignal  ( chw, entry ) ) return false;
    
    return true;
  }
};

void runCCQEAntiNu2DXSec()
{
  // Create the XSecLooper and tell it the input files
  // The file listing below is the latest official Eroica processing for the Minerva 1 playlist 
  XSecLooper loop("/minerva/data/users/cpatrick/mc_production/eroica20160202_v_1_93/grid/central_value/minerva/ana/v10r8p8/00/05/02/*/SIM_minerva_*_CCQEAntiNuTool_Ana_Tuple_*cpatrick.root");

  //  XSecLooper loop("/pnfs/minerva/persistent/users/drut1186/MMDCCQE_VtxStudy/MMDCCQE_ana_minerva1_vtx_150_mc/grid/central_value/minerva/ana/v10r8p8/00/01/02/*/SIM_*_MinModDepCCQE_Ana_Tuple_v10r8p8-drut1186.root"); 

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(-14);

  // Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
  loop.setNumUniv(0); 

  // Add the differential cross section dsigma/dpTdpZ
  // Flux-integrated over the range 0.0 to 100.0 GeV
  CCQEAntiNu2DXSec* ds_dpTdpZ = new CCQEAntiNu2DXSec("ds_dpTdpZ");
  double pt_edges[] = { 0.0, 0.15, 0.25, 0.4, 0.7, 1.0, 1.5 };
  int pt_nbins = 6; 
  double pz_edges[] = { 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0 };
  int pz_nbins = 11;
  
  ds_dpTdpZ->setBinEdges(pz_nbins, pz_edges, pt_nbins, pt_edges);
  ds_dpTdpZ->setVariable(XSec::kPZLep, XSec::kPTLep);
  ds_dpTdpZ->setIsFluxIntegrated(true);
  ds_dpTdpZ->setDimension(2);
  ds_dpTdpZ->setFluxIntLimits(0.0, 100.0);
  ds_dpTdpZ->setNormalizationType(XSec::kPerProton);  
  loop.addXSec(ds_dpTdpZ);
 
  CCQEAntiNu2DXSec* ds_dpT = new CCQEAntiNu2DXSec("ds_dpT");
  ds_dpT->setBinEdges(pt_nbins, pt_edges);
  ds_dpT->setVariable(XSec::kPTLep);
  ds_dpT->setIsFluxIntegrated(true);
  ds_dpT->setDimension(1);
  ds_dpT->setFluxIntLimits(0.0, 100.0);
  ds_dpT->setNormalizationType(XSec::kPerProton);  
  loop.addXSec(ds_dpT);

  loop.runLoop();

  // Get the output histograms and save them to file
  string geniefilename =  "/minerva/data/users/cpatrick/CCQEAntiNu2D/GENIEXSecExtractor/GENIE_CrossSection_pzpt_lowangleqelike_CCQEAntiNu2D.root";
  TFile fout(geniefilename.c_str(), "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i){
    if(loop.getXSecs()[i]->getDimension()==1) 
      loop.getXSecs()[i]->getXSecHist()->Write();
    else
      loop.getXSecs()[i]->get2DXSecHist()->Write();
  }
}

int main()
{
  runCCQEAntiNu2DXSec();
  return 0;
}

