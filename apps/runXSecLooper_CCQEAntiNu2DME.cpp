#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>
#include "TMath.h"
#include <cstdlib>
typedef unsigned int uint;

class CCQEAntiNu2DXSecME : public XSec
{
public:
  CCQEAntiNu2DXSecME(const char* name)
    : XSec(name)
    {};

bool isQELikeSignal( ChainWrapper& chw, int entry ) {


  int genie_n_muons         = 0;
  int genie_n_mesons        = 0;
  int genie_n_heavy_baryons_plus_pi0s = 0;
  int genie_n_photons       = 0;
  int genie_n_protons       = 0;

  int nparticles = (int)chw.GetValue("mc_nFSPart",entry);
  for(int i = 0; i < nparticles; ++i) {
     int pdg = (int)chw.GetValue("mc_FSPartPDG",entry,i);
     double energy = (double)chw.GetValue("mc_FSPartE",entry,i);

    if( abs(pdg) == 13 ) genie_n_muons++;
    else if( pdg == 22 && energy > 10 ) genie_n_photons++;
    else if( abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 323 || pdg == 111 || pdg == 130 || pdg == 310 || pdg == 311 || pdg == 313 ) genie_n_mesons++;
    else if( pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 4112 || pdg == 4122 || pdg == 4212 || pdg == 4222 || pdg == 411 || pdg == 421 || pdg == 111 ) genie_n_heavy_baryons_plus_pi0s++;
    else if( pdg == 2212 && energy > 1058.272 ) genie_n_protons++; //antinu
  }

  double theta              = 0.;
  double true_muon_px   = (double)chw.GetValue("mc_primFSLepton",entry,0);
  double true_muon_py   = (double)chw.GetValue("mc_primFSLepton",entry,1);
  double true_muon_pz   = (double)chw.GetValue("mc_primFSLepton",entry,2);
  double numi_beam_angle_rad = -0.05887;
  double pyprime = -1.0*sin(numi_beam_angle_rad)*true_muon_pz + cos(numi_beam_angle_rad)*true_muon_py;
  double pzprime =  1.0*cos(numi_beam_angle_rad)*true_muon_pz + sin(numi_beam_angle_rad)*true_muon_py;
  double pSquare = pow(true_muon_px,2) + pow(pyprime,2) + pow(pzprime,2);
  theta = acos( pzprime / sqrt(pSquare) );
  theta *= 180./3.14159;
    if( genie_n_muons         == 1 && 
	genie_n_mesons        == 0 && 
	genie_n_heavy_baryons_plus_pi0s == 0 && 
	genie_n_photons       == 0 &&
	genie_n_protons        == 0 &&
	theta<=20.0) return true;
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

void runCCQEAntiNu2DXSecME()
{
  // Create the XSecLooper and tell it the input files
  // The file listing below is the latest official Eroica processing for the Minerva 1 playlist 
  XSecLooper loop("/pnfs/minerva/persistent/users/drut1186/CCQENu_Anatuples/MuonKludge_MichelFit_BlobFit/MC_Merged/minervame6Gpass1//*.root");

  //  XSecLooper loop("/pnfs/minerva/persistent/users/drut1186/MMDCCQE_VtxStudy/MMDCCQE_ana_minerva1_vtx_150_mc/grid/central_value/minerva/ana/v10r8p8/00/01/02/*/SIM_*_MinModDepCCQE_Ana_Tuple_v10r8p8-drut1186.root"); 

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(-14);

  // Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
  loop.setNumUniv(0); 

  // Add the differential cross section dsigma/dpTdpZ
  // Flux-integrated over the range 0.0 to 100.0 GeV
  double pt_edges[] = {0.0,0.075,0.15,0.25,0.325,0.4,0.475,0.55,0.7,0.85,1.0,1.25,1.5,2.5};
  int pt_nbins = 13; 
  double pz_edges[] = {1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,15.0};
  int pz_nbins = 13;
  double q2_edges[] =  {0.0, 0.0125, 0.025,0.05, 0.1,  0.2, 0.4,  0.8, 1.2, 4.0};
  int q2_nbins = 9; 
  
  double enu_edges[]= {2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,15.0};
  int enu_nbins = 12;
  
    
  // Flux-integrated over the range 0.0 to 100.0 GeV
  CCQEAntiNu2DXSecME* ds_dpTdpZ = new CCQEAntiNu2DXSecME("ds_dpTdpZ");
  ds_dpTdpZ->setBinEdges(pz_nbins, pz_edges, pt_nbins, pt_edges);
  ds_dpTdpZ->setVariable(XSec::kPZLep, XSec::kPTLep);
  ds_dpTdpZ->setIsFluxIntegrated(true);
  ds_dpTdpZ->setDimension(2);
  ds_dpTdpZ->setFluxIntLimits(0.0, 100.0);
  ds_dpTdpZ->setNormalizationType(XSec::kPerNucleon);  
  loop.addXSec(ds_dpTdpZ);


  CCQEAntiNu2DXSecME* ds_dpTdQ2 = new CCQEAntiNu2DXSecME("ds_dpTdQ2");
  ds_dpTdQ2->setBinEdges(q2_nbins, q2_edges, pt_nbins, pt_edges);
  ds_dpTdQ2->setVariable(XSec::kQ2QE, XSec::kPTLep);
  ds_dpTdQ2->setIsFluxIntegrated(true);
  ds_dpTdQ2->setDimension(2);
  ds_dpTdQ2->setFluxIntLimits(0.0, 100.0);
  ds_dpTdQ2->setNormalizationType(XSec::kPerNucleon);  
  loop.addXSec(ds_dpTdQ2);
  
  CCQEAntiNu2DXSecME* ds_dQ2dEnu = new CCQEAntiNu2DXSecME("ds_dQ2dEnu");
  ds_dQ2dEnu->setBinEdges(enu_nbins,enu_edges,q2_nbins,q2_edges);
  ds_dQ2dEnu->setIsFluxIntegrated(true);
  ds_dQ2dEnu->setDimension(2);
  ds_dQ2dEnu->setFluxIntLimits(0.0, 100.0);
  ds_dQ2dEnu->setNormalizationType(XSec::kPerNucleon);  
  loop.addXSec(ds_dQ2dEnu);  
  

/* 
  MinModDepCCQEXSec* ds_dpT = new CCQEAntiNu2DXSecME("ds_dpT");
  ds_dpT->setBinEdges(pt_nbins, pt_edges);
  ds_dpT->setVariable(XSec::kPTLep);
  ds_dpT->setIsFluxIntegrated(true);
  ds_dpT->setDimension(1);
  ds_dpT->setFluxIntLimits(0.0, 100.0);
  ds_dpT->setNormalizationType(XSec::kPerNeutron);  
  loop.addXSec(ds_dpT);

  MinModDepCCQEXSec* ds_dpT2 = new CCQEAntiNu2DXSecME("ds_dpT2");
  ds_dpT2->setBinEdges(pt_nbins, pt_edges);
  ds_dpT2->setVariable(XSec::kPTLepSquare);
  ds_dpT2->setIsFluxIntegrated(true);
  ds_dpT2->setDimension(1);
  ds_dpT2->setFluxIntLimits(0.0, 10.0);
  ds_dpT2->setNormalizationType(XSec::kPerNeutron);  
  loop.addXSec(ds_dpT2);
*/

  CCQEAntiNu2DXSecME* ds_dQ2 = new CCQEAntiNu2DXSecME("ds_dQ2");
  ds_dQ2->setBinEdges(q2_nbins, q2_edges);
  ds_dQ2->setVariable(XSec::kQ2QE);  
  ds_dQ2->setIsFluxIntegrated(true);
  ds_dQ2->setDimension(1);
  ds_dQ2->setFluxIntLimits(0.0, 100.0);
  ds_dQ2->setNormalizationType(XSec::kPerNucleon);  
  loop.addXSec(ds_dQ2);
  
  
  loop.runLoop();

  // Get the output histograms and save them to file
  string geniefilename =  "GENIEXSECEXTRACT_CCQEAntiNu_me6A_NonueConstraint_ver2.root";
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
  TH1::AddDirectory(false);
  runCCQEAntiNu2DXSecME();
  return 0;
}

