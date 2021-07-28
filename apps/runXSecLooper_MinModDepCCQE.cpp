#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

#include <cstdlib>
typedef unsigned int uint;

class MinModDepCCQEXSec : public XSec
{
public:
  MinModDepCCQEXSec(const char* name)
    : XSec(name)
    {};

bool isQELikeSignal( ChainWrapper& chw, int entry ) {

  int genie_n_muons         = 0;
  int genie_n_mesons        = 0;
  int genie_n_heavy_baryons_plus_pi0s = 0;
  int genie_n_photons       = 0;

  int nparticles = (int)chw.GetValue("mc_nFSPart",entry);
  for(int i = 0; i < nparticles; ++i) {
     int pdg = (int)chw.GetValue("mc_FSPartPDG",entry,i);
     double energy = (double)chw.GetValue("mc_FSPartE",entry,i);

     if( abs(pdg) == 13 ) genie_n_muons++;
    else if( pdg == 22 && energy > 10 ) genie_n_photons++;
    else if( abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 323 || pdg == 111 || pdg == 130 || pdg == 310 || pdg == 311 || pdg == 313 ) genie_n_mesons++;
    else if( pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 4112 || pdg == 4122 || pdg == 4212 || pdg == 4222 || pdg == 411 || pdg == 421 || pdg == 111 ) genie_n_heavy_baryons_plus_pi0s++;
  }

  double theta              = 0.;
  double true_muon_px   = (double)chw.GetValue("mc_primFSLepton",entry,0)/1000;
  double true_muon_py   = (double)chw.GetValue("mc_primFSLepton",entry,1)/1000;
  double true_muon_pz   = (double)chw.GetValue("mc_primFSLepton",entry,2)/1000;
  double numi_beam_angle_rad = -0.05887;
  double pyprime = -1.0*sin(numi_beam_angle_rad)*true_muon_pz + cos(numi_beam_angle_rad)*true_muon_py;
  double pzprime =  1.0*cos(numi_beam_angle_rad)*true_muon_pz + sin(numi_beam_angle_rad)*true_muon_py;
  double pSquare = pow(true_muon_px,2) + pow(pyprime,2) + pow(pzprime,2);
  theta = acos( pzprime / sqrt(pSquare) );
  theta *= 180./3.14159;

  // int crazy_wgt             = 0;
  // double MFP_N_2  = (double)chw.GetValue("truth_genie_wgt_MFP_N",entry,2);
  // double MFP_N_4  = (double)chw.GetValue("truth_genie_wgt_MFP_N",entry,4);
  // double MFP_pi_2 = (double)chw.GetValue("truth_genie_wgt_MFP_pi",entry,2);
  // double MFP_pi_4 = (double)chw.GetValue("truth_genie_wgt_MFP_pi",entry,4);
   
  // if(MFP_N_2 > 1000. || MFP_N_4 > 1000. || MFP_pi_2 > 1000. || MFP_pi_4 > 1000.) crazy_wgt++; 
 

  //double theta = acos(pzprime/sqrt(pSquare))*180/3.14159;

  // CCQE-like: 1 muon (from neutrino) and no mesons/heavy baryons in final state. 
  // Any number of final state nucleons (protons or neutrons) allowed. 
  // Photons from nuclear de-excitation are kept. These tend to be < 10 MeV. Events with photons from other sources are excluded. 
  // GENIE simulates nuclear de-excitations only for Oxygen atoms at present. 
  if( genie_n_muons         == 1 &&
      genie_n_mesons        == 0 &&
      genie_n_heavy_baryons_plus_pi0s == 0 &&
      genie_n_photons       == 0 &&
      theta <=20.0 &&
      true_muon_pz >= 1.5) return true;
      
  return false;

}
  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(ChainWrapper& chw, int entry)
  {
    if((int)chw.GetValue("mc_incoming", entry)!=14) return false;
    if((int)chw.GetValue("mc_current", entry)!=1) return false;
    if(!isQELikeSignal  ( chw, entry ) ) return false;
    
    return true;
  }
};

void runMinModDepCCQEXSec()
{
  TH1::AddDirectory(false);
  // Create the XSecLooper and tell it the input files
  // Inputs should be the merged ntuples:
  //  XSecLooper loop("/pnfs/minerva/persistent/users/mateusc/CCQENu_v21r1p1_2019_minervame1M_MC_ntuples_ImprovedMichelUpdated_NEW_merged/CCQENu_mc_AnaTuple_*_Playlist.root");
  //  XSecLooper loop("/pnfs/minerva/persistent/users/drut1186/ME_CCQELike_ImprovedMichel_Merged/minervame1Bv1/*.root");
    XSecLooper loop("/pnfs/minerva/persistent/users/drut1186/HopefullyFinalMateusTuples/MC_Merged/minervame1Dpass1/*.root");

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(14);

  // Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
  loop.setNumUniv(0); 

  // Add the differential cross section dsigma/dpTdpZ

  double pt_edges[] = { 0.0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 1.25, 1.5, 2.5 };
  int pt_nbins = 13; 
  double q2_edges[] =  {0, 0.00625, 0.0125, 0.025, 0.0375, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 2.0, 4.0, 6.0, 8.0,10.0};
  int q2_nbins = 19;
  double pz_edges[] = { 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0 };
  int pz_nbins = 12;
 
  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dpTdpZ = new MinModDepCCQEXSec("ds_dpTdpZ");
  ds_dpTdpZ->setBinEdges(pz_nbins, pz_edges, pt_nbins, pt_edges);
  ds_dpTdpZ->setVariable(XSec::kPZLep, XSec::kPTLep);
  ds_dpTdpZ->setIsFluxIntegrated(true);
  ds_dpTdpZ->setDimension(2);
  ds_dpTdpZ->setFluxIntLimits(0.0, 100.0);
  ds_dpTdpZ->setNormalizationType(XSec::kPerNucleon);  
  loop.addXSec(ds_dpTdpZ);

  MinModDepCCQEXSec* ds_dpTdQ2 = new MinModDepCCQEXSec("ds_dpTdQ2");
  ds_dpTdQ2->setBinEdges(q2_nbins, q2_edges, pt_nbins, pt_edges);
  ds_dpTdQ2->setVariable(XSec::kQ2QE, XSec::kPTLep);
  ds_dpTdQ2->setIsFluxIntegrated(true);
  ds_dpTdQ2->setDimension(2);
  ds_dpTdQ2->setFluxIntLimits(0.0, 100.0);
  ds_dpTdQ2->setNormalizationType(XSec::kPerNucleon);  
  loop.addXSec(ds_dpTdQ2);

/* 
  MinModDepCCQEXSec* ds_dpT = new MinModDepCCQEXSec("ds_dpT");
  ds_dpT->setBinEdges(pt_nbins, pt_edges);
  ds_dpT->setVariable(XSec::kPTLep);
  ds_dpT->setIsFluxIntegrated(true);
  ds_dpT->setDimension(1);
  ds_dpT->setFluxIntLimits(0.0, 100.0);
  ds_dpT->setNormalizationType(XSec::kPerNucleon);  
  loop.addXSec(ds_dpT);

  MinModDepCCQEXSec* ds_dpT2 = new MinModDepCCQEXSec("ds_dpT2");
  ds_dpT2->setBinEdges(pt_nbins, pt_edges);
  ds_dpT2->setVariable(XSec::kPTLepSquare);
  ds_dpT2->setIsFluxIntegrated(true);
  ds_dpT2->setDimension(1);
  ds_dpT2->setFluxIntLimits(0.0, 10.0);
  ds_dpT2->setNormalizationType(XSec::kPerNucleon);  
  loop.addXSec(ds_dpT2);
*/

  MinModDepCCQEXSec* ds_dQ2 = new MinModDepCCQEXSec("ds_dQ2");
  ds_dQ2->setBinEdges(q2_nbins, q2_edges);
  ds_dQ2->setVariable(XSec::kQ2QE);
  ds_dQ2->setIsFluxIntegrated(true);
  ds_dQ2->setDimension(1);
  ds_dQ2->setFluxIntLimits(0.0, 100.0);
  ds_dQ2->setNormalizationType(XSec::kPerNucleon);  
  loop.addXSec(ds_dQ2);

/*
  vector<XSec::EVariable> vars;
  vars.push_back( XSec::kENu );
  varName = "Enu"
*/

  loop.runLoop();

  // Get the output histograms and save them to file
  string geniefilename =  "GENIEXSECEXTRACT_CCQENu_me1B.root";
  TFile fout(geniefilename.c_str(), "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i){
    if(loop.getXSecs()[i]->getDimension()==1){ 
      loop.getXSecs()[i]->getXSecHist()->Write();
      loop.getXSecs()[i]->getEvRateHist()->Write();}
    else
      loop.getXSecs()[i]->get2DXSecHist()->Write();
  }
}

int main()
{
  TH1::AddDirectory(false);
  runMinModDepCCQEXSec();
  return 0;
}

