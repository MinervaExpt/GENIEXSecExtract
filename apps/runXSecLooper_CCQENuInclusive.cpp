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
    :XSec(name)
  {
  };

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

  // CCQE-like: 1 muon (from neutrino) and no mesons/heavy baryons in final state. 
  // Any number of final state nucleons (protons or neutrons) allowed. 
  // Photons from nuclear de-excitation are kept. These tend to be < 10 MeV. Events with photons from other sources are excluded. 
  // GENIE simulates nuclear de-excitations only for Oxygen atoms at present. 
  if(!chw.GetValue("truth_is_fiducial",entry)) return false;
  // if( genie_n_muons         == 1 &&
  //     genie_n_mesons        == 0 &&
  //     genie_n_heavy_baryons_plus_pi0s == 0 &&
  //     genie_n_photons       == 0 &&
  //     theta <=20.0 &&
  //     pzprime >= 1.5) return true;
  if( pzprime >= 1.5 && theta <= 20.0 ) return true;
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
  //  XSecLooper loop("/pnfs/minerva/persistent/users/mateusc/CCQENu_v21r1p1_2019_minervame1D_MC_ntuples_ImprovedMichelUpdated_NEW_merged/CCQENu_mc_AnaTuple_*_Playlist.root");
    XSecLooper loop("/pnfs/minerva/persistent/users/drut1186/CCQENu_Anatuples/MuonKludge_SmallerMuonScaleUnc_ContainmentAndProtonVariables/MC_Merged/minervame1Lpass1/*.root");

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
 
  vector<double> recoil3Dbins;
  vector<double> pt3Dbins;
  vector<double> pz3Dbins;

  pt3Dbins.push_back(0.0);
  //  pt3Dbins.push_back(0.075);//added ME
  pt3Dbins.push_back(0.15);
  pt3Dbins.push_back(0.25);//added ME
  pt3Dbins.push_back(0.4);
  pt3Dbins.push_back(0.7);//added ME
  pt3Dbins.push_back(1.0);
  //  pt3Dbins.push_back(1.5);//added ME
  pt3Dbins.push_back(2.5);

  pz3Dbins.push_back(1.5);
  pz3Dbins.push_back(3.5);//added ME
  pz3Dbins.push_back(8.0);
  pz3Dbins.push_back(20.0);

  for(int i=0;i<10;i++)recoil3Dbins.push_back(i*40);
  for(int i=0;i<4;i++)recoil3Dbins.push_back(i*200+400);


  std::vector<std::vector<double> > full3D;
  full3D.push_back(recoil3Dbins);
  full3D.push_back(pt3Dbins);
  full3D.push_back(pz3Dbins);

  cout << pt3Dbins.size() << endl;
  int pzrecoil_nbins = (recoil3Dbins.size()+1)*(pz3Dbins.size()+1);
  double pzrecoil_edges[pzrecoil_nbins];
  double pt3D_edges[7];
  for(int i=0;i<pzrecoil_nbins;i++){
    pzrecoil_edges[i] = i;
  }
  for(int i=0;i<7;i++)pt3D_edges[i]=pt3Dbins[i];

  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dpTdpZ = new MinModDepCCQEXSec("ds_dpTdpZ");
  ds_dpTdpZ->setBinEdges(pz_nbins, pz_edges, pt_nbins, pt_edges);
  ds_dpTdpZ->setVariable(XSec::kPZLep, XSec::kPTLep);
  ds_dpTdpZ->setIsFluxIntegrated(true);
  ds_dpTdpZ->setDimension(2);
  ds_dpTdpZ->setFluxIntLimits(0.0, 100.0);
  ds_dpTdpZ->setNormalizationType(XSec::kPerNucleon);  
  loop.addXSec(ds_dpTdpZ);
  MinModDepCCQEXSec* ds_dpTdpZdRecoil = new MinModDepCCQEXSec("ds_dpTdpZdRecoil");
  ds_dpTdpZdRecoil->setBinEdges(pzrecoil_nbins-1, pzrecoil_edges, pt3Dbins.size()-1,&pt3Dbins[0]);
  ds_dpTdpZdRecoil->setVariable(XSec::kPZRecoil, XSec::kPTLep);
  ds_dpTdpZdRecoil->setIsFluxIntegrated(true);
  ds_dpTdpZdRecoil->setDimension(2);
  ds_dpTdpZdRecoil->setFluxIntLimits(0.0, 100.0);
  ds_dpTdpZdRecoil->setNormalizationType(XSec::kPerNucleon);  
  ds_dpTdpZdRecoil->setHyperDim(full3D,0);
  loop.addXSec(ds_dpTdpZdRecoil);
  

  loop.runLoop();

  // Get the output histograms and save them to file
  string geniefilename =  "GENIEXSECEXTRACT_CCQENuInclusive_me1L.root";
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

