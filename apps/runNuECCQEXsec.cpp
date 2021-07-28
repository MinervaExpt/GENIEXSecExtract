#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include <PlotUtils/MnvH1D.h>

#include <cstdlib>
#include <string>

typedef unsigned int uint;

class NuECCQELikeXSec : public XSec
{
public:
  NuECCQELikeXSec(const char* name)
    : XSec(name)
  {};

  virtual ~NuECCQELikeXSec() {};

  // Specify which events to include in this selection
  bool passesCuts(ChainWrapper& chw, int entry)
  {
    // cribbed from NuECCQE::tagTruth(), http://nusoft.fnal.gov/minerva/minervadat/software_doxygen/HEAD/MINERVA/classNuECCQE.html#aa700d96eead80c8484204c10a07eb7e
    
    // need a final-state electron or positron
    if (abs(int(chw.GetValue("mc_primaryLepton", entry))) != 11)
      return false;

    // loop over the remaining FS particles
    // and ensure that they are all nucleons
    int nparticles = int(chw.GetValue("mc_nFSPart", entry));
    for(int i = 0; i < nparticles; i++)
    {
      int pid = abs( int(chw.GetValue("mc_FSPartPDG", entry, i)) );
      if (
          pid != 11
          && pid != 2212
          && pid != 2112
          && !(pid == 22 && (chw.GetValue("mc_FSPartE", entry, i) < 15))  // de-excitation photons from oxygen...  can be ignored
          && pid < 1e9
         )
        return false;
    }

    // match the goofy analysis cuts
    if (getValue(chw, entry, kENu) > 10) // || getValue(chw, entry, kELep) < 1)
      return false;

    return true;
    
//    return abs(int(chw.GetValue("mc_incoming", entry))) == 12 && int(chw.GetValue("mc_intType", entry)) == 1 && int(chw.GetValue("mc_current", entry)) == 1;
  }
};

void runNuECCQELikeXSec()
{
  // Create the XSecLooper and tell it the input files
//  XSecLooper loop("/minerva/data/users/jwolcott/coherent_ntuples_resurrection_Minerva1+7+9+13.txt");
  char * ntupleListEnv = getenv("NTUPLE_LIST");
  std::string ntupleList;
  if (ntupleListEnv == NULL)
    ntupleList = "/minerva/data/users/jwolcott/ccinclusive_ntuples_resurrection_Minerva1+13C.txt";
  else
    ntupleList = ntupleListEnv;
  XSecLooper loop(ntupleList.c_str());

  // Tell the XSecLooper we're using both nu_es and antinu_es
  loop.addNuPDG(12);
  loop.addNuPDG(-12);

  //Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
   loop.setNumUniv(0); 

  // (1) dsigma/dtheta_e
  // flux-integrated over the range 1 to 10 GeV
  NuECCQELikeXSec* ds_dtheta = new NuECCQELikeXSec("ds_theta_e");
  double theta_edges[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 22, 30 };
  int theta_nbins = sizeof(theta_edges)/sizeof(theta_edges[0]) - 1;
  ds_dtheta->setBinEdges(theta_nbins,theta_edges);
  ds_dtheta->setVariable(XSec::kThetaLep);
  ds_dtheta->setIsFluxIntegrated(true);
  ds_dtheta->setFluxIntLimits(0, 10.0);
  ds_dtheta->setNormalizationType(XSec::kPerNeutron);
  loop.addXSec(ds_dtheta);
  
  // (2) dsigma/dE_e
  // flux-integrated over the range 1 to 10 GeV
  NuECCQELikeXSec* ds_dE = new NuECCQELikeXSec("ds_E_e");
  double Ee_edges[] = { 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7.5, 10 };
  int Ee_nbins = sizeof(Ee_edges)/sizeof(Ee_edges[0]) - 1;
  ds_dE->setBinEdges(Ee_nbins,Ee_edges);
  ds_dE->setVariable(XSec::kELep);
  ds_dE->setIsFluxIntegrated(true);
  ds_dE->setFluxIntLimits(0, 10.0);
  ds_dE->setNormalizationType(XSec::kPerNeutron);
  loop.addXSec(ds_dE);


  // (3) dsigma/dQ^2_QE
  // flux-integrated over the range 1 to 10 GeV
  NuECCQELikeXSec* ds_dQ2 = new NuECCQELikeXSec("ds_Q2_QE");
  double Q2_edges[] = { 0, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 1.0, 1.25, 2.0 };
  int Q2_nbins = sizeof(Q2_edges)/sizeof(Q2_edges[0]) - 1;
  ds_dQ2->setBinEdges(Q2_nbins,Q2_edges);
  ds_dQ2->setVariable(XSec::kQ2QE);
  ds_dQ2->setIsFluxIntegrated(true);
  ds_dQ2->setFluxIntLimits(0, 10.0);
  ds_dQ2->setNormalizationType(XSec::kPerNeutron);
  loop.addXSec(ds_dQ2);

  // (4) total cross-section vs energy (CCQE hypothesis)
  NuECCQELikeXSec* sigma = new NuECCQELikeXSec("sigma");
  double Enu_edges[] = { 0, 1, 2, 3, 4, 5, 7, 10 };
  int Enu_nbins = sizeof(Enu_edges)/sizeof(Enu_edges[0]) - 1;
  sigma->setBinEdges(Enu_nbins,Enu_edges);
  sigma->setVariable(XSec::kENuQE);
  sigma->setIsFluxIntegrated(false);
  sigma->setNormalizationType(XSec::kPerNeutron);
  loop.addXSec(sigma);

  // Once everything's set up, actually run the thing
  loop.runLoop();

  // Get the output histograms and save them to file
  char * outfilenameEnv = getenv("OUT_FILE");
  std::string outfilename;
  if (outfilenameEnv == NULL)
    outfilename =  "/minerva/data/users/jwolcott/nu_e/histos/GENIE/FHC/genie_crossSections.root";
  else
    outfilename = outfilenameEnv;
  TFile fout(outfilename.c_str(), "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i)
  {
    loop.getXSecs()[i]->getXSecHist()->Write();
    if (loop.getXSecs()[i]->getEvRateHist() != NULL)
      loop.getXSecs()[i]->getEvRateHist()->Write();
  }

  // write the calculated flux too, just in case
  loop.getFluxHist()->Write();
}

int main()
{
  runNuECCQELikeXSec();
  return 0;
}

