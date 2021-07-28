#include "GENIEXSecExtract/XSecLooper.h"


#include "TFile.h"
#include "TH1.h"
#include <PlotUtils/MnvH1D.h>

#include <TVector3.h>

typedef unsigned int uint;

class CCNpMpi0Transverse : public XSec
{
public:
  CCNpMpi0Transverse(const char* name)
  : XSec(name)
  {}

  // Override this method from the base class to decide what events to
  // include in this selection
  bool isCCNpNpi0(ChainWrapper& chw, int entry)
  {
    // Phase space cuts:
    const double mu_minMom = 1500.;
    const double mu_maxMom = 20000.;
    // const double mu_minTheta = 0.;
    const double mu_maxTheta = 20.;

    const double pr_minMom = 450.;
    // static const double mu_maxMom = 20000.;

    mu_mom.clear();
    pr_mom.clear();
    pi0_mom.clear();

    int n_mu = 0;
    int n_pr = 0;
    int n_ne = 0;
    int n_pi0 = 0;
    int n_other = 0;

    bool mu_ps_pass = false;
    bool pr_ps_pass = false;
    // bool pi0_ps_pass = false;

    int nparticles = (int)chw.GetValue("mc_nFSPart",entry);
    for(int i = 0; i < nparticles; i++){
      int pdg = (int)chw.GetValue("mc_FSPartPDG",entry,i);
      double energy = (double)chw.GetValue("mc_FSPartE",entry,i);

      double px = (double)chw.GetValue("mc_FSPartPx",entry,i);
      double py = (double)chw.GetValue("mc_FSPartPy",entry,i);
      double pz = (double)chw.GetValue("mc_FSPartPz",entry,i);      

      TVector3 pmom(px, py, pz);
      double theta = getTheta( px, py, pz );

      if(pdg == 13){
        n_mu++;
        mu_mom.push_back(pmom);
        if(theta < mu_maxTheta && mu_minMom <= pmom.Mag() && pmom.Mag() <= mu_maxMom) mu_ps_pass = true;
      }
      else if (pdg == 2212){
        n_pr++;
        if(pr_minMom <= pmom.Mag()) pr_ps_pass = true;
        pr_mom.push_back(pmom);
      }
      else if (pdg == 111){
        n_pi0++;
        pi0_mom.push_back(pmom);
      }
      else if (pdg == 2111) n_ne++;
      else n_other++;
    }

    SortParticles();

    return (n_mu > 0 && n_pr > 0 && n_pi0 > 0 && n_other == 0 && mu_ps_pass && pr_ps_pass);
  }

  virtual bool passesCuts(ChainWrapper& chw, int entry)
  {
    if((int)chw.GetValue("mc_incoming", entry)!=14) return false;
    if((int)chw.GetValue("mc_current", entry)!=1) return false;
    if(!isCCNpNpi0( chw, entry ) ) return false;
    return true;
  }

 };

 void runXSecLooper(const std::string &infilename)
 {

  XSecLooper loop(infilename.c_str());
  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(14);
  
  //Setting the number of Universes in the GENIE error band, default 100 universes put 0 if you do not want universes to be included
  // loop.setNumUniv(100); 
  loop.setNumUniv(0); 
  
   // Binning:
  const int nbins = 7;
  double binning[ nbins + 1 ] = {  -400., -200., -110., -30., 30., 110., 200., 400. };

  // Add the total CCQE cross section as a function of energy
  CCNpMpi0Transverse * ccNpMpi0Transverse = new CCNpMpi0Transverse("ccNpMpi0Transverse");
  ccNpMpi0Transverse->setBinEdges(nbins, binning);
  ccNpMpi0Transverse->setVariable(XSec::kTdpTT);
  // I guess this is in GeV:
  ccNpMpi0Transverse->setFluxIntLimits(0.0,50.);

  ccNpMpi0Transverse->setIsFluxIntegrated(true);
  ccNpMpi0Transverse->setNormalizationType(XSec::kPerNucleon);
  ccNpMpi0Transverse->setUniverses(0);//default value 100 universes, put 0 if you do not want universes to be included.

  loop.addXSec(ccNpMpi0Transverse);

  // Once everything's set up, actually run the thing
  loop.runLoop();

  // Get the output histograms and save them to file
  TFile fout("newCodeTest.root", "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i){
    loop.getXSecs()[i]->getXSecHist()->Write();
  }
}

int main(int argc, char *argv[])
{
  string ifilename = "/pnfs/minerva/persistent/users/dcoplowe/CCProtonPi0/MC/";
  bool rg2g = false;

  char cc;
    while((cc = getopt(argc, argv, "t::r:")) != -1){
      switch(cc){
        case 't':
        {
            ifilename += "run190618/minerva1_repro/grid/central_value/minerva/ana/v10r8p12/00/01/02/00/";
            ifilename += "SIM_minerva_00010200_Subruns_0001_CCProtonPi0_Ana_Tuple_v10r8p12.root";
            rg2g = true;
            break;
        }
        case 'r':{
          ifilename += std::string(optarg);
          ifilename =+ "/*/grid/central_value/minerva/ana/*/*/*/*/*/SIM*.root";
          rg2g = true;
          break;
        }
        default:
        {
          break;
        }
      }
    }

    if(rg2g){
      cout << "Producing MC XSec Result using file named:" << endl;
      cout << "    " << ifilename << endl;
      runXSecLooper(ifilename);
    }

  return 0;
}

