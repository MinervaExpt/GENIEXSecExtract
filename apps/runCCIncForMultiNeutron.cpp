//Copied by aolivier@ur.rochester.edu from code presumably written by Phil Rodriges, Rik Gran, Minerba Betancort, and maybe others.
//Cintex is only needed for older ROOT versions like the GPVMs.
//Let CMake decide whether it's needed.
#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "Math/AxisAngle.h"
#include "Math/Vector3D.h"

#include "TH1.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

#include <cstdlib>

typedef unsigned int uint;

template<class T> T sqr(T x) { return x*x; }

class CCIncXSec : public XSec
{
public:

  CCIncXSec(const char* name, const double maxEAvail)
    : XSec(name), m_max_E_avail(maxEAvail)
    {}

  bool passesMuKinematics(ChainWrapper& chw, int entry)
  {
    ROOT::Math::AxisAngle toBeamFrame(ROOT::Math::XYZVector(1., 0., 0.), -0.05887); //NuMI beam angle in mrad from PlotUtils
    ROOT::Math::XYZVector muon(chw.GetValue("mc_primFSLepton", entry, 0), //In MeV
                               chw.GetValue("mc_primFSLepton", entry, 1),
                               chw.GetValue("mc_primFSLepton", entry, 2));
    const double thetamu = (toBeamFrame * muon).theta();
    if(thetamu>=20.*M_PI/180.) return false; //This is a YAML parameter

    const double pMu = sqrt(muon.Mag2());
    if(pMu <= 2.0e3 || pMu >= 20.0e3) return false;

    //Has the <= and >= backwards compared to my analysis
    /*const double Emu=1e-3*chw.GetValue("mc_primFSLepton", entry, 3);
    // The Emu corresponding to pmu=1.5
    const double EmuMinCut=sqrt( sqr(2.0) + sqr(0.1056) ), //These are YAML parameters
                 EmuMaxCut=sqrt( sqr(20.0) + sqr(0.1056) );;
    if(Emu<EmuMinCut) return false;
    if(Emu>EmuMaxCut) return false;*/

    return true;
  }

  virtual bool passesCuts(ChainWrapper& chw, int entry)
  {
    //2+ neutrons above KE threshold
    const double neutronKEThreshold = 10; //MeV

    //Low Eavailable
    double E_avail = 0.; //In MeV
    const int nFS = chw.GetValue("mc_nFSPart", entry);
    int nNeutrons = 0;
    for(int whichFS = 0; whichFS < nFS; ++whichFS)
    {
      const int PDGCode = chw.GetValue("mc_FSPartPDG", entry, whichFS);
      const double energy = chw.GetValue("mc_FSPartE", entry, whichFS); //In MeV according to NucCCNeutrons code
      if(abs(PDGCode) == 211) E_avail += energy - 139.57;
      else if(PDGCode == 2212) E_avail += energy - 938.28;
      else if(PDGCode == 111) E_avail += energy;
      else if(PDGCode == 22) E_avail += energy;
      else if(PDGCode == 2112 && energy - 939.6 > neutronKEThreshold) ++nNeutrons;
      //Implicitly exclude neutrons, nuclei, kaons, and heavy baryons
    }

    if(nNeutrons < 2) return false;

    E_avail = std::max(0.0, E_avail);
                                                                        
    if(!(E_avail > -1 && E_avail < m_max_E_avail*1.e3)) return false;

    if(!passesMuKinematics(chw, entry)) return false;

    if((int)chw.GetValue("mc_current", entry)!=1) return false;

    return true;
  }

protected:
  double m_max_E_avail;
};

int runXSecLooper(const bool antinu, const double EAvailMax, const double Emin, const double Emax, const std::vector<const char*> fileNames)
{
  const char* fileName = "GENIEXSecExtract_NeutronMultiplicity_Tracker.root";
  auto outFile = TFile::Open(fileName, "CREATE");
  if(!outFile)
  {
    std::cerr << "GENIEXSecExtract failed to open a file for output named "
              << fileName << ".  Aborting the cross section loop...\n";
    return 1;
  }

  // Create the XSecLooper and tell it the input files
  if(fileNames.empty())
  {
    std::cerr << "No files for GENIEXSecExtract to read.  Aborting the cross section loop...\n";
    return 2;
  }

  XSecLooper loop(fileNames[0]);
  for(auto whichFile = fileNames.begin()+1; whichFile != fileNames.end(); ++whichFile)
  {
    loop.addFiles(*whichFile);
  }

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(antinu ? -14: 14);

  //Setting the number of Universes in the GENIE error band, default 100 universes put 0 if you do not want universes to be included
  loop.setNumUniv(0); 

  loop.setFiducial(5890, 8467);
  loop.setPlaylist(PlotUtils::FluxReweighter::minervame6A);

  const int nBins=22;
  double bins[nBins+1]={0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.5, 2.0, 2.5};

  CCIncXSec* CCMultiNeutrons=new CCIncXSec(TString::Format("cc_NeutronMultiplicity"), EAvailMax);
  CCMultiNeutrons->setVariable(XSec::kPTLep);
  CCMultiNeutrons->setDimension(1);
  CCMultiNeutrons->setBinEdges(nBins, bins);
  CCMultiNeutrons->setIsFluxIntegrated(true);
  CCMultiNeutrons->setFluxIntLimits(Emin, Emax);
  CCMultiNeutrons->setNormalizationType(XSec::kPerNucleon);
  CCMultiNeutrons->setUniverses(0);//default value, put 0 if you do not want universes to be included.

  loop.addXSec(CCMultiNeutrons);

  // Once everything's set up, actually run the thing
  loop.runLoop();

  // Get the output histograms and save them to file
  outFile->cd();
  for(uint i=0; i<loop.getXSecs().size(); ++i){
    loop.getXSecs()[i]->getXSecHist()->Write();
    loop.getXSecs()[i]->getEvRateHist()->Write();
  }

  return 0;
}

int main(int argc, char** argv)
{
  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable(); //Needed to look up dictionaries for PlotUtils classes like MnvH1D
  #endif

  TH1::AddDirectory(kFALSE); //Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).

  std::vector<const char*> fileNames(argv+1, argv+argc);

  double EAvailMax = 0.1; //GeV
  bool antinu=true;

  int Emin=0; //GeV
  int Emax=100; //GeV

  return runXSecLooper(antinu, EAvailMax, Emin, Emax, fileNames);
}
