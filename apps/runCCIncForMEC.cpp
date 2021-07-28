#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"

#include "TH1.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

#include <cstdlib>

typedef unsigned int uint;

template<class T> T sqr(T x) { return x*x; }

class CCIncXSec : public XSec
{
public:

  enum EEnergyVar
  {
    kTrueVisibleE,      // True visible energy: ~sum of charged hadrons
    kTrueq0,            // True q0

    kFracTrueVisibleE,  // Same as above, but divided by q3
    kFracTrueq0
  };

  CCIncXSec(const char* name, EEnergyVar energyVar, int intType, bool carbonOnly)
    : XSec(name), m_energyVar(energyVar), m_intType(intType), m_carbonOnly(carbonOnly)
    {}

  bool passesMuKinematics(ChainWrapper& chw, int entry)
  {
    const double thetamu=chw.GetValue("truth_muon_theta", entry);
    if(thetamu>20.*TMath::Pi()/180) return false;

    const double Emu=1e-3*chw.GetValue("mc_primFSLepton", entry, 3);
    // The Emu corresponding to pmu=1.5
    const double EmuCut=sqrt( sqr(1.5) + sqr(0.1056) );
    if(Emu<EmuCut) return false;

    return true;
  }

  virtual bool passesCuts(ChainWrapper& chw, int entry)
  {
    if(!passesMuKinematics(chw, entry)) return false;

    if((int)chw.GetValue("mc_current", entry)!=1) return false;
    const int mc_intType=(int)chw.GetValue("mc_intType", entry);

    if(m_carbonOnly){
      if((int)chw.GetValue("mc_targetZ", entry)!=6) return false;
    }

    switch(m_intType){
    case 0: // CC inclusive
    {
      return true;
    }
    case 1: // CC QE
    {
      return mc_intType==1;
    }
    case 2: // CC Delta
    {
      if(mc_intType!=2) return false;
      int mc_resID=(int)chw.GetValue("mc_resID", entry);
      return mc_resID==0;
    }
    case 3: // CC Res
    {
      return mc_intType==2;
    }
    case 4: // 2p2h nn
    {
      if(mc_intType!=8) return false;
      int targetNProtons=(int)chw.GetValue("mc_targetNucleon",
                                           entry)-2000000200;
      return targetNProtons==0;
    }
    case 5: // 2p2h np
    {
      if(mc_intType!=8) return false;
      int targetNProtons=(int)chw.GetValue("mc_targetNucleon",
                                           entry)-2000000200;
      return targetNProtons==1;
    }
    case 6: // 2p2h pp
    {
      if(mc_intType!=8) return false;
      int targetNProtons=(int)chw.GetValue("mc_targetNucleon",
                                           entry)-2000000200;
      return targetNProtons==2;
    }
    default:
      cout << "Unknown int type " << m_intType << endl;
      exit(1);
    }
  }

  std::vector<double> getVariableValues(ChainWrapper& chw, int entry)
  {

    std::vector<double> ret(2);
    double q3=1e-3*chw.GetValue("truth_q3", entry);
    ret[0]=q3;

    switch(m_energyVar){
    case kTrueVisibleE:
      ret[1]=1e-3*std::max(0., chw.GetValue("truth_visibleE", entry));
      break;
    case kTrueq0:
      ret[1]=1e-3*chw.GetValue("truth_nu", entry);
      break;
    case kFracTrueVisibleE:
      ret[1]=1e-3*std::max(0., chw.GetValue("truth_visibleE", entry))/q3;
      break;
    case kFracTrueq0:
      ret[1]=1e-3*chw.GetValue("truth_nu", entry)/q3;
      break;

    }

    // Stuff for debugging

    // int run=(int)chw.GetValue("mc_run", entry);
    // int subrun=(int)chw.GetValue("mc_subrun", entry);
    // int gate=(int)chw.GetValue("mc_nthEvtInFile", entry);

    // static int nCalls=0;
    // if(nCalls<100 && !m_plotTrueq0){
    //   printf("%d %d %d %.2f %.2f\n", run, subrun, gate, ret[0], ret[1]);
    //   ++nCalls;
    // }

    return ret;
  }

protected:
  EEnergyVar m_energyVar;
  int m_intType;
  bool m_carbonOnly;
};

class MyXSecLooper : public XSecLooper
{
public:
  MyXSecLooper(const char* inputFileGlob, bool doRPA, bool doPionWeight)
    : XSecLooper(inputFileGlob), m_doRPA(doRPA), m_doPionWeight(doPionWeight)
  {}

  virtual double getPionWeight(ChainWrapper& chw, int entry)
  {
    if(!m_doPionWeight) return 1;

    // Get the 1.5sigma weight
    double rvn1piWeight=1+1.5*(chw.GetValue("truth_genie_wgt_Rvn1pi", entry, 2)-1);

    int nPi=(int)chw.GetValue("truth_pionN", entry);
    double brandonWeight=1;
    if(nPi && chw.GetValue("mc_w", entry)<1800) brandonWeight=0.9;

    double cohWeight=1;
    if((int)chw.GetValue("mc_intType", entry)==4){
      double pionE=1e-3*chw.GetValue("truth_pionE", entry);
      if(pionE<0.450) cohWeight=0.5;
    }

    return rvn1piWeight*brandonWeight*cohWeight;
  }

  virtual double getRPAWeight(ChainWrapper& chw, int entry)
  {
    if(!m_doRPA) return 1;

    const int mc_intType=(int)chw.GetValue("mc_intType", entry);
    if(mc_intType!=1) return 1;

    //=====
    // This one is for using the 2D map in q3 (x) q0 (y)
    //
    static TH2D *hRPAratio=0;
    if(!hRPAratio){
      TFile f("outNievesRPA-Carbon3GeV.root");
      hRPAratio = (TH2D*)f.Get("hrelratio");
      hRPAratio->SetDirectory(0);
    }

    const double mc_q3=1e-3*chw.GetValue("truth_q3", entry);//getTrueq3(chw, entry);
    const double mc_q0=1e-3*chw.GetValue("truth_nu", entry); //getTrueq0(chw, entry);

    Int_t q3bin = Int_t(mc_q3*1000.);
    Int_t q0bin = Int_t(mc_q0*1000.);
    if(mc_q0 >= 1.20) q0bin = 1199;
    if(mc_q3 >= 1.20) q3bin = 1199;
    // Nieves does not calculate anything below binding energy.
    // I don't know what GENIE does, but lets be soft about this.
    if(mc_q0 < 0.017) q0bin = 17;
    Double_t thisrwtemp = hRPAratio->GetBinContent(q3bin,q0bin);

    // now trap bogus entries.  Not sure why they happen, but set to 1.0 not 0.0
    if(thisrwtemp <= 0.000001)thisrwtemp = 1.0;
    if(thisrwtemp >= 2.0)thisrwtemp = 2.0;
    if(!(thisrwtemp >= 0.000001 && thisrwtemp <= 2.0))thisrwtemp = 1.0;
    //cout << thisrwtemp << " rw " << q3bin << " " << q0bin << " " << endl;

    return thisrwtemp;
  }

  virtual double getSignalWeight(ChainWrapper& chw, int entry)
  {
    return getPionWeight(chw, entry)*getRPAWeight(chw, entry);
  }

protected:
  bool m_doRPA;
  bool m_doPionWeight;
};

void runXSecLooper(bool rpa, bool mec, bool pionWeight, bool antinu, bool carbonOnly, int Emin, int Emax)
{
  
  TString outfile("cc-inc-xsec");
  outfile+="-eroica";

  // First index is antinu, second index is mec
  TString fileGlobs[2][2];

  fileGlobs[false][false]="/pnfs/minerva/persistent/users/rodriges/mec/ana-2015-10-29/merged-lesstruth/mc/minerva*/merged_MECAnaTool*.root";
  fileGlobs[false][true]="/pnfs/minerva/persistent/users/rodriges/mec/ana-2015-10-29/merged-lesstruth/mecgenie/minerva13c/*.root";
  fileGlobs[true][false]="/pnfs/minerva/persistent/users/betan009/mec/MC/*.root";
  fileGlobs[true][true]="/pnfs/minerva/persistent/users/betan009/mec/ana-2017/mecgenie/minerva5/*.root";

  TString fileGlob=fileGlobs[antinu][mec];

  outfile+=rpa ? "-rpa" : "";
  outfile+=mec ? "-mec" : "";
  outfile+=pionWeight ? "" : "-nopionweight";
  outfile+=carbonOnly ? "-carbonOnly" : "";
  outfile+=antinu ? "-antinu" : "";
  outfile+=TString::Format("-Enu-%d-%d", Emin, Emax);
  outfile+=".root";

  // Create the XSecLooper and tell it the input files
  MyXSecLooper loop(fileGlob, rpa, pionWeight);
  if(mec){
    // The MEC files only contain MEC events, so the scaling to total
    // CC doesn't work. Instead we need to point it at a file with the
    // CC MEC spline, and normalize to that
    gSystem->Setenv("MPARAMFILESROOT", "/minerva/app/users/rodriges/cmtuser/mecgenie/MParamFiles");
    TString name(antinu ? "nu_mu_bar_C12/mec_cc_tensor" : "nu_mu_C12/mec_cc_tensor");
    loop.setRateSplineName(name);
  }

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(antinu ? -14: 14);
  
  //Setting the number of Universes in the GENIE error band, default 100 universes put 0 if you do not want universes to be included
  loop.setNumUniv(0); 

  const int nBinsq3=7;
  double binsq3[nBinsq3+1]={0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0};

  // Mangled bin edges so axis labels don't get drawn at the end in grids of plots
  const int nBinsq0=18;
  double binsq0[nBinsq0+1]={0, 0.02, 0.04, 0.06, 0.08, 0.1,
                            0.12, 0.14, 0.16, 0.2,
                            0.25, 0.299, 0.35, 0.399,
                            0.499, 0.599, 0.799, 0.999, 1.499};

  // Mangled bin edges so axis labels don't get drawn at the end in grids of plots
  const int nBinsq0_fine=40;
  double binsq0_fine[nBinsq0_fine+1]={0, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                                      0.10,  0.12, 0.14, 0.16, 0.18,
                                      0.20,  0.22, 0.24, 0.26, 0.28,
                                      0.30,  0.32, 0.34, 0.36, 0.38,
                                      0.399, 0.42, 0.44, 0.46, 0.48,
                                      0.499, 0.52, 0.54, 0.56, 0.58,
                                      0.599, 0.65, 0.70, 0.75, 0.80,
                                      0.9, 1.0};


  const int nBinsFracq0=20;
  double binsFracq0[nBinsFracq0+1]={0.  ,  0.05,  0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,
                                    0.45,  0.5 ,  0.55,  0.6 ,  0.65,  0.7 ,  0.75,  0.8 ,  0.85,
                                    0.9 ,  0.95, 1.00};


  for(int intType=0; intType<7; ++intType){
    TString tag;
    switch(intType){
    case 0: tag=""; break;
    case 1: tag="_qe"; break;
    case 2: tag="_delta"; break;
    case 3: tag="_res"; break;
    case 4: tag="_2p2h_nn"; break;
    case 5: tag="_2p2h_np"; break;
    case 6: tag="_2p2h_pp"; break;
    }
    // Add the CCQE differential cross section dsigma/dQ^2
    // flux-integrated over the range 1.5 to 12 GeV
    CCIncXSec* ccq0q3=new CCIncXSec(TString::Format("cc_q0q3%s", tag.Data()), CCIncXSec::kTrueq0, intType, carbonOnly);
    ccq0q3->setDimension(2);
    ccq0q3->setUniformBinning(80, 0, 2, 80, 0, 2);
    //ccq0q3->setVariable(XSec::kQ2);
    ccq0q3->setIsFluxIntegrated(true);
    ccq0q3->setFluxIntLimits(Emin, Emax);
    ccq0q3->setNormalizationType(XSec::kPerNucleon);
    if(carbonOnly) ccq0q3->setNormalizationValue(1);
    ccq0q3->setUniverses(0);//default value, put 0 if you do not want universes to be included.
  
    loop.addXSec(ccq0q3);

    CCIncXSec* ccVisEq3=new CCIncXSec(TString::Format("cc_visEq3%s", tag.Data()), CCIncXSec::kTrueVisibleE, intType, carbonOnly);
    ccVisEq3->setDimension(2);
    //ccVisEq3->setUniformBinning(80, 0, 2, 80, 0, 2);
    ccVisEq3->setBinEdges(nBinsq3, binsq3, nBinsq0, binsq0);

    ccVisEq3->setIsFluxIntegrated(true);
    ccVisEq3->setFluxIntLimits(Emin, Emax);
    ccVisEq3->setNormalizationType(XSec::kPerNucleon);
    if(carbonOnly) ccVisEq3->setNormalizationValue(1);
    ccVisEq3->setUniverses(0);
  
    loop.addXSec(ccVisEq3);

    CCIncXSec* ccVisEq3_fine=new CCIncXSec(TString::Format("cc_visEq3_fine%s", tag.Data()), CCIncXSec::kTrueVisibleE, intType, carbonOnly);
    ccVisEq3_fine->setDimension(2);
    ccVisEq3_fine->setBinEdges(nBinsq3, binsq3, nBinsq0_fine, binsq0_fine);
  
    ccVisEq3_fine->setIsFluxIntegrated(true);
    ccVisEq3_fine->setFluxIntLimits(Emin, Emax);
    ccVisEq3_fine->setNormalizationType(XSec::kPerNucleon);
    if(carbonOnly) ccVisEq3_fine->setNormalizationValue(1);
    ccVisEq3_fine->setUniverses(0);
  
    loop.addXSec(ccVisEq3_fine);

    CCIncXSec* ccFracq0q3=new CCIncXSec(TString::Format("cc_fracq0q3%s", tag.Data()), CCIncXSec::kFracTrueq0, intType, carbonOnly);
    ccFracq0q3->setDimension(2);
    ccFracq0q3->setUniformBinning(80, 0, 2, 80, 0, 1);
    ccFracq0q3->setIsFluxIntegrated(true);
    ccFracq0q3->setFluxIntLimits(Emin, Emax);
    ccFracq0q3->setNormalizationType(XSec::kPerNucleon);
    if(carbonOnly) ccFracq0q3->setNormalizationValue(1);
    ccFracq0q3->setUniverses(0);
  
    loop.addXSec(ccFracq0q3);

    CCIncXSec* ccFracVisEq3=new CCIncXSec(TString::Format("cc_fracvisEq3%s", tag.Data()), CCIncXSec::kFracTrueVisibleE, intType, carbonOnly);
    ccFracVisEq3->setDimension(2);
    ccFracVisEq3->setBinEdges(nBinsq3, binsq3, nBinsFracq0, binsFracq0);
    ccFracVisEq3->setIsFluxIntegrated(true);
    ccFracVisEq3->setFluxIntLimits(Emin, Emax);
    ccFracVisEq3->setNormalizationType(XSec::kPerNucleon);
    if(carbonOnly) ccFracVisEq3->setNormalizationValue(1);
    ccFracVisEq3->setUniverses(0);//default value, put 0 if you do not want universes to be included.
  
    loop.addXSec(ccFracVisEq3);
  }

  // Once everything's set up, actually run the thing
  loop.runLoop();

  // Get the output histograms and save them to file
  TFile fout(outfile, "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i){
    loop.getXSecs()[i]->get2DXSecHist()->Write();
  }
}

int main(int argc, char** argv)
{
  bool rpa=false;
  bool mec=false;
  bool pionWeight=true;
  bool antinu=false;
  bool carbonOnly=false;

  int Emin=2;
  int Emax=6;

  for(int i=1; i<argc; ++i){
    TString arg(argv[i]);
    arg.ToLower();
    if(arg=="rpa") rpa=true;
    if(arg=="mec") mec=true;
    if(arg=="antinu") antinu=true;
    if(arg=="carbon") carbonOnly=true;
    if(arg=="nopion" || arg=="nopionweight") pionWeight=false;
    if(arg=="highe"|| arg=="highenergy"){
      Emin=6; Emax=20;
    }
  }
  runXSecLooper(rpa, mec, pionWeight, antinu, carbonOnly, Emin, Emax);
  return 0;
}

