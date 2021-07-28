#include "GENIEXSecExtract/XSecLooper.h"

#include <cstdlib>
#include "TFile.h"
#include "TH1.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include <PlotUtils/MnvH1D.h>

//======================================
//! transverseCCQETwoTrack XSec
//======================================
class transverseCCQETwoTrackXSec : public XSec
{

public:
       transverseCCQETwoTrackXSec( const char* name, bool isQElike = true, 
                                  bool isQE = false, bool isRes = false, bool isDIS = false) : 
         XSec( name ),
         m_QElike( isQElike ),
         m_isQE( isQE ),
         m_isRes( isRes ),
         m_isDIS( isDIS ) {}

  static double getTheta( double x, double y, double z )
  {
    double pyprime = -1.0*sin(-0.05887)*z + cos(-0.05887)*y;
    double pzprime =  1.0*cos(-0.05887)*z + sin(-0.05887)*y;
    double den2    = x*x + pyprime*pyprime + pzprime*pzprime;

    if( den2 == 0 ) return -9999;
    else return acos( pzprime / sqrt(den2) )*180/TMath::Pi();
  }
  
       //! Return QE-like signal
       bool isQElikeSignal( ChainWrapper& chw, int entry ) {
         /*

a) muon theta angle (w.r.t. beam) < 20 deg
b) muon momentum 1.5-10 GeV
c) proton theta angle (w.r.t beam) < 70 deg
d) proton momentum 0.45-1.2 GeV

          */  
         bool pass_proton_fiducial = false;
         bool pass_muon_fiducial = false;

         int genie_n_muons         = 0;
         int genie_n_mesons        = 0;
         int genie_n_heavy_baryons = 0;
         int genie_n_photons       = 0;
         int genie_n_protons       = 0;
         int genie_n_neutrons      = 0;
         int genie_n_nucleons      = 0;
         int genie_n_pions         = 0;
         int genie_n_pi_zeros      = 0;
         int genie_n_charms        = 0;
         int genie_n_kaons         = 0;
         int genie_n_others        = 0;
         bool signal = false;

         const double muon_angle_MAXCut    =  20.0;
         const double muon_momentum_minCut    =  1500.0;
         const double muon_momentum_MAXCut    = 10000.0;
         const double proton_angle_MAXCut  =  70.0; 
         const double proton_momentum_minCut  =   450.0;
         const double proton_momentum_MAXCut  =  1200.0; 

         const int nparticles = (int)chw.GetValue("mc_nFSPart",entry);
         for(int i = 0; i < nparticles; i++) {
           //don't use abs here!
           const int pdg = (int)chw.GetValue("mc_FSPartPDG",entry,i);
           
           const double px = (double)chw.GetValue("mc_FSPartPx",entry,i);
           const double py = (double)chw.GetValue("mc_FSPartPy",entry,i);
           const double pz = (double)chw.GetValue("mc_FSPartPz",entry,i);
           const TVector3 pvec(px, py, pz);

           const double mom = pvec.Mag();
           const double theta = getTheta(px,py,pz);

           //muon- for neutrino
           if(pdg == 13){
             if( theta < muon_angle_MAXCut &&
                 mom > muon_momentum_minCut &&
                 mom < muon_momentum_MAXCut){

               pass_muon_fiducial = true;
             }
           }
           else if(pdg == 2212){
             if( theta < proton_angle_MAXCut &&
                 mom > proton_momentum_minCut &&
                 mom < proton_momentum_MAXCut){
            
               pass_proton_fiducial = true;
             }
           }

           //use fabs here!                                                                                                                                                                                                                      
           const int apd = fabs(pdg);
           if( apd == 13 ) genie_n_muons++;
           else if( apd == 22 || apd == 11 )  genie_n_photons++;
           else if( apd == 2212 ) { genie_n_protons++;  genie_n_nucleons++; }
           else if( apd == 2112 ) { genie_n_neutrons++; genie_n_nucleons++; }
           else if( apd == 211 )  { genie_n_pions++;    genie_n_mesons++; }
           else if( apd == 111 )  { genie_n_pi_zeros++; genie_n_photons++; }
           else if( apd == 411 || apd == 421 || apd == 431 ) { genie_n_charms++; genie_n_mesons++; }
           else if( apd == 321 || apd == 311 || apd == 310 || apd == 130 ) { genie_n_kaons++; genie_n_mesons++; }
           else if( apd > 3000 && apd < 5000 ) genie_n_heavy_baryons++;
           else genie_n_others++;

         }
         
         if( genie_n_muons         == 1 &&
             genie_n_mesons        == 0 &&
             genie_n_heavy_baryons == 0 &&
             genie_n_photons       == 0 &&
             genie_n_protons       != 0  ) {
           signal = true;
         }

         return pass_muon_fiducial && pass_proton_fiducial && signal;
      }

      //! Return QE signal
      bool isQElikeQE( ChainWrapper& chw, int entry ) {

         bool signal = false;

         if( isQElikeSignal(chw,entry) ) {
           if( (int)chw.GetValue("mc_intType",entry) == 1 && (int)chw.GetValue("mc_current",entry) == 1 ) signal = true;
         }

         return signal;
      }

      //! Return Resonant signal
      bool isQElikeResonant( ChainWrapper& chw, int entry ) {

         bool signal = false;

         if( isQElikeSignal(chw,entry) ) {
           if( (int)chw.GetValue("mc_intType",entry) == 2 && (int)chw.GetValue("mc_current",entry) == 1 ) signal = true;
         }

         return signal;
      }

      //! Return DIS signal
      bool isQElikeDIS( ChainWrapper& chw, int entry ) {

         bool signal = false;

         if( isQElikeSignal(chw,entry) ) {
           if( (int)chw.GetValue("mc_intType",entry) == 3 && (int)chw.GetValue("mc_current",entry) == 1 ) signal = true;
         }

         return signal;
      }

      //! Override this method from the base class to decide what events to include in this selection
      virtual bool passesCuts( ChainWrapper& chw, int entry ) {

         //! pass signal
         if( !m_QElike ) {
           printf("should not happen! only QElike is required! %d\n", m_QElike); exit(1);
         } 

         if( !isQElikeSignal(chw,entry) ) {
           return false;
         }

         if( m_isQE && !m_isRes && !m_isDIS) {
           return isQElikeQE(chw,entry);
         } 

         if( !m_isQE && m_isRes && !m_isDIS) {
           return isQElikeResonant(chw,entry);
         } 

         if( !m_isQE && !m_isRes && m_isDIS) {
           return isQElikeDIS(chw,entry);
         } 

         //total ccqe_like 
         if( !m_isQE && !m_isRes && !m_isDIS) {
           return true;
         }

         printf("passesCuts should not be here!! %d %d %d %d\n", m_QElike, m_isQE, m_isRes, m_isDIS); exit(1);

         return false;

      } //! end function passesCuts

      //! return the binning
      static void GetBins(const Int_t ivar, double *tmpbin, Int_t &nbins) {
        const double muonmomentumfullBins[] = {1, 1.250, 1.500, 1.750, 2.000, 2.250, 2.500, 2.750, 3.000, 3.250, 3.500, 3.750, 4.000, 4.250, 4.500, 4.750, 5.000, 5.250, 5.500, 5.750, 6.000, 6.250, 6.500, 6.750, 7.000, 7.250, 7.500, 7.750, 8.000, 8.250, 8.500, 8.750, 9.000, 9.250, 9.500, 9.750, 10.000, 10.250};//in GeV getMuonMomentumBins
        const double muonthetafullBins[] = {-0.1, 0.000, 1.000, 2.000, 3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.000, 11.000, 12.000, 13.000, 14.000, 15.000, 16.000, 17.000, 18.000, 19.000, 20.000, 21.000};//getMuonThetaBins;
        const double protonmomentumfullBins[] = {0.4, 0.425, 0.450, 0.600, 0.625, 0.650, 0.675, 0.700, 0.725, 0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 0.900, 0.925, 0.950, 0.975, 1.000, 1.025, 1.050, 1.075, 1.100, 1.125, 1.150, 1.175, 1.200, 1.225};//in GeV getProtonMomentumBins
        const double protonthetafullBins[] = {-0.1, 0.000, 2.500, 5.000, 7.500, 10.000, 12.500, 15.000, 17.500, 20.000, 22.500, 25.000, 27.500, 30.000, 32.500, 35.000, 37.500, 40.000, 42.500, 45.000, 47.500, 50.000, 52.500, 55.000, 57.500, 60.000, 62.500, 65.000, 67.500, 70.000, 72.500};//getProtonThetaBins
        const double dalphatfullBins[] = {-0.1, 0, 20, 40, 60, 80, 100, 120, 130, 140, 150, 160, 170, 180};//getDalphatBins
        const double dptfullBins[] = { -0.001, 0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.8, 1.0, 1.2, 2.0, 2.02};//in GeV getDptBins
        const double dphitfullBins[] = {-0.1, 0.00, 2.50, 5.00, 7.50, 10.00, 12.50, 15.00, 17.50, 20.00, 22.50, 25.00, 27.50, 30.00, 35.00, 40.00, 45.00, 50.00, 55, 60, 70, 85, 105, 130, 180};//getDphitBins
        const double neutronmomentumfullBins[] = { -0.001, 0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350, 0.375, 0.400, 0.425, 0.450, 0.475, 0.500, 0.525, 0.550, 0.575, 0.600, 0.650, 0.700, 0.8, 1.0, 1.2, 2.0, 2.02};//in GeV getNeutronMomentumBins
        
        switch(ivar){
        case(XSec::kTmuonmomentum): nbins = sizeof(muonmomentumfullBins)/sizeof(double); memcpy(tmpbin, muonmomentumfullBins, sizeof(muonmomentumfullBins)); break;
        case(XSec::kTmuontheta): nbins = sizeof(muonthetafullBins)/sizeof(double); memcpy(tmpbin, muonthetafullBins, sizeof(muonthetafullBins)); break;
        case(XSec::kTprotonmomentum): nbins = sizeof(protonmomentumfullBins)/sizeof(double); memcpy(tmpbin, protonmomentumfullBins, sizeof(protonmomentumfullBins)); break;
        case(XSec::kTprotontheta): nbins = sizeof(protonthetafullBins)/sizeof(double); memcpy(tmpbin, protonthetafullBins, sizeof(protonthetafullBins)); break;
        case(XSec::kTdalphat): nbins = sizeof(dalphatfullBins)/sizeof(double); memcpy(tmpbin, dalphatfullBins, sizeof(dalphatfullBins)); break;
        case(XSec::kTdpt): nbins = sizeof(dptfullBins)/sizeof(double); memcpy(tmpbin, dptfullBins, sizeof(dptfullBins)); break;
        case(XSec::kTdphit): nbins = sizeof(dphitfullBins)/sizeof(double); memcpy(tmpbin, dphitfullBins, sizeof(dphitfullBins)); break;
        case(XSec::kTneutronmomentum): nbins = sizeof(neutronmomentumfullBins)/sizeof(double); memcpy(tmpbin, neutronmomentumfullBins, sizeof(neutronmomentumfullBins)); break;
        default: printf("GetBins should not be hear %d\n", ivar); exit(1);
        }

        //nbins = boundary -1
        nbins--;
      }

      static string GetName(const Int_t ivar){
        switch(ivar){
        case XSec::kTmuonmomentum: return "muonmomentum";
        case XSec::kTmuontheta: return "muontheta";
        case XSec::kTprotonmomentum: return "protonmomentum";
        case XSec::kTprotontheta: return "protontheta";
        case XSec::kTdpt: return "dpt";
        case XSec::kTdphit: return "dphit";
        case XSec::kTdalphat: return "dalphat";
        case XSec::kTneutronmomentum: return "neutronmomentum";
        }
        printf("string GetName(const Int_t ivar) ivar error %d\n", ivar);
        exit(1);
      }
  
      //! input parameters
      bool m_QElike;
      bool m_isQE;
      bool m_isRes;
      bool m_isDIS;
};

//======================
//! main: 
//======================
void runXSecLooper_transverseCCQETwoTrack(const Int_t inopt, const int fsitarget)
{
  const string fsiname=XSecLooper::FSIName(fsitarget);
  printf("\nrunXSecLooper_transverseCCQETwoTrack inopt %d fsitarget %d %s\n", inopt, fsitarget, fsiname.c_str());

   //! get the ntuple directory names
   string directory   = "genie_xsection_ntuples";
   string output_name = "genie_xsections";

   output_name += "_plastic";
   output_name += "_minos_match";
   output_name += "_full_fiducial";
   output_name += Form("_histos%d.root", inopt);

   cout << "	creating a genie file with name = " << output_name << endl;

   //! Create the XSecLooper and tell it the input files
   string input;
   if(inopt==0){
     input = "/pnfs/minerva/persistent/users/xlu/NewProduction_protonPID_MichelTag/DataSet/mc_minerva1_20170827/grid/central_value/minerva/ana/v10r8p12/00/01/02/00/SIM_minerva_00010200_Subruns_0050_Ana_Tuple_v10r8p12.root";
   }
   else if(inopt==1){
     //can't use, contains 2p2h
     //input = "/pnfs/minerva/persistent/users/xlu/NewProduction_protonPID_MichelTag/DataSet/Ana20171115/mergedoutput/suboutputmc_minerva1/suboutputmc_minerva1_x05_mergeTreeOutput.root";
     input = "/pnfs/minerva/persistent/users/xlu/NewProduction_protonPID_MichelTag/DataSet/Ana20171115/mergedoutput/suboutputmc_minerva13C/suboutputmc_minerva13C_x08_mergeTreeOutput.root";
   }
   else if(inopt==2){
     input = "/pnfs/minerva/persistent/users/xlu/NewProduction_protonPID_MichelTag/DataSet/Ana20171115/mergedoutput/suboutputmc_minerva1/suboutputmc_minerva1_x00_mergeTreeOutput.root";
   }
   else if(inopt==3){
     input = "/pnfs/minerva/persistent/users/xlu/NewProduction_protonPID_MichelTag/DataSet/Ana20171115/mergedoutput/suboutputmc_minerva9/suboutputmc_minerva9_x00_mergeTreeOutput.root";
   }
   else if(inopt==4){
     input = "/pnfs/minerva/persistent/users/xlu/NewProduction_protonPID_MichelTag/DataSet/Ana20171115/mergedoutput/suboutputmc_minerva13C/suboutputmc_minerva13C_x00_mergeTreeOutput.root";
   }
   else if(inopt==5){
     input = "/pnfs/minerva/persistent/users/xlu/NewProduction_protonPID_MichelTag/DataSet/Ana20171115/mergedoutput/suboutputmc_minerva7/suboutputmc_minerva7_x00_mergeTreeOutput.root";
   }
   else if(inopt==101){
     input = "/pnfs/minerva/persistent/users/xlu/NewProduction_protonPID_MichelTag/DataSet/Ana20171115/finallist/mc_minerva1.txt";
   }
   else if(inopt==107){
     input = "/pnfs/minerva/persistent/users/xlu/NewProduction_protonPID_MichelTag/DataSet/Ana20171115/finallist/mc_minerva7.txt";
   }
   else if(inopt==109){
     input = "/pnfs/minerva/persistent/users/xlu/NewProduction_protonPID_MichelTag/DataSet/Ana20171115/finallist/mc_minerva9.txt";
   }
   else if(inopt==113){
     input = "/pnfs/minerva/persistent/users/xlu/NewProduction_protonPID_MichelTag/DataSet/Ana20171115/finallist/mc_minerva13C.txt";
   }
   else if(inopt==201){
     input = "/minerva/app/users/xlu/cmtuser/Minerva_v10r8p12/Ana/AnalysisFramework/External/GENIEXSecExtract/test/201_nofsi.txt";
   }
   else if(inopt==202){
     input = "/minerva/app/users/xlu/cmtuser/Minerva_v10r8p12/Ana/AnalysisFramework/External/GENIEXSecExtract/test/202_nofsi.txt";
   }
   else if(inopt==301){
     input = "/pnfs/minerva/persistent/users/xlu/SysXtalk_SpecialSample_only13Cnominal/DataSet/Ana20180123/finallist/mc_minerva13D.txt" ;
   }
   else if(inopt==302){
     input = "/pnfs/minerva/persistent/users/xlu/SysXtalk_SpecialSample_only13Cnominal/DataSet/Ana20180123/finallist/mc_minerva13E.txt";
   }

   printf("Input: %s\n", input.c_str());

   XSecLooper loop( input.c_str() );

   //! Tell the XSecLooper which neutrino type we're considering (mandatory)
   loop.setNuPDG(14);
  
   //! Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
   loop.setNumUniv(0); 

   //
   for(int ivar=XSec::kTransverseBegin+1; ivar<XSec::kTransverseEnd; ivar++){
     const TString xname = Form("ccqe_like_plastic_total_%s", transverseCCQETwoTrackXSec::GetName(ivar).c_str());

     printf("\n\n ivar %d %s\n", ivar, xname.Data());

     //! Add the CCQE differential cross section dsigma/dQ^2 
     transverseCCQETwoTrackXSec* xsec = new transverseCCQETwoTrackXSec(xname);

     //! container for bins
     int nbins = 0; 
     double tmpbin[100];
     transverseCCQETwoTrackXSec::GetBins(ivar, tmpbin, nbins);;

     //! set cross section data
     xsec->setBinEdges(nbins, tmpbin);
     xsec->setIsFluxIntegrated(true);
     xsec->setUniverses(0);

     xsec->setNormalizationType(XSec::kPerNucleon);

     xsec->setVariable((XSec::EVariable)ivar);

     //related to flux integral, has to be full range
     xsec->setFluxIntLimits(0.0,1000);

     //! add cross section
     loop.addXSec(xsec);

   } //! end loop over cross sections

   //! run
   loop.runLoop(0, fsitarget);

   TCanvas *c1 = new TCanvas;
   //! get the output histograms and save them to file
   string filename = Form("/minerva/data/users/%s/files/%s/genie/%s",getenv("USER"),getenv("MINERVA_RELEASE"),output_name.c_str());
   TFile f(filename.c_str(),"recreate");
   
   printf("Output: %s\n", filename.c_str());

   for(unsigned int i = 0; i < loop.getXSecs().size(); i++) {
     PlotUtils::MnvH1D * hh = loop.getXSecs().at(i)->getXSecHist();
     hh->Write();

     const TString hname = hh->GetName();
     gStyle->SetOptStat(0);

     printf("XSecLooper final %s bin 2: %e bin 10: %e\n", hh->GetName(), hh->GetBinContent(2), hh->GetBinContent(10));

     hh->Draw();
     c1->Print(Form("%s_%d_FSI%d%s.png", hname.Data(), inopt, abs(fsitarget), fsiname.c_str()));
     c1->Print(Form("%s_%d_FSI%d%s.cxx", hname.Data(), inopt, abs(fsitarget), fsiname.c_str()));
   }

   return;
}

int main(int argc, char* argv[])
{
  cout << "Enter running the GENIE Xsection Extraction for the 2track QE analysis" << endl;  

  if(argc!=3){
    printf("bad argc %d\n", argc); 
    return 1;
  }

  runXSecLooper_transverseCCQETwoTrack(atoi(argv[1]), atoi(argv[2]));
  
  cout << "completed genie xsection extraction" << endl;

  cout << "Exit running the GENIE XSection Extraction for the 2track QE analysis" << endl;
  return 0;
}
