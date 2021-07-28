#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include <PlotUtils/MnvH1D.h>

#include <cstdlib>
typedef unsigned int uint;

class CCProtonPi0XSec : public XSec
{
    public:
        CCProtonPi0XSec(const char* name) 
            : XSec(name)
        {};

        // Override this method from the base class to decide what events to
        // include in this selection
        virtual bool passesCuts(ChainWrapper& chw, int entry)
        {
            if((bool)chw.GetValue("truth_isSignal", entry)) return true;
            else return false;
        }

};

void runCCProtonPi0XSec()
{
    // Create the XSecLooper and tell it the input files
    XSecLooper loop("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_78/minerva1/grid/central_value/minerva/ana/v10r8p9/00/01/02/00/*.root");

    // Tell the XSecLooper which neutrino type we're considering (mandatory)
    loop.setNuPDG(14);

    //Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
    loop.setNumUniv(0); 

    // ------------------------------------------------------------------------
    // Muon Momentum
    // ------------------------------------------------------------------------
    CCProtonPi0XSec* muon_P = new CCProtonPi0XSec("muon_P");
    muon_P->setUniformBinning(10,0.0,10.0);
    muon_P->setVariable(XSec::kPLep);
    muon_P->setIsFluxIntegrated(true);
    muon_P->setFluxIntLimits(1.5, 20.0);
    muon_P->setNormalizationType(XSec::kPerNucleon);  
    loop.addXSec(muon_P);

    // ------------------------------------------------------------------------
    // Muon Theta
    // ------------------------------------------------------------------------
    CCProtonPi0XSec* muon_theta = new CCProtonPi0XSec("muon_theta");
    muon_theta->setUniformBinning(12,0.0,25.0);
    muon_theta->setVariable(XSec::kThetaLep);
    muon_theta->setIsFluxIntegrated(true);
    muon_theta->setFluxIntLimits(1.5, 20.0);
    muon_theta->setNormalizationType(XSec::kPerNucleon);  
    loop.addXSec(muon_theta);

    // ------------------------------------------------------------------------
    // Pi0 Momentum
    // ------------------------------------------------------------------------
    CCProtonPi0XSec* pi0_P = new CCProtonPi0XSec("pi0_P");
    pi0_P->setUniformBinning(17,0.0,1.7);
    pi0_P->setVariable(XSec::kPPi0);
    pi0_P->setIsFluxIntegrated(true);
    pi0_P->setFluxIntLimits(1.5, 20.0);
    pi0_P->setNormalizationType(XSec::kPerNucleon);  
    loop.addXSec(pi0_P);

    // ------------------------------------------------------------------------
    // Pi0 Kinetic Energy 
    // ------------------------------------------------------------------------
    CCProtonPi0XSec* pi0_KE = new CCProtonPi0XSec("pi0_KE");
    pi0_KE->setUniformBinning(17,0.0,1.7);
    pi0_KE->setVariable(XSec::kTPi0);
    pi0_KE->setIsFluxIntegrated(true);
    pi0_KE->setFluxIntLimits(1.5, 20.0);
    pi0_KE->setNormalizationType(XSec::kPerNucleon);  
    loop.addXSec(pi0_KE);

    // ------------------------------------------------------------------------
    // Pi0 Theta 
    // ------------------------------------------------------------------------
    CCProtonPi0XSec* pi0_theta = new CCProtonPi0XSec("pi0_theta");
    pi0_theta->setUniformBinning(18,0.0,180.0);
    pi0_theta->setVariable(XSec::kThetaPi0);
    pi0_theta->setIsFluxIntegrated(true);
    pi0_theta->setFluxIntLimits(1.5, 20.0);
    pi0_theta->setNormalizationType(XSec::kPerNucleon);  
    loop.addXSec(pi0_theta);

    // ------------------------------------------------------------------------
    // QSq
    // ------------------------------------------------------------------------
    CCProtonPi0XSec* QSq = new CCProtonPi0XSec("QSq");
    QSq->setUniformBinning(40,0.0,4.0);
    QSq->setVariable(XSec::kQ2);
    QSq->setIsFluxIntegrated(true);
    QSq->setFluxIntLimits(1.5, 20.0);
    QSq->setNormalizationType(XSec::kPerNucleon);  
    loop.addXSec(QSq);

    // Once everything's set up, actually run the thing
    loop.runLoop();

    // Get the output histograms and save them to file
    string rootDir =  "/minerva/data/users/oaltinok/NTupleAnalysis/GENIEXSec/GENIEXSec.root";
    TFile f_out(rootDir.c_str(), "RECREATE");
    for(uint i=0; i<loop.getXSecs().size(); ++i){
        loop.getXSecs()[i]->getXSecHist()->Write();
    }
}

int main()
{
    runCCProtonPi0XSec();
    return 0;
}

