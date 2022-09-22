#include "GENIEXSecExtract/XSecLooper.h"

#include "GENIEXSecExtract/ChainWrapper.h"
#include "GENIEXSecExtract/XSec.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"

#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>
#include <PlotUtils/MnvApplication.h>
#include <PlotUtils/FluxReweighter.h>


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <sstream>

#include <algorithm>

using namespace std;
using namespace PlotUtils;

bool kprint=false;

//======================================================================
XSecLooper::XSecLooper(const char* inputFileGlob)
: m_chw(new ChainWrapper("Truth")),
m_fluxHist(0), m_nUniv(0), m_rateSplineName(""),
m_fiducialZMin(5990), m_fiducialZMax(8340), m_apothem(850),
m_useFiducialBranch(true), m_playlist(FluxReweighter::minervame1D)
{
    PlotUtils::Initialize();
    addFiles(inputFileGlob);
    
    //setup CV weighting infrastructure
    char *mparalocation = std::getenv("MPARAMFILESROOT");
    string directory_data = string(mparalocation)+"/data/Reweight/";
    weight_cv_2p2h= new weight_2p2h(directory_data+"/fit-mec-2d-noScaleDown-penalty00300-best-fit");
    weight_cv_and_var_RPA = new weightRPA(directory_data+"/outNievesRPAratio-nu12C-20GeV-20170202.root");
}

//======================================================================
void XSecLooper::addFiles(const char* inputFileGlob)
{
    const TString inlist(inputFileGlob);
    if(inlist.Contains(".root")){
        m_chw->Add(inputFileGlob);
    }
    else{
        ifstream infile;
        infile.open(inputFileGlob);
        if( !infile ) { cout << "filename = " << inputFileGlob  << ", could not be open! exiting!" << endl; exit(-1); }
        
        cout<<"XSecLooper::addFiles list: "<<inputFileGlob<<endl;
        
        int nfile=0;
        string line;
        while( !infile.eof() ) {
            infile >> line;
            if(line==""){ continue;}
            
            TString tmpl(line);
            if(tmpl.Contains("#")){
                continue;
            }
            
            m_chw->Add(line.c_str());
            nfile++;
            line="";
        }
        infile.close();
        printf("XSecLooper::addFiles nfile %d\n", nfile);
    }
}

//======================================================================
XSecLooper::~XSecLooper()
{
    delete m_chw;
}


//======================================================================
double XSecLooper::GetTrueQ3(const int entry){
    double px = m_chw->GetValue("mc_primFSLepton",entry,0)-m_chw->GetValue("mc_incomingPartVec",entry,0);
    double py = m_chw->GetValue("mc_primFSLepton",entry,1)-m_chw->GetValue("mc_incomingPartVec",entry,1);
    double pz = m_chw->GetValue("mc_primFSLepton",entry,2)-m_chw->GetValue("mc_incomingPartVec",entry,2);
    double q3 = sqrt(px*px+py*py+pz*pz);
    return q3/1000.0;
}
//======================================================================
double XSecLooper::GetTrueQ0(const int entry){
    double q0 = -1;
    q0 = m_chw->GetValue("mc_incomingPartVec",entry,3)-m_chw->GetValue("mc_primFSLepton",entry,3);
    return q0/1000.0;
}


//======================================================================
double XSecLooper::Get2p2hWeight(const int entry){
    if(m_chw->GetValue("mc_intType",entry)!=8){
        return 1.0; //if given invalid variations or the user doesn't want to apply this modification or the MC event is not CCQE or MEC
    }
    double q3=XSecLooper::GetTrueQ3( entry );//GeV
    double q0=XSecLooper::GetTrueQ0( entry );//GeV
    double ret = 1.0;
    ret=weight_cv_2p2h->getWeight(q0,q3);
    return ret;
}
//======================================================================
double XSecLooper::GetRPAWeight(const int entry){
    double ret = 1.0;
    if(m_chw->GetValue("mc_intType",entry)!=1) return 1.0;
    if(m_chw->GetValue("mc_targetZ",entry)<6) return 1.0;
    double q3=GetTrueQ3( entry );//GeV
    double q0=GetTrueQ0( entry );//GeV
    ret=weight_cv_and_var_RPA->getWeight(q0,q3);
    return ret;
}




//====================================================================== 
void XSecLooper::addXSec(XSec* xsec)
{
    m_xsecs.push_back(xsec);
}

//====================================================================== 
//http://cdcvs.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/AnalysisFramework/Ana/CCPionInc/ana/CCNuPionInc_Ana/src/TruthTuple.cxx?rev=1.9.2.3;content-type=text%2Fplain;cvsroot=mnvsoft
void XSecLooper::GetMotherFSIParticles( std::vector<unsigned int>& particles, const int entry ) const {
    
    const int mc_er_nPart = m_chw->GetValue("mc_er_nPart", entry);
    
    particles.clear();
    for (unsigned int i1 = 0; i1 != (unsigned int)mc_er_nPart; ++i1) {
        const int mc_er_status = m_chw->GetValue("mc_er_status", entry, i1);
        if (mc_er_status == 14) {
            unsigned int index = i1;
            while (true) {
                const int mc_er_mother = m_chw->GetValue("mc_er_mother", entry, index);
                int mother = mc_er_mother;
                if (mother < 0) break;
                const int mc_er_status = m_chw->GetValue("mc_er_status", entry, mother);
                if (mc_er_status != 14) break;
                index = (unsigned int)mother;
            }
            
            if ( std::find(particles.begin(), particles.end(), index) == particles.end() ) particles.push_back(index);
        }
    }
}

void XSecLooper::GetFSParticles( int& pionindex, int & protonid, const int entry ) const {
    
    pionindex = -1;
    unsigned int n_FS = 0;
    
    const int mc_er_nPart = m_chw->GetValue("mc_er_nPart", entry);
    
    for (unsigned int i1 = 0; i1 != (unsigned int)mc_er_nPart; ++i1) {
        const int mc_er_status = m_chw->GetValue("mc_er_status", entry, i1);
        if (mc_er_status == 1) {
            n_FS++;
            
            const int mc_er_ID = m_chw->GetValue("mc_er_ID", entry, i1);
            if ( abs(mc_er_ID) == 211 ) pionindex = i1;
            else if( mc_er_ID == 2212 ) protonid = i1;
        }
    }
    
    const int mc_nFSPart = m_chw->GetValue("mc_nFSPart", entry);
    if (n_FS != (unsigned int)mc_nFSPart) std::cout<<"Problem with number of FS particles!"<<std::endl;
    //if (pionindex == -1) std::cout<<"Could not find final state pion!"<<std::endl;
}

bool XSecLooper::IsMuonCollinear(const TVector3 selectedp, const int entry) const
{
    //have to go to transverse plan! logitudinal angular difference is big!
    const double neutrinopx = m_chw->GetValue("mc_incomingPartVec", entry, 0);
    const double neutrinopy = m_chw->GetValue("mc_incomingPartVec", entry, 1);
    const double neutrinopz = m_chw->GetValue("mc_incomingPartVec", entry, 2);
    TVector3 tmpneutrino(neutrinopx, neutrinopy, neutrinopz);
    const TVector3 unitneutrino = tmpneutrino.Unit();
    
    const double muonpx = m_chw->GetValue("mc_primFSLepton", entry, 0);
    const double muonpy = m_chw->GetValue("mc_primFSLepton", entry, 1);
    const double muonpz = m_chw->GetValue("mc_primFSLepton", entry, 2);
    const TVector3 muv(muonpx, muonpy, muonpz);
    
    const TVector3 plmuon   = unitneutrino *   muv.Dot(unitneutrino);
    const TVector3 plproton = unitneutrino * selectedp.Dot(unitneutrino);
    
    const TVector3 pTmuon   = muv   - plmuon;
    const TVector3 pTproton = selectedp - plproton;
    
    const double costheta = pTproton.Dot(-pTmuon)/pTproton.Mag()/pTmuon.Mag();
    
    const double theta = TMath::ACos(costheta)*TMath::RadToDeg();
    
    const double epsilon = 2.5;
    
    return fabs(theta)< epsilon;
}

bool XSecLooper::IsScaled(const TVector3 selectedp, const TVector3 primaryproton) const
{
    const double epsilon = 1e-6;
    
    const double costheta = selectedp.Dot(primaryproton)/selectedp.Mag()/primaryproton.Mag();
    
    return fabs(fabs(costheta)-1)<epsilon;
}

string XSecLooper::FSIName(const int imode)
{
    if(imode<0){
        return "Default";
    }
    
    switch(imode){
        case ProtonNonInteracting: return "ProtonNonInteracting";
        case ProtonAcceleration: return "ProtonAcceleration";
        case ProtonDeceleration: return "ProtonDeceleration";
        case PionAbsorption: return "PionAbsorption";
        default:
            printf("not known FSI name!\n"); exit(1);
    }
    
    return "NOFSI";
}

void XSecLooper::PrintProton(const int fsitype, int & kprint, const int entry, const int protonindex, const int it_primary_proton) const
{
    if(kprint<5){
        printf("\nPrintProton eventType %.0f FSI %d %s: selected %d %f, %f, %f, %f -- primary %d %f, %f, %f, %f\n", m_chw->GetValue("mc_intType",entry), fsitype, FSIName(fsitype).c_str(), protonindex, m_chw->GetValue("mc_FSPartE",  entry, protonindex), m_chw->GetValue("mc_FSPartPx", entry, protonindex), m_chw->GetValue("mc_FSPartPy", entry, protonindex), m_chw->GetValue("mc_FSPartPz", entry, protonindex), it_primary_proton, m_chw->GetValue("mc_er_E",  entry, it_primary_proton), m_chw->GetValue("mc_er_Px", entry, it_primary_proton), m_chw->GetValue("mc_er_Py", entry, it_primary_proton), m_chw->GetValue("mc_er_Pz", entry, it_primary_proton));
        kprint++;
    }
}

int XSecLooper::FSIType(const int entry, const int protonindex) const {
    
    std::vector<unsigned int> particles;
    GetMotherFSIParticles(particles, entry);
    bool pi0 = false, pim = false, pip = false;
    unsigned int n_pions = 0;
    int it_primary_pion = -1;
    int it_primary_proton = -1;
    for (std::vector<unsigned int>::iterator itP = particles.begin(); itP != particles.end(); ++itP) {
        const int mc_er_ID = m_chw->GetValue("mc_er_ID", entry, *itP);
        if ( mc_er_ID == 211 ) {
            it_primary_pion = *itP;
            pip = true; n_pions++;
        }
        else if (mc_er_ID == -211 ) {
            pim = true; n_pions++;
        }
        else if (mc_er_ID == 111) {
            pi0 = true; n_pions++;
        }
        else if (mc_er_ID == 2212) {
            it_primary_proton = *itP;
        }
    }
    
    int it_fs_pion = -999;
    int it_fs_proton = -999;
    GetFSParticles(it_fs_pion, it_fs_proton, entry);
    if(it_fs_pion>=0){
        printf("XSecLooper::FSIType found pion!!! %d\n", it_fs_pion); exit(1);
    }
    
    if(n_pions==0) {
        //with FSI, proton may be secondary, e.g. from pion absorption
        if(it_primary_proton<0){
            printf("no primary proton found in FSI particles!\n");
            exit(1);
        }
        
        const TVector3 selproton(m_chw->GetValue("mc_FSPartPx", entry, protonindex), m_chw->GetValue("mc_FSPartPy", entry, protonindex), m_chw->GetValue("mc_FSPartPz", entry, protonindex));
        const TVector3 priproton(m_chw->GetValue("mc_er_Px", entry, it_primary_proton), m_chw->GetValue("mc_er_Py", entry, it_primary_proton), m_chw->GetValue("mc_er_Pz", entry, it_primary_proton));
        /*
         //expected to have difference due to fiducial and highest momentum requirement
         const TVector3 fsproton(m_chw->GetValue("mc_er_Px", entry, it_fs_proton), m_chw->GetValue("mc_er_Py", entry, it_fs_proton), m_chw->GetValue("mc_er_Pz", entry, it_fs_proton));
         if(fabs(selproton.Mag()/fsproton.Mag()-1)>1E-6){
         printf("\nselected proton not equal to fs proton %f %f %f\n",selproton.Mag(), fsproton.Mag(), selproton.Mag()/fsproton.Mag()-1);
         printf("\nselected:\n");
         selproton.Print();
         printf("\nFS:\n");
         fsproton.Print();
         }
         */
        
        const double bindingE = 25;//in MeV
        const double deltaE = m_chw->GetValue("mc_FSPartE", entry, protonindex)-m_chw->GetValue("mc_er_E", entry, it_primary_proton)+ bindingE;
        
        /*
         const TVector3 delp=selproton-priproton;
         const double deltaT  = delp.Mag();
         const double deltaPx = fabs(delp.X());
         const double deltaPy = fabs(delp.Y());
         const double deltaPz = fabs(delp.Z());
         */
        /*
         const double deltaT  = delp.Mag()/priproton.Mag();
         const double deltaPx = fabs(delp.X()/priproton.X());
         const double deltaPy = fabs(delp.Y()/priproton.Y());
         const double deltaPz = fabs(delp.Z()/priproton.Z());
         */
        
        const double epsilon = 1.e-6;
        
        const bool iselastic = fabs(deltaE)<epsilon;
        //const bool isscaled = IsScaled(selproton, priproton);
        if(iselastic){
            static int kprintProtonNonInteracting = 0; PrintProton(ProtonNonInteracting, kprintProtonNonInteracting, entry, protonindex, it_primary_proton);
            return ProtonNonInteracting;
        }
        else if(deltaE>=epsilon){
            static int kprintProtonAcceleration = 0; PrintProton(ProtonAcceleration, kprintProtonAcceleration, entry, protonindex, it_primary_proton);
            return ProtonAcceleration;
        }
        else{
            static int kprintProtonDeceleration = 0; PrintProton(ProtonDeceleration, kprintProtonDeceleration, entry, protonindex, it_primary_proton);
            return ProtonDeceleration;
        }
    }
    //n_pions!=0
    else {
        return PionAbsorption;
    }
    
    std::cout<<"WARNING: Is not any FSI type!"<<std::endl;
    exit(1);
    
    return NOFSI;
}

//======================================================================
void XSecLooper::runLoop(int entries, const int fsitarget)
{
    if(m_nuPDGs.size()==0){
        cerr << "XSecLooper::runLoop: neutrino PDG not set. Use XSecLooper::setNuPDG(). Bailing" << endl;
        exit(1);
    }
    if(m_nuPDGs.size() != 1) std::cerr << "You are using " << m_nuPDGs.size() << " neutrino PDGs, but the flux will only be reweighted for PDG = " << *m_nuPDGs.begin() << ".  I'm not sure how to calculate a flux for multiple neutrino PDGs.\n";
    if(m_nUniv==0){
        cout<<"Not using the multi-verse GENIE error band :)"<<endl;
    }
    else{
        cout<<"Number of universes that will be used for GENIE error band: "<<m_nUniv<<endl;
    }
    
    cout << "Getting entries..." << flush;
    if(entries <= 0) entries = m_chw->GetEntries();
    cout << "done" << endl;
    cout << entries << " total entries" << endl;
    
    int iSignalEvent=0;
    int iCCEvent = 0;
    
    //This fluxreweighter matters which one you pick. There needs to be some sort of external way to decide whether you want a wiggle and what playlist you want... --Anne
    bool useNuEConstraint = true; //The old default behavior
    /*if(*m_nuPDGs.begin() < 0)
    {
      std::cerr << "Disabling the nu-e flux constraint because it wasn't officially ready for antineutrino analyses as of May 2021.  Remove line ~379-383 of XSecLooper.cxx when it's ready!\n";
      useNuEConstraint = false;
    }*/
    FluxReweighter* fluxReweighter = new FluxReweighter( *m_nuPDGs.begin(), useNuEConstraint, m_playlist, FluxReweighter::gen2thin, FluxReweighter::g4numiv6 );
    
    //std::ofstream fout("/scratch/minerva/jwolcott/geniexsecextract_vals.txt");
    for(int i=0; i<entries; ++i){
//        if(i % 100000 == 0) cout << (i/1000) << "k " << flush;
        
        const int mc_intType = m_chw->GetValue("mc_intType", i);
        
        int nuPDG=(int)m_chw->GetValue("mc_incoming", i);
        if(std::find(m_nuPDGs.begin(), m_nuPDGs.end(), nuPDG) == m_nuPDGs.end())
            continue;
        
        double Enu=1e-3*m_chw->GetValue("mc_incomingE", i);
        
        //using weights
        double wgt_val = fluxReweighter->GetFluxCVWeight(Enu, nuPDG );
        
        //Additional weights
        bool isGenieNonRes1pi = m_chw->GetValue("truth_genie_wgt_Rvn1pi",i,2) < 1.0 ||
        m_chw->GetValue("truth_genie_wgt_Rvp1pi",i,2) < 1.0;
        bool isGenieCCRes = m_chw->GetValue("mc_current",i) == 1 && m_chw->GetValue("mc_intType",i) == 2;
        double deuteriumMaRes          = 0.94;
        double maresval = m_chw->GetValue("truth_genie_wgt_MaRES",i,2);
        double genieMaRes1sig          = 0.2 * 1.12;
        double genieMaRes              = 1.12;
        bool doshift = false;
        // Reduced MvRES error from electroproduction data fit
        if(doshift){
            if(isGenieCCRes){
                wgt_val *= 1.15;
                wgt_val *= 1.0 + fabs( deuteriumMaRes - genieMaRes ) * ( maresval - 1.0 ) / genieMaRes1sig;
            }
        }
        //if NonResPion
	if(isGenieNonRes1pi) wgt_val *= 0.43;
        
        //2p2h and rpa weights
	wgt_val*=XSecLooper::Get2p2hWeight(i);
	wgt_val*=XSecLooper::GetRPAWeight(i);
        
        
        int run = (int)m_chw->GetValue("mc_run",i);
        int subrun = (int)m_chw->GetValue("mc_subrun",i);
        int gate = (int)m_chw->GetValue("mc_nthEvtInFile",i);
        
        
        // Only use CC events on carbon in tracker for normalization hist
        if(isCCRateEvent(*m_chw, i))
        {
            double cc_rate_wgt = fluxReweighter->GetFluxCVWeight(Enu, nuPDG );
            if (m_ccRateHists.find(nuPDG) == m_ccRateHists.end())
                m_ccRateHists[nuPDG] = new MnvH1D( Form("ccRateHist_%d", rand()), "", 1200, 0., 120. );
            m_ccRateHists[nuPDG]->Fill(Enu,cc_rate_wgt);
            //test printf("rfill %d %f %f\n", iCCEvent, Enu,wgt_val);
            iCCEvent++;
        }
        
        // Only use seleted fiducial events for the extracted xsec
	if(!isFiducial(*m_chw, i)) continue;
        
        //get signal weight (usually 1 unless reweighting to different xsec model)
        const double wgt_signal = getSignalWeight(*m_chw,i);
        
        //Creating the vector of GENIE weights
        vector<double> w_gen;
        if(m_nUniv!=0) w_gen = getGenV(*m_chw,i);
        
        for(unsigned int j=0; j<m_xsecs.size(); ++j){
            kprint = (j==0);
            
            XSec& xsec=*m_xsecs[j];
            
            if(!xsec.passesCuts(*m_chw, i))
                continue;

            bool isInEnergyRange = true;
	    if (xsec.isFluxIntegrated())
	      isInEnergyRange = Enu > xsec.getFluxIntEmin() && Enu < xsec.getFluxIntEmax();
	    else
	      isInEnergyRange = Enu > xsec.getXSecHist()->GetBinLowEdge(1) && Enu < (xsec.getXSecHist()->GetBinLowEdge(xsec.getXSecHist()->fN-1) + xsec.getXSecHist()->GetBinWidth(xsec.getXSecHist()->fN-1));
            
            if(!isInEnergyRange)
                continue;
            
            if(fsitarget!=-999){
                const int protonindex = XSec::GetProtonIndex(*m_chw,i);
                const int fsitype = FSIType(i, protonindex);
                if(fsitype!=fsitarget)
                    continue;
            }
            
            //fout << xsec.getName() << " ";
            if(xsec.getDimension()==1){
                const double val = xsec.getVariableValue(*m_chw, i);
                
                //fout << val << "\n";
                xsec.getXSecHist()->Fill(val,wgt_val*wgt_signal);
                
                //test if(j==XSec::kTdalphat){ printf("\niSignalEvent %d %f %f %f\n", iSignalEvent, m_chw->GetValue("mc_primFSLepton", i, 0), val,wgt_val*wgt_signal );}
                if(j==0) iSignalEvent++;
                if(m_nUniv!=0)xsec.getXSecHist()->FillVertErrorBand("GENIE",val,w_gen,wgt_val*wgt_signal);
            }
            else if(xsec.getDimension()==2){
                vector<double> var = xsec.getVariableValues(*m_chw, i);
                //fout << var.at(0) << " " << var.at(1) << "\n";
                xsec.get2DXSecHist()->Fill(var.at(0),var.at(1),wgt_val*wgt_signal);
		if(j==0) iSignalEvent++;
                if(m_nUniv!=0)xsec.get2DXSecHist()->FillVertErrorBand("GENIE",var.at(0),var.at(1),w_gen,wgt_val*wgt_signal);
            }
            else{
                cerr << "XSecLooper::runLoop: Couldn't determine X-Section dimension(1D,2D...). Use setDimension() for that. Bailing" << endl;
                exit(1);
            }
        } // for (j)
        //should not put things here because if cut fails, j loop is broken!
    } // for (i)
    cout << "\nTotal CC events filled: "<<iCCEvent<<" and total signale events filled: "<<iSignalEvent<<endl;
    
//    m_ccRateHists[14]->SaveAs("CCRateHist.root");
    
    //fout.close();
    
    // Divide CC rate histograms through by the bin width
    for (std::map<int, PlotUtils::MnvH1D*>::iterator it_rateHist = m_ccRateHists.begin(); it_rateHist != m_ccRateHists.end(); ++it_rateHist){
        divideBinsByWidth( it_rateHist->second );
	//        it_rateHist->second->Write();
    }
    
    // Normalize all xsec histograms
    for(unsigned int j=0; j<m_xsecs.size(); ++j){
        kprint = (j==0);
        
        XSec& xsec=*m_xsecs[j];
        const double normFactor=getNormFactor(xsec);
        //const double normFactor = 1.0;
        MnvH1D *xsecHist = NULL;
        MnvH2D *xsecHist2D = NULL;
        if(xsec.getDimension()==1)
        {
            xsecHist=xsec.getXSecHist();
            
            // save the event rate first, in case the user needs to look at it
            xsec.setEvRateHist( dynamic_cast<PlotUtils::MnvH1D*>(xsecHist->Clone( (std::string(xsecHist->GetName()) + "_evRate").c_str() )) );
        }
        else if( xsec.getDimension()==2)
	  {
            xsecHist2D=xsec.get2DXSecHist();
	    xsecHist2D->SaveAs(Form("%s_2Devtrate.root",xsecHist2D->GetName()));
	  }
        else{
            //Just in case
            cerr << "XSecLooper::runLoop: Couldn't determine X-Section dimension(1D,2D...). Use setDimension() for that. Bailing" << endl;
            exit(1);
        }
        
        //printf("XSecLooper raw %s bin 2: %e bin 10: %e\n", xsecHist->GetName(), xsecHist->GetBinContent(2), xsecHist->GetBinContent(10));
        
        if(xsec.getDimension()==1){
            xsecHist->Scale(normFactor);//xsec had carbon in denorminator, so the normFactor is direct multiplied.
//            if(kprint)
                printf("normFactor %f 1/normFactor %f\n", normFactor, 1/normFactor);
        }
        else if(xsec.getDimension()==2)
            xsecHist2D->Scale(normFactor);
        
        if(xsec.isFluxIntegrated()){
            // get_flux_hist() returns a histogram that is essentially
            // the flux*pot*(number of C12), and the flux-integrated
            // cross section is (in Q^2, say):
            //
            // N(Q^2) / ( (number of targets)*pot*(flux integral)*(bin width) )
            //
            // so we just need to scale by the ratio of number of
            // targets to number of C12, the flux histogram integral,
            // and the bin width
            if(xsec.getDimension()==1)
                divideBinsByWidth(xsecHist);
            else if(xsec.getDimension()==2)
                divideBinsByWidth(xsecHist2D);
            
            const double Emin=xsec.getFluxIntEmin();
            const double Emax=xsec.getFluxIntEmax();
            const double fluxTimesNC12=getFluxHist()->Integral(m_ccRateHists.begin()->second->FindFixBin(Emin),
                                                               m_ccRateHists.begin()->second->FindFixBin(Emax),
                                                               "width");
            
            MnvH1D* flux = (MnvH1D*)getFluxHist(); //Get the CCRate Hist Flux
            

            
//            if(kprint)
                printf("fluxTimesNC12 %e Emin %e Emax %e\n", fluxTimesNC12, Emin, Emax);
            
            if(xsec.getDimension()==1){
              MnvH1D* Integrated_flux = fluxReweighter->GetIntegratedFluxReweighted_FromInputFlux(flux, xsecHist, Emin, Emax);
              xsecHist->Divide(xsecHist,Integrated_flux, 1., 1.);
            }else if( xsec.getDimension()==2){
                xsecHist2D->Scale(1.0/fluxTimesNC12); //Haven't dealt with any sort of 2D reweighter...
            }
            else{
                //Just in case
                cerr << "XSecLooper::runLoop: Couldn't determine X-Section dimension(1D,2D...). Use setDimension() for that. Bailing" << endl;
                exit(1);
            }
        }
        
        else{
            if(xsec.getDimension()==1)
                normalizeEnuHist(xsecHist);
            else if(xsec.getDimension()==2){
                //!@todo
                cerr << "XSecLooper::runLoop: Flux-unfolded 2D Cross sections method is not implemented yet. Bailing" << endl;
                exit(1);
            }
        }
    }
    
    
}

//======================================================================
void XSecLooper::divideBinsByWidth(PlotUtils::MnvH1D* h)
{
    //todo: maybe respect bin normalization width
    h->Scale( 1., "width" );
}

//======================================================================
void XSecLooper::divideBinsByWidth(PlotUtils::MnvH2D* h)
{
    //todo: maybe respect bin normalization width
    h->Scale( 1., "width" );
}

//======================================================================
double XSecLooper::getNormFactor(const XSec& xsec)
{
    const XSec::ENormType normType = xsec.getNormalizationType();
    
    // if the XSec is storing its own normalization return it
    if( XSec::kSelfNorm == normType )
//        cout<<"Anne!!! Normalization "<<xsec.getNormalizationValue()<<endl;
        return xsec.getNormalizationValue();
    
    // if normalization doesn't matter return 1.
    if(normType==XSec::kNoNorm) return 1;
    
    
    //otherwise the fiducial volume is the same that was used for the normalization.
    //note that the norm factor will be wrong if you implement your own isFiducial
    //  unless you use self normalization.
    
    
    // Hold details of the nuclei that make up the Minerva tracker and
    // their mass and number fractions. Numbers are taken from the
    // spreadsheet at
    // https://docs.google.com/spreadsheet/ccc?key=0AtkEbkjYmJo6dE92Ql9ZSlF6OFhBZV9oeEJWMmNCM2c#gid=0
    // which I used to calculate the mass and number fractions
    struct Nucleus
    {
        const char* name;
        int Z;
        double A;
        double massFraction;
        double numberFraction;
    };
    
    //note: set number fraction right after construction of Nuclei
    //      this makes it easier to change the mass fractions
    
    //=== Using mass fractions from material assay
    /*
     Nucleus nuclei[kNNuclei]={
     {"C12",   6, 12.011, 0.8760, 0. },
     {"H1",    1,  1.008, 0.0747, 0. },
     {"O16",   8, 15.999, 0.0315, 0. },
     {"Ti48", 22, 47.867, 0.0069, 0. },
     {"Al27", 13, 26.982, 0.0025, 0. },
     {"Si28", 14, 28.085, 0.0030, 0. },
     {"Cl35", 17, 35.450, 0.0054, 0. }
     }; //==== END mass assay fractions
     */
    
    /*
     //==== Using mass fractions from Titan (recommended)
     Nucleus nuclei[kNNuclei]={
     {"C12",   6, 12.011, 0.875264, 0. },
     {"H1",    1,  1.008, 0.0824587, 0. },
     {"O16",   8, 15.999, 0.029275, 0. },
     {"Ti48", 22, 47.867, 0.0043074, 0. },
     {"Al27", 13, 26.982, 0.00176433, 0. },
     {"Si28", 14, 28.085, 0.00183692, 0. },
     {"Cl35", 17, 35.450, 0.00509393, 0. }
     }; //==== END titan MC mass fractions
     */
    //use Eroica info 20171212 Xianguo
    /*
    Nucleus nuclei[]={
        {"C12",   6, EroicaAtomicMass::C,  EroicaMassFractionMC::C,  0. },
        {"H1",    1, EroicaAtomicMass::H,  EroicaMassFractionMC::H,  0. },
        {"O16",   8, EroicaAtomicMass::O,  EroicaMassFractionMC::O,  0. },
        {"Ti48", 22, EroicaAtomicMass::Ti, EroicaMassFractionMC::Ti, 0. },
        {"Al27", 13, EroicaAtomicMass::Al, EroicaMassFractionMC::Al, 0. },
        {"Si28", 14, EroicaAtomicMass::Si, EroicaMassFractionMC::Si, 0. },
        {"Cl35", 17, EroicaAtomicMass::Cl, EroicaMassFractionMC::Cl, 0. }
    };
    */

    Nucleus nuclei[]={
        {"C12",   6, NXAtomicMass::C,  NXMassFractionMC::C,  0. },
        {"H1",    1, NXAtomicMass::H,  NXMassFractionMC::H,  0. },
        {"O16",   8, NXAtomicMass::O,  NXMassFractionMC::O,  0. },
        {"Ti48", 22, NXAtomicMass::Ti, NXMassFractionMC::Ti, 0. },
        {"Al27", 13, NXAtomicMass::Al, NXMassFractionMC::Al, 0. },
        {"Si28", 14, NXAtomicMass::Si, NXMassFractionMC::Si, 0. },
        {"Cl35", 17, NXAtomicMass::Cl, NXMassFractionMC::Cl, 0. }
    };
 //==== END titan MC mass fractions
    
    const int kNNuclei=sizeof(nuclei)/sizeof(Nucleus);
    if(kprint) printf("kNNuclei %d\n", kNNuclei);
    
    double totMassF=0;
    //set the numberFraction
    {
        //one loop to get total atoms in gram of tracker
        double totalAtomsPerG = 0.;
        for(int i=0; i<kNNuclei; ++i){
            const Nucleus& nucleus=nuclei[i];
            const double atomsPerG = nucleus.massFraction / nucleus.A;
            totMassF += nucleus.massFraction;
            totalAtomsPerG += atomsPerG;
        }
        //another loop to set numberFraction
        for(int i=0; i<kNNuclei; ++i){
            Nucleus& nucleus=nuclei[i];
            const double atomsPerG = nucleus.massFraction / nucleus.A;
            nucleus.numberFraction = atomsPerG / totalAtomsPerG;
            
            if(kprint){
                std::cout << " Nucleus: " << nucleus.name << ", massFrac = " << nucleus.massFraction << ", numberFrac = " << nucleus.numberFraction << std::endl;
            }
        }
    }
    if(kprint) printf("Total Mass Fraction -1 = %e\n", totMassF-1);
    
    double ret=0;
    
    for(int i=0; i<kNNuclei; ++i){
        const Nucleus& nucleus=nuclei[i];
        
        switch(normType){
            case XSec::kPerNucleon:
                ret+=nucleus.A*nucleus.numberFraction;
                break;
            case XSec::kPerNeutron:
                ret+=(nucleus.A-nucleus.Z)*nucleus.numberFraction;
                break;
            case XSec::kPerProton:
                ret+=nucleus.Z*nucleus.numberFraction;
                break;
            default:
                assert(0 && "Unknown cross section normalization");
        }
    }
    
    // normalizeEnuHist normalizes the histogram to the carbon cross
    // section, so we need the number of {neutrons, protons, nucleons}
    // per carbon atom (nuclei[0] is C12)
    const double normFactor = nuclei[0].numberFraction/ret;
    
    if(kprint) printf("getNormFactor: nuclei[0].numberFraction %e ret %e normFactor %e\n", nuclei[0].numberFraction, ret, normFactor);
    return normFactor;
}

//======================================================================
TH1* XSecLooper::getFluxHist()
{
    if(!m_fluxHist)
    {
        //only 1 flux from the CV
        m_fluxHist = dynamic_cast<TH1D*>(m_ccRateHists.begin()->second->GetCVHistoWithStatError().Clone("hFlux"));
        m_fluxHist->Reset();
        
        std::map<int, GraphSpline*> grCC;
        for (std::set<int>::const_iterator it_pdg = m_nuPDGs.begin();
             it_pdg != m_nuPDGs.end();
             ++it_pdg)
            grCC.insert( std::make_pair(*it_pdg, new GraphSpline(getCCSpline(*it_pdg))) );  // default is for nu
        
        std::map<int, TF1*> grCCSpline;
        std::stringstream ss;
        for (std::map<int, GraphSpline*>::iterator it_graph = grCC.begin();
             it_graph != grCC.end();
             ++it_graph)
        {
            ss.str("");
            ss << it_graph->first;
            grCCSpline.insert( std::make_pair(it_graph->first, new TF1( ("f1gcc" + ss.str()).c_str(), it_graph->second, 0, 120, 0, "GraphSpline" )) );
        }
        
        double totRate = 0, totSpInt = 0, totFlux=0;
        for(int i=1; i<m_fluxHist->GetNbinsX()+1; ++i)
        {
            const double binLowEdge=m_fluxHist->GetBinLowEdge(i);
            const double binWidth=m_fluxHist->GetBinWidth(i);
            double flux = 0;
            double err = 0;
            for (std::map<int, PlotUtils::MnvH1D*>::iterator it_rateHist = m_ccRateHists.begin();
                 it_rateHist != m_ccRateHists.end();
                 ++it_rateHist)
            {
                // avg XSec = integral(xsec)/width
                int pdg = it_rateHist->first;
                PlotUtils::MnvH1D * h_rate = it_rateHist->second;
                TF1 * spline = grCCSpline[pdg];
                double xsecInt = spline->Integral(binLowEdge, binLowEdge+binWidth);
                double xsec = xsecInt/binWidth;
                
                if( xsec < 1E-6 )
                    std::cout << " Warning, xsec is 0 for bin from " << binLowEdge << " - " <<  binLowEdge+binWidth << " for spline for PDG " << pdg << std::endl;
                else
                {
                    //(Note that the more widespread cross section units, 10−38cm2, are used when the cross section data are exported to a ROOT format for inclusion in user analysis code.
                    const double splXsecUnit = 1E-38;
                    xsec *= splXsecUnit;
                    
                    flux += h_rate->GetBinContent(i) / xsec;
                    err += h_rate->GetBinError(i) / xsec;  // I think this is right, since division is distributive and the propagation of an additive error is additive?
                    
                    totRate += h_rate->GetBinContent(i);
                    totSpInt += xsec;
                    totFlux += h_rate->GetBinContent(i)/xsec;
                }
                //test static int ii=0; printf("ii %d binLowEdge %f binWidth %f xsecInt %e ", ii, binLowEdge, binWidth, xsecInt); ii++;
            }
            
            //test printf("flux %e\n", flux);
            
            m_fluxHist->SetBinContent(i, flux);
            m_fluxHist->SetBinError(i, err);
        } // for (i)
        printf("totRate %f totSpInt %e totFlux %e\n", totRate, totSpInt, totFlux);
    }  // if (!m_fluxHist)
    
//    m_fluxHist->SaveAs("flux_histogram_from_looper.root");
    
    return m_fluxHist;
}

//======================================================================
//void XSecLooper::normalizeEnuHist(PlotUtils::MnvH1D* h) //Old version... didn't use the same rebinner as the fluxreweighter
//{
//    HistSpline* fluxHistSpline=new HistSpline(getFluxHist());
//    TF1 spline("f1sdf", fluxHistSpline, 0, 120, 0, "HistSpline");
//
//    //use the same bins as the histogram we are normalizing
//    fluxAvg = new TH1D( *dynamic_cast<TH1D*>(h) );
//
//    fluxAvg->SetName("fluxAvg");
//    fluxAvg->Reset();
//
////    SplineAvg = new TH1D( *dynamic_cast<TH1D*>(h) );
//
//
//    for(int i=1; i<h->GetNbinsX()+1; ++i){
//        const double binLowEdge=h->GetBinLowEdge(i);
//        const double binWidth=h->GetBinWidth(i);
//        const double fluxIntegral = spline.Integral(binLowEdge, binLowEdge+binWidth);
//
//        /*cout << "Flux Average: " << endl;
//         cout << i << " " << fluxIntegral << endl;*/
//        fluxAvg->SetBinContent( i, fluxIntegral );
//
//        //assume that flux contributes 0 error.
//        //this is not correct when the xsec and normalization histograms use different events.
//        //this error should be stat error on noramlization sample,
//        //  but I'm not sure how to get it.
//        //The error will be small and cancel perfectly in xsec ratios so I don't care about it.
//        fluxAvg->SetBinError( i, 0. );
//    }
//    //  TCanvas c;
//    //  fluxAvg->Draw();
//    //  c.Print("flux.png");
//    //
//    //  h->Draw();
//    //  c.Print("rate.png");
//    fluxAvg->SaveAs("flux_for_normalization.root");
//
//    h->DivideSingle( h, fluxAvg );
//}

//======================================================================
void XSecLooper::normalizeEnuHist(PlotUtils::MnvH1D* h)
{

   

    FluxReweighter* fluxReweighter = new FluxReweighter( 14, false, m_playlist, FluxReweighter::gen2thin, FluxReweighter::g4numiv6 );
  

    TH1D* flux = (TH1D*)getFluxHist();

    MnvH1D* hFluxrebin = (MnvH1D*)fluxReweighter->GetRebinnedFluxReweighted_FromInputFlux(flux, h);
   // hFluxrebin->SaveAs("Rebinned_Flux_XSecLooper_FluxReweighter.root");
//    hFluxrebin->SaveAs("flux_for_normalization.root");
    h->DivideSingle( h, hFluxrebin );
}




//======================================================================
TGraph* XSecLooper::getCCSpline(int nuPDG)
{
    TString splineFileName(gSystem->Getenv("MPARAMFILESROOT"));
    splineFileName+="/data/GENIE/spline_files/gxspl-nuclear-MINERVA_Full_v2126.root";
    TFile splineFile(splineFileName);
    if(splineFile.IsZombie()){
        cerr << "Can't find GENIE spline file at:" << endl;
        cerr << splineFileName << endl;
        cerr << "Make sure $MPARAMFILESROOT points at a nightly/release dated after 7 Dec 2012" << endl;
        cerr << "Exiting" << endl;
        exit(1);
    }
    
    TString splineName(m_rateSplineName);
    //if no spline name is set, form the usual as total CC on C12
    if( 0 == splineName.Length() )
    {
        splineName = "nu_";
        
        switch(abs(nuPDG)){
            case 14:
                splineName+="mu";
                break;
            case 12:
                splineName+="e";
                break;
            default:
                cerr << "XSecLooper::getCCSpline(): nu PDG is set to " << nuPDG
                << ", which we don't have splines for. Bailing" << endl;
                exit(1);
        }
        if(nuPDG<0) splineName+="_bar";
        
        splineName+="_C12/tot_cc";
    } // if (splineName.
    
    TGraph* gr=dynamic_cast<TGraph*>(splineFile.Get(splineName));
    //    cout << " got spline: " << splineName << endl;
    if(!gr){
        cerr << "Can't find spline " << splineName << " in file " << splineFileName << ". Bailing." << endl;
        exit(1);
    }
    return gr;
}

//======================================================================
bool XSecLooper::isFiducialTracker(ChainWrapper& chw, int entry)
{

  if(m_useFiducialBranch) return (bool)chw.GetValue("truth_is_fiducial",entry);

  //else  
  // fiducual volume is mod 27-79 inclusive inside 850mm apothem
  double z=chw.GetValue("mc_vtx", entry, 2);
  //  int true_module   = (int)chw.GetValue("truth_vtx_module",entry); //Dan, this is the hacky bit where I'd want to get the right volume for the scattering center... --Anne
  //    if( true_module <27 || true_module > 79 ) return false;
  if(z<m_fiducialZMin || z>m_fiducialZMax) return false;
  double x=fabs(chw.GetValue("mc_vtx", entry, 0));
  double y=fabs(chw.GetValue("mc_vtx", entry, 1));
  
  if(x*x + y*y < m_apothem*m_apothem) return true;
  
  double lenOfSide = m_apothem * ( 2 / sqrt(3) );
  
  if( x > m_apothem )
      return false;
  
  if( y < lenOfSide/2.0 )
      return true;
  
  double slope = (lenOfSide / 2.0) / m_apothem;
  if( y < lenOfSide - x*slope )
      return true;
  
  return false;
}

//=====================================================================
bool XSecLooper::isCCRateEvent(ChainWrapper& chw, int entry)
{
    //requires a CC event...
    int current=(int)chw.GetValue("mc_current", entry);
    if( 1 != current ) return false;
    
    //... on carbon ...
    int target=(int)chw.GetValue("mc_targetNucleus", entry);
    if( 1000060120 != target ) return false;
    
    int type=(int)chw.GetValue("mc_intType",entry);
    
    //..... not 2p2h for now....
    if(type==8) {
        //This is hacky and verbose, but I don't want people to forget... --Anne
      //        std::cout<<entry<<" I'm skipping 2p2h events until the GENIE Splines are fixed. If you have fixed the GENIE splines, fix the isCCRateEvent function"<<std::endl;
        return false;
    }
    //... in the tracker fiducial region.
    if( !isFiducialTracker(chw, entry) ) return false;
    
    return true;
}

//=====================================================================
bool XSecLooper::isFiducial(ChainWrapper& chw, int entry)
{
    return isFiducialTracker( chw, entry );
}
//===================================================================
vector<double> XSecLooper::getGenV(ChainWrapper& chw, int entry)
{
    vector<double> genie_wgts;
    for(int i=0; i<m_nUniv; ++i){
        double univ =(double)chw.GetValue("mc_wgt_GENIE", entry, i);
        genie_wgts.push_back(univ);
    }
    return genie_wgts;
}
//======================================================================
double XSecLooper::getSignalWeight(ChainWrapper& chw, int entry)
{
    //this next line is just an example to avoid compiler warnings, not used
    (double)chw.GetValue("mc_cvweight_totalXsec", entry );
    
    //remember weight is relative to the existing CV weight
    return 1.;
}
