#include "GENIEXSecExtract/XSec.h"

#include "TH1.h"
#include "TVector3.h"
#include "TMath.h"

#include <iostream>
#include <cmath>
#include <cassert>

#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

#include <algorithm>


using namespace std;
using namespace PlotUtils;

bool Sort(TVector3 const &lhs, TVector3 const &rhs)
{
  return (lhs.Mag() > rhs.Mag());
}

// Convenience function for squaring things
template<class T>
T sqr(T x) { return x*x; }

//======================================================================
    XSec::XSec(const char* name)
: m_name(name), m_xsecHist(0), m_evRateHist(NULL), m_2DxsecHist(0), m_variable_x(kNullVariable), m_variable_y(kNullVariable), m_isFluxIntegrated(false),
    m_Emin(-1), m_Emax(200), m_Univ(100), m_normType(kNoNorm), m_normValue(-1), m_dimension(1)
{

}

//======================================================================
void XSec::setUniformBinning(int nbins, double xmin, double xmax)
{
    if(m_xsecHist){
        cout << "XSec::setUniformBinning " << m_name << " deleting existing histogram. This may or may not be what you want" << endl;
        delete m_xsecHist;
    }
    m_xsecHist=new MnvH1D(Form("%s_xsec", m_name.c_str()), "", nbins, xmin, xmax);
    if(m_Univ!=0)m_xsecHist->AddVertErrorBand("GENIE",m_Univ);
}

//======================================================================
void XSec::setUniformBinning(int x_nbins, double xmin, double xmax, int y_nbins, double ymin, double ymax)
{
    if(m_2DxsecHist){
        cout << "XSec::setUniformBinning " << m_name << " deleting existing histogram. This may or may not be what you want" << endl;
        delete m_2DxsecHist;
    }
    m_2DxsecHist=new MnvH2D(Form("%s_xsec", m_name.c_str()), "", x_nbins, xmin, xmax, y_nbins, ymin, ymax);
    if(m_Univ!=0)m_2DxsecHist->AddVertErrorBand("GENIE",m_Univ);
}

//======================================================================
void XSec::setBinEdges(int nbins, double* binEdges)
{
    if(m_xsecHist){
        cout << "XSec::setBinEdges " << m_name << " deleting existing histogram. This may or may not be what you want" << endl;
        delete m_xsecHist;
    }
    m_xsecHist=new MnvH1D(Form("%s_xsec", m_name.c_str()), "", nbins, binEdges);
    if(m_Univ!=0)m_xsecHist->AddVertErrorBand("GENIE",m_Univ);
}

//======================================================================
void XSec::setBinEdges(int x_nbins, double* x_binEdges, int y_nbins, double* y_binEdges)
{
    if(m_2DxsecHist){
        cout << "XSec::setBinEdges " << m_name << " deleting existing histogram. This may or may not be what you want" << endl;
        delete m_2DxsecHist;
    }
    m_2DxsecHist=new MnvH2D(Form("%s_xsec", m_name.c_str()), "", x_nbins, x_binEdges, y_nbins, y_binEdges);
    if(m_Univ!=0)m_2DxsecHist->AddVertErrorBand("GENIE",m_Univ);
}

//======================================================================
vector<int> XSec::getFSPDGs(ChainWrapper& chw, int entry)
{
    const int nParticles=(int)chw.GetValue("mc_nFSPart", entry);
    vector<int> ret;
    for(int i=0; i<nParticles; ++i){
        int pdg=(int)chw.GetValue("mc_FSPartPDG", entry, i);
        ret.push_back(pdg);
    }
    return ret;
}

//======================================================================
double XSec::getVariableValue(ChainWrapper& chw, int entry)
{
    return getValue(chw, entry, m_variable_x);
}

//======================================================================
std::vector<double> XSec::getVariableValues(ChainWrapper& chw, int entry)
{
    std::vector<double> vec;
    vec.push_back( getValue( chw, entry, m_variable_x ) );
    vec.push_back( getValue( chw, entry, m_variable_y ) );

    return vec;
}

//======================================================================
int XSec::GetProtonIndex(ChainWrapper& chw, const int entry)
{
  TVector3 protonfullp;
  int protonindex = -999;
  const int nparticles = (int)chw.GetValue("mc_nFSPart",entry);
  for(int tmpipar = 0; tmpipar < nparticles; tmpipar++) {
    const int pdg = (int)chw.GetValue("mc_FSPartPDG",entry,tmpipar);
    if(pdg!=2212){
      continue;
    }

    const double px = (double)chw.GetValue("mc_FSPartPx",entry,tmpipar);
    const double py = (double)chw.GetValue("mc_FSPartPy",entry,tmpipar);
    const double pz = (double)chw.GetValue("mc_FSPartPz",entry,tmpipar);
    const TVector3 pvec(px, py, pz);

    const double tmpmom = pvec.Mag();
    const double tmptheta = getTheta(px,py,pz);

    if( tmptheta < proton_angle_MAXCut() &&
        tmpmom > proton_momentum_minCut() &&
        tmpmom < proton_momentum_MAXCut()){

      if(tmpmom>protonfullp.Mag()){
        protonindex = tmpipar;
        protonfullp = pvec;
      }
    }
  }
  return protonindex;
}

double XSec::getTheta( double x, double y, double z )
{
  double pyprime = -1.0*TMath::Sin(-0.05887)*z + TMath::Cos(-0.05887)*y;
  double pzprime =  1.0*TMath::Cos(-0.05887)*z + TMath::Sin(-0.05887)*y;
  double den2    = x*x + pyprime*pyprime + pzprime*pzprime;

  if( den2 == 0 ) return -9999;
  else return TMath::ACos( pzprime / sqrt(den2) )*180/TMath::Pi();
}

double XSec::getValue(ChainWrapper& chw, int entry, XSec::EVariable var)
{
    // Some particle masses in GeV
    static const double Me=0.00051100;
    static const double Mmu=0.1056583;
    static const double Mpi0=0.135;
    static const double Mpip=0.13957;

    const double Mlep = abs(chw.GetValue("mc_incoming", entry)) == 14 ? Mmu : Me;

    // Which sort of pion are we dealing with (if we're dealing with a pion)?
    int piPDG;
    double Mpi;
    switch(var){
        case kPPi0:
        case kTPi0:
        case kThetaPi0:
        case kCosPi0:
            piPDG=111;
            Mpi=Mpi0;
            break;
        case kPPiPlus:
        case kTPiPlus:
        case kCosPiPlus:
        case kThetaPiPlus:
            piPDG=211;
            Mpi=Mpip;
            break;
        case kThetaChargedPion:
        case kTChargedPion:
            Mpi=Mpip;
            piPDG=211;
            break;
        default:
            piPDG=-1;
            Mpi=-1;
    }

    int piIndex=-1;
    // Find the pion in the list of final state particles, if we're looking at a pion variable
    if(piPDG>0){
        vector<int> fsPDGs=getFSPDGs(chw, entry);
        double leading_E = -1.0;
        for(unsigned int i=0; i<fsPDGs.size(); ++i){
            if (var == XSec::kThetaChargedPion || var == XSec::kTChargedPion) {
                if (abs(fsPDGs[i]) == piPDG) {
                    //select leading charged pion by energy
                    double E = (double)chw.GetValue("mc_FSPartE", entry, i);
                    if (E > leading_E) {	
                        piIndex=i;
                        if (leading_E >= 0.0) cout<<"More than one charged pion in entry "<<entry<<endl;
                        leading_E = E;

                    }
                }
            }
            else if(fsPDGs[i]==piPDG){
                piIndex=i;
                break;
            }
        }
    }

    // Make sure we found the pion if we're trying to get a pion variable
    assert(piIndex!=-1 || piPDG==-1);

    if (piPDG > 0 && piIndex==-1) cout<<"assert failed! Pion index: "<<piIndex<<" PDG: "<<piPDG<<endl;

    if(var>kTransverseBegin && var<kTransverseEnd){
      //==================================== reading muon ================================================
      const double muonpx = chw.GetValue("mc_primFSLepton", entry, 0);
      const double muonpy = chw.GetValue("mc_primFSLepton", entry, 1);
      const double muonpz = chw.GetValue("mc_primFSLepton", entry, 2);

      TVector3 muonfullp(muonpx, muonpy, muonpz);
      muonfullp *= 1E-3;//MeV to GeV
      if(var== kTmuonmomentum){
        return muonfullp.Mag();
      }

      if(var== kTmuontheta){
        return getTheta(muonpx,muonpy,muonpz);//in Deg
      }

      //==================================== reading proton ================================================
      const int protonindex = GetProtonIndex(chw, entry);
      const double px = (double)chw.GetValue("mc_FSPartPx",entry,protonindex);
      const double py = (double)chw.GetValue("mc_FSPartPy",entry,protonindex);
      const double pz = (double)chw.GetValue("mc_FSPartPz",entry,protonindex);
      TVector3 protonfullp(px, py, pz);
      const double protontheta=getTheta(px,py,pz);

      protonfullp *= 1E-3;//MeV to GeV
      if(var== kTprotonmomentum){
        return protonfullp.Mag();
      }

      if(var== kTprotontheta){
        return protontheta;
      }
      //==================================== calculate single-transverse ================================================
      const double neutrinopx = chw.GetValue("mc_incomingPartVec", entry, 0);
      const double neutrinopy = chw.GetValue("mc_incomingPartVec", entry, 1);
      const double neutrinopz = chw.GetValue("mc_incomingPartVec", entry, 2);
      TVector3 tmpneutrino(neutrinopx, neutrinopy, neutrinopz);
      const TVector3 unitneutrino = tmpneutrino.Unit();

      const TVector3 plmuon   = unitneutrino *   muonfullp.Dot(unitneutrino);
      const TVector3 plproton = unitneutrino * protonfullp.Dot(unitneutrino);

      const TVector3 pTmuon   = muonfullp   - plmuon;
      const TVector3 pTproton = protonfullp - plproton;
      const TVector3 vdPt     = pTmuon + pTproton;

      if(var== kTdpt){
        return vdPt.Mag();//in GeV
      }

      const TVector3 unitqt = -pTmuon.Unit();
      if(var== kTdphit){
        return TMath::ACos(pTproton.Dot(unitqt)/pTproton.Mag())*TMath::RadToDeg(); //in Deg                                                                                                                                                                      
      }

      if(var== kTdalphat){
        return TMath::ACos(    vdPt.Dot(unitqt)/vdPt.Mag()    )*TMath::RadToDeg(); //in Deg
      }

      //==================================== calculate neutron momentum ================================================
      const Double_t muonmass = MuonMass();//GeV

      //modified from original codes by J. Sobczyk 14 Nov 2017
      //https://journals.aps.org/prc/abstract/10.1103/PhysRevC.95.065501
      //Eq. 11
      //all in GeV

      const double Mp = ProtonMass();//GeV                                                                                                                                                                                                                 
      const double Mn = NeutronMass();//GeV                                                                                                                                                                                                                
      const double MA = 6*Mn + 6*Mp - 92.162/1E3;//GeV
      const double Bin=27.13/1E3;//GeV //average excitation energy instead of sampling it as in the paper
      const double MAstar = MA - Mn + Bin; //GeV
      const double Epprim = sqrt(Mp*Mp + protonfullp.Mag2());
      const double kprimL = plmuon.Mag();
      const double pprimL = plproton.Mag();
      const double Eprim = sqrt(muonmass*muonmass+muonfullp.Mag2());
      const double factor = MA - Eprim - Epprim + kprimL + pprimL;
      const double pT = vdPt.Mag();
      const double pL = -(MAstar*MAstar + pT*pT-factor*factor)/2.0/factor;

      if(var== kTneutronmomentum){
        return sqrt(pL*pL + pT*pT);
      }
    }

    if(var == kTdpTT){
      if(mu_mom.size() > 0 && pr_mom.size() > 0 && pi0_mom.size() > 0){        
        const double nu_px = chw.GetValue("mc_incomingPartVec", entry, 0);
        const double nu_py = chw.GetValue("mc_incomingPartVec", entry, 1);
        const double nu_pz = chw.GetValue("mc_incomingPartVec", entry, 2);
        TVector3 neutrino(nu_px, nu_py, nu_pz);
        return GetdpTT(neutrino, mu_mom[0], pr_mom[0], pi0_mom[0]);
      }
      cout << "WARNING: Bad TVector3 in dpTT Calculation. :-(" << endl;
      return -9999.;
    }

    switch(var){
        case kx:
            return chw.GetValue("mc_Bjorkenx", entry );
        case kxExp:
            {
                static double Mn = 1.0;
                //neutron
                if( (int)chw.GetValue("mc_targetNucleon", entry) == 2112)
                    Mn = 0.93956;
                //proton
                else if( (int)chw.GetValue("mc_targetNucleon", entry) == 2212)
                    Mn = 0.938272;
                else
                    Mn = (1.5*0.93956 + 0.938272 )/2.5; //weighted average because xsec is bigger on n
                const double Q2 = getValue(chw, entry, kExpQ2 );
                const double Ehad = 1e-3*( chw.GetValue("mc_incomingE", entry) - chw.GetValue("mc_primFSLepton", entry, 3) );
                return Q2 / (2*Mn*Ehad);
            }  
        case ky:
            return chw.GetValue("mc_Bjorkeny", entry );
        case kW:
            return 1e-3*chw.GetValue("mc_w", entry );
        case kExpW:
            {
                double Mn = 1.0;
                //neutron
                if( (int)chw.GetValue("mc_targetNucleon", entry) == 2112)
                    Mn = 0.93956;
                //proton
                else if( (int)chw.GetValue("mc_targetNucleon", entry) == 2212)
                    Mn = 0.938272;
                else
                    Mn = (1.5*0.93956 + 0.938272 )/2.5; //weighted average because xsec is bigger on n  
                const double Q2 = getValue(chw, entry, kExpQ2 );
                const double Ehad = 1e-3*( chw.GetValue("mc_incomingE", entry) - chw.GetValue("mc_primFSLepton", entry, 3) );
                double W2 = pow(Mn, 2) + 2*Mn*Ehad - Q2;
                if( W2 > 0)	
                    return sqrt(W2);
                else
                    return -1.0;  

            }
        case kExpQ2:
            {
                double theta = chw.GetValue("truth_muon_theta", entry);	
                const double Q2 = 1e-6*4*( chw.GetValue("mc_incomingE", entry)*chw.GetValue("mc_primFSLepton", entry, 3) )*pow(sin(theta/2), 2);	 
                return Q2;  

            }  
       case kq3:
          {
	    double Q2 = 1e-6*chw.GetValue("mc_Q2", entry);
	    double Ehad = 1e-3*( chw.GetValue("mc_incomingE", entry) - chw.GetValue("mc_primFSLepton", entry, 3) );
	    return sqrt(Q2 + (Ehad*Ehad));
          }
        case kQ2:
            return 1e-6*chw.GetValue("mc_Q2", entry);
        case kQ2QE:
            {

                const double EnuQE=getValue(chw, entry, kENuQE);
                const double Elep=1e-3*chw.GetValue("mc_primFSLepton", entry, 3);
                const double coslep=getValue(chw, entry, kCosLep);

                return 2*EnuQE*(Elep - sqrt(sqr(Elep) - sqr(Mlep))*coslep) - sqr(Mlep);
            }
        case kQ2QEProton:
            {
                double max_p = 0;
                const int nParticles = (int)chw.GetValue("mc_nFSPart",entry);
                for(int i = 0; i < nParticles; ++i) {
                    if( (int)chw.GetValue("mc_FSPartPDG",entry,i) == 2212 ) {
                        const double px = chw.GetValue("mc_FSPartPx",entry,i);
                        const double py = chw.GetValue("mc_FSPartPy",entry,i);
                        const double pz = chw.GetValue("mc_FSPartPz",entry,i);
                        double p = sqrt( px*px + py*py + pz*pz );
                        if( p > max_p ) max_p = p; 
                    }
                } 

                max_p /= 1000;

                static const double Mp = 0.938272;
                const double ekin = sqrt( max_p*max_p + Mp*Mp ) - Mp;
                static const double Mn_prime = 0.93956-0.034;

                return Mn_prime*Mn_prime - Mp*Mp + 2*Mn_prime*(ekin + Mp - Mn_prime); 
            }
        case kENuQE:
            {
                // NOTE: this is for the neutrino version.
                // for the antineutrino version, the binding energy is different (30 MeV)
                // and the neutron and proton masses need to be swapped!

                // Neutron mass less the binding energy.
                static const double Mn=0.93956-0.034;

                // Proton mass
                static const double Mp=0.938272;

                const double coslep=getValue(chw, entry, kCosLep);
                const double Elep=getValue(chw, entry, kELep);

                const double numerator=2*Mn*Elep - ( sqr(Mn) + sqr(Mlep) - sqr(Mp) );
                const double denominator=2*(Mn - Elep + sqrt(sqr(Elep)-sqr(Mlep))*coslep);
                //        cout << numerator << " " << denominator << endl;

                return numerator/denominator;
            }
        case kPPi0:
        case kPPiPlus:
            {
                const double px=chw.GetValue("mc_FSPartPx", entry, piIndex);
                const double py=chw.GetValue("mc_FSPartPy", entry, piIndex);
                const double pz=chw.GetValue("mc_FSPartPz", entry, piIndex);

                double p = sqrt(px*px + py*py + pz*pz);
                return p*1e-3; // GeV
                //const double Epi=chw.GetValue("mc_FSPartE", entry, piIndex);
                //return sqrt(Epi*Epi-Mpi*Mpi);
            }
        case kCosPi0:
        case kCosPiPlus:
            {
                const double px=chw.GetValue("mc_FSPartPx", entry, piIndex);
                const double py=chw.GetValue("mc_FSPartPy", entry, piIndex);
                const double pz=chw.GetValue("mc_FSPartPz", entry, piIndex);
                return pz/sqrt(px*px+py*py+pz*pz);
            }
        case kThetaPi0:
        case kThetaChargedPion:
        case kThetaPiPlus:
            {
                const double px=chw.GetValue("mc_FSPartPx", entry, piIndex);
                const double py=chw.GetValue("mc_FSPartPy", entry, piIndex);
                const double pz=chw.GetValue("mc_FSPartPz", entry, piIndex);

                const double numi_beam_angle_rad = -0.05887; // Taken from MinervaPhysicalConstants.h
                double pyp = -1.0 *sin( numi_beam_angle_rad )*pz + cos( numi_beam_angle_rad )*py;
                double pzp = cos( numi_beam_angle_rad )*pz + sin( numi_beam_angle_rad )*py;
                double denom2 = pow(px,2) + pow(pyp,2) + pow(pzp,2);
                if( 0. == denom2 )
                    return -999;
                else 
                    return (180.0/3.141592654)*acos(pzp / sqrt(denom2) );
            }
        case kThetaLep:
            {

                const double px=chw.GetValue("mc_primFSLepton", entry, 0);
                const double py=chw.GetValue("mc_primFSLepton", entry, 1);
                const double pz=chw.GetValue("mc_primFSLepton", entry, 2);

                TVector3 pl(px, py, pz);

                // Copied (and modified) from MinervaCoordSysTool::thetaWRTBeam
                const double numi_beam_angle_rad = -0.05887; // Taken from MinervaPhysicalConstants.h
                double pyp = -1.0 *sin( numi_beam_angle_rad )*pz + cos( numi_beam_angle_rad )*py;
                double pzp = cos( numi_beam_angle_rad )*pz + sin( numi_beam_angle_rad )*py;
                assert(pzp != 0);
                return TMath::RadToDeg()*atan2(sqrt( sqr(px) + sqr(pyp)), pzp);

            }

        case kCosLep:
            return cos( TMath::DegToRad() * getValue(chw, entry, kThetaLep ) );
        case kENu:
            return 1e-3*chw.GetValue("mc_incomingE", entry);
        case kEHad:
            return 1e-3*( chw.GetValue("mc_incomingE", entry) - chw.GetValue("mc_primFSLepton", entry, 3) );
        case kEAvail:
          {
	    double recoil = 0;
	    int n_parts = chw.GetValue("mc_nFSPart", entry);
	    double mass_pion = 139.57;
	    double mass_proton = 938.27;
	    for(int i=0;i<n_parts;i++){
	      int pdg = chw.GetValue("mc_FSPartPDG",entry,i);
	      if(pdg == 22) {
          recoil+=chw.GetValue("mc_FSPartE",entry,i);//total energy
        } else if (pdg == 211 || pdg == -211) {
          recoil+=chw.GetValue("mc_FSPartE",entry,i)-mass_pion;//kinetic 
        } else if (pdg == 111) {
          recoil+=chw.GetValue("mc_FSPartE",entry,i);//total energy             
        } else if(pdg == 2212) {
          recoil+=chw.GetValue("mc_FSPartE",entry,i)-mass_proton;//kinetic  
        } else if (pdg>1000000000) {// neucleons and GENIE pseudo particles
          // do nothing, assume negligible energy for nucleons
          // save me the trouble of asking what that nucleon's mass was.
        } else if (pdg>=2000) {//primarily strange baryons 3122s
          recoil+=chw.GetValue("mc_FSPartE",entry,i)-mass_proton;
        } else if (pdg<-2000) {//primarily anti-baryons -2112 and -2212s,
          recoil+=chw.GetValue("mc_FSPartE",entry,i)+mass_proton;
        } else if (abs(pdg)==11 || abs(pdg)==13) {
          // do nothing. Don't include charged lepton
        } else {//primarily kaons and eta mesons.
          recoil+=chw.GetValue("mc_FSPartE",entry,i);
        }
	    }
	    return 1e-3*recoil;
          }
        case kTLep:
            {
                const double Elep=1e-3*chw.GetValue("mc_primFSLepton", entry, 3);
                return Elep-Mlep;
            }
        case kELep:
            return 1e-3*chw.GetValue("mc_primFSLepton", entry, 3);
        case kTPi0:
        case kTPiPlus:
            {
                const double Epi=1e-3*chw.GetValue("mc_FSPartE", entry, piIndex);
                return Epi-Mpi;
            }
        case kTChargedPion:
            {
                const double Epi=1e-3*chw.GetValue("mc_FSPartE", entry, piIndex);
                return 1.0e3*(Epi-Mpi);
            }
        case kPLep:
            {
                const double px = chw.GetValue("mc_primFSLepton", entry, 0);
                const double py = chw.GetValue("mc_primFSLepton", entry, 1);
                const double pz = chw.GetValue("mc_primFSLepton", entry, 2);
                double P = sqrt( px*px + py*py + pz*pz );

                return 1.0e-3*P; //GeV
            }
        case kPTLep:
            {
                const double px=chw.GetValue("mc_primFSLepton", entry, 0);
                const double py=chw.GetValue("mc_primFSLepton", entry, 1);
                const double pz=chw.GetValue("mc_primFSLepton", entry, 2);
                // Copied (and modified) from MinervaCoordSysTool::thetaWRTBeam
                const double numi_beam_angle_rad = -0.05887; // Taken from MinervaPhysicalConstants.h
                double pyp = -1.0 *sin( numi_beam_angle_rad )*pz + cos( numi_beam_angle_rad )*py;
                double pt   = sqrt( pow(px,2) +pow(pyp,2) );

                return 1.0e-3*pt; //GeV
            }
        case kPZLep:
            {
                const double pylep=chw.GetValue("mc_primFSLepton", entry, 1);
                const double pzlep=chw.GetValue("mc_primFSLepton", entry, 2);
                // Copied (and modified) from MinervaCoordSysTool::thetaWRTBeam
                const double numi_beam_angle_rad = -0.05887; // Taken from MinervaPhysicalConstants.h
                double pzp = cos( numi_beam_angle_rad )*pzlep + sin( numi_beam_angle_rad )*pylep; 

                return 1.0e-3*pzp; //GeV
            }
        case kPZRecoil:
	    {
                const double pxlep=chw.GetValue("mc_primFSLepton", entry, 0);
                const double pylep=chw.GetValue("mc_primFSLepton", entry, 1);
                const double pzlep=chw.GetValue("mc_primFSLepton", entry, 2);
                // Copied (and modified) from MinervaCoordSysTool::thetaWRTBeam
                const double numi_beam_angle_rad = -0.05887; // Taken from MinervaPhysicalConstants.h
                double pzp = cos( numi_beam_angle_rad )*pzlep + sin( numi_beam_angle_rad )*pylep; 
                double pyp = -1.0 *sin( numi_beam_angle_rad )*pzlep + cos( numi_beam_angle_rad )*pylep;
                double pt   = sqrt( pow(pxlep,2) +pow(pyp,2) );

		
		//recoil value
		double recoil = 0;
		int n_parts = chw.GetValue("mc_nFSPart", entry);
		double mass_pion = 139.57;
		double mass_proton = 938.27;
		for(int i=0;i<n_parts;i++){
		  int pdg = chw.GetValue("mc_FSPartPDG",entry,i);
		  if(pdg == 22) recoil+=chw.GetValue("mc_FSPartE",entry,i);//total energy
		  if(pdg == 211 || pdg == -211) recoil+=chw.GetValue("mc_FSPartE",entry,i)-mass_pion;//kinetic
		  if(pdg == 111) recoil+=chw.GetValue("mc_FSPartE",entry,i);//total energy
		  if(pdg == 2212) recoil+=chw.GetValue("mc_FSPartE",entry,i)-mass_proton;//kinetic
		}
		
		//call hyperdim
		vector<double> val_vect;
		val_vect.push_back(recoil);
		val_vect.push_back(pt/1000.);
		val_vect.push_back(pzp/1000.);
		int binnum = m_HyperDim->GetBin(val_vect).first;
		//return bin number
                return binnum+0.5; //Bin number(int+0.5)
	    }	        
        case kPTLepSquare:
            {
                const double px=chw.GetValue("mc_primFSLepton", entry, 0);
                const double py=chw.GetValue("mc_primFSLepton", entry, 1);
                const double pz=chw.GetValue("mc_primFSLepton", entry, 2);
                // Copied (and modified) from MinervaCoordSysTool::thetaWRTBeam
                const double numi_beam_angle_rad = -0.05887; // Taken from MinervaPhysicalConstants.h
                double pyp = -1.0 *sin( numi_beam_angle_rad )*pz + cos( numi_beam_angle_rad )*py;
                double pt   = sqrt( pow(px,2) +pow(pyp,2) );

                return 1.0e-6*pt*pt; //GeV^2
            }
        default:
            assert(0 && "Unknown variable");
    }


    //this is never returned.
    //if variable isn't known then assert fails
    return -666.;
}

// David Coplowe: Added 26/07/18
void XSec::SortParticles()
{
  std::sort(mu_mom.begin(), mu_mom.end(), Sort);
  std::sort(pr_mom.begin(), pr_mom.end(), Sort);
  std::sort(pi0_mom.begin(), pi0_mom.end(), Sort);
}

double XSec::GetdpTT(const TVector3 nu_dir, const TVector3 lep, const TVector3 had1, const TVector3 had2)
{
  TVector3 cross = nu_dir.Cross(lep);
  TVector3 had_sum = had1 + had2;
  return (cross.Unit().Dot(had_sum));
}
// David Coplowe: End 26/07/18

