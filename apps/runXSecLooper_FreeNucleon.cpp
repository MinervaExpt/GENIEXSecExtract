#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include <PlotUtils/MnvH1D.h>
using namespace std;

#include "NukeCC_bins.cc"

namespace
{
  //neutrino energy from 2-20 or 5-20 GeV
  const double MIN_E = 5.;
  const double MIN_E_INEL = 5.0;
  const double MAX_E = 50.;
	const double MIN_Q2 = 1.0;
	const double MAX_Q2 = 1000.0;
	const double MIN_W = 2.0;
	const double MAX_W = 1000.0;
  //muon theta from 0-17 deg
  const double MIN_COSMU = cos( TMath::DegToRad() * 17. );
  const double MAX_COSMU = 1.1;

  //pdg
  const int PROTON  = 2212;
  const int NEUTRON = 2112;
  const int NUMU    = 14;
  const int NUMUBAR = -14;
}


typedef pair<double,double> KineCut;
typedef map<XSec::EVariable, KineCut> KineMap;
typedef pair<int,int> IPair;


//=================================
// FreeNucleon XSec
//  requires CC and kinematic cuts
//=================================
class FreeNucleon_XSec : public XSec
{
  public:
    explicit FreeNucleon_XSec(const char* name, int targetNucleon = 0, bool inel = false, int channel = -1) :
      XSec(name),
      targetNucleon_(targetNucleon),
      inel_(inel),
      channel_(channel)
  {};

    // Override this method from the base class to decide what events to
    // include in this selection
    virtual bool passesCuts(ChainWrapper& chw, int entry)
    {
      //must be CC
      if((int)chw.GetValue("mc_current", entry)!=1) 
        return false;
	
      if(channel_ != -1 && (int)chw.GetValue("mc_intType", entry) != channel_ )
        return false;
      
      //No charm in CCQE
      if(channel_ == 1 && (int)chw.GetValue("mc_charm", entry) == 1 )
        return false;
      
      //must be inelastic (skip QE )
      if( inel_ && ( (int)chw.GetValue("mc_intType", entry)==1 || (int)chw.GetValue("mc_intType", entry)==2 ) )
        return false;

      //require a target nucleon?
      if( 0 != targetNucleon_ && ( (int)chw.GetValue("mc_targetNucleon", entry) != targetNucleon_ ) )
        return false;

      //must pass kinematic cuts
      for( KineMap::const_iterator iKin = kineMap_.begin(); iKin != kineMap_.end(); ++iKin )
      {
        EVariable var = iKin->first;
        double lo = iKin->second.first;
        double hi = iKin->second.second;
        double val = getValue(chw, entry, var );
        if( val < lo || hi < val ) 
          return false;
      }

      return true;
    }

    void AddKineCut( EVariable var, double lo, double hi )
    {
      kineMap_[var] = KineCut(lo,hi);
    }

    void AddKineCut( EVariable var, KineCut x )
    {
      kineMap_[var] = x;
    }

    //use these to impose a kinematic region
    KineMap kineMap_;

    int targetNucleon_;
    bool inel_;
    int channel_;

};



//====================================
// XSecLooper for FreeNucleon
//====================================
class FreeNucleon_XSecLooper : public XSecLooper
{
  public:
    FreeNucleon_XSecLooper(const char* inputFileGlob, int nucleon, int neutrino) :   
      XSecLooper(inputFileGlob),
      m_nucleon(nucleon)
  {
    setNuPDG( neutrino );
    TString splineName = "nu_";

    if (m_nuPDGs.size() != 1)
    {
      cerr << "FreeNucleon_XSecLooper: only single-neutrino cross-sections are currently implemented." << endl;
      exit(1);
    }

    int pdg = *(m_nuPDGs.begin());
    switch(abs(pdg))
    {
      case 14:
        splineName+="mu";
        break;
      case 12:
        splineName+="e";
        break;
      default:
        cerr << "FreeNucleon_XSecLooper: nu PDG is set to " << pdg
          << ", which we don't have splines for. Bailing" << endl;
        exit(1);
    }
    if(pdg<0) splineName+="_bar";

    if( m_nucleon == NEUTRON )
      splineName+="_n/tot_cc";
    else if( m_nucleon == PROTON )
      splineName+="_H1/tot_cc";
    else
    {
      cerr << "FreeNucleon_XSecLooper: nucleon PDG set to " << m_nucleon
        << ", which we don't have splines for.  Bailing" << endl;
      exit(1);
    }

    setRateSplineName(splineName);
  }

    //override event selections
    bool isCCRateEvent(ChainWrapper& chw, int entry)
    {
      //requires a cc event...
      int current=(int)chw.GetValue("mc_current", entry);
      if( 1 != current ) return false;

      //... on the right nucleon
      int target=(int)chw.GetValue("mc_targetNucleon", entry);
      if( m_nucleon != target ) return false;

      return true;
    }

    // free nucleon samples generated on target mix not a detector, so all is fiducial.
    bool isFiducial(ChainWrapper& chw, int entry)
    {
      return true;
    }

  private:
    int m_nucleon;

};

//============================
// Crete a new XSec object
//============================
FreeNucleon_XSec* GetNewXSec( XSec::EVariable var, KineMap kineCuts, const string tag, int targetNucl, bool inel = false, int channel = -1,  bool fine = false )
{

  string varName = "";
  switch(var)
  {
    case XSec::kx:
      varName = "xGen";
      break;
    case XSec::kxExp:
      varName = "x";
      break;  
    case XSec::ky:
      varName = "y";
      break;
    case XSec::kW:
      varName = "W";
      break;
    case XSec::kQ2:
      varName = "Q2";
      break;
    case XSec::kTLep:
      varName = "Tmu";
      break;
    case XSec::kELep:
      varName = "Emu";
      break;
    case XSec::kENu:
      varName = "Enu";
      break;
    case XSec::kEHad:
      varName = "Ehad";
      break;
    case XSec::kThetaLep:
      varName = "ThetaMu";
      break;
    case XSec::kCosLep:
      varName = "CosMu";
      break;
    default:
      varName = "UNKNOWN";
      break;
  }

  string nuclStr = "";
  switch(targetNucl)
  {
    case 0:
      break;
    case PROTON:
      nuclStr = "freep";
      break;
    case NEUTRON:
      nuclStr = "freen";
      break;
    default:
      nuclStr = "UNKNOWN";
      break;
  }

  if( fine )
    varName += "_fine";

  vector<double> bins;
  GetBins( var, bins, fine, inel );

  TString name = Form( "%s_%s_%s", nuclStr.c_str(), varName.c_str(), tag.c_str() );

  cout << "  Creating new XSec with name: " << name << endl;

  FreeNucleon_XSec* xsec = new FreeNucleon_XSec(name.Data(), targetNucl, inel, channel);

  xsec->setBinEdges( bins.size()-1, &bins.front() );
  xsec->setVariable(var);
  xsec->setUniverses(0);//default value 100 universes, put 0 if you do not want universes to be included
  if( XSec::kENu == var )
    xsec->setIsFluxIntegrated(false);
  else
    xsec->setIsFluxIntegrated(true);

  //don't limit the E range unless the E cut is applied
  KineMap::const_iterator itEnuCut = kineCuts.find(XSec::kENu);
  if( XSec::kENu != var && itEnuCut != kineCuts.end() )
  {
    KineCut enuCut = itEnuCut->second;
    xsec->setFluxIntLimits(enuCut.first, enuCut.second);
  }

  xsec->setNormalizationType( XSec::kNoNorm );

  //apply all kine cuts that apply to variables except the plotted variable
  for( KineMap::const_iterator i = kineCuts.begin(); i != kineCuts.end(); ++i )
  {
    bool apply = true;

    //do not apply kine cut to plotted variable
    if( var == i->first )
      apply = false;
    //special check for thetamu/cosThetaMu.  we apply the cut in cos theta mu (for some reason)
    if( var == XSec::kThetaLep && i->first == XSec::kCosLep )
      apply = false;

    //if(apply)
      xsec->AddKineCut( i->first, i->second );
  }


  return xsec;
}


//=======================================
// MAIN
//=======================================
void runXSecLooper_FreeNucleon()
{
  // Create the XSecLooper and tell it the input files
  FreeNucleon_XSecLooper loop_n("/minerva/app/users/drimal/cmtuser/Minerva_v10r8/IsoCor/anaFiles/freen.root", NEUTRON, NUMU);
  FreeNucleon_XSecLooper loop_p("/minerva/app/users/drimal/cmtuser/Minerva_v10r8/IsoCor/anaFiles/freep.root", PROTON,  NUMU);

  //Setting the number of Universes in the GENIE error band
  loop_n.setNumUniv(0); 
  loop_p.setNumUniv(0); 

  // define common kinematic regions
  KineMap noCuts;
  KineMap stdCuts;
  stdCuts[XSec::kCosLep] = KineCut( MIN_COSMU, MAX_COSMU );
  stdCuts[XSec::kENu]   = KineCut( MIN_E, MAX_E );
  KineMap stdCuts_inel;
	
  stdCuts_inel[XSec::kCosLep] = KineCut( MIN_COSMU, MAX_COSMU );
  stdCuts_inel[XSec::kENu]   = KineCut( MIN_E_INEL, MAX_E );

	KineMap stdCuts_dis;
 	stdCuts_dis[XSec::kExpQ2]  = KineCut( MIN_Q2, MAX_Q2 );
 	stdCuts_dis[XSec::kExpW]   = KineCut( MIN_W, MAX_W );
	stdCuts_dis[XSec::kCosLep] = KineCut( MIN_COSMU, MAX_COSMU );
 	stdCuts_dis[XSec::kENu]    = KineCut( MIN_E, MAX_E );
	

  //reject QE/RES events
  const bool inelOnly = true;

  //use finer bins
  const bool fine = true;

  vector<XSec::EVariable> vars;
  vars.push_back( XSec::kENu );
//  vars.push_back( XSec::kEHad );
//  vars.push_back( XSec::kTMu );
//  vars.push_back( XSec::kCosMu );
  vars.push_back( XSec::kx );
  vars.push_back( XSec::kxExp );
//  vars.push_back( XSec::ky );
//  vars.push_back( XSec::kW );
//  vars.push_back( XSec::kQ2 );


  //with no kinematic cuts.
  {
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var )
    {
      loop_n.addXSec( GetNewXSec( *var, noCuts, "full", NEUTRON) );
      loop_p.addXSec( GetNewXSec( *var, noCuts, "full", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, noCuts, "full", NEUTRON, !inelOnly, fine) );
      loop_p.addXSec( GetNewXSec( *var, noCuts, "full", PROTON, !inelOnly, fine ) );
    }
  }

  //with std kinematic cuts...
  {
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, stdCuts, "std", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, stdCuts, "std", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, stdCuts, "std", NEUTRON, !inelOnly, fine ) );
      loop_p.addXSec( GetNewXSec( *var, stdCuts, "std", PROTON, !inelOnly, fine ) );
    }
  }

  //with no kinematic cuts and inelastic
  {
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var )
    {
      loop_n.addXSec( GetNewXSec( *var, noCuts, "full_inel", NEUTRON, inelOnly ) );
      loop_p.addXSec( GetNewXSec( *var, noCuts, "full_inel", PROTON, inelOnly ) );
      loop_n.addXSec( GetNewXSec( *var, noCuts, "full_inel", NEUTRON, inelOnly, fine ) );
      loop_p.addXSec( GetNewXSec( *var, noCuts, "full_inel", PROTON, inelOnly, fine ) );
    }
  }


  //with std kinematic cuts and inelastic
  {
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var )
    {
      loop_n.addXSec( GetNewXSec( *var, stdCuts_inel, "std_inel", NEUTRON, inelOnly ) );
      loop_p.addXSec( GetNewXSec( *var, stdCuts_inel, "std_inel", PROTON, inelOnly ) );
      loop_n.addXSec( GetNewXSec( *var, stdCuts_inel, "std_inel", NEUTRON, inelOnly, fine ) );
      loop_p.addXSec( GetNewXSec( *var, stdCuts_inel, "std_inel", PROTON, inelOnly, fine ) );
    }
  }


  //with std kinematic cuts QE only...
  {
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, stdCuts, "stdQE", NEUTRON, !inelOnly, 1) );
      loop_p.addXSec( GetNewXSec( *var, stdCuts, "stdQE", PROTON, !inelOnly, 1 ) );
      loop_n.addXSec( GetNewXSec( *var, stdCuts, "stdQE", NEUTRON, !inelOnly, 1, fine ) );
      loop_p.addXSec( GetNewXSec( *var, stdCuts, "stdQE", PROTON, !inelOnly, 1, fine ) );
    }
  }

  //with std kinematic cuts resonance only...
  {
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, stdCuts, "stdRES", NEUTRON, !inelOnly, 2) );
      loop_p.addXSec( GetNewXSec( *var, stdCuts, "stdRES", PROTON, !inelOnly, 2 ) );
      loop_n.addXSec( GetNewXSec( *var, stdCuts, "stdRES", NEUTRON, !inelOnly, 2, fine ) );
      loop_p.addXSec( GetNewXSec( *var, stdCuts, "stdRES", PROTON, !inelOnly, 2, fine ) );
    }
  }
 
  //with std kinematic cuts delta only...
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(0, 1.35);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "stdDRES", NEUTRON, inelOnly, -1) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "stdDRES", PROTON, inelOnly, -1 ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "stdDRES", NEUTRON, inelOnly, -1, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "stdDRES", PROTON, inelOnly, -1, fine ) );
    }
  }  

  //with std kinematic cuts dis only...
  {
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, stdCuts_dis, "stdDIS", NEUTRON, !inelOnly, 3) );
      loop_p.addXSec( GetNewXSec( *var, stdCuts_dis, "stdDIS", PROTON, !inelOnly, 3 ) );
      loop_n.addXSec( GetNewXSec( *var, stdCuts_dis, "stdDIS", NEUTRON, !inelOnly, 3, fine ) );
      loop_p.addXSec( GetNewXSec( *var, stdCuts_dis, "stdDIS", PROTON, !inelOnly, 3, fine ) );
    }
  }  

/*
  //==================
  // DIS type cuts
  //==================
  //... and w>1.4
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(1.4, 1.E9);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "w_14", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "w_14", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "w_14", NEUTRON,!inelOnly, -1, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "w_14", PROTON, !inelOnly, -1, fine ) );
    }
  }
  
   //... and w>1.7
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(1.7, 1.E9);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "w_17", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "w_17", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "w_17", NEUTRON,!inelOnly, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "w_17", PROTON, !inelOnly, fine ) );
    }
  }
  
  
  //... and 1.7 < W <2
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(1.7, 2.0);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "w_172", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "w_172", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "w_172", NEUTRON,!inelOnly, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "w_172", PROTON, !inelOnly, fine ) );
    }
  }


  //... and w>2.
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(2., 1.E9);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var )
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "w_2", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "w_2", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "w_2", NEUTRON,!inelOnly, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "w_2", PROTON,!inelOnly, fine ) );
    }
  }


  //... and q2>1 all W
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kQ2] = KineCut(1., 1.E9);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2_1", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2_1", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2_1", NEUTRON,!inelOnly, -1, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2_1", PROTON, !inelOnly, -1, fine ) );
    }
  }
  
  //... and w>1.4 and q2 > 1
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(1.4, 1.E9);
    wCuts[XSec::kQ2] = KineCut(1., 1.E9);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_14", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_14", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_14", NEUTRON,!inelOnly, -1, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_14", PROTON, !inelOnly, -1, fine ) );
    }
  }



 
  //... and w>1.7 and q2 > 1
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(1.7, 1.E9);
    wCuts[XSec::kQ2] = KineCut(1., 1.E9);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var )
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_17", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_17", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_17", NEUTRON,!inelOnly, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_17", PROTON,!inelOnly, fine ) );
    }
  }

    //... and 1.7 < W < 2.0 and q2 > 1
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(1.7, 2);
    wCuts[XSec::kQ2] = KineCut(1., 1.E9);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var )
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_172", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_172", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_172", NEUTRON,!inelOnly, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_172", PROTON,!inelOnly, fine ) );
    }
  }

  //... and w>2. and q2 > 1
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(2., 1.E9);
    wCuts[XSec::kQ2] = KineCut(1., 1.E9);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var )
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_2", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_2", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_2", NEUTRON,!inelOnly, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2_1_w_2", PROTON,!inelOnly, fine ) );
    }
  }
  
    //... and q2<1 all W
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kQ2] = KineCut(0., 1.);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2lt1", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2lt1", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2lt1", NEUTRON,!inelOnly, -1, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2lt1", PROTON, !inelOnly, -1, fine ) );
    }
  }
  
  //... and w>1.4 and q2 < 1
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(1.4, 1.E9);
    wCuts[XSec::kQ2] = KineCut(0., 1.0);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_14", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_14", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_14", NEUTRON,!inelOnly, -1, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_14", PROTON, !inelOnly, -1, fine ) );
    }
  }

  //... and w>1.7 and q2 < 1
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(1.7, 1.E9);
    wCuts[XSec::kQ2] = KineCut(0., 1.0);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_17", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_17", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_17", NEUTRON,!inelOnly, -1, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_17", PROTON, !inelOnly, -1, fine ) );
    }
  }  

  //... and w>2 and q2 < 1
  {
    KineMap wCuts = stdCuts;
    wCuts[XSec::kW] = KineCut(2.0, 1.E9);
    wCuts[XSec::kQ2] = KineCut(0., 1.0);
    for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
    {
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_2", NEUTRON ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_2", PROTON ) );
      loop_n.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_2", NEUTRON,!inelOnly, -1, fine ) );
      loop_p.addXSec( GetNewXSec( *var, wCuts, "q2lt1_w_2", PROTON, !inelOnly, -1, fine ) );
    }
  } */

  //================================================
  // Once everything's set up, actually run the thing
  loop_n.runLoop(  );
  loop_p.runLoop(  );

  // Get the output histograms and save them to file
  string fname = "FreeNucleon_xsec.root";
  TFile fout(fname.c_str(), "RECREATE");

  for(uint i=0; i<loop_n.getXSecs().size(); ++i)
    loop_n.getXSecs()[i]->getXSecHist()->Write();

  for(uint i=0; i<loop_p.getXSecs().size(); ++i)
    loop_p.getXSecs()[i]->getXSecHist()->Write();

  //write the flux too
  loop_p.getFluxHist()->Write();

}


int main()
{
  runXSecLooper_FreeNucleon();
  return 0;
}
