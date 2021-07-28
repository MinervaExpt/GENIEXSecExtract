#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"

#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/NuclModUtils.h>
#include <PlotUtils/MnvNuclearModelWeight.h>

#include <iostream>
#include <vector>
#include <map>
#include <cassert>
using namespace std;
using namespace PlotUtils;

#include "NukeCC_bins.cc"


namespace
{
  int nLoop      = 0; //0 for all
/*  bool run_all   = true;
  bool run_cv    = true;
  bool run_by03  = true;
  bool run_by13  = true;
  bool run_noMod = true;
*/
  //neutrino energy from 2-20 or 5-20 GeV
  const double MIN_E = 2.;
  const double MIN_E_INEL = 2.;
  const double MAX_E = 20.;

  //muon theta from 0-17 deg
  const double MIN_COSMU = cos( TMath::DegToRad() * 17. );
  const double MAX_COSMU = 1.1;

  //pdg
  const int PROTON  = 2212;
  const int NEUTRON = 2112;


  namespace NuclModVariation
  {
    const std::string BY03 = "BY03";
    const std::string BY03_Bug = "BY03_bug";
    const std::string BY13 = "BY13";
    const std::string KPNoIso  = "KulaginNoIso";
    const std::string KPWithTwist  = "KulaginTwistOne";
    const std::string KPNoTwist  = "KulaginNukeCorr";
    const std::string NoMod = "NoModel";
    const std::string cv = BY03_Bug;
  }
}


typedef pair<double,double> KineCut;
typedef map<XSec::EVariable, KineCut> KineMap;
typedef pair<int,int> IPair;


// number of nucleons in standard fid volume / number in target
double GetNormValue( int targetZ, int targetNucleon = 0 )
{
  const double trackerAtomsC = 2.22311e+27 * 92.;

  double passiveNucleons = -999.;
  if( 0 == targetNucleon ) //total nucleons
  {
    if( 6 == targetZ )  //tgt 3
      passiveNucleons = 1.00e+29;
    else if( 26 == targetZ ) //tgt 1 + 2 + 3 + 5
      passiveNucleons =  /*1.95e+29 +*/ 1.95e+29 + 1.02e+29 + 9.86e+28;
    else if( 82 == targetZ ) //tgt 1 + 2 + 3 + 4 + 5
      passiveNucleons = /*1.59e+29 +*/ 1.59e+29 + 7.34e+28 + 1.37e+29 + 8.02e+28;
    else if( 0 == targetZ )
      passiveNucleons = 3.05063e+28 * 92.; //92 planes
    else if( 1 == targetZ )
      passiveNucleons = 2.51551e+27 * 92.; //92 planes
  }
  else if( PROTON == targetNucleon )
  {
    if( 6 == targetZ )  //tgt 3
      passiveNucleons = 5.00e+28;
    else if( 26 == targetZ ) //tgt 1 + 2 + 3 + 5
      passiveNucleons =  /*9.09e+28 +*/ 9.09e+28 + 4.74e+28 + 4.59e+28;
    else if( 82 == targetZ ) //tgt 1 + 2 + 3 + 4 + 5
      passiveNucleons = /*6.29e+28 +*/ 6.29e+28 + 2.91e+28 + 5.42e+28 + 3.18e+28;
    else if( 0 == targetZ )
      passiveNucleons = 1.64697e+28 * 92.; //92 planes
    else if( 1 == targetZ )
      passiveNucleons = 2.49570e+27 * 92.; //92 planes    
  }
  else if( NEUTRON == targetNucleon )
  {
    if( 6 == targetZ )  //tgt 3
      passiveNucleons = 5.00e+28;
    else if( 26 == targetZ ) //tgt 1 + 2 + 3 + 5
      passiveNucleons = /*1.04e+29 +*/ 1.04e+29 + 5.44e+28 + 5.27e+28;
    else if( 82 == targetZ ) //tgt 1 + 2 + 3 + 4 + 5
      passiveNucleons = /*9.60e+28 +*/ 9.60e+28 + 4.44e+28 + 8.27e+28 + 4.85e+28;
    else if( 0 == targetZ )
      passiveNucleons = 1.33624e+28 * 92.; //92 planes
    else if( 1 == targetZ )
      passiveNucleons = 1.98158e+25 * 92.; //92 planes 
  }
  else
  {
    assert( false && "Target nucleons can be all, proton or neutron only" );
  }

  if( passiveNucleons < 0 )
    assert( false && "Normalizations only known for Z = 1,6,26,82 and 0(scint)" );

  return trackerAtomsC / passiveNucleons;
}

//=================================
// NukeCC XSec
//  requires CC and kinematic cuts
//=================================
class NukeCC_XSec : public XSec
{
  public:
    explicit NukeCC_XSec(const char* name, int targetZ = 0, int targetNucleon = 0, bool inel = false) :
      XSec(name),
      targetZ_(targetZ),
      targetNucleon_(targetNucleon),
      inel_(inel)
  {};

    // Override this method from the base class to decide what events to
    // include in this selection
    virtual bool passesCuts(ChainWrapper& chw, int entry)
    {
      //must be CC
      if((int)chw.GetValue("mc_current", entry)!=1) 
        return false;

      //must be inelastic (skip QE and Res)
      if( inel_ && ( (int)chw.GetValue("mc_intType", entry)==1 ||(int)chw.GetValue("mc_intType", entry)==2 ) )
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


      const double tracker_min = 6.117;
      const double tracker_max = 8.193;

      //pass fiducial cuts
      double z=1e-3*fabs(chw.GetValue("mc_vtx", entry, 2));

      if( 0 != targetZ_ )
      {
        //require the right targetZ
        int target=(int)chw.GetValue("mc_targetZ", entry);
        if( targetZ_ != target)
          return false;

        //z must be in tgt2,3,4,5 ? 
        //bool tgt1 = ( 4.46 < z && z < 4.495 );
        bool tgt2 = ( 4.68 < z && z < 4.72 );
        bool tgt3 = ( 4.899 < z && z < 4.985 );
        bool tgt4 = ( 5.637 < z && z < 5.65 );
        bool tgt5 = ( 5.76 < z && z < 5.79 );

        if( 1 == targetZ_ )
        {
          if( z < tracker_min || tracker_max < z )
            return false;
        }
        else if( ! (tgt2 || tgt3 || tgt4 || tgt5 ) )
          return false;
      }
      else
      {
        if( z < tracker_min || tracker_max < z )
          return false;
      }

      //note hexagon cut already made by XSecLooper

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

    int targetZ_;
    int targetNucleon_;
    bool inel_;

};



//====================================
// XSecLooper for NukeCC
// can slect events based on targetZ
//====================================
class NukeCC_XSecLooper : public XSecLooper
{
  private:
    std::string _nuclMod;
    MnvNuclearModelWeight *mnvNuke;

  public:
    explicit NukeCC_XSecLooper(const char* inputFileGlob, const std::string nukeMod) :   
      XSecLooper(inputFileGlob),
      _nuclMod(nukeMod),
      mnvNuke(0)
  {mnvNuke = new MnvNuclearModelWeight(_nuclMod);}

    std::string getNuclMod( ) { return _nuclMod; }
    void setNuclMod( const std::string& mod ) { _nuclMod = mod; }

    bool isFiducial(ChainWrapper& chw, int entry)
    {
      //just make sure it's in the hexagon
      double x=1e-3*fabs(chw.GetValue("mc_vtx", entry, 0));
      double y=1e-3*fabs(chw.GetValue("mc_vtx", entry, 1));

      const double apothem=0.85;
      if(x*x + y*y < apothem*apothem) return true;

      double lenOfSide = apothem * ( 2 / sqrt(3) );

      if( x > apothem )
        return false;

      if( y < lenOfSide/2.0 )
        return true;

      double slope = (lenOfSide / 2.0) / apothem;
      if( y < lenOfSide - x*slope )
        return true;

      return false;
    }

    double getSignalWeight(ChainWrapper& chw, int entry)
    {
      
      //do not apply to charm events
      if( 1 == (int)chw.GetValue( "mc_charm", entry) )
        return 1.;

      int intType = (int)chw.GetValue( "mc_intType", entry );
      double w = (double)chw.GetValue( "mc_w", entry ) * 1e-3;

      //weight for DIS and possibly res xsec
      bool modifyRes  = true;
      double res_lowW = 1.3;  //<-- only modify res with W>this

      //modify dis events
      bool modify = intType == 3;
      //and resonance events if desired and W is high enough
      modify = modify || ( modifyRes && intType == 2 && res_lowW < w );

      //use weight 1 if no modification
      if(!modify)
        return 1.;

      const double x  = (double)chw.GetValue( "mc_Bjorkenx", entry );
      const double y  = (double)chw.GetValue( "mc_Bjorkeny", entry );
      const double q2 = (double)chw.GetValue( "mc_Q2", entry ) * 1e-6;
      const double Enu  = (double)chw.GetValue( "mc_incomingE", entry ) * 1e-3;
      const int    tgtNucleon = (int)chw.GetValue("mc_targetNucleon", entry);
      const int    A  = (int)chw.GetValue("mc_targetA", entry);      
      double rw =  mnvNuke->CalculateWeight(q2, x, y, Enu, tgtNucleon, A);
      return rw;

    }
};


//============================
// Crete a new XSec object
//============================
NukeCC_XSec* GetNewXSec( XSec::EVariable var, KineMap kineCuts, const string tag, IPair targetPair, bool inel = false, bool fine = false )
{
  const int targetZ = targetPair.first;
  const int targetNucl = targetPair.second;

  string varName = "";
  switch(var)
  {
    case XSec::kx:
      varName = "x";
      break;
    case XSec::kxExp:
      varName = "xExp";
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
    case XSec::kCosLep:
      varName = "CosMu";
      break;
    case XSec::kThetaLep:
      varName = "ThetaMu";
      break;
    default:
      varName = "UNKNOWN";
      break;
  }

  string tgtStr = "";
  switch(targetZ)
  {
    case 0:
      tgtStr = "tracker";
      break;
    case 1:
      tgtStr = "p";
      break;
    case 6:
      tgtStr = "carbon";
      break;
    case 26:
      tgtStr = "iron";
      break;
    case 82:
      tgtStr = "lead";
      break;
    default:
      tgtStr = "UNKNOWN";
      break;
  }

  string nuclStr = "";
  switch(targetNucl)
  {
    case 0:
      break;
    case PROTON:
      nuclStr = "_proton";
      break;
    case NEUTRON:
      nuclStr = "_neutron";
      break;
    default:
      nuclStr = "UNKNOWN";
      break;
  }

  if( fine )
    varName += "_fine";


  vector<double> bins;
  GetBins( var, bins, fine, inel );

  TString name = Form( "%s%s_%s_%s", tgtStr.c_str(), nuclStr.c_str(), varName.c_str(), tag.c_str() );
  NukeCC_XSec* xsec = new NukeCC_XSec(name.Data(), targetZ, targetNucl, inel);

  xsec->setBinEdges( bins.size()-1, &bins.front() );
  xsec->setVariable(var);
  xsec->setUniverses(0);
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

  xsec->setNormalizationValue( GetNormValue( xsec->targetZ_, xsec->targetNucleon_ ) );

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

    if(apply)
      xsec->AddKineCut( i->first, i->second );
  }

  return xsec;
}


//=======================================
// MAIN
//=======================================
void runXSecLooper_NuclMod()
{
  // Create the XSecLooper and tell it the input files

  NukeCC_XSecLooper loop_by13("/minerva/data/users/tice/NukeTuples/v19/NukeTruth_AnaTuple_minerva*root", NuclModVariation::BY13);
  NukeCC_XSecLooper loop_kp("/minerva/data/users/tice/NukeTuples/v19/NukeTruth_AnaTuple_minerva*root", NuclModVariation::KPNoTwist);


  // Tell the XSecLooper which neutrino type we're considering (mandatory)

  loop_by13.setNuPDG(14);
  loop_kp.setNuPDG(14);


  //Set the number of Universes to be used in GENIE error band (default 100)
  //Put setNumUniv(0) if you do not want to use any universe. 

  loop_by13.setNumUniv(0);
  loop_kp.setNumUniv(0);

  // define common kinematic regions
  KineMap noCuts;
  KineMap stdCuts;
  stdCuts[XSec::kCosLep] = KineCut( MIN_COSMU, MAX_COSMU );
  stdCuts[XSec::kENu]   = KineCut( MIN_E, MAX_E );
  KineMap stdCuts_inel;
  stdCuts_inel[XSec::kCosLep] = KineCut( MIN_COSMU, MAX_COSMU );
  stdCuts_inel[XSec::kENu]   = KineCut( MIN_E_INEL, MAX_E );

  //some reject QE events
  const bool inelOnly = true;

  //use finer bins
  const bool fine = true;

  vector< IPair > targets;
  targets.push_back(  IPair(0,0) );
  targets.push_back(  IPair(6,0) );
  targets.push_back( IPair(26,0) );
  targets.push_back( IPair(82,0) );


  vector<XSec::EVariable> vars;
  vars.push_back( XSec::kENu );
//  vars.push_back( XSec::kx );
  vars.push_back( XSec::kxExp );
//  vars.push_back( XSec::kThetaMu );
//  vars.push_back( XSec::kQ2 );

  // loop over kinematic limits and add XSec objects
  for( vector<IPair>::iterator targetZ = targets.begin(); targetZ != targets.end(); ++targetZ )
  {
    //with no kinematic cuts.
    {
      for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var )
      {

       loop_by13.addXSec( GetNewXSec( *var, noCuts, "full_by13", *targetZ ) );
       loop_kp.addXSec( GetNewXSec( *var, noCuts, "full_kp", *targetZ ) );
        loop_by13.addXSec( GetNewXSec( *var, noCuts, "full_by13", *targetZ, !inelOnly, fine ) );
        loop_kp.addXSec( GetNewXSec( *var, noCuts, "full_kp", *targetZ, !inelOnly, fine ) );

      }
    }

    //with std kinematic cuts...
    {
      for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ) 
      {

        loop_by13.addXSec( GetNewXSec( *var, stdCuts, "std_by13", *targetZ ) );
        loop_kp.addXSec( GetNewXSec( *var, stdCuts, "std_kp", *targetZ ) );
        loop_by13.addXSec( GetNewXSec( *var, stdCuts, "std_by13", *targetZ, !inelOnly, fine ) );
        loop_kp.addXSec( GetNewXSec( *var, stdCuts, "std_kp", *targetZ, !inelOnly, fine ) );

      }
    }

    //with no kinematic cuts and inelastic
    {
      for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var )
      {
        loop_by13.addXSec( GetNewXSec( *var, noCuts, "full_inel_by13", *targetZ, inelOnly ) );
        loop_kp.addXSec( GetNewXSec( *var, noCuts, "full_inel_kp", *targetZ, inelOnly ) );
        loop_by13.addXSec( GetNewXSec( *var, noCuts, "full_inel_by13", *targetZ, inelOnly, fine ) );
        loop_kp.addXSec( GetNewXSec( *var, noCuts, "full_inel_kp", *targetZ, inelOnly, fine ) );
      }
    }


    //with std kinematic cuts and inelastic
    {
      for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var )
      {
        loop_by13.addXSec( GetNewXSec( *var, stdCuts_inel, "std_inel_by13", *targetZ, inelOnly ) );
        loop_kp.addXSec( GetNewXSec( *var, stdCuts_inel, "std_inel_kp", *targetZ, inelOnly ) );
        loop_by13.addXSec( GetNewXSec( *var, stdCuts_inel, "std_inel_by13", *targetZ, inelOnly, fine ) );
        loop_kp.addXSec( GetNewXSec( *var, stdCuts_inel, "std_inel_kp", *targetZ, inelOnly, fine ) );
     //   loop_noMod.addXSec( GetNewXSec( *var, stdCuts_inel, "std_inel_noMod", *targetZ, inelOnly, fine ) );
      }
    }
  }

  // Once everything's set up, actually run the thing
  /*if( run_cv || run_all )
    loop.runLoop(      nLoop );
  if( run_by03 || run_all )
    loop_by03.runLoop( nLoop );
  if( run_by13 || run_all )*/
//    loop.runLoop( nLoop );  
    loop_kp.runLoop( nLoop );
    loop_by13.runLoop( nLoop );
  //if( run_noMod || run_all )
//    loop_noMod.runLoop( nLoop );

  // Get the output histograms and save them to file
  string fname = "NukeCC_NuclMod_KPNukeEffNothing.root";

  TFile fout(fname.c_str(), "RECREATE");

  for(uint i=0; i<loop_kp.getXSecs().size(); ++i){
  //  loop.getXSecs()[i]->getXSecHist()->Write();
     loop_kp.getXSecs()[i]->getXSecHist()->Write();
     loop_by13.getXSecs()[i]->getXSecHist()->Write();
  //  loop_noMod.getXSecs()[i]->getXSecHist()->Write();
    
    /*if( run_cv || run_all )
      loop.getXSecs()[i]->getXSecHist()->Write();
    if( run_by03 || run_all )
      loop_by03.getXSecs()[i]->getXSecHist()->Write();
    if( run_by13 || run_all )
      loop_by13.getXSecs()[i]->getXSecHist()->Write();
    if( run_noMod || run_all )
      loop_noMod.getXSecs()[i]->getXSecHist()->Write();*/
  }

  fout.Close();
}


int main()
{
  runXSecLooper_NuclMod();
  return 0;
}
