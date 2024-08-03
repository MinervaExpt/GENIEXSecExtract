//  by cnguyen@FNAL.GOV
// Helium GENIEXSecExtract Tool
//Oct 24 2022

#include "GENIEXSecExtract/XSecLooper.h"
#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif


#include "TFile.h"
#include "TH1.h"
#include <PlotUtils/MnvH1D.h>
//#include "TDatabasePDG.h"
#include "TMath.h"
#include <assert.h>
#include <PlotUtils/FluxReweighter.h>
#include <PlotUtils/TargetUtils.h>
//TDatabasePDG *pdg_DATABASEobject = TDatabasePDG::Instance();// Could not get this to work , Calculated particle mass by kinmatic equation

// varibles cases kPZLep, kPTLep, kELep, and kThetaLep only looking at kPTLep

///////////////////////////////////////////////////////////////////////////////
///// Helium Cut parameter
///////////////////////////////////////////////////////////////////////////////

namespace HeliumCUTConsts {
  const int NeutrinoPDGType = 14;
  const int HeliumZ= 2;
  const double Maximum_MuonAngle = 12.0;
  const double Maximum_SigCurvature = -5.0;
  const int MC_currentType = 1; //1 CC 2 NC
  const int Neutrino_Muon_Chargetype = -1;
  const int GreaterthanNTracks= 1;
  const double Muon_Minimum_Energy = 2000; // MeV
  const double Muon_Maximum_Energy = 50000; // MeV
  const double VextexMiniumDisToCryoTankInnerTank = 50.0; // mm For Cut Value


  const double Truth_pion_Minimum_Energy = .060; // GeV
  const double Truth_proton_Minimum_Energy = .105; // GeV
  const double Truth_dimuon_Maximum_Energy = .060; // GeV
  const double Truth_HardonAngle_wrtb_Maximum = 55; //deg
  const double Truth_EFractionAreTrueDigits = 0.95; //

  const double Maximum_secTrkwrtb_CUT = 55;
  const double Maximum_secTrkwrtb_GreatestKE = 55; // For Truth Section of index of second track with Greatest KE

  const double Maximum_Longtrack_VertexChiSqrt = 40.0; //25 // 40
  const double Maximum_OtherLongtrack_VertexChiSqrt = 40; //20 // 40
  const double Maximum_Short_VertexChiSqrt = 25.0; //10 // 25


}


///////////////////////////////////////////////////////////////////////////////
///// CryoTank parameters
///////////////////////////////////////////////////////////////////////////////


namespace CryoTankConsts {

////////////////////////////
// All units are in mm /////
////////////////////////////

const double CryoDishedHeadInnerRadius = 1476.0;
const double CryoDishedHeadInnerRadiusMin = 1271.24977876;
const double CryoDishedHeadInnerTheta = 0.533012304503;
const double CryoDishedHeadOuterRadius = 1483.62;
const double CryoDishedHeadTheta = 0.534952705856;
const double CryoHeatShieldAngle = 1.57079632679;
const double CryoHeatShieldCapThickness = 1.5875;
const double CryoHeatShieldInnerRadius = 838.0;
const double CryoHeatShieldLength = 1711.16;
const double CryoHeatShieldOuterRadius = 839.5875;
const double CryoHeatShieldThickness = 1.5875;
const double CryoInnerVesselDishedHeadThickness = 7.62;
const double CryoInnerVesselInnerDiameter = 1500.0;
const double CryoInnerVesselInnerRadius = 750.0;
const double CryoInnerVesselLength = 1187.46;
const double CryoInnerVesselOuterRadius = 756.35;
const double CryoInnerVesselThickness = 6.35;
const double CryoInnerVesselTotalLength = 3727.04417087;
const double CryoInnerVesselZPosInVacuumVessel = -206.0;
const double CryoVacuumVesselCylinderLength = 2262.7;
const double CryoVacuumVesselDishedHeadInnerRadius = 1933.0;
const double CryoVacuumVesselDishedHeadOuterRadius = 1940.62;
const double CryoVacuumVesselDishedHeadTheta = 0.490639019349;
const double CryoVacuumVesselDishedHeadThickness = 7.62;
const double CryoVacuumVesselEndHemiInnerRadius = 1066.8;
const double CryoVacuumVesselEndHemiOuterRadius = 1071.1942;
const double CryoVacuumVesselEndHemiTheta = 1.02969680084;
const double CryoVacuumVesselEndHemiThickness = 4.3942;
const double CryoVacuumVesselInnerRadius = 906.526;
const double CryoVacuumVesselOuterRadius = 914.4;
const double CryoVacuumVesselThickness = 7.874;
const double DefaultOffset = 3093.91;
const double FBackOfCryoInnerVessel = 3687.44726506;
const double FBackOfCryoVessel = 4225.26;
const double FrontOfCryoInnerVessel = 2088.37273494;
const double FrontOfCryoVessel = 1733.6288223;
const double ZZCryoInnerVesselLength = 1596.96044248;
const double ZZCryoVacuumVesselBackToCryoInnerVessel = 537.812734936;
const double ZZCryoVacuumVesselFrontToCryoInnerVessel = 354.743912632;
const double ZZCryoVacuumVesselLength = 2491.6311777;

/////////////////////////////////////////////////////////////
////Length of Upstream and DownStream Bulges of CyroTank/////
const double Lenghtofbulge =   abs(CryoDishedHeadInnerRadiusMin - CryoDishedHeadInnerRadius);
/////////////////////////////////////////////////////////////
// Z position of the DOWNstream
////////////////////////////////////////////////////////////
const double endbulgepoint = FBackOfCryoInnerVessel - CryoDishedHeadInnerRadius;
/////////////////////////////////////////////////////////////
// Z position of the END of upstream Bulge
////////////////////////////////////////////////////////////
const double Zpostion_EndOfupstreamBulge = FrontOfCryoInnerVessel + Lenghtofbulge;
/////////////////////////////////////////////////////////////
// Z position of the START of Downstream Bulge
////////////////////////////////////////////////////////////
const double Zpostion_StartOfdownstreamendBulge = FBackOfCryoInnerVessel - Lenghtofbulge;
const double inverse_CryoDishedHeadInnerRadius   = 1.0 / CryoDishedHeadInnerRadius;
const double FrontCryotank = CryoDishedHeadInnerRadius + FrontOfCryoInnerVessel;

const double startofradiusUpstreamVector = FrontOfCryoInnerVessel + CryoDishedHeadInnerRadius;
const double startofradiusDownstreamVector = FBackOfCryoInnerVessel - CryoDishedHeadInnerRadius;

const double Lenght_of_CryotankBarrel = Zpostion_StartOfdownstreamendBulge - Zpostion_EndOfupstreamBulge;
const double Zpostion_MidCryoTank = FrontOfCryoInnerVessel + Lenghtofbulge + (Lenght_of_CryotankBarrel / 2.0);

}

///////////////////////////////////////////////////////////////////////////////
///// PDG
///////////////////////////////////////////////////////////////////////////////

namespace Helium_PDG {

const int pdg_Pi0 = 111;
const int pdg_neutron = 2112;
const int pdg_antineutron = -2112;
const int pdg_Genie_bindingE = 2000000101;
const int pdg_Sigma0 =3212;
const int pdg_antiSigma0 = -3212;
const int pdg_Lambda0 =3122;
const int pdg_antiLambda0 = -3122;
const int pdg_Nu_e = 12;
const int pdg_Nu_mu = 14;
const int pdg_Proton = 2212;
const int pdg_Pion_neg = -211;
const int pdg_Pion_pos = 211;
const int pdg_Photon = 22;

}

///////////////////////////////////////////////////////////////////////////////
///// Helium GENIEXSecExtract Class
///////////////////////////////////////////////////////////////////////////////
typedef unsigned int uint;

class CCInclusive_Helium1DXSec : public XSec
{
public:
   CCInclusive_Helium1DXSec(const char* name) : XSec(name)
    {


    }

  // Override this method from the base class to decide what events to
  // include in this selection
  ///////////////////////////////////////////////////////////////////////////////
  ///// Cut Function
  ///////////////////////////////////////////////////////////////////////////////
  virtual bool passesCuts(ChainWrapper& chw, int entry)
  {

    //std::cout<<"inside Helium PassesCuts "<< std::endl;
    if(!IsNeutrino( chw,  entry)){return false;}
    if(!IsCCInteraction( chw,  entry)){return false;}
    if(!IsMuonEnergyGood( chw,  entry)){return false;}
    if(!IsMuonAngleGood( chw,  entry)){return false;}
    if(!TargetisHelium( chw,  entry)){return false;}
    if(!TRUTH_GreaterthanONEFS(chw, entry)){return false;}
    if(!IsInFiducalVolumeFromtheInnerEdgeTRUTH_new( chw,  entry, HeliumCUTConsts::VextexMiniumDisToCryoTankInnerTank)){return false;}
    //if(!TRUTH_Is2ndTrk_maxiumAngle_threshold_No_Neutral(chw,  entry, HeliumCUTConsts::Maximum_secTrkwrtb_CUT)){return false;}
    if(!TRUTH_Is2ndTrk_maxiumAngle_threshold_No_Neutral_WITHProtonAndPion_thresholds( chw,  entry ,HeliumCUTConsts::Maximum_secTrkwrtb_CUT,HeliumCUTConsts::Truth_pion_Minimum_Energy,HeliumCUTConsts::Truth_proton_Minimum_Energy)){return false; }// includes 2ndtrkenergy threshold cut

    return true;
  }
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////
  bool IsNeutrino(ChainWrapper& chw, int entry){
    if ((int)chw.GetValue("mc_incoming",entry) != HeliumCUTConsts::NeutrinoPDGType)
    {return false;}
    else{return true;}
  }
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////
  bool  IsCCInteraction(ChainWrapper& chw, int entry){
    if ((int)chw.GetValue("mc_current",entry) != HeliumCUTConsts::MC_currentType)
    {return false;}
    else{return true;}
  }
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////

  bool IsMuonEnergyGood(ChainWrapper& chw, int entry){
    double TrueMuonEnergy = (double)chw.GetValue("mc_primFSLepton",entry,3);
    if (TrueMuonEnergy <= HeliumCUTConsts:: Muon_Minimum_Energy || TrueMuonEnergy >= HeliumCUTConsts::Muon_Maximum_Energy )
    {return false;}
    else{return true;}
  }
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////
  bool  IsMuonAngleGood(ChainWrapper& chw, int entry){
    double TrueMuonAngle = (double)chw.GetValue("truth_true_muon_theta",entry)* TMath::RadToDeg();
    if (TrueMuonAngle > HeliumCUTConsts::Maximum_MuonAngle ||  TrueMuonAngle < 0.0 )
    {return false;}
    else{return true;}
  }
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////
  bool TargetisHelium(ChainWrapper& chw, int entry){
    if((int)chw.GetValue("mc_targetZ",entry) != HeliumCUTConsts::HeliumZ)
    {return false;}
    else{return true;}
  }
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////
  bool TRUTH_GreaterthanONEFS(ChainWrapper& chw, int entry){
    int nparticles = (int)chw.GetValue("mc_nFSPart", entry);
    if (nparticles == 1 || nparticles == 0)
    {return false;}
    else{return true;}
  }
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////
  double Distance_to_innerTank(double Zinput, double Rinput ){
    double delta_to_shell_Sphereical_Cap_upstream = 99999;
    double delta_to_shell_Sphereical_Cap_downstream = 99999;
    double delta_to_shell_Barrel = 99999;
    double Rsmall_upstream = 99999;
    double Rsmall_downstream = 99999;
    //double Rsmall_Center = 99999;

    //////////////
    //upstream Cap
    //////////////
    if((Zinput < CryoTankConsts::Zpostion_EndOfupstreamBulge)){
      double lenghtupstream_sphereical  = abs(CryoTankConsts::startofradiusUpstreamVector  -  Zinput);
      Rsmall_upstream = sqrt( pow(lenghtupstream_sphereical,2) + pow(Rinput,2) );
      delta_to_shell_Sphereical_Cap_upstream =  CryoTankConsts::CryoDishedHeadInnerRadius  - Rsmall_upstream;
      delta_to_shell_Barrel  = CryoTankConsts::CryoInnerVesselInnerRadius - Rinput;

      if(delta_to_shell_Sphereical_Cap_upstream < delta_to_shell_Barrel ){return delta_to_shell_Sphereical_Cap_upstream;}
      else{return delta_to_shell_Barrel; }

    }

    //////////////
    //  Barrel
    //////////////
    else if ((CryoTankConsts::Zpostion_EndOfupstreamBulge <= Zinput) && (Zinput <= CryoTankConsts::Zpostion_StartOfdownstreamendBulge)){
      //delta_to_shell  = (AddRegionOutside + CryoTankConsts::CryoInnerVesselInnerRadius) - Rinput;
      ///////
      // if the vertex is outside of the barrel then the closest postion would be to the Cylindrical outter wall and the caps are not considered
      /////

      ///////////
      /// IF the vertex point is "inside" the Cyrotank's Inner walls then the distances
      //to the  Barrel's wall and upstream and downstream Caps' are checked and the smallest distance is considered
      //////////////////////////////////////
      //// Find Distance to barrel wall
      /////////////////////////////////////
      delta_to_shell_Barrel  = CryoTankConsts::CryoInnerVesselInnerRadius - Rinput;
      //////////////////////////////////////
      //// Find Distance to upstream Cap
      /////////////////////////////////////
      if(CryoTankConsts::Zpostion_MidCryoTank > Zinput){
        // MoreUpstream towards the beam
        double lenghtupstream_sphereical  = abs(CryoTankConsts::startofradiusUpstreamVector  -  Zinput);
        Rsmall_upstream = sqrt( pow(lenghtupstream_sphereical,2) + pow(Rinput,2) );
        delta_to_shell_Sphereical_Cap_upstream =  CryoTankConsts::CryoDishedHeadInnerRadius  - Rsmall_upstream;
      }
      //////////////////////////////////////
      //// Find Distance to Downstream Cap
      /////////////////////////////////////
      else if(CryoTankConsts::Zpostion_MidCryoTank < Zinput){
        // more downstream towawds Minerva
        double lenghtdownstream_sphereical  = abs(Zinput - CryoTankConsts::startofradiusDownstreamVector);
        Rsmall_downstream = sqrt( pow(lenghtdownstream_sphereical,2) + pow(Rinput,2) );
        delta_to_shell_Sphereical_Cap_downstream =  CryoTankConsts::CryoDishedHeadInnerRadius - Rsmall_downstream;
      }
      //////////////////////////////////////
      ////Pick the shortest distance tank's inner wall to be the return distance
      /////////////////////////////////////
      if (delta_to_shell_Barrel <= delta_to_shell_Sphereical_Cap_upstream && delta_to_shell_Barrel <= delta_to_shell_Sphereical_Cap_downstream)
      {return delta_to_shell_Barrel;}
      else if (delta_to_shell_Sphereical_Cap_upstream <= delta_to_shell_Barrel && delta_to_shell_Sphereical_Cap_upstream <= delta_to_shell_Sphereical_Cap_downstream)
      {return delta_to_shell_Sphereical_Cap_upstream;}
      else{return delta_to_shell_Sphereical_Cap_downstream;}
    }

    //////////////
    // End Barrel
    //////////////
    else if(Zinput >CryoTankConsts::Zpostion_StartOfdownstreamendBulge){
      double lenghtdownstream_sphereical  = abs(Zinput - CryoTankConsts::startofradiusDownstreamVector);
      Rsmall_downstream = sqrt( pow(lenghtdownstream_sphereical,2) + pow(Rinput,2) );
      delta_to_shell_Sphereical_Cap_downstream  =  CryoTankConsts::CryoDishedHeadInnerRadius - Rsmall_downstream;
      delta_to_shell_Barrel  = CryoTankConsts::CryoInnerVesselInnerRadius - Rinput;

      if(delta_to_shell_Sphereical_Cap_downstream < delta_to_shell_Barrel ){return delta_to_shell_Sphereical_Cap_downstream;}
      else{return delta_to_shell_Barrel; }

      //return delta_to_shell_Sphereical_Cap_downstream;
    }

    else{
      std::cout<<"INside:Return_RECODistance_to_innerTank: a Case failed to correctly calculate the distance to the innerVessels edge "<<std::endl;
      assert(false);
      return -9999;
    }



  }//end of function
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////
  double CenterVolume(double DistancefromEdge){
    double radius = (CryoTankConsts::CryoInnerVesselInnerRadius - DistancefromEdge);
    double Volume = TMath::Pi()*radius*radius*CryoTankConsts::CryoInnerVesselLength ;
    return Volume;

  }
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////
  double CapVolume(double DistancefromEdge ){
    double h = CryoTankConsts::Lenghtofbulge - DistancefromEdge;
    double a = CryoTankConsts::CryoInnerVesselInnerRadius - DistancefromEdge;
    double volume = (TMath::Pi()/6.0) * h *  (3 * a * a + h * h );
    return volume;
  }
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////
  double CryoTankInnerVolume(double DistancefromEdge){
    double volume = CenterVolume(DistancefromEdge) + 2 * CapVolume( DistancefromEdge );
    return volume;
  }
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////
  double CryoTankInnerVolume_metercube(double DistancefromEdge){
    double volume  = CryoTankInnerVolume( DistancefromEdge) * pow(10,-9);
    return volume;

  }

///////////////////////////////////////////////////////////////////////////////
///// Number of Helium atoms
///////////////////////////////////////////////////////////////////////////////
double GetNormNumberHeliumTargets(double  DistancefromEdge){

    double Fiduical_Volume_metercubed = CryoTankInnerVolume_metercube(DistancefromEdge);
    double Helium_Density =  121.403; // kg / m^3 Averged measured denisty taked from DocDB: 9639 slide 21 converted from 7.579 lb/ft^3
    double mass = Fiduical_Volume_metercubed*Helium_Density; // Kg Mass of Helium
    std::cout<< " Helium Mass [kg] = " << mass<<std::endl;
    double helium_molar_mass = 4.002602 / 1000.0; // g/mol * kg/1000g
    double moles = mass / helium_molar_mass;
    std::cout<< "Helium moles = " << moles <<std::endl;
    double NHelium = 6.02214 * moles * pow(10,23); // Moles to number of helium atoms using Avogadro constant

    std::cout<<"NHelium =  " << NHelium << std::endl;

    return NHelium;

  }
  ///////////////////////////////////////////////////////////////////////////////
  /////
  ///////////////////////////////////////////////////////////////////////////////
bool IsInFiducalVolumeFromtheInnerEdgeTRUTH(ChainWrapper& chw, int entry, double Min_distance_toShell ){

  double Xinput = (double)chw.GetValue("mc_vtx",entry,0);
  double Yinput = (double)chw.GetValue("mc_vtx",entry,1);
  double Zinput = (double)chw.GetValue("mc_vtx",entry,2);

  double Rinput = sqrt( pow(Xinput,2) + pow(Yinput,2));



  double Rsmall, delta_to_shell;

  if( (CryoTankConsts::FrontOfCryoInnerVessel <= Zinput) && (Zinput < CryoTankConsts::Zpostion_EndOfupstreamBulge) ){
    double lenghtupstream  = abs(CryoTankConsts::startofradiusUpstreamVector  -  Zinput);
    Rsmall = sqrt( pow(lenghtupstream,2) + pow(Rinput,2) );
    delta_to_shell = CryoTankConsts::CryoDishedHeadInnerRadius - Rsmall;
  }

  else if ((CryoTankConsts::Zpostion_EndOfupstreamBulge <= Zinput) && (Zinput <= CryoTankConsts::Zpostion_StartOfdownstreamendBulge)){
    delta_to_shell  = CryoTankConsts::CryoInnerVesselInnerRadius - Rinput;
  }

  else if((Zinput > CryoTankConsts::Zpostion_StartOfdownstreamendBulge) && (CryoTankConsts::FBackOfCryoInnerVessel >= Zinput) ){
    double lenghtdownstream  = abs(Zinput - CryoTankConsts::startofradiusDownstreamVector);
    Rsmall = sqrt( pow(lenghtdownstream,2) + pow(Rinput,2) );
    delta_to_shell = CryoTankConsts::CryoDishedHeadInnerRadius - Rsmall;
  }

  else{Rsmall = -9999.00;
    return false;
  }

  if (Min_distance_toShell <= delta_to_shell && delta_to_shell  > 0.0 ) {
    return true;

  }

  else{
    return false;

  }

}
///////////////////////////////////////////////////////////////////////////////
/////
///////////////////////////////////////////////////////////////////////////////
bool IsInFiducalVolumeFromtheInnerEdgeTRUTH_new(ChainWrapper& chw, int entry, double Min_distance_toShell )  {

  double Xinput = (double)chw.GetValue("mc_vtx",entry,0);
  double Yinput = (double)chw.GetValue("mc_vtx",entry,1);
  double Zinput = (double)chw.GetValue("mc_vtx",entry,2);

  double Rinput = sqrt( pow(Xinput,2) + pow(Yinput,2));

  double Distance_toEdge = Distance_to_innerTank(Zinput, Rinput );
  if(Distance_toEdge > Min_distance_toShell){return true; }
  else{return false;}
}
///////////////////////////////////////////////////////////////////////////////
/////
///////////////////////////////////////////////////////////////////////////////
double ThetaWRTBeam(double x , double y, double z, double bias) const {
  double beamconstant = -0.05887;
  double pyp = -1.0 *TMath::Sin( beamconstant + bias )*z + TMath::Cos( beamconstant + bias )*y;
  double pzp = TMath::Cos(beamconstant + bias )*z + TMath::Sin( beamconstant + bias )*y;
  double denom2 = pow(x,2) + pow(pyp,2) + pow(pzp,2);
  if( 0. == denom2 || denom2 < 0.0){
    return -999;
  }
  else{
    return acos(pzp / sqrt(denom2) );
  }

}
///////////////////////////////////////////////////////////////////////////////
/////
///////////////////////////////////////////////////////////////////////////////
bool  TRUTH_Is2ndTrk_maxiumAngle_threshold_No_Neutral_WITHProtonAndPion_thresholds(ChainWrapper& chw, int entry, double Max_deg, double Pion_Energy , double Proton_Energy)
{
  int nparticles = (int)chw.GetValue("mc_nFSPart", entry);
  std::vector<double> Angle_trklist;
  std::vector<double> Energy_trklist;
  std::vector<int> PDG_trklist;

  double mass = -9999;
  double KE = -999;


  for(int i = 0; i < nparticles; ++i)
{
    int pdg = (int)chw.GetValue("mc_FSPartPDG",entry,i);
    double px = (double)chw.GetValue("mc_FSPartPx",entry,i);
    double py = (double)chw.GetValue("mc_FSPartPy",entry,i);
    double pz = (double)chw.GetValue("mc_FSPartPz",entry,i);
    double totalE = (double)chw.GetValue("mc_FSPartE",entry,i);
    double rad_angle = ThetaWRTBeam(px, py, pz, 0.0);

            if (pdg == 1000130270)   {mass = 26.90*931.49;} //27Al //MeV
       else if (pdg == 1000020040)   {mass = 3.72*931.49;} //4He //MeV
       else if (pdg == 1000020030)   {mass = 3.02*931.49;} //3He //MeV
       else if (pdg == 1000010020)   {mass = 2.01*931.49;} //2 deuterium //MeV
       else if (pdg == 1000010060)   {mass = 15.99*931.49;} //16Oxygen) //MeV
       else if (pdg == 1000030060)   {mass = 6.02*931.49;} //Li6 //MeV
       else if (pdg == 2000000101)   {mass=0;}
       else{
        // auto Particle_type =  pdg_DATABASEobject->GetParticle(pdg);
        // mass = Particle_type->Mass()*1000; //// Must convert the Mass GeV to MeV
           mass = sqrt(pow(totalE,2) - (px*px + py*py + pz*pz)); // should be in MeV  can't get this to work "TDatabasePDG.h" so i have to solve for mass
  //        std::cout<<"PDG = "<< pdg  <<"Mass [MeV]= "<< mass <<std::endl;
          //KE = Total_E_particles.at(i) - mass;
       }


       totalE = round( totalE * 100.0 ) / 100.0;
       mass = round( mass * 100.0 ) / 100.0; // did this before I got DataBASEobject working

       KE = (totalE - mass)*.001; // Solving for KE and Convert MeV to GeV

       if(mass == -9999 || KE==-999){assert(false && "FAILED TO FIND mass and calculate KE" );}
       if(KE < 0){std::cout<<"This has a Neg KE ,Thats Weird!! the pdg = "<< pdg<<std::endl;assert(false); }

       Angle_trklist.push_back(rad_angle * TMath::RadToDeg());
       PDG_trklist.push_back(pdg);
       Energy_trklist.push_back(KE);

  } // End of loop

  if(PDG_trklist.size()==0) return false; // No FS Particle return false
  for(unsigned int i = 1; i < PDG_trklist.size(); ++i )
     {
       double TrueHardonangle_wrtb = Angle_trklist.at(i);

       if (TrueHardonangle_wrtb < Max_deg &&
        Helium_PDG::pdg_Pi0 != PDG_trklist.at(i) &&
        Helium_PDG::pdg_neutron != PDG_trklist.at(i) &&
        Helium_PDG::pdg_antineutron!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_Genie_bindingE!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_Sigma0!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_antiSigma0!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_Nu_e!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_Nu_mu!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_Lambda0!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_antiLambda0 != PDG_trklist.at(i) &&
        Helium_PDG::pdg_Photon != PDG_trklist.at(i)
        )
        {
          if(PDG_trklist.at(i)== Helium_PDG::pdg_Proton && Energy_trklist.at(i) > Proton_Energy ){return true;}
          else if( (PDG_trklist.at(i) == Helium_PDG::pdg_Pion_pos
                                      ||
                  PDG_trklist.at(i) == Helium_PDG::pdg_Pion_neg) &&
                    Energy_trklist.at(i) > Pion_Energy)
                  {return true;}

          else if(PDG_trklist.at(i) != Helium_PDG::pdg_Pion_pos  &&
                  PDG_trklist.at(i) != Helium_PDG::pdg_Pion_neg  &&
                  PDG_trklist.at(i) != Helium_PDG::pdg_Proton)
                  {return true;}
        }// End of if statement checkin if theres a neutral particle

  } // End of loop
  return false;
}// End of function
///////////////////////////////////////////////////////////////////////////////
/////
///////////////////////////////////////////////////////////////////////////////
bool  TRUTH_Is2ndTrk_maxiumAngle_threshold_No_Neutral(ChainWrapper& chw, int entry, double Max_deg){

   int nparticles = (int)chw.GetValue("mc_nFSPart", entry);
   std::vector<double> Angle_trklist;
   std::vector<int> PDG_trklist;

  for(int i = 0; i < nparticles; ++i)
   {
     int pdg = (int)chw.GetValue("mc_FSPartPDG",entry,i);
     double px = (double)chw.GetValue("mc_FSPartPx",entry,i);
     double py = (double)chw.GetValue("mc_FSPartPy",entry,i);
     double pz = (double)chw.GetValue("mc_FSPartPz",entry,i);

     double rad_angle = ThetaWRTBeam( px , py, pz, 0.0);
     Angle_trklist.push_back(rad_angle * TMath::RadToDeg());
     PDG_trklist.push_back(pdg);
    } // End of loop


  if(PDG_trklist.size()==0) return false;
  for(unsigned int i = 1; i < PDG_trklist.size(); ++i ){
    double TrueHardonangle_wrtb = Angle_trklist.at(i);
    if (TrueHardonangle_wrtb < Max_deg &&
        Helium_PDG::pdg_Pi0 != PDG_trklist.at(i) &&
        Helium_PDG::pdg_neutron != PDG_trklist.at(i) &&
        Helium_PDG::pdg_antineutron!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_Genie_bindingE!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_Sigma0!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_antiSigma0!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_Nu_e!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_Nu_mu!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_Lambda0!= PDG_trklist.at(i) &&
        Helium_PDG::pdg_antiLambda0 != PDG_trklist.at(i) &&
        Helium_PDG::pdg_Photon != PDG_trklist.at(i)
      )
    {return true;}

  }

  return false;
}
}; // End of Class

///////////////////////////////////////////////////////////////////////////////
///// Helium XSec Loop
///////////////////////////////////////////////////////////////////////////////
void runXSecLooper()
{
  // Create the XSecLooper and tell it the input files


  PlotUtils::TargetUtils Target_Tool;
  XSecLooper loop("/pnfs/minerva/persistent/users/cnguyen/ME_playlist_1L_Merg/*.root");
  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setPlaylist(PlotUtils::FluxReweighter::minervame1L);
  loop.setNuPDG(14);
  //Setting the number of Universes in the GENIE error band, default 100 universes put 0 if you do not want universes to be included
  loop.setFiducial(false);
  loop.CheckFiducialStatus();
  loop.setNumUniv(0);
  // Add the total CCQE cross section as a function of energy
  CCInclusive_Helium1DXSec* ccheliumPT=new CCInclusive_Helium1DXSec("ccheliumPT");

  double PT_edges[] ={0.0, 0.075, 0.15, 0.25, 0.325, 0.4 ,0.475, 0.55, 0.7 ,0.85, 1.0, 1.25, 1.5, 2.5};

  int PT_nbins = 13;
  ccheliumPT->setBinEdges(PT_nbins, PT_edges);
  ccheliumPT->setVariable(XSec::kPTLep);
  ccheliumPT->setIsFluxIntegrated(true);
  ccheliumPT->setFluxIntLimits(0.,120.);
  ccheliumPT->setDimension(1);
  ccheliumPT->setUniverses(0);
  double NHelium = ccheliumPT->GetNormNumberHeliumTargets(HeliumCUTConsts::VextexMiniumDisToCryoTankInnerTank);

  double trackerAtomsC = Target_Tool.GetTrackerNCarbonAtoms( 92, true );
  double TrackeratomsH = Target_Tool.GetTrackerElementNAtoms( 1, 92, true);
  double trackerCH = TrackeratomsH + trackerAtomsC ;
  std::cout<<"trackerAtomsCH = "<< trackerCH << std::endl;

  ccheliumPT->setNormalizationType(XSec::kSelfNorm);
  ccheliumPT->setNormalizationValue(NHelium/trackerCH);
  loop.addXSec(ccheliumPT);


  // Add the CCQE differential cross section dsigma/dQ^2
  // flux-integrated over the range 1.5 to 12 GeV
  /*
  CCInclusive_Helium1DXSec* ccheliumPZ=new CCInclusive_Helium1DXSec("ccheliumPZ");

  double PZ_edges[] ={1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0};
  int PZ_nbins = 12;
  ccheliumPZ->setBinEdges(PZ_nbins,PZ_edges);
  ccheliumPZ->setVariable(XSec::kPZLep);
  ccheliumPZ->setIsFluxIntegrated(false);
  ccheliumPZ->setFluxIntLimits(0.,120.);
  ccheliumPZ->setUniverses(0);
  ccheliumPZ->setNormalizationType(XSec::kSelfNorm);
  ccheliumPZ->setNormalizationValue(NHelium/trackerCH);
  loop.addXSec(ccheliumPZ);

  CCInclusive_Helium1DXSec* ccheliumMuonE=new CCInclusive_Helium1DXSec("ccheliumMuonE");

  double E_edges[] ={0.0, 2.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 26.0, 32.0, 42.0, 50.0};
  int E_nbins = 13;
  ccheliumMuonE->setBinEdges(E_nbins,E_edges);
  ccheliumMuonE->setVariable(XSec::kELep);
  ccheliumMuonE->setIsFluxIntegrated(false);
  ccheliumMuonE->setFluxIntLimits(0.,120.);
  ccheliumMuonE->setUniverses(0);
  ccheliumMuonE->setNormalizationValue(NHelium/trackerCH);
  ccheliumMuonE->setNormalizationType(XSec::kSelfNorm);
  loop.addXSec(ccheliumMuonE);

  CCInclusive_Helium1DXSec* ccheliumMuonTheta=new CCInclusive_Helium1DXSec("ccheliumMuonTheta");

  double theta_edges[] ={0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0};
  int theta_nbins = 12;
  ccheliumMuonTheta->setBinEdges(theta_nbins,theta_edges);
  ccheliumMuonTheta->setVariable(XSec::kThetaLep);
  ccheliumMuonTheta->setIsFluxIntegrated(false);
  ccheliumMuonTheta->setFluxIntLimits(0.,120.);
  ccheliumMuonTheta->setUniverses(100);
  ccheliumMuonTheta->setNormalizationType(XSec::kSelfNorm);
  ccheliumMuonTheta->setNormalizationValue(NHelium/trackerCH);
  loop.addXSec(ccheliumMuonTheta);
*/
  // Once everything's set up, Lets actually run the thing ;)
  loop.runLoop();

  // Get the output histograms and save them to file
  TFile fout("GENIEXSecExtract_Helium_CCSemi_ME1L.root", "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i)
     {
       loop.getXSecs()[i]->getXSecHist()->Write();
       loop.getXSecs()[i]->getEvRateHist()->Write();
     }// end of loop

  fout.Close();

std::cout<< "Finished Running runLoop for Helium "<< std::endl;


}// End of rrunXSecLooper

///////////////////////////////////////////////////////////////////////////////
///// Main
///////////////////////////////////////////////////////////////////////////////

int main()
{

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable(); //Needed to look up dictionaries for PlotUtils classes like MnvH1D
  #endif


  TH1::AddDirectory(false);
  runXSecLooper();
  std::cout<<"Finished all scripts in main()"<<std::endl;

 return 0;
}
