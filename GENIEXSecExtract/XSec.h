#ifndef XSEC_H
#define XSEC_H

#include "PlotUtils/ChainWrapper.h"

#include <string>
#include <vector>
#include <TVector3.h>
#include "PlotUtils/HyperDimLinearizer.h"
namespace PlotUtils
{
  class MnvH1D;
  class MnvH2D;
  class HyperDimLinearizer;
}

// A class representing a particular (possibly differential) cross
// section. Needs to know the selection of events, the variable which
// we're plotting, and a few more details. To use, subclass and
// implement the pure virtual function passesCuts, which defines
// whether a given event contributes to this cross section.
class XSec
{
public:
  // First, some types...

  // What are we normalizing the cross section to?
  enum ENormType { kNoNorm, kPerNucleon, kPerNeutron, kPerProton, kSelfNorm, kNNormTypes };

  // Variables whose distributions we might want to produce. If you
  // want to draw a variable that's not included here, either add it
  // to the list (and implement it in getValue() ), and check it in,
  // or override getVariableValue() in your derived class
  enum EVariable
  {
    kNullVariable,
    kPPi0,             // pi0 momentum
    kCosPi0,           // pi0 costheta
    kThetaPi0,         // pi0 theta
    kTPi0,             // pi0 kinetic energy
    kPPiPlus,          // pi+ momentum
    kTPiPlus,          // pi+ kinetic energy
    kTChargedPion,     // charged pion kinetic energy
    kThetaPiPlus,      // pi+ theta WRTBeam
    kThetaChargedPion, // charged pion theta WRTBeam
    kCosPiPlus,        // pi+ costheta
    kx,                // true generated x
    kxExp,             // true experimental x
    ky,                // true y
    kW,                // true W
    kExpW,             // true experimental W
    kq3,
    kq3LowRecoil,      // uses a different Q^2 calculation
    kQ2,               // true Q^2
    kExpQ2,            // true experimental Q^2
    kQ2QE,             // Q^2 based on lepton kinematics, CCQE assumption
    kQ2QEProton,       // Q^2 based on the final state proton, CCQE assumption
    kENuQE,            // Enu based on lepton kinematics, CCQE assumption
    kThetaLep,          // Muon theta (deg)
    kCosLep,            // Muon costheta
    kTLep,              // Muon kinetic energy
    kELep,              // Muon energy
    kPLep,              // Muon momentum 
    kENu,              // True neutrino energy - use this for total cross sections vs Enu
    kEHad,             // True neutrino energy minus true lepton energy (our definition of Ehad)
    kEAvail,
    kPTLep,             // Muon Transverse Momentum
    kPZLep,             // Muon Z Momentum
    kPTLepSquare,       // Muon Transverse Momentum Square
    kTransverseBegin,    //for transverse analysis ->
    kTmuonmomentum,
    kTmuontheta,
    kTprotonmomentum,
    kTprotontheta,
    kTdpt,
    kTdphit,
    kTdalphat,
    kTneutronmomentum,
    kTransverseEnd,      //for transverse analysis
    kTdpTT,              //Double transverse momentum for CC Np Mpi0
    kPZRecoil,             // Muon Z Momentum convoluted with recoil
    kNVars
  };

  // Construct an XSec object. Make sure to call this from derived classes
  XSec(const char* name);
  virtual ~XSec() {};
  
  // Name of the cross section, used for naming the output histogram
  void setName(const char* name) { m_name=name; }
  std::string getName() const { return m_name; }

  // Does the event pass the selection to contribute to this cross
  // section? Override in your derived class
  virtual bool passesCuts(PlotUtils::ChainWrapper& chw, int entry)=0;

  // Set the binning for the output histogram to be nbins equal-sized
  // bins between xmin and xmax
  void setUniformBinning(int nbins, double xmin, double xmax);
  void setUniformBinning(int x_nbins, double xmin, double xmax, int y_nbins, double ymin, double ymax);
  // Set the binning for the output histogram to be nbins bins with
  // edges defined by the binEdges array (of size nbins+1)
  void setBinEdges(int nbins, double* binEdges);
  void setBinEdges(int x_nbins, double* x_binEdges, int y_nbins, double* y_binEdges);

  // Set the variable that this XSec is for
  void setVariable(EVariable var) { m_variable_x=var; }
  void setVariable(EVariable var_x, EVariable var_y) { m_variable_x=var_x; m_variable_y=var_y;}

  // Get the value of the variable that we're using for the given
  // entry in the chain(wrapper). Override this if you want to use a
  // variable that isn't defined in EVariable above
  virtual double getVariableValue(PlotUtils::ChainWrapper& chw, int entry);
  // If working in 2 dimensions, return a vector with both values
  virtual std::vector<double> getVariableValues(PlotUtils::ChainWrapper& chw, int entry);


  // Get the histogram of the cross section. Will be wrong if called
  // before XSecLooper::runLoop has run
  PlotUtils::MnvH1D* getXSecHist() const { return m_xsecHist; }

  // Get or set the histogram of the event rate. Will be wrong if called
  // before XSecLooper::runLoop has run
  PlotUtils::MnvH1D* getEvRateHist() const { return m_evRateHist; }
  void setEvRateHist(PlotUtils::MnvH1D * h) { m_evRateHist = h; }

  // Get the histogram of the cross section. Will be wrong if called
  // before XSecLooper::runLoop has run
  PlotUtils::MnvH2D* get2DXSecHist() const { return m_2DxsecHist; }

  // Is the cross section integrated over flux? (Total cross sections
  // as a function of Enu aren't, differential cross sections are)
  void setIsFluxIntegrated(bool isFluxInt) { m_isFluxIntegrated=isFluxInt; }
  bool isFluxIntegrated() const { return m_isFluxIntegrated; }

  //Set the number of universes to be used in the GENIE error band
  void setUniverses(int nUniv){m_Univ = nUniv;};
  int getUniverses(){return m_Univ;}

  // If we're flux-integrating the cross section, set the neutrino
  // energy limits over which we're integrating the flux
  void setFluxIntLimits(double Emin, double Emax) { m_Emin=Emin; m_Emax=Emax; }
  // Get the minimum and maximum Enu set by setFluxIntLimits
  double getFluxIntEmin() { return m_Emin; }
  double getFluxIntEmax() { return m_Emax; }

  // How do we want the cross section normalized? (per proton, neutron, etc)
  void setNormalizationType(ENormType normType) { m_normType=normType; }
  ENormType getNormalizationType() const { return m_normType; }

  // Set a custom normalization value for this cross section
  void setNormalizationValue( double norm ) { m_normType=kSelfNorm; m_normValue = norm; }
  double getNormalizationValue() const { return m_normValue; }
  
  // Set Dimension number
  void setDimension( int nDim ){ m_dimension = nDim; }
  int  getDimension(){return m_dimension;}

  static Double_t ProtonMass(){ return 938.272/1e3;}//in GeV //google = wiki
  static Double_t NeutronMass(){ return 939.565/1e3;}//in GeV //wiki
  static Double_t MuonMass(){ return 105.65837/1e3; }//in GeV //google = wiki

  static int GetProtonIndex(PlotUtils::ChainWrapper& chw, const int entry);
  //get theta angle wrt beam, where x, y, z in detector coordinate
  static double getTheta( double x, double y, double z );

  static double proton_angle_MAXCut(){ return  70.0;}
  static double proton_momentum_minCut(){ return   450.0;}
  static double proton_momentum_MAXCut(){ return  1200.0;}
  // Set HyperDim
  void setHyperDim(std::vector< std::vector <double> > vec, int type){m_HyperDim = new PlotUtils::HyperDimLinearizer(vec,type);}
protected:

    // DC: Added 26/07/18  
  std::vector<TVector3> mu_mom;
  std::vector<TVector3> pr_mom;
  std::vector<TVector3> pi0_mom;
  void SortParticles();

  static double GetdpTT(const TVector3 nu_dir, const TVector3 lep, const TVector3 had1, const TVector3 had2);

  // Get a list of PDG codes for all the particles in the final state,
  // for entry number 'entry' in the chain(wrapper).
  std::vector<int> getFSPDGs(PlotUtils::ChainWrapper& chw, int entry);

  // Get the value of variable 'var' in entry in the chain
  double getValue(PlotUtils::ChainWrapper& chw, int entry, XSec::EVariable var);

  std::string m_name;       // Name of the cross section
  PlotUtils::MnvH1D* m_xsecHist;          // The output cross section histogram
  PlotUtils::MnvH1D* m_evRateHist;        // the event rate which was used to make the XSec
  PlotUtils::MnvH2D* m_2DxsecHist;          // The output cross section histogram
  EVariable m_variable_x;     // The variable we're making this cross section in (x axis)
  EVariable m_variable_y;     // The variable we're making this cross section in (y axis)
  bool m_isFluxIntegrated;  // Are we making a flux-integrated cross section?
  double m_Emin;            // If flux-integrated, the minimum neutrino energy used
  double m_Emax;            // If flux-integrated, the maximum neutrino energy used
  int m_Univ;               // Number of universes to be included
  ENormType m_normType;     // The normalization type (per proton, nucleon, etc)
  double m_normValue;       // The normalization value for self normalization
  int m_dimension;         // X-Section dimension (1D, 2D...)
  PlotUtils::HyperDimLinearizer*  m_HyperDim;
};

#endif
