#include "GENIEXSecExtract/XSecLooper.h"

#include <cstdlib>
#include "TFile.h"
#include "TH1.h"
#include <PlotUtils/MnvH1D.h>

namespace {
   const double TRACKER_ZMIN = 6117; 
   const double TRACKER_ZMAX = 8193; 
}

//======================================
//! CCQETwoTrack XSec
//======================================
class CCQETwoTrackXSec : public XSec
{

    public:

      CCQETwoTrackXSec( const char* name, int targetZ = 0, bool isQElike = false, 
                        bool isQE = false, bool isRes = false, bool isDIS = false, string proton_threshold = "no_cut_proton" ) : 
        XSec( name ),
        m_nucleus( targetZ ),
        m_QElike( isQElike ),
        m_isQE( isQE ),
        m_isRes( isRes ),
        m_isDIS( isDIS ),
        m_proton_threshold( proton_threshold ) {};
      
      //! Return QE-like signal
      bool isQElikeSignal( ChainWrapper& chw, int entry ) {

         bool signal = false;
         bool pass_proton_energy_threshold = false;
         double proton_momentum_threshold  = 450.0;

         int n_muons         = 0;
         int n_mesons        = 0;
         int n_heavy_baryons = 0;
         int n_em_particles  = 0;
         int n_protons       = 0;

         int nparticles = (int)chw.GetValue("mc_nFSPart",entry);
         for(int i = 0; i < nparticles; i++) {
           int pdg = fabs( (int)chw.GetValue("mc_FSPartPDG",entry,i) );
          
           double px = (double)chw.GetValue("mc_FSPartPx",entry,i);
           double py = (double)chw.GetValue("mc_FSPartPy",entry,i);
           double pz = (double)chw.GetValue("mc_FSPartPz",entry,i);
           double p  = sqrt( px*px + py*py + pz*pz );

           if( pdg == 2212 && p >= proton_momentum_threshold ) pass_proton_energy_threshold = true;

           if( pdg == 13 ) n_muons++;
           else if( pdg == 22   || pdg == 111 || pdg == 11 ) n_em_particles++;
           else if( pdg == 2212 ) { n_protons++;  }
           else if( pdg > 3000  && pdg <  5000 ) n_heavy_baryons++;
           else if( pdg == 211  || pdg == 130 || pdg == 321 || pdg == 311 || pdg == 310 ||
                    pdg == 411  || pdg == 421 || pdg == 431 ) n_mesons++; 
         }


         if( m_proton_threshold == "no_cut_proton" ) {
           if( n_muons == 1 && n_mesons == 0 && n_heavy_baryons == 0 && n_em_particles == 0 && n_protons != 0 ) 
             signal = true;
         } else if( m_proton_threshold == "exactly_one_proton_450MeV" ) {
           if( n_muons == 1 && n_mesons == 0 && n_heavy_baryons == 0 && n_em_particles == 0 && n_protons == 1 ) {
             if( pass_proton_energy_threshold ) signal = true;
           }
         } else if( m_proton_threshold == "atleast_one_proton_450MeV" ) {
           if( n_muons == 1 && n_mesons == 0 && n_heavy_baryons == 0 && n_em_particles == 0 && n_protons != 0 ) {
             if( pass_proton_energy_threshold ) signal = true; 
           }
         }

         return signal;
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

      //! Return the target number of the nucleus 
      int getTargetNumber( double z ) {
        
         double centers[5]   = { 4481.05, 4702.12, 4940.82, 5644.74, 5777.37 };
         double thickness[5] = { 25.75,   25.75,   76.30,   8.0,     13.0    };
        
         for(int i = 0; i < 5; i++) {
           if( fabs(z-centers[i]) <= thickness[i]/2. ) return i+1;
         }

         return 0;
      } 

         
      //! Override this method from the base class to decide what events to include in this selection
      virtual bool passesCuts( ChainWrapper& chw, int entry ) {
        
         //! pass signal
         if( !m_QElike ) {
           if( (int)chw.GetValue("mc_current",entry) != 1 ) return false;
           if( (int)chw.GetValue("mc_intType",entry) != 1 ) return false;
         } else if( m_QElike ) {
           if( !isQElikeSignal(chw,entry) )     return false;

           if( m_isQE ) {
             if( !isQElikeQE(chw,entry) )       return false;
           } else if( m_isRes && !m_isDIS ) {
             if( !isQElikeResonant(chw,entry) ) return false;
           } else if( m_isDIS && !m_isRes ) {
             if( !isQElikeDIS(chw,entry) )      return false;
           } else if( m_isRes &&  m_isDIS ) {
             bool pass = false;
             if( isQElikeResonant(chw,entry) || isQElikeDIS(chw,entry) ) pass = true;
             if( !pass ) return false; 
           }
         }

         //! pass z range
         double z = chw.GetValue("mc_vtx",entry,2);
 
         if( m_nucleus != 0 ) {
           if( (int)chw.GetValue("mc_targetZ",entry) != m_nucleus ) return false;
           else { if( getTargetNumber(z) == 0 ) return false; }
         } else {
           if( z < TRACKER_ZMIN ) return false;
           if( z > TRACKER_ZMAX ) return false;
         } 
        
         return true;

      } //! end function passesCuts

      //! return the binning
      std::vector< double > get2TrackMuonQSquaredBins() {
          std::vector< double > tmp;
          tmp.push_back( 0. );
          tmp.push_back( 0.2 );
          tmp.push_back( 0.4 );
          tmp.push_back( 0.8 );
          tmp.push_back( 1.2 );
          tmp.push_back( 2.0 );
          
          return tmp;
      }

      std::vector< double > get2TrackProtonQSquaredBins() {
          std::vector< double > tmp;
          tmp.push_back( 0. );
          tmp.push_back( 0.15 );
          tmp.push_back( 0.29 );
          tmp.push_back( 0.36 );
          tmp.push_back( 0.46 );
          tmp.push_back( 0.59 );
          tmp.push_back( 0.83 );
          tmp.push_back( 1.33 );
          tmp.push_back( 2.0 );

          return tmp;
      }

      //! input parameters
      int  m_nucleus;
      bool m_QElike;
      bool m_isQE;
      bool m_isRes;
      bool m_isDIS;
      string m_proton_threshold;

};

//======================
//! main: 
//======================
void runXSecLooper_CCQETwoTrack( string& run_type, string& target_type, string& muon_type, string& proton_threshold )
{
   //! initialize
   bool use_one_track_style = false;
   bool use_minos_match     = false;

   //! get the ntuple directory names
   string directory   = "genie_xsection_ntuples";
   string output_name = "genie_xsections";

   if( run_type.find("no_fsi") != string::npos ) {
     directory    = "genie_xsection_no_fsi_ntuples";
     output_name += "_no_fsi";
   } else if( run_type.find("hN") != string::npos ) {
     directory = "genie_xsections_hN_ntuples";
     output_name += "_hN";
   }

   output_name += "_" + target_type;
   if( muon_type.find("one_track") != string::npos ) {
     use_one_track_style = true;
     output_name += "_one_track";
   }

   if( muon_type.find("minos_match") != string::npos ) {
     use_minos_match = true;
     output_name += "_minos_match";
   } 

   if( proton_threshold.find("no_cut_proton") == string::npos ) {
     output_name += "_" + proton_threshold;
   }

   output_name += "_histos.root";

   cout << "	creating a genie file with name = " << output_name << endl;

   //! Create the XSecLooper and tell it the input files
   string input = Form("/minerva/data/users/%s/%s/*.root",getenv("USER"),directory.c_str());
   XSecLooper loop( input.c_str() );

   //! Tell the XSecLooper which neutrino type we're considering (mandatory)
   loop.setNuPDG(14);
  
   //! Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
   loop.setNumUniv(100); 

   //! store cross section names
   const char* channel  = use_one_track_style ? "ccqe" : "ccqe-like";

   std::vector< string > targets;
   if( target_type == "nuclei" ) {
     targets.push_back( "carbon" );
     targets.push_back( "iron" );
     targets.push_back( "lead" );
   } else { targets.push_back( "plastic" ); }


   std::vector< string > particles;
   if(  use_minos_match )     particles.push_back( "muon" );
   if( !use_one_track_style ) particles.push_back( "proton" );

   string processes[5] = { "total", "quasielastic", "resonant", "dis", "resonant_dis" };

   std::vector< string > names;
   for(unsigned int t = 0; t < targets.size(); t++) {
     for(unsigned int p = 0; p < particles.size(); p++) {
       for(int i = 0; i < 5; i++) {
         names.push_back( Form("%s_%s_%s_%s",channel,particles[p].c_str(),targets[t].c_str(),processes[i].c_str()) );
       }
     }
   }

   //! loop over the cross section extraction
   for(unsigned int i = 0; i < names.size(); i++) {

     //! target
     int target = 0;
     if( string(names[i]).find("carbon")    != string::npos ) target = 6;
     else if( string(names[i]).find("iron") != string::npos ) target = 26;
     else if( string(names[i]).find("lead") != string::npos ) target = 82;

     //! is this a QE-like extracted cross section
     bool isQElike = string(names[i]).find("ccqe-like")    != string::npos ? true : false;
     bool isQE     = string(names[i]).find("quasielastic") != string::npos ? true : false;
     bool isRes    = string(names[i]).find("resonant")     != string::npos ? true : false;
     bool isDIS    = string(names[i]).find("dis")          != string::npos ? true : false;

     //! Add the CCQE differential cross section dsigma/dQ^2 
     CCQETwoTrackXSec* xsec = new CCQETwoTrackXSec( names[i].c_str(), target, isQElike, isQE, isRes, isDIS, proton_threshold );

     //! container for bins
     std::vector< double > bins;

     //! get the binning
     if( string(names[i]).find("muon") != string::npos ) {
       bins = xsec->get2TrackMuonQSquaredBins();
     } else if( string(names[i]).find("proton") != string::npos ) {
       bins = xsec->get2TrackProtonQSquaredBins();
     }

     const int size   = (int)bins.size();
     double* q2_bins  = new double[size];
     int     q2_nbins = (int)bins.size() - 1;     

     for(unsigned int j = 0; j < bins.size(); j++) q2_bins[j] = bins[j]; 
     
     //! set cross section data
     xsec->setBinEdges(q2_nbins,q2_bins);
     xsec->setIsFluxIntegrated(true);
     xsec->setUniverses(0);

     if( isQElike ) xsec->setNormalizationType(XSec::kPerNucleon);
     else xsec->setNormalizationType(XSec::kPerNeutron);

     if( string(names[i]).find("muon") != string::npos ) { 
       xsec->setVariable(XSec::kQ2QE);
     } else if( string(names[i]).find("proton") != string::npos ) {
       xsec->setVariable(XSec::kQ2QEProton);
     }

     if( use_minos_match ) {
       xsec->setFluxIntLimits(1.5,10);
     } else { 
       xsec->setFluxIntLimits(0.,100.);
     }

     //! add cross section
     loop.addXSec(xsec);

   } //! end loop over cross sections

   //! run
   loop.runLoop();

   //! get the output histograms and save them to file
   string filename = Form("/minerva/data/users/%s/files/%s/genie/%s",getenv("USER"),getenv("MINERVA_RELEASE"),output_name.c_str());
   TFile f(filename.c_str(),"recreate");
   
   for(unsigned int i = 0; i < loop.getXSecs().size(); i++) {
     loop.getXSecs().at(i)->getXSecHist()->Write();
   }

   return;
}

int main( int argc, char *argv[] )
{
  cout << "Enter running the GENIE Xsection Extraction for the 2track QE analysis" << endl;  

  if( argc == 1 ) {
    cout << "============================================================================================================" << endl;
    cout << "	unable to run the xsection extraction, must specify correct arguments" << endl;
    cout << "	  nuclei, plastic ( option: default = plastic )" << endl;
    cout << "	  1track_style, minos_match, all_muons ( option: default = all_muons )" << endl;
    cout << "	  no_cut_proton, exactly_one_proton_450MeV, atleast_one_proton_450MeV ( option: default = no_cut_proton )" << endl;
    cout << "	  run_xsection, run_no_fsi_xsection, run_hN_xsection ( option: run_xsection )" << endl;
    cout << "	" << endl;
    cout << "=============================================================================================================" << endl;
  } else {

       string target_type      = "plastic";
       string muon_type        = "all_muons";
       string proton_threshold = "no_cut_proton";
       string run_type         = "";

       for(int i = 1; i < argc; i++) {
         if( string(argv[i]).find("no_fsi")      != string::npos ) run_type    = "no_fsi";
         else if( string(argv[i]).find("hN")     != string::npos ) run_type    = "hN";
         else if( string(argv[i]).find("1track") != string::npos ) muon_type   = "one_track_minos_match";
         else if( string(argv[i]).find("minos")  != string::npos ) muon_type   = "minos_match";
         else if( string(argv[i]).find("nuclei") != string::npos ) target_type = "nuclei";
         else if( string(argv[i]) == "exactly_one_proton_450MeV" || string(argv[i]) == "atleast_one_proton_450MeV" ) 
            proton_threshold = string(argv[i]);
       }

       cout << "running genie xsection extraction" << endl;

       runXSecLooper_CCQETwoTrack( run_type, target_type, muon_type, proton_threshold );

       cout << "completed genie xsection extraction" << endl;
       
  } 

  cout << "Exit running the GENIE XSection Extraction for the 2track QE analysis" << endl;
  return 0;
}
