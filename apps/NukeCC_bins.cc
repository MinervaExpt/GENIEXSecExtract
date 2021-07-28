#include <vector>
#include "GENIEXSecExtract/XSec.h"
#include "TMath.h"

using namespace std;

//==============================
// Binning utilities
//==============================
void AddBins( std::vector<double>& binsLowEdge, const int nBins, const double binWidth )
{
  // If the bins have no defined low entry, assume it should be 0
  double x = 0.;
  if( binsLowEdge.size() > 0 )
    x = binsLowEdge[ binsLowEdge.size() - 1 ];
  else
    binsLowEdge.push_back( x );

  // add nBins spaced by this width
  for( int iBin = 0; iBin < nBins; ++iBin )
  {
    x += binWidth;
    binsLowEdge.push_back( x );
  }
}

void GetBins( XSec::EVariable var, std::vector<double>& bins, bool fine = false, bool inel = false )
{
  bins.clear();

  switch(var)
  {
    case XSec::kx:
      if(fine)
      {
        //0-.3 in .01 bins
        AddBins( bins, 11, .1 );
        //.3-1.25 in .025 bins
        //AddBins( bins, 48, .025 );
      }
      else
      {
        /*bins.push_back(0.);

        // 0-.1 (shadowing)
        bins.push_back(0.01);
        bins.push_back(0.1);

        //.1-.3 (antishadowing)
        bins.push_back(0.2);
        bins.push_back(0.3);

        //.3-.7 (EMC)
        bins.push_back(0.5);
        bins.push_back(0.75);

        //elasticish; overflow
        bins.push_back(0.9);
        bins.push_back(1.1);
        bins.push_back(1.5);*/
	
       //New bins
    bins.push_back(0.);

    // 0-.1 (shadowing)
    bins.push_back(0.005);
    bins.push_back(0.05);
    bins.push_back(0.1);

    //.1-.3 (antishadowing)
    bins.push_back(0.2);
    bins.push_back(0.3);

    //.3-.7 (EMC)
    bins.push_back(0.5);
    bins.push_back(0.8);

    //elasticish; overflow
    //bins.push_back(0.9);
    //bins.push_back(1.1);
    bins.push_back(1.3);
    
   /*    bins.push_back(0.0);
       bins.push_back(0.1);
       bins.push_back(0.2);
       bins.push_back(0.3);
       bins.push_back(0.4);
       bins.push_back(0.75);
       
       //overflow
       bins.push_back(10.0); */
      }

      break;
        case XSec::kxExp:
      if(fine)
      {
        //0-.3 in .01 bins
        AddBins( bins, 33, .1 );
        //.3-1.25 in .025 bins
        //AddBins( bins, 48, .025 );
      }
      else
      {
       //New bins
   // bins.push_back(0.);

    // 0-.1 (shadowing)
    bins.push_back(0.005);
    bins.push_back(0.05);
    bins.push_back(0.1);

    //.1-.3 (antishadowing)
    bins.push_back(0.2);
    bins.push_back(0.3);

    //.3-.7 (EMC)
    bins.push_back(0.5);
    bins.push_back(0.8);

    //elasticish; overflow
    //bins.push_back(0.9);
    //bins.push_back(1.3);
    //bins.push_back(1.5);
   /* 
       //New bins
       bins.push_back(0.0);
       bins.push_back(0.1);
       bins.push_back(0.2);
       bins.push_back(0.3);
       bins.push_back(0.4);
       bins.push_back(0.75);
       bins.push_back(1.0);
       
       //overflow
       bins.push_back(10.0); */
      }

      break;  
      
    case XSec::ky:
      if( fine )
      {
        AddBins(bins, 50, .02);
      }
      else
      {
        double tmpBins[] = { 0, 0.0195312, 0.0539453, 0.102092, 0.159776, 0.229381, 0.313886, 0.430302, 0.612201, 1.};
        bins.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
      }
      break;
    case XSec::kW:
      if(fine)
      {
        bins.push_back(0.);
        bins.push_back(.5);
        bins.push_back(.925);
        bins.push_back(.95);
        bins.push_back(1.);
        //up to 3GeV
        AddBins( bins, 40, .05 );
        //up to 8Gev
        AddBins( bins, 50, .1 );
      }
      else
      {
        double tmpBins[]  = {0, 2., 2.5, 3.5, 5., 8., 10., 20.};
        //double tmpBins[] = {0, 1.0, 1.3, 1.7, 2., 2.5, 3., 3.6, 4.2, 5., 6., 10., 25., 100. };
        bins.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
      }
      break;
    case XSec::kExpW:
      if(fine)
      {
        bins.push_back(0.);
        bins.push_back(.5);
        bins.push_back(.925);
        bins.push_back(.95);
        bins.push_back(1.);
        //up to 3GeV
        AddBins( bins, 40, .05 );
        //up to 8Gev
        AddBins( bins, 50, .1 );
      }
      else
      {
        double tmpBins[]  = {0, 2., 2.5, 3.5, 5., 8., 10., 20.};
        //double tmpBins[] = {0, 1.0, 1.3, 1.7, 2., 2.5, 3., 3.6, 4.2, 5., 6., 10., 25., 100. };
        bins.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
      }
      break;  
      
    case XSec::kQ2:
      if(fine)
      {
        //0-.1
        AddBins( bins, 20, .005 );
        //.1-.25
        AddBins( bins, 15, .01 );
        //.25-1
        AddBins( bins, 30, .025 );
        //1-2
        AddBins( bins, 20, .05 );
        //2-5
        AddBins( bins, 30, .1 );
        //5-20
        AddBins( bins, 30, .5);
      }
      else
      {
        double tmpBins[] = { 0, 1.0, 1.25, 3.0, 5.0, 8.0, 30., 100. };
        //double tmpBins[] = { 0, 0.015, 0.07, 0.2, 0.5, 1.0, 1.75, 3., 6., 10., 20., 50., 100. };
        bins.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
      }
      break;
    case XSec::kExpQ2:
      if(fine)
      {
        //0-.1
        AddBins( bins, 20, .005 );
        //.1-.25
        AddBins( bins, 15, .01 );
        //.25-1
        AddBins( bins, 30, .025 );
        //1-2
        AddBins( bins, 20, .05 );
        //2-5
        AddBins( bins, 30, .1 );
        //5-20
        AddBins( bins, 30, .5);
      }
      else
      {
        double tmpBins[] = { 0, 1.0, 1.25, 3.0, 5.0, 8.0, 30., 100. };
        //double tmpBins[] = { 0, 0.015, 0.07, 0.2, 0.5, 1.0, 1.75, 3., 6., 10., 20., 50., 100. };
        bins.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
      }
      break;  
      
    case XSec::kENu:
      //.25GeV bins up to 100
      if(fine)
      {
        AddBins( bins, 400, .25 );
      }
      else
      {
        //inelastic sample has very poor efficiency below 6 GeV, use bigger bins
//        if( inel )
            if( true )
        {
            double tmpBins[]  = { 5., 7.5, 10., 12.5, 15., 20., 30., 40., 50.};//, 60., 70.,80.,90.,100.,110.,120. };
//            double tmpBins[]  = { 5000., 7500., 10000., 12500, 15000., 20000., 30000., 40000., 50000. };
        bins.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
         /*AddBins( bins, 6, 1. );
         //2 GeV bins up to 10
         AddBins( bins, 2, 2. );
         //5GeV bins up to 30GeV
         AddBins( bins, 4, 5. );
         //10GeV Bins up to 50 GeV
         AddBins( bins, 2, 10 );
         AddBins( bins, 1, 50. ); //50-100*/
	 
	 //1 GeV bins up to 5 GeV
         //These bins are largely outside the DIS acceptance range
         //But making them 1 GeV wide menas the framework will automatically bin width normalize to 1 GeV...
         /*AddBins( bins, 5, 1. );
         //5GeV bins up to 30GeV
         AddBins( bins, 6, 2.5 );
         //10GeV Bins up to 50 GeV
         AddBins( bins, 2, 5. );
         AddBins( bins, 2, 10. ); //50-100
         //Overflow
         AddBins( bins, 1, 50.);
         AddBins( bins, 1, 100.);
	 */
        }
        
	else
        {
        double tmpBins[]  = { 5., 7.5, 10., 12.5, 15., 20., 30., 40., 50. };
        bins.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
         /*AddBins( bins, 6, 1. );
         //2 GeV bins up to 10
         AddBins( bins, 2, 2. );
         //5GeV bins up to 30GeV
         AddBins( bins, 4, 5. );
         //10GeV Bins up to 50 GeV
         AddBins( bins, 2, 10 );
         AddBins( bins, 1, 50. ); //50-100*/
	 
	 //1 GeV bins up to 5 GeV
         //These bins are largely outside the DIS acceptance range
         //But making them 1 GeV wide menas the framework will automatically bin width normalize to 1 GeV...
         //AddBins( bins, 5, 1. );
         //5GeV bins up to 30GeV
         //AddBins( bins, 5, 5. );
         //10GeV Bins up to 50 GeV
         //AddBins( bins, 2, 10 );
         //AddBins( bins, 1, 50. ); //50-100
         //Overflow
         //AddBins( bins, 1, 100);
       }
      }
      break;
    case XSec::kTLep:
    {
    // 1 deg bins up to 8
    //AddBins( binsLowEdge, 10, 1.5 );
    AddBins( bins, 8, 1.0);
    // 2 deg bins up to 14, 
    //AddBins( binsLowEdge, 3, 2.0 );
    // 3 deg bins up to 26
    AddBins( bins, 4, 3. );
    //--- overflow --//
    // 5 deg bins up to 51
    AddBins( bins, 5, 5. );
    // 2 20 deg bins up to 91
    AddBins( bins, 2, 20. );
    }
    case XSec::kELep:
    {
        //for DIS: one underflow bin (0 - 2 GeV) same otherwise
        double tmpBins[]  = { 2., 3., 4., 5., 6., 8., 10., 13., 16., 20., 25., 35., 50., 120. };
        bins.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
      //same bins as Enu (but never with inelastic binning)
      //GetBins(XSec::kENu, bins, fine, false);
      break;

    case XSec::kThetaLep:
      if( fine )
      {
        //0-25
        AddBins( bins, 100 , .25 );
        //25-45
        AddBins( bins, 20, 1. );
      }
      else
      {
       // 1 deg bins up to 8
       //AddBins( binsLowEdge, 10, 1.5 );
       AddBins( bins, 8, 1.0);
       // 2 deg bins up to 14, 
       //AddBins( binsLowEdge, 3, 2.0 );
       // 3 deg bins up to 26
       AddBins( bins, 4, 3. );
       //--- overflow --//
       // 5 deg bins up to 51
       AddBins( bins, 5, 5. );
       // 2 20 deg bins up to 91
       AddBins( bins, 2, 20. );
      }
      /*{
        double tmpBins[]  = { 0, 2., 3., 4., 5., 6., 8., 10., 13., 16., 20., 25., 35., 50 };
        bins.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
        // 1.5 deg bins up to 15
        //AddBins( bins, 10, 1.5 );
        // 2 deg bins up to 17, our theta cutoff
        //AddBins( bins, 1, 2. );
        // 3 deg bins up to 25
        //AddBins( bins, 3, 3. );
        //--- overflow --//
        // 5 deg bins up to 50
        //AddBins( bins, 5, 5. );
        // 2 20 deg bins up to 90
        //AddBins( bins, 2, 20. );
      }*/
      break;

    case XSec::kCosLep:
      GetBins(XSec::kThetaLep, bins, fine);
      for( vector<double>::iterator i = bins.begin(); i != bins.end(); ++i )
          *i = cos( (TMath::Pi()/180.0) * (*i) );
//        *i = cos( TMath::DegToRad() * (*i) );
      reverse( bins.begin(), bins.end() );
      break;
    case XSec::kEHad:
      if( fine)
      {
        //0-1
        AddBins( bins, 40, .025 );
        //1-5
        AddBins( bins, 40, .1 );
        //5-10
        AddBins( bins, 20, .25 );
        //10-50
        AddBins( bins, 20, 2. );
      }
      else
      {
        bins.push_back(0.);
        bins.push_back(0.25);
        bins.push_back(1);
        bins.push_back(2.5);
        bins.push_back(5);
        bins.push_back(10);
        bins.push_back(15);
        bins.push_back(20.);
        bins.push_back(30.);
        bins.push_back(50.);
        bins.push_back(100.);
      }
      break;

    default:
      break;
  }
}
