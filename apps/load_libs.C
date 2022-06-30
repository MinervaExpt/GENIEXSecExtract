#include "TSystem.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include <iostream>
#include <string>

void load_libs()
{
  gSystem->SetAclicMode(TSystem::kDebug);

  // hide override warnings
  TString makeSharedLib(gSystem->GetMakeSharedLib());
  makeSharedLib.ReplaceAll("-Woverloaded-virtual", "-Wno-overloaded-virtual");
  gSystem->SetMakeSharedLib(makeSharedLib);

  //add GENIEXSecExtact stuff
  {
    gInterpreter->AddIncludePath( gSystem->ExpandPathName("$GENIEXSECEXTRACTROOT") );
    std::string newpath = std::string(gROOT->GetMacroPath()) + ":" + std::string(gSystem->ExpandPathName("$GENIEXSECEXTRACTROOT"));
    gROOT->SetMacroPath( newpath.c_str() );
    gSystem->Load( gSystem->ExpandPathName("$GENIEXSECEXTRACTROOT/libGENIEXSecExtract.so") );
  }

  //add MAT/PlotUtils stuff
  {
    gInterpreter->AddIncludePath( gSystem->ExpandPathName("$PLOTUTILSROOT") );
    std::string newpath = std::string(gROOT->GetMacroPath()) + ":" + std::string(gSystem->ExpandPathName("$PLOTUTILSROOT"));
    gROOT->SetMacroPath( newpath.c_str() );
    gSystem->Load( gSystem->ExpandPathName("$PLOTUTILSROOT/libMAT.so") );
  }

  // for ROOT5
  //gSystem->Load( "libCintex.so" );  // needed to process the dictionaries for the objects
  //Cintex::Enable();
}
