void load_libs()
{
  gSystem->SetAclicMode(TSystem::kDebug);

  //add GENIEXSecExtact stuff
  {
  gInterpreter->AddIncludePath( gSystem->ExpandPathName("$GENIEXSECEXTRACTROOT") );
  string newpath = string(gROOT->GetMacroPath()) + ":" + string(gSystem->ExpandPathName("$GENIEXSECEXTRACTROOT")) + "/GENIEXSecExtract";
  gROOT->SetMacroPath( newpath.c_str() );
  gSystem->Load( gSystem->ExpandPathName("$GENIEXSECEXTRACTROOT/$CMTCONFIG/libGENIEXSecExtract.so") );
  }

  //add PlotUtils stuff
  {
    gInterpreter->AddIncludePath( gSystem->ExpandPathName("$PLOTUTILSROOT") );
    string newpath = string(gROOT->GetMacroPath()) + ":" + string(gSystem->ExpandPathName("$PLOTUTILSROOT")) + "/PlotUtils";
    gROOT->SetMacroPath( newpath.c_str() );
    gSystem->Load( gSystem->ExpandPathName("$PLOTUTILSROOT/$CMTCONFIG/libplotutils.so") );
  }

  gSystem->Load( "libCintex.so" );  // needed to process the dictionaries for the objects
  Cintex::Enable();
}
