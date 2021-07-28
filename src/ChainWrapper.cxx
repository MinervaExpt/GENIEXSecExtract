#include "GENIEXSecExtract/ChainWrapper.h"

#include <fstream>
#include <iostream>
#include "TChain.h"

#ifndef __CINT__
#include "glob.h"
#endif

ChainWrapper::ChainWrapper(const char* name)
  : TreeWrapper(new TChain(name))
{

  wrappingChain=true;
}

//===========================================================================

int ChainWrapper::Add(const char* globStr)
{
  // Urgh, downcast
  TChain* ch=(TChain*)tree;

  // if globStr ends with .txt or .dat, it's a playlist
  std::string globAsStr = globStr;
  if(globAsStr.find(".txt") != std::string::npos || globAsStr.find(".dat") != std::string::npos)
  {
    ifstream f(globStr);
    if(!f) {
      cerr << "unable to open playlist \"" << globStr << "\"" << endl;
      exit(1);
    }

    int n = 0;
    for(string s; getline(f, s); ) {
      if(!s.length()) continue; // ignore empty lines
      if(s[0] == '#') continue; // ignore comments

      ch->Add(s.c_str()); n++;
    }

    return n;
  }

  glob_t g;
  glob(globStr, 0, 0, &g);
  for(int i=0; i<(int)g.gl_pathc; ++i){
    ch->Add(g.gl_pathv[i]);
  }
  int ret=g.gl_pathc;

  globfree(&g);

  return ret;
}
