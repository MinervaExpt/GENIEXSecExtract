# GENIEXSecExtract
Closure test programs for MINERvA analyses.  This is how you check that your event loop is self-consistent enough for unfolding.

`` git clone https://github.com/MinervaExpt/GENIEXSecExtract.git  ``

In your working area opt directory:
1. Make a build dir for this package
``mkdir buildGENIEXSecExtract ``
2. Enter build dir
``cd buildGENIEXSecExtract ``
3. Run cmake
``cmake ../../GENIEXSecExtract -DCMAKE_INSTALL_PREFIX=`pwd`/.. -DCMAKE_BUILD_TYPE=Release ``
4. Install the build
``make install``


This builds the GENIE extraction package to be used by your analysis. 
