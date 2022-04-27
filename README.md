# GENIEXSecExtract
Closure test programs for MINERvA analyses.  This is how you check that your event loop is self-consistent enough for unfolding.

In your working area opt directory:

`mkdir buildGENIEXSecExtract
cd buildGENIEXSecExtract
cmake ../../GENIEXSecExtract -DCMAKE_INSTALL_PREFIX=`pwd`/.. -DCMAKE_BUILD_TYPE=Release`

This builds the GENIE extraction package to be used by your analysis. 
