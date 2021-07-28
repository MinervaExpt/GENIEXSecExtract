#!/usr/bin/env python
##############################################
# Look for applications in the apps dir
#   that end in cpp.
# Print lines that should go in the cmt
#   requirements file to build them.
#
###############################################
import os, glob, sys

if not os.environ["GENIEXSECEXTRACTROOT"]:
  sys.exit( "You must setup the package first." )

for f in glob.glob( "%s/apps/*cpp" % os.environ["GENIEXSECEXTRACTROOT"] ):
  app, ext = os.path.splitext( os.path.basename(f) )
  print "application %(app)s  ${GENIEXSECEXTRACTROOT}/apps/%(app)s.cpp" % locals()
