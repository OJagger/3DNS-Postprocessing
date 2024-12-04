#! /usr/bin/python
# Script to calculate a blade using Mises.
# --
# Import Python modules
import MisesFunctions
import math
import os
# Define constants
Nrun=3 # Number of times to run ises
KSPEC=1 # Value of SPEC for polar - (1) gives Sinl
WriteIdats=False # Write Idat files for each incidence in a loss loop (False by default)
ext='mises' # File extension used for mises files

# Say hello
print '  Running Mises. Reqd files are blade.%s, stream.%s and ises.%s.' % (ext,ext,ext)
print '  Add export PATH="/home/hot20/bin:$PATH" to ~/.bashrc if programs not found.' 

# Read ises file to get inputs
Ises=MisesFunctions.ReadIsesFile(ext=ext)

# Write key inputs to screen
print '  - Reynolds Number = %.3fE+06' % (Ises['reyn']/1000000)
print '  - Turbulence Level = %.1f%%' % (Ises['turb'])
print u'  - Inlow Angle = %.2f\u00B0' % (Ises['binl'])

# Run iset and ises to initialise the case
MisesFunctions.RunIset(command='iset', ext=ext, Sinl=Ises['sinl'])
MisesFunctions.RunIses(command='ises', ext=ext, Nrun=Nrun)
# Setup for polar: generate appropriate SPEC
if KSPEC==1:
	SPEC=[Ises['sinl']]   
	print u'  - Polar will run at a SINGLE INFLOW ANGLE of %.2f\u00B0' % (((180/math.pi)*math.atan(float(SPEC[0]))))
# Setup for polar: generate the spec file
MisesFunctions.GenerateSpecFile(ext=ext, WriteIdats=WriteIdats, KSPEC=KSPEC, SPEC=SPEC, Positive=True)
# Run polar
MisesFunctions.RunPolar(command='polar',ext=ext, ClearFiles=True)

# Remove files
os.remove('input.ises.mises')
os.remove('input.iset.mises')
