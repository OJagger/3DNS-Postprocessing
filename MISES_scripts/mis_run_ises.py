#! /usr/bin/python
# Script to calculate operating point of a section using Mises.
# --

# Import Python modules
import MisesFunctions
import math
import os

# Define constants
Nrun=3 # Number of times to run ises
ext='mises' # File extension used for mises files

# Say hello
print '  Running Mises. Reqd files are blade.%s, stream.%s and ises.%s.' % (ext,ext,ext)

# Read ises file to get inputs
Ises=MisesFunctions.ReadIsesFile(ext=ext)

# Write key inputs to screen
print '  - Reynolds Number = %.3fE+06' % (Ises['reyn']/1000000)
print '  - Turbulence Level = %.1f%%' % (Ises['turb'])
print '  - Inlow Angle = %.2f deg' % (Ises['binl'])

# Run iset at a single incidence
MisesFunctions.RunIses(command='ises', ext=ext, Nrun=Nrun)

# Remove files
os.remove('input.ises.mises')
