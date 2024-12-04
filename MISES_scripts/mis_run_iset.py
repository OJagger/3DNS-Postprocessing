#! /usr/bin/python
# Script to setup a MISES calculation with iset
# --

# Import Python modules
import MisesFunctions
import math
import os

# Define constants
ext='mises' # File extension used for mises files

# Say hello
print '  Running Mises. Reqd files are blade.%s, stream.%s and ises.%s.' % (ext,ext,ext)

# Read ises file to get inputs
Ises=MisesFunctions.ReadIsesFile(ext=ext)

# Run iset to initialise the case
MisesFunctions.RunIset(command='iset', ext=ext, Sinl=Ises['sinl'])

# Remove files
os.remove('input.iset.mises')
