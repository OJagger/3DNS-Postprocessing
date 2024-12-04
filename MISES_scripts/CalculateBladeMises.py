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
SinglePolar=False # Run polar at a single value of SPEC (True), or run at many different values of SPEC, i.e. for a loss loop (False)
WriteIdats=True # Write Idat files for each incidence in a loss loop (False by default)
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

# Run iset, ises and polar for sinle or multiple incidence
if SinglePolar: # Run for just one value of SPEC
	# Run iset and ises to initialise the case
	MisesFunctions.RunIset(command='iset', ext=ext, Sinl=Ises['sinl'])
	#os.system('iset mises')
	MisesFunctions.RunIses(command='ises', ext=ext, Nrun=Nrun)
	# Setup for polar: generate appropriate SPEC
	if KSPEC==1:
		SPEC=[Ises['sinl']]
		print u'  - Polar will run at a SINGLE INFLOW ANGLE of %.2f\u00B0' % (((180/math.pi)*math.atan(float(SPEC[0]))))
	# Setup for polar: generate the spec file
	MisesFunctions.GenerateSpecFile(ext=ext, WriteIdats=WriteIdats, KSPEC=KSPEC, SPEC=SPEC, Positive=True)
	# Run polar
	MisesFunctions.RunPolar(command='polar',ext=ext, ClearFiles=True)
else: # Run over multiple values of SPEC
	# Remove files beginning with 'idat'
	path = os.getcwd()
	files = os.listdir(path)
	for ii in range(len(files)):
		if files[ii][0:4] == 'idat':
			if 'idatfiles'	in locals():		
				idatfiles.extend([files[ii]])
			else:
				idatfiles = [files[ii]]
	if 'idatfiles'	in locals():
		for ii in range(len(idatfiles)):
			os.remove(idatfiles[ii])
	# Run iset and ises to initialise the case
	#MisesFunctions.RunIset(command='iset', ext=ext, Sinl=Ises['sinl'])
	os.system('iset mises')		
	os.system('ises mises')	
	#MisesFunctions.RunIses(command='ises', ext=ext, Nrun=Nrun)
	# Setup for polar: generate vector of SPEC values (corresponding to KSPEC)
	if KSPEC==1: # Generate incidences. This requires polar to run twice, once for positive angles and again for negative
		#MaxPosAngle=8 # Positive range 
		#MaxNegAngle=-6 # Negative range
		MaxPosAngle=10 # Positive range
		MaxNegAngle=-8 # Negative range
		Resolution=2 # Number of points per degree
		SPECpos=[x/float(Resolution) for x in range(0,(MaxPosAngle*Resolution)+1,1)] # Create +ve pertubation list (currently 0 to 12 in 0.5 degree steps)
		SPECneg=[x/float(Resolution) for x in range(0,(MaxNegAngle*Resolution)-1,-1)] # Create -ve pertubation list (currently 0 to -12 in 0.5 degree steps)
		SPECpos=[(x+Ises['binl']) for x in SPECpos] # Add design incidence to get list of flow angles
		SPECneg=[(x+Ises['binl']) for x in SPECneg] 
		SPECpos=[math.tan(x*(math.pi/180)) for x in SPECpos] # Convert to tangent of angle in radians
		SPECneg=[math.tan(x*(math.pi/180)) for x in SPECneg] 
		# Work out max and min angles and write to screen
		SPECmin=((180/math.pi)*math.atan(float(SPECneg[-1])))
		SPECmax=((180/math.pi)*math.atan(float(SPECpos[-1])))
		print u'  - Polar will run a loss loop with %i INFLOW ANGLES from %.2f to %.2f\u00B0' % ((len(SPECpos)+len(SPECneg)),SPECmin,SPECmax)
		print '  - Remember to check variables and constraints in ises.%s' % (ext)
		# Setup for polar: generate the spec file for positive changes in SPEC
		MisesFunctions.GenerateSpecFile(ext=ext, WriteIdats=WriteIdats, KSPEC=KSPEC, SPEC=SPECpos, Positive=True)
		# Run polar
		MisesFunctions.RunPolar(command='polar', ext=ext, ClearFiles=True)
		# Rename files
		os.rename('output.polar.%s' % (ext) ,'output.polar1.%s' % (ext))
		os.rename('spec.%s' % (ext) ,'spec1.%s' % (ext))
		# Polar sometimes fails to converge on the first point when polar is re-run (something to do with polar.xxx file). Uncomment next line if this happens
		#os.rename('polar.%s' % (ext) ,'polar.%s1' % (ext))
		# Run iset and ises to initialise the case again ready for negative changes
		#MisesFunctions.RunIset(command='iset', ext=ext, Sinl=Ises['sinl'])
		#MisesFunctions.RunIses(command='ises', ext=ext, Nrun=Nrun)
		# Setup for polar: generate the spec file for negative changes in SPEC
		os.system('iset mises')		
		os.system('ises mises')	
		MisesFunctions.GenerateSpecFile(ext=ext, WriteIdats=WriteIdats, KSPEC=KSPEC, SPEC=SPECneg, Positive=False)
		# Run polar
		#MisesFunctions.RunPolar(command='polar', ext=ext, ClearFiles=False)
		os.system('polar mises')
		# Rename files
		#os.rename('output.polar.%s' % (ext) ,'output.polar5.%s' % (ext))
		#os.rename('spec.%s' % (ext) ,'spec5.%s' % (ext))

# Remove files
os.remove('input.ises.mises')
os.remove('input.iset.mises')

# Read in polarx
# MisesFunctions.ReadPolarx(ext=ext);

# Plot results
