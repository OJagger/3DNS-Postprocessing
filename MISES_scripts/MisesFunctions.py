# Module containing function definitions for running Mises
# To use, add "import MisesFunctions" to your python script.
# To call, use (for example) MisesFunctions.RunIset(...)

###############################################################################
# FUNCTIONS FOR RUNNING PROGRAMS IN THE MISES SUITE                           #
###############################################################################

# Function to run iset
# -----------------------------------------------------------------------------
def RunIset(command='iset', ext='mises', Sinl=0.89526):
	"""Function to run iset. Inputs:
	1. command: Command to run iset
	2. ext: File extention for mises imputs """
	# Import os for operating system interaction and sys for exiting
	import os
	import sys
	print '  Writing %s input file and running %s...' % (command,command)
	# Generate file names
	IsetInfileName='input.%s.%s' % (command,ext)
	IsetOutfileName='output.%s.%s' % (command,ext)
	# Create iset input file
	IsetInfile=open(IsetInfileName,'w')
	#IsetInfile.write('%-0.0f \n \n%-0.0f \n \n%-0.0f \n%-0.0f \n%-0.0f \n' % (1,2,3,4,0))
	IsetInfile.write('%-0.0f \n%-0.0f \n%-0.0f \n \n%-0.0f \n%-0.0f \n%-0.0f \n' % (1,-999,2,3,4,0))
	#IsetInfile.write('%-0.0f \n \n%-0.0f \n%s \n \n \n%-0.0f \n%-0.0f \n%-0.0f \n' % (1,2,'D',3,4,0))
	#IsetInfile.write('%-0.0f \n%-0.0f \n%-0.0f \n%s \n \n \n%-0.0f \n%-0.0f \n%-0.0f \n' % (1,-999,2,'D',3,4,0))
	IsetInfile.close()
	# Run iset using the input file	
	os.system("%s %s < %s > %s" % (command,ext,IsetInfileName,IsetOutfileName))
	# Check output
	GridOK=0
	for line in open(IsetOutfileName):
		if "Grid not initialized" in line:
			sys.exit("  Error: Grid not initialized!")
		if "Number of streamlines" in line:
			print '  Grid check ok.'
			GridOK=1
	if GridOK==0:
		sys.exit("  Error: Grid generation failed!")
	return	

# Function to run ises
# -----------------------------------------------------------------------------
def RunIses(command='ises', ext='mises', Nrun=3):
	"""Function to run ises. Inputs:
	1. command: Command to run ises
	2. ext: File extention for mises imputs """
	# Import os for operating system interaction and sys for exiting
	import os
	import sys
	print '  Writing %s input file and running %s %0.0f times...' % (command,command,Nrun)
	# Generate file names
	IsesInfileName='input.%s.%s' % (command,ext)
	IsesOutfileName='output.%s.%s' % (command,ext)
	# Write ises input file
	IsesInfile=open(IsesInfileName,'w')
	for ii in range(Nrun):
		IsesInfile.write('%-0.0f \n' % (15))
	# Finish off file
	IsesInfile.write('%-0.0f \n' % (0))
	IsesInfile.close()
	# Run iset using the input file
	os.system("%s %s < %s > %s" % (command,ext,IsesInfileName,IsesOutfileName))
	# Check output
	ConvOK=0
	storeline = []
	for line in open(IsesOutfileName):
		storeline.append(line)
	#for line in storeline[len(storeline)-5:len(storeline)]
	if "Converged on tolerance" in storeline[len(storeline)-3]:
		print '  Convergence check ok.'
		ConvOK=1
	if ConvOK==0:
		sys.exit("  Error: Convergence check failed!")
	return

# Function to run polar
# -----------------------------------------------------------------------------
def RunPolar(command='polar', ext='mises', ClearFiles=True):
	"""Function to run polar. Inputs:
	1. command: Command to run polar
	2. ext: File extention for mises imputs """
	# import os for operating system interaction
	import os
	# Generate Filenames
	PolarOutfileName='output.%s.%s' % (command,ext)
	if ClearFiles: # Clear polar and polarx files
		print '  Clearing polar.%s and polarx.%s files...' % (ext,ext)
		open('polar.%s' % (ext),'w')
		open('polarx.%s' % (ext),'w')
	# Run polar
	print '  Running %s...' % (command)
	os.system("%s %s > %s" % (command,ext,PolarOutfileName))
	return

# Function to run mises boundary layer solver independently
# -----------------------------------------------------------------------------
def RunMblrun(command='mblrun', ext='mises'):
	"""Function to run mblrun. Inputs:
	1. command: Command to run mblrun
	2. ext: File extention for mises imputs """
    # Import os for operating system interaction and sys for exiting
	import os
	import sys
	# Tell code what the input and output files are called
	BLInfileName='input.%s.%s' % (command,ext)
	BLDataName='%s.%s' % (command,ext)
	BLOutfileName='output.%s.%s' % (command,ext)
	print '  Running %s...  Interupt with ctrl+c if nothing happens within a few seconds.' % (command)
	# Run boundary layer solver
	os.system("%s %s %s > %s" % (command,BLInfileName,BLDataName,BLOutfileName))
	# Check output
	ConvOK=0
	for line in open(BLOutfileName):
		if "Convergence failed" in line:
			ConvOK=ConvOK+1
	if ConvOK>=4:
		sys.exit("  Error: Over four convergence errors!")
	return

###############################################################################
# FUNCTIONS FOR READING AND WRITING FILES FOR USE BY THE MISES CODES          #
###############################################################################

# Function to read an ises.xxx file
# ----------------------------------------------------------------------------
def ReadIsesFile(ext='mises'):
	"""Function to read in an ises.xxx file. Inputs:
	1. ext: File extention for mises imputs """
	# Imports to gain maths functions
	import math 
	# Initialise output dictionary (like a struct in Matlab)
	Ises={}
	# Define filename and state progress
	IsesFileName='ises.%s' % (ext)
	print '  Reading %s file...' % (IsesFileName)
	# Open file, read and split into lines
	IsesFile=open(IsesFileName,'r')
	data=IsesFile.read()
	data=data.splitlines()
	# Extract global variables (1st line), ignoring comments
	t=data[0].split()
	Ises['gvar']=[]
	for ii in range(len(t)):
		try:
			tt=int(t[ii])
			Ises['gvar'].append(tt)			
		except ValueError:
			pass
	# Extract global constraints (2nd line), ignoring comments
	t=data[1].split()
	Ises['gcon']=[]
	for ii in range(len(t)):
		try:
			tt=int(t[ii])
			Ises['gcon'].append(tt)
		except ValueError:
			pass
	# Extract inflow boundary conditions
	t=data[2].split()
	Ises['minl']=float(t[0])
	Ises['p1pt']=float(t[1])
	Ises['sinl']=float(t[2])
	Ises['xinl']=float(t[3])
	try:
		Ises['v1at']=float(t[4]);
	except ValueError:
		pass	
	# Extract outflow boundary conditions
	t=data[3].split()
	Ises['mout']=float(t[0])
	Ises['p2pt']=float(t[1])
	Ises['sout']=float(t[2])
	Ises['xout']=float(t[3])
	try:
		Ises['v2at']=float(t[4]);
	except ValueError:
		pass
	# Extract splitter mass frac and non-adiabatic wall (could also get gamma here too)
	t=data[4].split()
	Ises['mfr']=float(t[0])
	Ises['hwrat']=float(t[1])	
	# Extract reynolds no and ncrit (could also get ncrit2 and tsrat here too)
	t=data[5].split()
	Ises['reyn']=float(t[0])
	Ises['ncrit']=float(t[1])
	# Extract forced transition point
	t=data[6].split()
	Ises['strp']=[float(t[0]), float(t[1])]
	# Extract algoritm controls
	t=data[7].split()
	Ises['ismom']=int(t[0])
	Ises['mcrit']=float(t[1])
	Ises['mucon']=float(t[2])
	# Extract stream-tube thickness mode amplitudes if present
	if (len(data)-1)>=8:
		t=data[8].split()
        Ises['bvr']=[float(t[0]), float(t[1])]
	# Extract geometry movement, scaling and rotation mode amplitudes if present
	if (len(data)-1)>=9:
		t=data[9].split()
		Ises['movx']=float(t[0])
		Ises['movy']=float(t[1])
		Ises['scal']=float(t[2])
		Ises['rota']=float(t[3])
	# Extract geometry shape mode amplitudes if present
	if (len(data)-1)>=10:
		t=data[10].split()
		Ises['kmod']=float(t[0])
		Ises['gmod']=float(t[1])
	# Close the ises.xxx file to free resources
	IsesFile.close()
	# Deal with the turbulence intensity parameters
	if Ises['ncrit']>=0: # Just use value for e^n transition model
		Ises['ncrit']=Ises['ncrit']; 
		TuPrime=(100*(math.exp((8.43+Ises['ncrit'])/(-2.4))));
		Ises['turb']=((2.7/2)*(math.log((1+(TuPrime/2.7))/(1-(TuPrime/2.7)))));
	elif Ises['ncrit']<0: # Calculate the value of ncrit used by mises (see manuals)
		Ises['turb']=-Ises['ncrit'];
		TuPrime=(2.7*math.tanh((Ises['turb'])/2.7));
		Ises['ncrit']=(-8.43-(2.4*math.log(TuPrime/100)));
	# Generate human readable angles
	Ises['binl']=((180/math.pi)*math.atan(Ises['sinl']))
	Ises['bout']=((180/math.pi)*math.atan(Ises['sout']))
	# Return the Ises dictionary for further use
	return Ises

# Function to read data from a polarx file
# -----------------------------------------------------------------------------
def ReadPolarx(ext='mises'):
    """Function to read in polarx files.
    """
    pass
    return

# Function for generating the spec.xxx file for polar
# -----------------------------------------------------------------------------
def GenerateSpecFile(ext='mises', WriteIdats=False, KSPEC=1, SPEC=[0.90], Positive=False):
	"""Function to generate a spec file for polar. Inputs:
	1. ext: File extension for mises inputs
	2. WriteIdats: Switch to write all the idat files
	3. KSPEC: Type of BC in spec file (see Mises manual)
	4. SPEC: Vector of KSPEC variable """
	# Define filename and state progress
	SpecFileName='spec.%s' % (ext)
	print '  Writing %s file...' % (SpecFileName)
	# Open file (if existing, will be erased)
	SpecFile=open(SpecFileName,'w')
	# Write Header
	SpecFile.write('%-2.0f \n' % (KSPEC))	
	# Write File
	for ii in range(len(SPEC)):
		# Write every other idat file
		if WriteIdats and ii%2==0:
			if Positive:
				SpecFile.write('%-10.6f %-2.0f \n' % (SPEC[ii],ii/2+1));
			else:
				SpecFile.write('%-10.6f %-2.0f \n' % (SPEC[ii],ii/2+51));
		else:
			SpecFile.write('%-10.6f \n' % (SPEC[ii]));
	# Close the spec.xxx file to free resources
	SpecFile.close()
	return

###############################################################################
# FUNCTIONS FOR PLOTTING THE RESULTS OF MISES CALCULATIONS                    #
###############################################################################

# Function to plot data
# -----------------------------------------------------------------------------
def PlotMisesData():
	pass

# Function to plot data
# -----------------------------------------------------------------------------
def PlotMblrunData():
    pass

