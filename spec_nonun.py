import os
import gc
import array
import scipy.io
import numpy as np
from scipy import spatial

import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

#from turb_postproc import read3DNSorder

def smallCubeSpec(u,v,w,boxL):

    #Takes inlet turbulence box and resamples onto cube the size of smallest dimension of box
    #to calculate 3D spectra of turbulence box.
    #u,v,w: velocity components (arrays of floats of size ni,nj,nk)
    #boxL: array of three floats, holding box dimensions in real space

    lx = boxL[0]
    ly = boxL[1]
    lz = boxL[2]

    lenv = np.array([lx,ly,lz])
    #minimum index
    minidx = np.argmin(lenv)
    #minimum edge length
    minl = lenv[minidx]
    #print('minimum edge length ',minl)

    #number of points in minimum edge length direction
    n = np.shape(u)[minidx]


    [uc, vc, wc] = doRegInterp(u,v,w,boxL,[minl,minl,minl],[n,n,n])

    # 3D spec
    tke_spectrum, wave_numbers, cutoff = calcSpec(uc, vc, wc, n, minl)
    
    plt.figure()
    plt.loglog(wave_numbers[1:cutoff[0][0]],tke_spectrum[1:cutoff[0][0]])
    plt.show()

    #write to file if desired
    scipy.io.savemat('3dspec-ax-box-c.mat',dict(E=tke_spectrum[0:cutoff[0][0]],k=wave_numbers[0:cutoff[0][0]]))
    
    return wave_numbers[1:cutoff[0][0]],tke_spectrum[1:cutoff[0][0]]

def floBox(grid,flo,boxPnts,nOff,dOff,minl):
    
    #This function takes a single chunk of flow from a given file, interpolates it on a 
    #regular cubic grid, and then calculates its 3D spectra.
    #Use this on your flo and nod files that are 3DNS outputs. 
    # You will need to read them in and concatenate them correctly yourself, before feeding them to this function.

    #grid: x,y,z: gridpoints in physical space
    #flo: u,v,w: velocity fluctuation components
    #boxPnts: ni,nj,nk: number of gridpoints in flow file
    #nOff: niO,njO,nkO: gridpoint offset 
    #dOff: real space offset of the cube from the chunk. This is a fudge factor. When first using this function, go to doBaryInterp, and uncomment the 
    # figure command, to calibrate your dOff values.
    #minl: minimum length in domain, from which (and of which size) cube is taken (this code defaults to assuming it's the span)
    
    u = flo[:,:,:,0] 
    v = flo[:,:,:,1]
    w = flo[:,:,:,2]
    
    x = grid[:,:,:,0] 
    y = grid[:,:,:,1]
    z = grid[:,:,:,2]
    
    ni = boxPnts[0] 
    nj = boxPnts[1]
    nk = boxPnts[2] 
    
    #offset
    [nio,njo] = nOff
    
    [offx,offy] = dOff

    [uc, vc, wc ] = doBaryInterp(u,v,w,x,y,z,ni,nj,nk,nio,njo,minl,offx,offy)
    
    print('Calculating spectra')
    
    # 3D spec
    tke_spectrum, wave_numbers, cutoff = calcSpec(uc, vc, wc, n, minl)

    plt.figure()
    plt.loglog(wave_numbers[1:cutoff[0][0]],tke_spectrum[1:cutoff[0][0]])
    plt.show()

    #scipy.io.savemat('../../DNS/3dspec-flo-1-6MIt2.mat',dict(E=tke_spectrum[0:cutoff[0][0]],k=wave_numbers[0:cutoff[0][0]]))

    return tke_spectrum[0:cutoff[0][0]], wave_numbers[0:cutoff[0][0]]

def floBoxN(grid,flo,boxPnts,nBoxes,nOff,dOff,addPnt,minl,wdw):
    
    #This function takes several chunks of flow from a given file, interpolates them on a 
    #regular cubic grid, and then calculates their 3D spectra, and averages them.

    #It works in the same way as floBox (so if you're confused about this, use floBox first, then move onto here.)

    #grid: x,y,z: gridpoints in physical space
    #flo: u,v,w: velocity fluctuation components
    #boxPnts: ni,nj,nk: number of gridpoints taken up by chunk in flow file
    #nBoxes: number of flow boxes taken.
    #nOff: nio,njo,nko: vectors holding (number-of-) gridpoint offsets to be taken within the flow chunks (i.e. defining the bottom
    #left-hand corner for each chunk)
    #dOff: real space offset of the cube from the chunk. This is a fudge factor. When first using this function, go to doBaryInterp, and uncomment the 
    # figure command, to calibrate your dOff values.
    #addPnts: any extra points that need to be added to nj in specific cases where the grid gets very skewed(basically another fudge factor)
    #minl: minimum length in domain, from which (and of which size) cube is taken (the code assumes this is the span)
    #wdw: bool, indicating whether to use a 3D Hamming window (True or False)

    #span = 0.009465
    #ni = 101#101#  100 150 
    #nj = 101 #101 #1000 150
    #nk = 96  #144 #
    #offset
    #nio = 96  #64#122 #0 #64
    #number j-offset
    #njoCb = [31,128,225] # [57,230,403] #
    #offxCb = [0.0005,0.0015,0.0020]
    #offyCb = [0.00075,0.00075,0.00075]
    #addPnt = [0,0,15] #[0,20,20] #[20,20,20]#
    #njo1 = 128
    #njo2 = 225

    x = grid[:,:,:,0]
    y = grid[:,:,:,1]
    z = grid[:,:,:,2]
    
    u = flo[:,:,:,0]
    v = flo[:,:,:,1]
    w = flo[:,:,:,2]

    [ni,nj,nk] = boxPnts
    
    [nio,njo] = nOff
    
    [offxCb,offxCb,offxCb] = dOff

    TKE = 0
    n = nk 

    if wdw == True:
        #Putting together a 3D Hamming window
        hammx = np.zeros((n,1,1))
        hammx[0:n,0,0] = np.hamming(n)
        hammx = np.tile(hammx,[1,n,n])
        hammy = np.zeros((1,n,1))
        hammy[0,0:n,0] = np.hamming(n)
        hammy = np.tile(hammy,[n,1,n])
        hammz = np.zeros((1,1,n))
        hammz[0,0,0:n] = np.hamming(n)
        hammz = np.tile(hammz,[n,n,1])
        hamm3D = hammx*hammy*hammz

    urms = 0
    vrms = 0
    wrms = 0


    for cb in range (0,nBoxes):        
        
        nio_loop = nio[cb]
        njo_loop = njo[cb]
        offx_loop = offx[cb]
        offy_loop = offy[cb]
        #adding extra range in the j-direction as grid can get very skewed
        nj = nj+addPnt[cb]
        
        #interpolate
        uc,vc,wc = doBaryInterp(u,v,w,x,y,z,ni,nj,nk,nio_loop,njo_loop,minl,offx_loop,offy_loop)

        #scipy.io.savemat('../../DNS/3dspec-flo-wd-1-6MIt.mat',dict(w=wc,wMeanOff=wc2))
        nanidxuc = np.argwhere(np.isnan(uc)) #Again, they're either all nan, or none nan
        print('Nans at ',nanidxuc)

        urms = urms+np.mean(uc.flatten()**2)
        vrms = vrms+np.mean(vc.flatten()**2)
        wrms = wrms+np.mean(wc.flatten()**2)


        print('Calculating spectra')

        #window where applicable
        if wdw == True:
            uc = uc*hamm3D
            vc = vc*hamm3D
            wc = wc*hamm3D
            
        # Calculate 3D spec
        tke_spectrum, wave_numbers, cutoff = calcSpec(uc, vc, wc, n, minl)
        TKE = TKE + tke_spectrum

        

    urms = urms/nBoxes
    vrms = vrms/nBoxes
    wrms = wrms/nBoxes
    urms2 = np.sqrt(1/3*(urms+vrms+wrms))
    print (urms2)

    #Averaging TKE spectrum.
    TKE = TKE/nBoxes
    #n.b. wavenumbers will be the same for each of the n sub-cubes, as they depend on the sampling of the sub-cube,
    #which in turn depends on the sampling in the span-wise direction, which should be constant for all parts of the flow.
    
    #plt.figure()
    #plt.loglog(wave_numbers[1:cutoff[0][0]],TKE[1:cutoff[0][0]])#tke_spectrum[1:cutoff[0][0]])
    #plt.show()

    #scipy.io.savemat('../../DNS/3dspec-flo-1-6MIt2Test-window.mat',dict(E=tke_spectrum[0:cutoff[0][0]],k=wave_numbers[0:cutoff[0][0]]))
    #scipy.io.savemat('../../DNS/3dspec-flo-hfRe-tf2.mat',dict(E=tke_spectrum[0:cutoff[0][0]],k=wave_numbers[0:cutoff[0][0]]))

    return tke_spectrum[0:cutoff[0][0]], wave_numbers[0:cutoff[0][0]]

def calcSpec(u,v,w,n,L):

    #You might not need to call this function directly.

    #calculate 3D spectra on a uniform cubic grid.
    #u,v,w: velocity components.
    #n: number of grid points in cube (per edge)
    #L: length of cube edges in real space

    #Code to calculate the power spectrum adapted from McDermott, R. J. (2005), “Toward One-Dimensional Turbulence Subgrid
    #Closure for Large-Eddy Simulation” 
    # https://github.com/firemodels/fds-smv    
    
    nt = n**3 

    #Do fft in all three dimensions, and normalise by total number of pints
    uh = np.fft.fftn(u) / nt
    vh = np.fft.fftn(v) / nt
    wh = np.fft.fftn(w) / nt

    #calc power spectral density
    tkeh = 0.5 * (uh * np.conj(uh) + vh * np.conj(vh) + wh * np.conj(wh)).real

    #Find minimum wave number.
    k0x = 2.0 * np.pi / L 
    k0y = 2.0 * np.pi / L 
    k0z = 2.0 * np.pi / L 

    knorm = (k0x + k0y + k0z) / 3.0
    #As this only really works on a cube, knorm = k0x=k0y=k0z.
    print('knorm = ', knorm)

    #kxmax = nx / 2
    #kymax = ny / 2
    #kzmax = nz / 2

    # dk = (knorm - kmax)/n
    # wn = knorm + 0.5 * dk + arange(0, nmodes) * dk

    wave_numbers = knorm * np.arange(0, n)

    tke_spectrum = np.zeros(len(wave_numbers))

    for kx in range(-n//2, n//2-1):
        for ky in range(-n//2, n//2-1):
            for kz in range(-n//2, n//2-1):
                #radial wave number
                rk = np.sqrt(kx**2 + ky**2 + kz**2)
                #round it up to "bin" it for the histogram
                k = int(np.round(rk))
                tke_spectrum[k] += tkeh[kx, ky, kz]
    
    tke_spectrum = tke_spectrum / knorm

    #  tke_spectrum = tke_spectrum[1:]
    #  wave_numbers = wave_numbers[1:]
    #if smooth:
    #    tkespecsmooth = movingaverage(tke_spectrum, 5)  # smooth the spectrum
    #    tkespecsmooth[0:4] = tke_spectrum[0:4]  # get the first 4 values from the original data
    #    tke_spectrum = tkespecsmooth

    knyquist = knorm * n/2 

    cutoff = np.where(wave_numbers>=knyquist)

    return [tke_spectrum, wave_numbers, cutoff]

def doRegInterp(u,v,w,inL,outL,outN):

    #This funcion interpolates from one regular grid onto another.
    #It uses the scipy function RegularGridInterpolator

    nxi = np.shape(u)[0] 
    nyi = np.shape(u)[1] 
    nzi = np.shape(u)[2] 

    #in grid
    xi = np.linspace(0,inL[0],nxi)
    yi = np.linspace(0,inL[1],nyi)
    zi = np.linspace(0,inL[2],nzi)

    #desired out grid
    gxo = np.linspace(0,outL[0],outN[0])
    gyo = np.linspace(0,outL[1],outN[1])
    gzo = np.linspace(0,outL[2],outN[2])

    #desired out indices
    i = np.linspace(0,outN[0]-1,outN[0])
    j = np.linspace(0,outN[1]-1,outN[1])
    k = np.linspace(0,outN[2]-1,outN[2])

    #interpolator functions
    uint = RegularGridInterpolator((xi, yi, zi), u)
    vint = RegularGridInterpolator((xi, yi, zi), v)
    wint = RegularGridInterpolator((xi, yi, zi), w)

    #empty arrays to save out field
    uo = np.zeros([outN[0],outN[1],outN[2]])
    vo = np.zeros([outN[0],outN[1],outN[2]])
    wo = np.zeros([outN[0],outN[1],outN[2]])

    #numpy 3D meshgrid to hold the out grid
    #for cube xo, yo and zo hold exactly the same values, but it's less confusing this way
    xo, yo, zo = np.meshgrid(gxo, gyo, gzo, indexing='ij')

    #numpy 3D meshgrid containing the cube grid indices
    ic, jc, kc = np.meshgrid(i, j, k, indexing='ij')

    #flattened array containing all cube grid point indicies
    ijk = np.transpose(np.array([ic.flatten(),jc.flatten(),kc.flatten()],dtype=np.int32))

    #vector to index the ijk array
    idx = np.linspace(0,outN[0]*outN[1]*outN[2]-1,num=outN[0]*outN[1]*outN[2],dtype=np.int32)

    #This is a vectorisation of the loop below. The code is mess, but the speed-up is ridiculous.
    #points at which we want to find the field (spatial position within the cube)
    pnts = np.array([xo[ijk[idx,0],ijk[idx,1],ijk[idx,2]], \
        yo[ijk[idx,0],ijk[idx,1],ijk[idx,2]],zo[ijk[idx,0],ijk[idx,1],ijk[idx,2]]])

    #for i in range(0,n):
    #    for j in range (0,n):
    #        for k in range (0,n):
    #            pnts = [xc[i,j,k],yc[i,j,k],zc[i,j,k]]
    #            uc[i,j,k] = ui(pnts)
    #            vc[i,j,k] = vi(pnts)
    #            wc[i,j,k] = wi(pnts)

    #print(np.shape(pnts))
    
    pnts = np.transpose(pnts)
    uo[ijk[idx,0],ijk[idx,1],ijk[idx,2]] = uint(pnts)
    vo[ijk[idx,0],ijk[idx,1],ijk[idx,2]] = vint(pnts)
    wo[ijk[idx,0],ijk[idx,1],ijk[idx,2]] = wint(pnts)


    return uo,vo,wo

def doBaryInterp(u,v,w,x,y,z,ni,nj,nk,nio,njo,minl,offx,offy):

    from scipy import spatial

    #This function takes a ni,nj,nk-large chunk out of a given flow file, and does Barycentric interpolation in 3D
    #onto a regular cubic grid, to then calculate spectra. 
    #It assumes a structured grid in x,y, and a uniform structured grid in z,
    #in line with Andy's axial cascade case grids.
    #It also assumes that the smallest number of points is in the z (k) direction.
    
    #u,v,w; x,y,z velocity components / grid points in real space of flow fed in
    #ni,nj,nk: size of the chunk of flow being made into a cube (as grid point integers). Assumes nk<=ni,nj.
    #nio, njo grid point offset at which bottom left hand corner of chunk is to be positioned
    #minl: lenght of minimum dimension of flow (in real space). This code assumes that it is the span.
    #offx, offy: offsets in x aand y (real space), at which the cube grid will be taken. These might or might not be needed,
    #and suitable values potentially need some trial and error. Useful for curvilinear starting grids.
    
    n = nk

    #to be flattened version of x,y,z,u,v,w
    #the bit of the grid being taken out is ni, nj points wide.
    #offsetted into the i and j directions by nio, njo.
    xf = x[nio:ni+nio,njo:nj+njo,:]
    yf = y[nio:ni+nio,njo:nj+njo,:]
    zf = z[nio:ni+nio,njo:nj+njo,:]
    uf = u[nio:ni+nio,njo:nj+njo,:]
    vf = v[nio:ni+nio,njo:nj+njo,:]
    wf = w[nio:ni+nio,njo:nj+njo,:]

    #offsetting the grid in case its minimum x or y are negative
    #This caused problems with the barycentric interpolation beforehand.
    xf = xf+np.abs(np.min(x[nio:ni+nio,njo:nj+njo,:].flatten()))
    yf = yf+np.abs(np.min(y[nio:ni+nio,njo:nj+njo,:].flatten()))

    #Placeholder cube velocity components
    uc = np.zeros([n,n,n])
    vc = np.zeros([n,n,n])
    wc = np.zeros([n,n,n])

    #x-offset and y-offset to make sure the cubic grid points lie within the
    #offsetted grid. This is hack-y but it works. Value might have to be adjusted.
    #offx = 0.0005
    #offy = 0.00075
    xo = np.abs(np.min(xf.flatten()))+offx
    yo = np.abs(np.min(yf.flatten()))+offy

    #desired cube grid
    gcx = np.linspace(xo,minl+xo,n)
    gcy = np.linspace(yo,minl+yo,n)
    gcz = np.linspace(0,minl,n)

    #desired cube indices
    i = np.linspace(0,n-1,n)
    j = np.linspace(0,n-1,n)
    k = np.linspace(0,n-1,n)

    #numpy 3D meshgrid to hold the cube grid
    #xc, yc and zc hold exactly the same values, but it's less confusing this way
    xc, yc, zc = np.meshgrid(gcx, gcy, gcz, indexing='ij')

    #Uncomment the figure command below and check that all the blue dots are within the "frame" made by the orange x's,
    #to check that your offx, offy are okay.
    '''
    plt.figure()
    plt.scatter(xc[:,:,0], yc[:,:,0], marker='o')
    plt.scatter(xf[:,:,0], yf[:,:,0], marker='x')
    plt.show()
    '''
    
    #numpy 3D meshgrid containing the cube grid indices
    ic, jc, kc = np.meshgrid(i, j, k, indexing='ij')

    #flattened array containing all cube grid point indicies
    ijk = np.transpose(np.array([ic.flatten(),jc.flatten(),kc.flatten()],dtype=np.int32))

    #vector to index the ijk array
    idx = np.linspace(0,n*n*n-1,n*n*n,dtype=np.int32)
    print('Interpolating')

    #This is a vectorisation of the query point assignment. The code is mess, but the speed-up is ridiculous.
    #points at which we want to find the field (spatial position within the cube)
    Qpnts = np.array([xc[ijk[idx,0],ijk[idx,1],ijk[idx,2]], \
        yc[ijk[idx,0],ijk[idx,1],ijk[idx,2]], zc[ijk[idx,0],ijk[idx,1],ijk[idx,2]]])
    Qpnts = np.transpose(Qpnts)

    #use kdTree to find closest point.
    #constructing a KDTree 
    #construct index arrays to find index of point closest to query
    xidx = np.zeros((ni,1,1))
    xidx[0:ni,0,0] = np.linspace(0,ni-1,ni,dtype= np.int)
    xidx = np.tile(xidx,[1,nj,nk])
    yidx = np.zeros((1,nj,1))
    yidx[0,0:nj,0] = np.linspace(0,nj-1,nj,dtype= np.int)
    yidx = np.tile(yidx,[ni,1,nk])
    zidx = np.zeros((1,1,nk))
    zidx[0,0,0:nk] = np.linspace(0,nk-1,nk,dtype= np.int)
    zidx = np.tile(zidx,[ni,nj,1])

    #index array
    Idx = np.transpose(np.array([np.reshape(xidx,ni*nj*nk), \
        np.reshape(yidx,ni*nj*nk),np.reshape(zidx,ni*nj*nk)]))

    #construct KDTree
    Grd = np.transpose(np.array([np.reshape(xf,ni*nj*nk), \
        np.reshape(yf,ni*nj*nk),np.reshape(zf,ni*nj*nk)]))

    print('Constructing KDTree')
    tree = spatial.KDTree(Grd)
    print('Querying')

    fracTol = 1e-4
    span = minl 
    zSpac = span/(nk-1)

    #find distance to closest grid point and its idx in flattened grid (distance at clst[0], index at clst[1])
    clst = np.array(tree.query(Qpnts,k=1,eps=0))
    clst = np.transpose(clst)
    #find its idx in 3D grid
    clstIdx = np.array(Idx[clst[:,1].astype(int)].astype(int))

    #find point at this idx in 3D grid 
    clstPnt = np.array([xf[clstIdx[:,0],clstIdx[:,1],clstIdx[:,2]], \
        yf[clstIdx[:,0],clstIdx[:,1],clstIdx[:,2]],zf[clstIdx[:,0],clstIdx[:,1],clstIdx[:,2]]])
    clstPnt = np.transpose(clstPnt)
    
    #calculate factors for linear interpolation in regular direction along grid (z-direction)
    #These should usually all be equal to zero, as we expect no interpolation in the regular z-Direction
    zDist = Qpnts[:,2]-clstPnt[:,2]

    #interpolate in 3D using barycentric interpolation in plane (x,y), and linear in z-direction
    #Usually the box should be aligned with the cube in z, so all zDist should be zero.
    for i in range(0,n**3):

        if zDist[i] != 0:
            zDir = np.sign(zDist[i])
            zFact = np.abs(zDist[i]/zSpac)
        else:
            zDir = 0
            zFact = 1
        
        [ud,vd,wd] = interpBar3D(Qpnts[i],clstIdx[i],xf[:,:,0],yf[:,:,0],zDir,zFact,[uf,vf,wf],fracTol) 

        #test whether returned as not in any triangle
        if np.isnan(ud): # only testing for one of them vd or wd == np.nan:
            #print('Queried points ' ,Qpnts[i])
            itmp = clstIdx[i][0]
            jtmp = clstIdx[i][1]
            ktmp = clstIdx[i][2]
            found = False
            #iterate over all eight neighbours in x-y (this is done if the point in question is not within
            # any of the four surrounding triangles, which can be the case in very skewed grids.)
            for ii in range(-1,2):
                for jj in range (-1,2):
                    #print('ii,jj ',ii,jj)
                    [ud,vd,wd] = interpBar3D(Qpnts[i],[itmp+ii,jtmp+jj,ktmp],xf[:,:,0],yf[:,:,0],zDir,zFact,[uf,vf,wf],fracTol)
                    if not np.isnan(ud): #and vd and wd #(only need to check one, as they are either all nan or none nan)
                        uc[ijk[i,0],ijk[i,1],ijk[i,2]] = ud
                        vc[ijk[i,0],ijk[i,1],ijk[i,2]] = vd
                        wc[ijk[i,0],ijk[i,1],ijk[i,2]] = wd
                        found = True
                        break 
                #Use "found" flag to break out of multiple loops at once
                if found == True:
                    break
            if found == False:
                raise Exception("Point doesn't seem to be in any triangles close by")
        else:
            uc[ijk[i,0],ijk[i,1],ijk[i,2]] = ud
            vc[ijk[i,0],ijk[i,1],ijk[i,2]] = vd
            wc[ijk[i,0],ijk[i,1],ijk[i,2]] = wd


    #nanidxuc = np.argwhere(np.isnan(uc)) #Again, they're either all nan, or none nan
    #print('Nans at ',nanidxuc)

    return uc, vc, wc

def interpBar3D(spacCoords,clstIdx,xx,yy,zDir,zFact,tIntrp,fracTol):

    #barycentric interpolation of tIntrp onto a new point found at spacCoords.
    #the point closest to the new one was found using a kdTree. 
    #In this function the two other points which make up the triangle around the point of interest are found, and
    #barycentric interpolation is performed.
    #Similar to interpBar2D, but in this case taking in a variable which tells us in which direction to finish the linear interpolation in 3D.
    #(But interpBar2D isn't included in this pack...). Also this isn't always stable, so might need the occasional fudging.    
    #get the four grd pnts around and repeat the first one to complete the "circle"
    barPnts = np.tile(clstIdx,[5,1])
    barPnts = barPnts+np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0],[1,0,0]])
    #closest point is the centre (vertex 0)
    x0 = xx[clstIdx[0],clstIdx[1]]
    y0 = yy[clstIdx[0],clstIdx[1]]
    x = spacCoords[0]
    y = spacCoords[1]
    lmbd = np.zeros([4,3])
    sln = np.zeros([4])

    ni,nj = np.shape(xx)

    intrpOut = np.zeros([len(tIntrp)])
    
    #calculate the barycentric coefficients for each triangle
    for i in range(0,4):
        
        #if statement to skip edge points which are outside of the grid.
        if (0<=barPnts[i,0]<ni) and (0<=barPnts[i,1]<nj) and (0<=barPnts[i+1,0]<ni) and (0<=barPnts[i+1,1]<nj):
            #vertex 1
            x1 = xx[barPnts[i,0],barPnts[i,1]]
            y1 = yy[barPnts[i,0],barPnts[i,1]]
            #vertex 2
            x2 = xx[barPnts[i+1,0],barPnts[i+1,1]]
            y2 = yy[barPnts[i+1,0],barPnts[i+1,1]]
        else:
            continue
        #third vertex is always the centrepoint. 
        #set up linear equation to solve for lambdas:
        #matrix
        m = np.array([[x0,x1,x2],[y0,y1,y2],[1,1,1]])
        try:
            #vector
            v = np.array([x,y,1])
            #solve for lambdas
            lmbd[i,:] = np.linalg.solve(m,v)
        except np.linalg.LinAlgError:
            print(m)   
            raise
        #check whether lambdas are all +ve and sum to one, in which case the point lies within the triangle
        #fracTol in there bc if point is on a line, one of the three vales will be zero, and if it lies on a point it'll return two zeros
        #slightly negative values will push it over a line or over the point, and so checking they are above a certian tiny negative threshold extends the capture of a triangle.
        #this can be helpful when the point lies on a line (edge case) - otherwise numerical precision errors might lead to the program not finding a point
        #when it lies on a line
        if all(lmbd[i,:]>= -fracTol):
            sln[i] = 1
    tst = np.sum(sln)
    if tst == 0:
      #print(clstIdx[0],clstIdx[1])  
      #print(x,y)  
      #print(lmbd)
      #print('Projection not in either triangle, exiting')
      return np.nan*np.ones((len(tIntrp))) 
    elif tst == 1:
      #print('Single triangle')
      #print(lmbd)
      #find the non-zero line of sln
      solTri = np.nonzero(sln)[0][0]
      #lambdas of the line in solTri corresponding to the triangle
      #print(np.array(lmbd[solTri,:]))
      for l in range(0,len(tIntrp)):
        #find values of tIntrp on vertices to interpolate with
        intrpVx0 = tIntrp[l][clstIdx[0],clstIdx[1],clstIdx[2]]
        intrpVx1 = tIntrp[l][barPnts[solTri,0],barPnts[solTri,1],barPnts[solTri,2]]
        intrpVx2 = tIntrp[l][barPnts[solTri+1,0],barPnts[solTri+1,1],barPnts[solTri+1,2]]
        #find values of tIntrp on vertices to interpolate with (z-shift)
        intrpVx0z = tIntrp[l][clstIdx[0],clstIdx[1],np.int32(clstIdx[2]+zDir)]
        intrpVx1z = tIntrp[l][barPnts[solTri,0],barPnts[solTri,1],np.int32(barPnts[solTri,2]+zDir)]
        intrpVx2z = tIntrp[l][barPnts[solTri+1,0],barPnts[solTri+1,1],np.int32(barPnts[solTri+1,2]+zDir)]
        #interpolating between the centre, and the vertecies, which will be the points corresponding to barPnts of solTri and solTri+1
        ur1  = lmbd[solTri,:]*[intrpVx0,intrpVx1,intrpVx2]
        ur2  = lmbd[solTri,:]*[intrpVx0z,intrpVx1z,intrpVx2z]
        #print(lmbd[solTri,:],[intrpVx0,intrpVx1,intrpVx2])
        #print(lmbd[solTri,:],[intrpVx0z,intrpVx1z,intrpVx2z])
        #linear interpolation in z-direction
        intrpOut[l] = np.sum((1-zFact)*ur1+zFact*ur2)
      return intrpOut
    elif tst == 2:
      #print('On triangle boundary, defaulting to first triangle')
      solTri = np.nonzero(sln)[0][0]
      #print('sol Triangle number', solTri) 
      #print ('lmd', lmbd)
      #print(np.array(lmbd[solTri,:]))
      for l in range(0,len(tIntrp)):
        intrpVx0 = tIntrp[l][clstIdx[0],clstIdx[1],clstIdx[2]]
        intrpVx1 = tIntrp[l][barPnts[solTri,0],barPnts[solTri,1],barPnts[solTri,2]]
        intrpVx2 = tIntrp[l][barPnts[solTri+1,0],barPnts[solTri+1,1],barPnts[solTri+1,2]]
        #find values of tIntrp on vertices to interpolate with (z-shift)
        intrpVx0z = tIntrp[l][clstIdx[0],clstIdx[1],np.int32(clstIdx[2]+zDir)]
        intrpVx1z = tIntrp[l][barPnts[solTri,0],barPnts[solTri,1],np.int32(barPnts[solTri,2]+zDir)]
        intrpVx2z = tIntrp[l][barPnts[solTri+1,0],barPnts[solTri+1,1],np.int32(barPnts[solTri+1,2]+zDir)]
        #interpolating between the centre, and the vertecies, which will be the points corresponding to barPnts of solTri and solTri+1
        ur1  = lmbd[solTri,:]*[intrpVx0,intrpVx1,intrpVx2]
        ur2  = lmbd[solTri,:]*[intrpVx0z,intrpVx1z,intrpVx2z]
        #print(lmbd[solTri,:],[intrpVx0,intrpVx1,intrpVx2])
        #print(lmbd[solTri,:],[intrpVx0z,intrpVx1z,intrpVx2z])
        #linear interpolation in z-direction
        intrpOut[l] = np.sum((1-zFact)*ur1+zFact*ur2)
        #print('interpolated value', intrpOut[l])
      return intrpOut
    elif tst == 4:
      #print('On centre')
      for l in range(0,len(tIntrp)):
        #in this case only the closest point (vertex 0) and its neighbour in the z-direction are needed
        intrpVx0 = tIntrp[l][clstIdx[0],clstIdx[1],clstIdx[2]]
        #find values of tIntrp on vertices to interpolate with (z-shift)
        intrpVx0z = tIntrp[l][clstIdx[0],clstIdx[1],np.int32(clstIdx[2]+zDir)]
        #no in-plane interpolation needed
        ur1  = intrpVx0
        ur2  = intrpVx0z
        #print(intrpVx0, intrpVx0z)
        #linear interpolation in z-direction
        intrpOut[l] = np.sum((1-zFact)*ur1+zFact*ur2)
      return intrpOut
    else:
        print(tst)
        print('Something else went wrong')
        return     

def doRegInterp2(u,v,w,xi,yi,zi,xo,yo,zo):

    #This funcion interpolates from one regular grid onto another.
    #It uses the scipy function RegularGridInterpolator
    # Does the same thing as doRegInterp, just takes the input arguments explicitly

    #print(np.min(np.min(np.min(xi))),np.max(np.max(np.max(xi))),np.min(np.min(np.min(yi))),np.max(np.max(np.max(yi))),np.min(np.min(np.min(zi))),np.max(np.max(np.max(zi))))
    #print(np.min(np.min(np.min(xo))),np.max(np.max(np.max(xo))),np.min(np.min(np.min(yo))),np.max(np.max(np.max(yo))),np.min(np.min(np.min(zo))),np.max(np.max(np.max(zo))))
    #print(np.min(np.min(xo[-1,:,:])))

    [ni,nj,nk]=np.shape(xo)
    #desired out indices
    i = np.linspace(0,ni-1,ni)
    j = np.linspace(0,nj-1,nj)
    k = np.linspace(0,nk-1,nk)

    #interpolator functions
    uint = RegularGridInterpolator((xi[:,0,0], yi[0,:,0], zi[0,0,:]), u)
    vint = RegularGridInterpolator((xi[:,0,0], yi[0,:,0], zi[0,0,:]), v)
    wint = RegularGridInterpolator((xi[:,0,0], yi[0,:,0], zi[0,0,:]), w)

    #empty arrays to save out field
    uo = np.zeros([ni,nj,nk])
    vo = np.zeros([ni,nj,nk])
    wo = np.zeros([ni,nj,nk])

    #numpy 3D meshgrid to hold the out grid
    #for cube xo, yo and zo hold exactly the same values, but it's less confusing this way
    #xo, yo, zo = np.meshgrid(gxo, gyo, gzo, indexing='ij')

    #numpy 3D meshgrid containing the slice grid indices
    ic, jc, kc = np.meshgrid(i, j, k, indexing='ij')

    #flattened array containing all slice grid point indicies
    ijk = np.transpose(np.array([ic.flatten(),jc.flatten(),kc.flatten()],dtype=np.int32))

    #vector to index the ijk array
    idx = np.linspace(0,ni*nj*nk-1,num=ni*nj*nk,dtype=np.int32)

    #This is a vectorisation of the loop below. The code is mess, but the speed-up is ridiculous.
    #points at which we want to find the field (spatial position within the cube)
    pnts = np.array([xo[ijk[idx,0],ijk[idx,1],ijk[idx,2]], \
        yo[ijk[idx,0],ijk[idx,1],ijk[idx,2]],zo[ijk[idx,0],ijk[idx,1],ijk[idx,2]]])

    #for i in range(0,n):
    #    for j in range (0,n):
    #        for k in range (0,n):
    #            pnts = [xc[i,j,k],yc[i,j,k],zc[i,j,k]]
    #            uc[i,j,k] = ui(pnts)
    #            vc[i,j,k] = vi(pnts)
    #            wc[i,j,k] = wi(pnts)

    #print(np.shape(pnts))
    
    pnts = np.transpose(pnts)
    uo[ijk[idx,0],ijk[idx,1],ijk[idx,2]] = uint(pnts)
    vo[ijk[idx,0],ijk[idx,1],ijk[idx,2]] = vint(pnts)
    wo[ijk[idx,0],ijk[idx,1],ijk[idx,2]] = wint(pnts)


    return uo,vo,wo


if __name__ == '__main__':

    fC = '../Data/turbin_grid_ax_c.dat'#trbTst2.dat' #nC = '../../DNS/in_dns_Re2.dat'

    c = array.array('f')
    c.fromfile(open(fC, 'rb'), os.path.getsize(fC) // c.itemsize)

    C = np.frombuffer(c, dtype=np.float64)  
    C = np.reshape(C,(584,415,96,3))

    u = C[:,:,:,0]
    v = C[:,:,:,1]
    w = C[:,:,:,2]

    u = np.swapaxes(u,1,2)
    v = np.swapaxes(v,1,2)
    w = np.swapaxes(w,1,2)

    tke_spectrum, wave_numbers = smallCubeSpec(u, v, w, [0.081618141935311,0.0408,0.009533835077005])
    


