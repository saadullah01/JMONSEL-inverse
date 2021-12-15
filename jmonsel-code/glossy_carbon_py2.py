# This script determines the yield for electrons incident on 1 or more
# trapezoidal (with top corner radii) resist lines on a 3-layer substrate.

# The script is for understanding the distribution of
# the number of SE per incident electron. This script doesn't add randomness of
# source particle therefore no source shot noise included.
# The number of sub-acquisitions n = 100 (or in this script, acquisitionNum = 100).

import gov.nist.microanalysis.EPQLibrary as epq
#from gov.nist.microanalysis import EPQLibrary as epq
import gov.nist.microanalysis.EPQTools as ept
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.nanoscalemetrology.JMONSEL as mon
import gov.nist.microanalysis.Utility as nmu
import java.io as jio
import java.util as jutil
import java.lang as jl
import jarray
import java.nio.charset as cs
import random
#import numpy as np

hnm = 15.  #resist line height in nm, 0 means that there is only substrace
# determine where to save the results
dest = r"C:\Users\vaibh\Documents\503-project\JMONSEL\JMONSEL\Examples\glC_dataset";
jio.File(dest).mkdirs()
# filename = dest+PathSep+"LinesOnLayers_results_100acqusition_Height"+str(int(hnm))+"SameMaterial1.txt"
# file = open(filename,'w')
#print ("Output will be to file: " + filename)

# Make a record of the random seed we use so we can exactly repeat this calculation
# (same random number sequence) if necessary.
seed = nmu.Math2.rgen.nextLong() # Pick a random seed

# To exactly repeat a previous calculation (e.g., for bug fix) uncomment the next line
# and replace the question marks with the seed that was recorded in the previous calculation's
# output file.

#seed = -1411824173088723636L
nmu.Math2.initializeRandom(seed)
#print >>file,"Random number seed: ",seed
#print "Random number seed: ",seed
for i in range(0,10):
    r = nmu.Math2.rgen.nextDouble()
#    print >>file, r
#    print r

# Model parameters
acquisitionNum = 1
nTrajectories = 1   #5000 is the default value

# Shape parameters.
pitchnm = 5. 	# Distance between line centers in nm
nlines = 1.  	#number of lines
#hnm = 10.  #resist line height in nm, 0 means that there is only substrace
wnm = 4.	#resist line bottom width in nm
linelengthnm = 4. # resist line length in nm
# Note that sidewall angles are specified with respect to vertical,
# so 0. is vertical, positive angles have bottom wider than top, and
# negative angles are the reverse (undercut).
thetardeg = 0.	#resist line right sidewall angle in degrees
thetaldeg = 0.	#resist line left sidewall angle in degrees
radrnm = 0.		#resist line top right corner radius in nm
radlnm = 0.		#resist line top left corner radius in nm
layer1thicknessnm = 80.	# Thickness in nm of the 1st layer (immediately below the lines)
layer2thicknessnm = 200.# Thickness in nm of the 2nd layer
# We'll make the substrate infinitely thick.
beamEeVvals = [500.] # Beam energies in eV
beamsizenm = 0.5	# beam size in nm
# The following parameter needs a bit of explaining. In the model we'll build below, the infinitely deep
# layer 3 will be artificially divided into two parts, a skin layer that is close to the surface and a
# deep part that is farther from the surface. The definition of "deep" is set by
# the deepnm parameter on the next line. Both the skin region and the deep region will contain the same
# material (Si), but in the deep region we'll make a model in which electrons with energies less than
# 50 eV are dropped from the simulation. This can save lots of time (particularly if beam energies are large)
# because there are lots of secondary electrons with energies < 50 eV, and lots of simulation time
# must be devoted to tracking them. We can't afford to drop them when they are generated near the surface,
# because they might escape and be detected. I.e., they're important there. However, low energy electrons
# that are deep inside the sample can't escape, so there is no harm done in not tracking them. Thus, the
# parameter below should be set to several times the typical escape depth (so there's little chance of
# dropping an electron that would have escaped). This is only important for high beam energies, because
# only then will the electrons have sufficient range to reach the deeper layer, but in such cases there
# can be a significant time savings.
deepnm = 15. # Depth below which to use the "deep model."

trajImg = 1
trajImgMaxTraj = 50
trajImgSize = 100.e-9

VRML = 0
VRMLImgMaxTraj = 0 # Include no trajectories in VRML (Show sample only.) Leaving trajectories
# out makes a VRML that displays easily. It's good for checking the sample. I find that adding
# trajectories significantly slows down the display. If you want to try it, keep the number of
# displayed trajectories small (20 is a reasonable number) and turn off collision detection in
# your VRML viewer.


# Make materials

#	A Secondary Electron vaccum
vacuum = mon.SEmaterial()
vacuum.setName("SE vacuum")
vacuumBarrier = mon.ExpQMBarrierSM(vacuum)
vacuumMSM = mon.MONSEL_MaterialScatterModel(vacuum)
vacuumMSM.setBarrierSM(vacuumBarrier)

# PMMA: Scattering tables for the DFT model of PMMA are not yet available.
# Instead the code below implements a backup model based on the FittedInelSM class,
# supplemented with a charge trapping model. Each of these has two free parameters.
# These are chosen to provide the best fit to measured SE yield vs. energy.
# There is no guarantee that a model so constructed will match topographic yield
# (yield vs. angle of incidence) since all of the data used to determine the
# parameters were measured at normal incidence. Nevertheless, it's the best we
# can do for now.

breakE = epq.ToSI.eV(45.)
density = 1190.
workfun = 5.5
bandgap = 5. # width of band gap in eV, based on TPP. Svorcik, Lyutakov, Huttel get about 5.3
EFermi = -bandgap # This puts the Fermi level at the top of the valence band.
potU = -workfun-EFermi # Gives 0.5 for this material for now. This is based on Sayyah et al.,
# Int. J. Polymeric Mat. 54, p 505 (2005). It's about mid-range for the values they report for
# different forms of PMMA.

# Material defined in terms of its constituent elements and their weight fractions
# Elemental Constituents
C = epq.Element.C
Ox = epq.Element.O
H = epq.Element.H
PMMAcomp = epq.Composition()
PMMAcomp.defineByMoleFraction([C,Ox,H],[5,2,8])
PMMA = mon.SEmaterial(PMMAcomp,density)
PMMA.setName("PMMA")
PMMAWorkfunction=epq.ToSI.eV(workfun)
PMMA.setWorkfunction(PMMAWorkfunction)
PMMA.setBandgap(epq.ToSI.eV(bandgap))
PMMA.setEnergyCBbottom(epq.ToSI.eV(potU))

# Create scatter mechanisms
PMMANISTMott = mon.SelectableElasticSM(PMMA,mon.NISTMottRS.Factory)
PMMACSD = mon.JoyLuoNieminenCSD(PMMA,breakE)
PMMAfittedInel = mon.FittedInelSM(PMMA,epq.ToSI.eV(65.4),PMMACSD)
# Parameters in the next line are from my fits. (See PMMAOptimization.nb)
PMMApolaron = mon.GanachaudMokraniPolaronTrapSM(2.e7,1./epq.ToSI.eV(4.))
#PMMAClassicalBarrier = mon.ExpQMBarrierSM(PMMA)
# Make a material scatter model
# to be used in thin layer
PMMAMSM = mon.MONSEL_MaterialScatterModel(PMMA)
PMMAMSM.addScatterMechanism(PMMANISTMott)
PMMAMSM.addScatterMechanism(PMMAfittedInel)
PMMAMSM.addScatterMechanism(PMMApolaron)
PMMAMSM.setCSD(PMMACSD)
#PMMAMSM.setBarrierSM(PMMAClassicalBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
PMMAMSMDeep = mon.MONSEL_MaterialScatterModel(PMMA)
PMMAMSMDeep.addScatterMechanism(PMMANISTMott)
PMMAMSMDeep.addScatterMechanism(PMMAfittedInel)
PMMAMSMDeep.addScatterMechanism(PMMApolaron)
PMMAMSMDeep.setCSD(PMMACSD)
#PMMAMSMDeep.setBarrierSM(PMMAClassicalBarrier)
PMMAMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))

# TODO: Generate an ARC model.

# BEGIN TEMPORARY
# Replace the following lines with an ARC model when available.
# I'm replacing the ARC with PMMA during this test phase.
ARCMSM = PMMAMSM
# END TEMPORARY

# #	glassy carbon, set up for TabulatedInelasticSM mode 3 with energy levels as follows:
# #	CB bottom at -25.4 eV relative to vacuum = 0 eV.
# #	Interface barrier is gradual.

density = 1800.
workfun = 5.0
bandgap = 0. # width of band gap in eV
EFermi = 20.4 #
potU = -workfun-EFermi
glC = mon.SEmaterial([epq.Element.C],[1.],density,"glassy Carbon")
glCWorkfunction=epq.ToSI.eV(workfun)
glC.setWorkfunction(glCWorkfunction)
glC.setEnergyCBbottom(epq.ToSI.eV(potU))
glC.setBandgap(epq.ToSI.eV(bandgap))
glC.setCoreEnergy([epq.ToSI.eV(284.2)])
# Edit the string below so it is the path to the folder where you have stored the glassy carbon
# scattering tables that I provide
tablePath = r"C:\Users\vaibh\Documents\ScatteringTables\glassyCTables"+PathSep
glCTables = [tablePath +"IIMFPPennInterpglassyCSI.csv", tablePath +"interpNUSimReducedDeltaEglassyCSI.csv", tablePath +"interpsimTableThetaNUglassyCSI.csv", tablePath +"interpSimESE0NUglassyCSI.csv"]
# Create scatter mechanisms
glCNISTMott = mon.SelectableElasticSM(glC,mon.NISTMottRS.Factory)
glCDS = mon.TabulatedInelasticSM(glC,3,glCTables)
#glCAbruptBarrier = mon.ExpQMBarrierSM(glC,0.)
# Make a material scatter model
# MSM to be used in thin layer (includes SE generation)
glCMSM = mon.MONSEL_MaterialScatterModel(glC)
glCMSM.addScatterMechanism(glCNISTMott)
glCMSM.addScatterMechanism(glCDS)
#glCMSM.setBarrierSM(glCAbruptBarrier)
# MSM to be used deep inside (drops electrons with E<50 eV)
glCMSMDeep = mon.MONSEL_MaterialScatterModel(glC)
glCMSMDeep.addScatterMechanism(glCNISTMott)
glCMSMDeep.addScatterMechanism(glCDS)
#glCMSMDeep.setBarrierSM(glCAbruptBarrier)
glCMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))

# # Si
# phononE = 0.063 # I've seen the number reported as 510 cm^-1. this is conversion of that to eV.
# phononStrength = 3. # Turner & Inkson dispersion curves appear to show 3 LO phonon modes converging to the same
# # energy at the Gamma point.
# density = 2330.
# workfun = 4.85
# bandgap = 1.1 # width of band gap in eV
# EFermi = -bandgap # This puts the Fermi level at the top of the valence band.
# potU = -workfun-EFermi
# Si = mon.SEmaterial([epq.Element.Si],[1.],density,"Silicon")
# SiWorkfunction=epq.ToSI.eV(workfun)
# Si.setWorkfunction(SiWorkfunction)
# Si.setEnergyCBbottom(epq.ToSI.eV(potU))
# Si.setBandgap(epq.ToSI.eV(bandgap))
# Si.setCoreEnergy([epq.ToSI.eV(99.2),epq.ToSI.eV(99.8),epq.ToSI.eV(149.7),epq.ToSI.eV(1839.)])
# # Edit the string below so it is the path to the folder where you have stored the silicon scattering tables
# # that I provide
# #tablePath = "C:\Users\mxpeng\Researches\JMONSEL_distribution\JMONSEL\ScatteringTables\SiTables"+PathSep
# tablePath=r"C:\Users\vaibh\Documents\ScatteringTables\SiTables"+PathSep
# SiTables = [tablePath +"IIMFPFullPennInterpSiSI.csv",
# tablePath +"interpNUSimReducedDeltaEFullPennSiSI.csv",
# tablePath +"interpNUThetaFullPennSiBGSI.csv",
# tablePath +"interpSimESE0NUSiBGSI.csv"]
# # Create scatter mechanisms
# SiNISTMott = mon.SelectableElasticSM(Si,mon.NISTMottRS.Factory)
# SiDS = mon.TabulatedInelasticSM(Si,3,SiTables,epq.ToSI.eV(13.54))
# # The eps0 value is n^2, where n=3.4155 is taken from Palik. epsInfinity is from Palik's 2000 eV value of
# # n = 0.9999048
# Siphonon = mon.GanachaudMokraniPhononInelasticSM(phononStrength,epq.ToSI.eV(phononE),300.,11.7,1.)
# #SiAbruptBarrier = mon.ExpQMBarrierSM(Si,0.)
# # Make a material scatter model
# # MSM to be used in thin layer (includes SE generation)
# SiMSM = mon.MONSEL_MaterialScatterModel(Si)
# SiMSM.addScatterMechanism(SiNISTMott)
# SiMSM.addScatterMechanism(SiDS)
# SiMSM.addScatterMechanism(Siphonon)
# #SiMSM.setBarrierSM(SiAbruptBarrier) # Omitting this line causes the barrier to default to gradual/classical
# # MSM to be used deep inside (drops electrons with E<50 eV)
# SiMSMDeep = mon.MONSEL_MaterialScatterModel(Si)
# SiMSMDeep.addScatterMechanism(SiNISTMott)
# SiMSMDeep.addScatterMechanism(SiDS)
# SiMSMDeep.addScatterMechanism(Siphonon)
# #SiMSMDeep.setBarrierSM(SiAbruptBarrier) # Omitting this line causes the barrier to default to gradual/classical
# SiMSMDeep.setMinEforTracking(epq.ToSI.eV(50.))	#default

# Conversions of shape parameters to SI units.
# Shape parameters.
meterspernm = 1.e-9	# conversion from nanometers to meters
pitch = pitchnm*meterspernm
h = hnm*meterspernm
w = wnm*meterspernm
linelength = linelengthnm*meterspernm
# Note that sidewall angles are specified with respect to vertical,
# so 0. is vertical, positive angles have bottom wider than top, and
# negative angles are the reverse (undercut).
radperdeg = jl.Math.PI/180.	# conversion from degrees to radians
thetar = thetardeg*radperdeg
thetal = thetaldeg*radperdeg
radr = radrnm*meterspernm
radl = radlnm*meterspernm
layer1thickness = layer1thicknessnm*meterspernm
layer2thickness = layer2thicknessnm*meterspernm
beamsize = beamsizenm*meterspernm
deep = deepnm*meterspernm

# create an instance of the model
monte=nm.MonteCarloSS() #creates an instance of the model with all default characteristics
eg = nm.GaussianBeam(beamsize) # makes electron gun, Gaussian with standard deviation = beamsize
monte.setElectronGun(eg) # This gun is "attached" to the model.

# SAMPLE DESCRIPTION

# NISTMonte provides us with a "chamber" in the form of a 0.1 m sphere inside of which we build
# out sample. Replace the default vacuum in the chamber with SEvacuum. (SEmaterials define additional
# properties, such as work function, that are needed by JMONSEL.)
chamber = monte.getChamber()
chamber.updateMaterial(chamber.getScatterModel(),vacuumMSM)

# Generate the sample. The GaussianBeam electron gun has this pecularity: it defines the +z axis to
# be in the direction of travel of the electrons. When we describe the sample in this coordinate
# system, it is inverted along the z direction.

# Make sample component shapes
normalvector = [0.,0.,-1.]

# First we make the layers. The simplest way to do this is to define each as a MultiPlaneShape with a single
# plane, each nested inside the previous one.
layer1 = mon.NormalMultiPlaneShape()
layer1.addPlane(normalvector,[0.,0.,0.]) #layer 1 is now the half space of everything above the x-y plane
# This region has shape defined by layer1, scattering properties defined for ARC, and is a subregion of (is
# wholly contained within) the chamber.
layer1Region = monte.addSubRegion(chamber,ARCMSM,layer1)

layer2 = mon.NormalMultiPlaneShape()
layer2.addPlane(normalvector,[0.,0.,layer1thickness]) #layer 2 starts layer1thickness farther up.
# We give it the properties of carbon, and make it a subregion of layer1Region. At this point, layer1Region
# extends only for 0<=z<=layer1thickness. //Here I've changed the Carbon to be Si.
layer2Region = monte.addSubRegion(layer1Region,glCMSM,layer2)

layer3 = mon.NormalMultiPlaneShape()
layer3.addPlane(normalvector,[0.,0.,layer1thickness+layer2thickness]) #layer 3 starts
# yet another layer2thickness farther up.
# We give it the properties of Si, and make it a subregion of layer2Region. At this point, layer2Region
# extends only for layer1thickness<=z<=layer1thickness+layer2thickness.
layer3Region = monte.addSubRegion(layer2Region,glCMSM,layer3)

layer4 = mon.NormalMultiPlaneShape()
layer4.addPlane(normalvector,[0.,0.,layer1thickness+layer2thickness+deep]) #layer 4 starts even farther
# from the surface.
# We make it a subregion of layer3Region. At this point, layer3Region
# extends only for layer1thickness+layer2thickness<=z<=layer1thickness+layer2thickness+deep.
deepRegion = monte.addSubRegion(layer3Region,glCMSMDeep,layer4)

# Make the array of lines. The integer divide in (nlines/2) truncates fractions.
# The result always has one line centered at (x,y)=(0,0). If nlines is odd
# the remaining lines are placed symmetrically left and right of this one.
# If nlines is even, there will be one more line on the right side than on the
# left side.
leftmostLineCenterx = -pitch*(nlines/2)
# for i in range(nlines):
# 	xcenter = leftmostLineCenterx+i*pitch
# 	line = mon.NShapes.createLine(-h,w,linelength,thetal,thetar,radl,radr)
# 	line.translate([xcenter,0.,0.])
# 	lineRegion = monte.addSubRegion(chamber,PMMAMSM,line)

for coordx in range(-10,10,4):
    for coordy in range(10,-10,-4):
        filename = dest + PathSep + "glC" + str(int(coordx))+"," + str(int(coordy)) + "data.txt"
        file = open(filename, 'w')
        print filename
        line = mon.NShapes.createLine(-h,w,linelength,thetal,thetar,radl,radr)
        line.translate([coordx,coordy,0.])
        lineRegion = monte.addSubRegion(chamber,PMMAMSM,line)

        # Scan parameters

        #yvals = [0.] #A single line scan at y = 0 (center of lines)
        yvals = []
        deltay = 1.
        ybottom = 0.
        ystart = ybottom - 10.5
        ystop = ybottom + 10.5
        y = ystart
        while y<ystop:
            yvals.append(y)
            y += deltay
        yvals.append(ystop)
        # The following parameters are set for a 201 nm scan centered on the position where the right
        # edge of the center line intersects the substrate. We could sample x at equal intervals (e.g., every nm)
        # but fine sampling is mainly important where the topography changes rapidly, so the code here samples
        # every 5 nm except within the zone that starts 25 nm to the left of the top corner and ends 25 nm to the
        # right of the bottom corner.
        xbottom = wnm/2.
        xtop = wnm/2.-hnm*jl.Math.tan(thetar)
        xstart = xbottom - 10.5
        xstop = xbottom + 10.5
        xfinestart = xtop-20.5
        if thetar<0.:	# undercut line
            xfinestop = xtop+20.5
        else: # normal line
            xfinestop = wnm/2.+20.5	# 20.5 nm to the right of the bottom corner

        xvals = []
        deltax = 1.
        x = xstart
        #while x<xfinestart:
        #	xvals.append(x)
        #	x += deltax
        #x = xfinestart
        #deltax = 1.
        #while x<xfinestop:
        #	xvals.append(x)
        #	x += deltax
        #x = xfinestop
        #deltax = 5.
        while x<xstop:
            xvals.append(x)
            x += deltax
        xvals.append(xstop)

        # Print simulation parameters to window and output file.
        #print >>file, "# Trajectories at each landing position: ",nTrajectories
        #print "# Trajectories at each landing position: ",nTrajectories
        # print >>file, "The only material is Si"
        # print "The only material is Si"
        # print >>file, "# Pitch of lines (nm): ",pitchnm
        # print "# Pitch of lines (nm): ",pitchnm
        # print >>file, "# lines: ",nlines
        # print "# lines: ",nlines
        # print >>file, "Line height (nm): ",hnm
        # print "Line height: ",hnm
        # print >>file, "Line bottom width (nm): ",wnm
        # print "Line bottom width (nm): ",wnm
        # print >>file, "Line length (nm): ",linelengthnm
        # print "Line length (nm): ",linelengthnm
        # print >>file, "Left and right sidewall angles (deg): ",thetaldeg,thetardeg
        # print "Left and right sidewall angles (deg): ",thetaldeg,thetardeg
        # print >>file, "Left and right top corner radii (nm): ",radlnm,radrnm
        # print "Left and right top corner radii (nm): ",radlnm,radrnm
        # print >>file, "Thicknesses of 1st and second layers (nm): ",layer1thicknessnm,layer2thicknessnm
        # print "Thicknesses of 1st and second layers (nm): ",layer1thicknessnm,layer2thicknessnm
        # print >>file, "Beam landing energies (eV): ",beamEeVvals
        # print "Beam landing energies (eV): ",beamEeVvals
        # print >>file, "Beam size (standard deviation, in nm): ",beamsizenm
        # print "Beam size (standard deviation, in nm): ",beamsizenm
        #
        # print >>file	# Blank line before start of calculation results
        # print
        #
        # print >>file, "beamE (eV)	x(nm)	y (nm)	SE count"
        # print "beamE (eV)	x(nm)	y (nm)	SE count"
        #print >>file, "beamE (eV)	x(nm)	y (nm)	BSE yield	SE yield    SE count"
        #print "beamE (eV)	x(nm)	y (nm)	BSE yield	SE yield	SE count"

        binSizeEV = 10.	# Width (in eV) of bins in energy histogram

        for beamEeV in beamEeVvals:
            beamE = epq.ToSI.eV(beamEeV)
            monte.setBeamEnergy(beamE) # sets this model's beam energy
            for ynm in yvals:
                y = ynm*meterspernm
                for xnm in xvals:
                    x = xnm*meterspernm
                    eg.setCenter([x,y,-h-20.*meterspernm]) # Aims the gun at x,y.

                    ## Define our backscatter detector.
                    #back=nm.BackscatterStats(monte)
                    #nbins = int(beamEeV/binSizeEV)
                    #monte.addActionListener(back)
                    #back.setEnergyBinCount(nbins)

                    # Add a trajectory image
                    if trajImg:  # output the trajectory image
                        #img=nm.TrajectoryImage(2048,2048,trajImgSize)
                        img=nm.TrajectoryImage(2048,2048,trajImgSize)
                        img.setMaxTrajectories(trajImgMaxTraj)
                        img.setYRange(-h-trajImgSize/10,.9*trajImgSize-h)
                        img.setXRange(x-trajImgSize/2.,x+trajImgSize/2.)
                        monte.addActionListener(img)

                    # Add a vrml
                    if VRML:  # output the trajectory image
                        #fos=jio.FileOutputStream("%s/angle - %lg deg.wrl" % (dest, phi))
                        fos=jio.FileOutputStream( "%s.wrl" % (dest))
                        tw=jio.OutputStreamWriter(fos,cs.Charset.forName("UTF-8"))
                        vrml=nm.TrajectoryVRML(monte,tw)
                        vrml.setMaxTrajectories(VRMLImgMaxTraj)
                        vrml.setTrajectoryWidth(0.25e-9)
                        vrml.setDisplayBackscatter(1)
                        vrml.addView("Gun",[0.0,0.0,-5.0e-7],[0.0,0.0,0.0])
                        vrml.addView("X-Axis",[1.0e-6,0.0,0.0,0.0],[0.0,0.0,0.0])
                        vrml.addView("Y-Axis",[0.0,1.0e-6,0.0,0.0],[0.0,0.0,0.0])
                        vrml.addView("Close perspective",[-110.e-8,100.e-8,-100.e-8],[-100.e-9,0.,0.])
                        vrml.renderSample()
                        monte.addActionListener(vrml)

                    print >>file, "%8.1f \t%8.1f \t%8.1f" % (beamEeV,xnm,ynm),
                    print "%8.1f \t%8.1f \t%8.1f" % (beamEeV,xnm,ynm),

                    # Run the simulation for acquisitionNum times
                    for ii in range(acquisitionNum):
                        # Define our backscatter detector.
                        back=nm.BackscatterStats(monte)
                        nbins = int(beamEeV/binSizeEV)
                        monte.addActionListener(back)
                        back.setEnergyBinCount(nbins)

                        monte.runMultipleTrajectories(nTrajectories)

                        hist = back.backscatterEnergyHistogram()
                        #fhist = back.forwardscatterEnergyHistogram()
                        energyperbineV = beamEeV/hist.binCount()
                        maxSEbin = 50./energyperbineV	# bin number of the one with 50 eV
                        totalSE = 0

                        for j in range(0,int(maxSEbin)):
                            totalSE = totalSE+hist.counts(j)
                            #print j,hist.counts(j)
                            #totalSE = totalSE+hist.counts(j)+fhist.counts(j)

                        #print >>file, "The totalSE count is %3.1f\n" %(float(totalSE))
                        #print "The totalSE count is %3.1f\n" %(float(totalSE))
                        SEcount = float(totalSE)
                        SEf = float(totalSE)/nTrajectories
                        bsf = back.backscatterFraction()-SEf#+float(totalFwd)/nTrajectories

                        print >>file, "\t%3.1f" % (SEcount),
                        print "\t%3.1f" % (SEcount),
                        #print >>file, "%8.1f \t%8.1f \t%8.1f \t%3.3f \t%3.3f \t%8.1f" % (beamEeV,xnm,ynm,bsf,SEf,SEcount),
                        #print "%8.1f \t%8.1f \t%8.1f \t%3.3f \t%3.3f \t%8.1f" % (beamEeV,xnm,ynm,bsf,SEf,SEcount),
                        #for j in range(0,int(maxSEbin)):
                        #	print >>file, "\t%8.1f	\t%8.1f  \t%8.1f" %(j,hist.counts(j), hist.maxValue(j)),
                        #	print "\t%8.1f	\t%8.1f  \t%8.1f" %(j,hist.counts(j), hist.maxValue(j)),

                        back.dump(jio.FileOutputStream(dest+PathSep+"backscatter.prn"))
                        monte.removeActionListener(back)
                    print " "
                    print >>file, " "

                    #back.dump(jio.FileOutputStream(dest+PathSep+"backscatter.prn"))
                    #monte.removeActionListener(back)

                    if trajImg:  # output the trajectory image
                        img.dumpToFile(dest)
                        monte.removeActionListener(img)

                    if VRML:
                        tw.flush()
                        fos.close()
                        monte.removeActionListener(vrml)

        file.close()