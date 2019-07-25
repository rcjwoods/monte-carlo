#!/APSshare/anaconda/x86_64/bin/python2.7

#!/APSshare/epd/rh6-x86_64/bin/python2.7

# Descrpition: A Monte Carlo simulation of the Germanium Strip Detector
#			   to investigate charge splitting.
#
# Authors: R.Woods
#    Date: 05/03/15
#		   08/03/15
#		   08/27/15
#          10/21/15 - Update: Removing Low-field approximation from diffusion calc.

import numpy
import time
import math
import pylab
from pylab import Rectangle
import sys
import NIST_data

#-----------------------------------------------
# PHYSICAL CONSTANTS
KBOLTZMAN = 1.38E-23		# J/K
ELEMCHARGE = 1.602E-19		# C

# BEAM CONSTANTS
ENERGY = 20.0				# keV
MEANBEAMX = 0				# microns
MEANBEAMY = 0				# microns
SIGMABEAMX = 1				# microns
SIGMABEAMY = 1				# microns

# DETECTOR CONSTANTS
TEMP = 300.0							# Kelvin -- NOTE: For Germanium, use 77K (only temp that works)
BIASVOLTAGE = 300.0						# Volts  -- NOTE: For Germanium, use 300 V (E=1000 V/cm)
BIASHEIGHT = 3000.0						# microns
EFIELD = BIASVOLTAGE/BIASHEIGHT*10000	# V/cm
DRIFTMAX = 2500.0						# microns
STRIPWIDTH = 125.0						# microns
STRIPHEIGHT = 5000.0					# microns
GAP = 30.0								# microns - Strip Spacing
PROXIMITY = GAP/2 						#6.3 #12.6, microns - "c" in Charge Sharing Effect: 2nd Order

# SENSOR PROPERTIES						# NOTE: set 'atomicNum':14 for Silicon, or 32 for Germanium, other values are initially set to zero
propertiesDict = {'atomicNum':32, 'fano':0, 'energyPerHole':0, 'massDensity':0, 'massAttenuation':0, 'electronMobility':0, 'holeMobility':0, 'driftVelocity':0}
NIST_data.getElementProperties(propertiesDict, ENERGY, TEMP, EFIELD)
FANO = propertiesDict['fano']
ENERGYPERHOLE =  propertiesDict['energyPerHole']
MASSDENSITY = propertiesDict['massDensity']
MASSATTENUATION = propertiesDict['massAttenuation']
ELECTRONMOBILITY = propertiesDict['electronMobility']
HOLEMOBILITY = propertiesDict['holeMobility']
DRIFTVELOCITY = propertiesDict['driftVelocity']			# cm/s

# EVENTS
numEvents = 5000			# Number of Photons Generated
numCounts = 0				# Number of Detector Counts
plotEvent = False
numLostPhotons = 0
spectrumArray1= []
spectrumArray2= []
xPositionArray = []
zPositionArray = []
tStripArray = []
event = 0

#-----------------------------------------------
# MAIN

for event in range(0, numEvents):
#while(1):
	#event += 1
	print '---------------------'
	print 'Event = ' + str(event)


	#-------------------------
	# Sample the xray hit position (x, y) or just x
	#x = numpy.random.normal(MEANBEAMX, SIGMABEAMX, size=1)
	x = numpy.random.uniform(0, STRIPWIDTH/2 + GAP/2)		# microns
	#x = 0
	y = 0
	xPositionArray.append(x)
	print 'x = ' + str(x)


	#-------------------------
	# Sample the absorption length using PDF (see notes) = 1/scale * exp(-x/scale)
	# NOTE: Photon depth == Hole Drift Distance!
	photonDepth = numpy.random.exponential(scale = 1/(MASSDENSITY * MASSATTENUATION))	# units: cm
	z = photonDepth*10000								# microns
	#z = math.floor(photonDepth*10000)					# microns
	#z = 2500											# microns
	print 'z = ' + str(z) + ' microns'
	if (z > DRIFTMAX) or (z <= 0):
		print 'Photon not absorbed by detector!'
		numLostPhotons += 1
		continue										# NOTE: this cut pulls the average away from the expectation value (1/mu = 0.929mm)

	zPositionArray.append(z)
	if(len(zPositionArray)) == 5000: break


	#-------------------------
	# Sample the number of electron/hole pairs for Energy, from G(Energy, Fano) 
	numElectrons = math.floor(numpy.random.normal(ENERGY/ENERGYPERHOLE*1000, math.sqrt(FANO*ENERGYPERHOLE*ENERGY*1000)/ENERGYPERHOLE))
	print 'numElectrons = ' + str(numElectrons)


	#-------------------------
	# Calculate the cloud width due to transverse diffusion
	sigmaCloud = math.sqrt(2 * (KBOLTZMAN * TEMP / ELEMCHARGE) * HOLEMOBILITY * z / DRIFTVELOCITY * 10000)
	# Low-field Approximation:
	#sigmaCloud = math.sqrt(2 * (KBOLTZMAN * TEMP / ELEMCHARGE) * (BIASHEIGHT/BIASVOLTAGE) * z)


	#-------------------------
	# Sample a list of electron positions by sampling the Gauss(x, sigmaCloud), Gauss(y, sigmaCloud)
	positionArrayX = numpy.random.normal(x, sigmaCloud, size=numElectrons)
	positionArrayY = numpy.random.normal(y, sigmaCloud, size=numElectrons)


	#-------------------------
	# Calculate Drift Time = tDrift 
	#tDrift = z / (HOLEMOBILITY * BIASVOLTAGE / BIASHEIGHT) * 10E-8		# units: s
	tDrift = z / DRIFTVELOCITY *10000									# units: s
	tStripArray.append(tDrift)
	print 'vDrift = ' + str(DRIFTVELOCITY)
	print 'tDrift = ' + str(tDrift)


	#-------------------------
	# Hit Positions in each Strip Region
	stripArray0 = positionArrayX[ (positionArrayX < (-STRIPWIDTH/2-GAP)) ]
	proxArray0  = positionArrayX[ (positionArrayX > (-STRIPWIDTH/2-GAP)) & (positionArrayX < (-STRIPWIDTH/2-GAP+PROXIMITY)) ]
	gapArray1   = positionArrayX[ (positionArrayX > (-STRIPWIDTH/2-GAP+PROXIMITY)) & (positionArrayX < (-STRIPWIDTH/2-PROXIMITY)) ]
	proxArray1  = positionArrayX[ (positionArrayX > (-STRIPWIDTH/2-PROXIMITY)) & (positionArrayX < -STRIPWIDTH/2) ]
	stripArray1 = positionArrayX[ (positionArrayX >  -STRIPWIDTH/2) & (positionArrayX < STRIPWIDTH/2) ]
	proxArray2  = positionArrayX[ (positionArrayX >   STRIPWIDTH/2) & (positionArrayX < (STRIPWIDTH/2+PROXIMITY)) ]
	gapArray2   = positionArrayX[ (positionArrayX >  (STRIPWIDTH/2+PROXIMITY)) & (positionArrayX < (STRIPWIDTH/2+GAP-PROXIMITY)) ]
	proxArray3  = positionArrayX[ (positionArrayX >  (STRIPWIDTH/2+GAP-PROXIMITY)) & (positionArrayX < (STRIPWIDTH/2+GAP)) ]
	stripArray2 = positionArrayX[ (positionArrayX >  (STRIPWIDTH/2+GAP)) & (positionArrayX < (STRIPWIDTH/2+GAP+STRIPWIDTH)) ]

	print 'hits on strip0 = ' + str(numpy.size(stripArray0))
	print 'hits in prox0  = ' + str(numpy.size(proxArray0))
	print 'misses in gap1 = ' + str(numpy.size(gapArray1))
	print 'hits in prox1  = ' + str(numpy.size(proxArray1))
	print 'hits on strip1 = ' + str(numpy.size(stripArray1))
	print 'hits in prox2  = ' + str(numpy.size(proxArray2))
	print 'misses in gap2 = ' + str(numpy.size(gapArray2))
	print 'hits in prox3  = ' + str(numpy.size(proxArray3))
	print 'hits on strip2 = ' + str(numpy.size(stripArray2))
	print 'total holes collected = ' + str(numpy.size(stripArray0) +numpy.size(proxArray0) +numpy.size(gapArray1) +numpy.size(proxArray1) +numpy.size(stripArray1) +numpy.size(proxArray2) +numpy.size(gapArray2) +numpy.size(proxArray3) +numpy.size(stripArray2))

	# Strip 0 holes
	chargeCarriers0 = numpy.size(stripArray0) + numpy.size(proxArray0)

	# Strip 1 holes
	chargeCarriers1 = numpy.size(proxArray1) + numpy.size(stripArray1) + numpy.size(proxArray2)
	
	# Strip 2 holes
	chargeCarriers2 = numpy.size(proxArray3) + numpy.size(stripArray2)
	
	# Event Energy
	spectrumArray1.append(chargeCarriers1)
	spectrumArray2.append(chargeCarriers2)
	print 'chargeCarriers0 = ' + str(chargeCarriers0)
	print 'chargeCarriers1 = ' + str(chargeCarriers1)
	print 'chargeCarriers2 = ' + str(chargeCarriers2)


	if (plotEvent == True):
		# Plot Electron hit/miss positions
		pylab.figure(1)
		pylab.scatter(positionArrayX, positionArrayY, s=5, facecolor='blue', alpha=0.15)		# s=size, alpha=Transparency
		pylab.title('Hole Cloud Projection on Readout Plane, ' + str(ENERGY) + ' keV')
		pylab.xlabel('x position (um)')
		pylab.ylabel('y position (um)')
		pylab.grid(True)
		
		# Draw strips, gaps, and proximity regions
		currentAxis = pylab.gca()
		currentAxis.add_patch(Rectangle((-3*STRIPWIDTH/2-GAP, -2500), STRIPWIDTH, 5000, facecolor='blue', alpha=0.10))
		currentAxis.add_patch(Rectangle((-STRIPWIDTH/2-GAP, -2500), GAP, 5000, facecolor='yellow', alpha=0.10))
		currentAxis.add_patch(Rectangle((-STRIPWIDTH/2, -2500), STRIPWIDTH, 5000, facecolor='blue', alpha=0.10))
		currentAxis.add_patch(Rectangle(( STRIPWIDTH/2, -2500), GAP, 5000, facecolor='yellow', alpha=0.10))
		currentAxis.add_patch(Rectangle(( STRIPWIDTH/2+GAP, -2500), STRIPWIDTH, 5000, facecolor='blue', alpha=0.10))
		currentAxis.add_patch(Rectangle((-STRIPWIDTH/2-GAP, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))
		currentAxis.add_patch(Rectangle((-STRIPWIDTH/2-PROXIMITY, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))
		currentAxis.add_patch(Rectangle(( STRIPWIDTH/2, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))
		currentAxis.add_patch(Rectangle(( STRIPWIDTH/2+GAP-PROXIMITY, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))

		# Set the Display Region
		pylab.xlim([-3*STRIPWIDTH/2-GAP, 3*STRIPWIDTH/2+GAP])
		pylab.ylim([-100, 100])
		pylab.show()


#--------------------------------------------
# Analysis

print '-------------------------------------'
print 'Simulation finished'
print 'Beam Energy = ' + str(ENERGY) + ' keV'
print 'Mass Desnsity = ' + str(MASSDENSITY)
print 'Mass Attenuation = ' + str(MASSATTENUATION)
print 'Number of X-rays = ' + str(numEvents)
print 'X-rays not absorbed = ' + str(numLostPhotons)
print 'X-rays absorbed = ' + str(len(zPositionArray))
print 'Average z = ' + str(numpy.sum(zPositionArray)/len(zPositionArray))

# Histogram spectrumArray
pylab.figure(3)
counts, bins, ignored = pylab.hist(spectrumArray1, 200, normed=False, facecolor='green', alpha=0.4, log=True)
pylab.title('Spectrum, ' + str(ENERGY) + ' keV')
pylab.xlabel('Holes Collected by Strip_0')
pylab.ylabel('Counts')
pylab.grid(True)

# Histogram spectrumArray2
pylab.figure(6)
counts, bins, ignored = pylab.hist(spectrumArray2, 200, normed=False, facecolor='green', alpha=0.4, log=True)
pylab.title('Spectrum, ' + str(ENERGY) + ' keV')
pylab.xlabel('Holes Collected by Strip_2')
pylab.ylabel('Counts')
pylab.grid(True)



#----------------------------------------
# Histogram xPositionArray
pylab.figure(4)
counts, bins, ignored = pylab.hist(xPositionArray, 25, normed=False, facecolor='red', alpha=0.4)
pylab.title('X-ray beam positions')
pylab.xlabel('x (microns)')
pylab.ylabel('Counts')
pylab.grid(True)

# Histogram zPositionArray
pylab.figure(5)
binWidth = 50
pylab.xlim([0-(binWidth/2), 2500+(binWidth/2)])
pylab.ylim(0,300)
#pylab.hist(zPositionArray, bins=numpy.arange((0-binWidth/2), (2501+binWidth/2), binWidth), normed=False, facecolor='blue', alpha=0.4)
pylab.hist(zPositionArray, bins=numpy.arange((0), (2501), binWidth), normed=False, facecolor='blue', alpha=0.4, log=False, label='Energy = ' +str(ENERGY)+'keV')
pylab.title('Hole Cluster Drift Distance, ' + str(ENERGY) + ' keV')
pylab.xlabel('z (microns)')
pylab.ylabel('Counts')
pylab.grid(True)
pylab.label = 'extra stuff'
#pylab.legend(loc='upper right', frameon=True, fancybox=True, shadow=False, mode='none', title='Stats' )
pylab.show()
