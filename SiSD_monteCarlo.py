#!/usr/bin/python

#!/APSshare/anaconda/x86_64/bin/python2.7

# Descrpition: A Monte Carlo of semiconducter detector physics
#              to investigate charge splitting.
#
# Authors: R.Woods, D.Siddons
#    Date: 05/03/15
#          08/03/15
#          08/27/15
#          05/15/18
#          06/08/18 - added ENC to each pixel

import numpy
import time
import math
import pylab
from pylab import Rectangle
import sys
import NIST_data

#-----------------------------------------------
# PHYSICAL CONSTANTS
KBOLTZMAN = 1.38E-23    # J/K
ELEMCHARGE = 1.602E-19  # C

# BEAM CONSTANTS
ENERGY = 6.0            # keV
MEANBEAMX = 0           # microns
MEANBEAMY = 0           # microns
SIGMABEAMX = 1          # microns
SIGMABEAMY = 1          # microns

# DETECTOR CONSTANTS
NPIXELS = 5             # number of pixels
SIGMAENC  = 15          # RMS quanta added by readout electronics (for now, each pixel will have same ENC RMS)
TEMP = 270              # Kelvin
BIASVOLTAGE = 150       # volts
BIASHEIGHT = 500        # microns
DRIFTMAX = 500          # microns
STRIPWIDTH = 15.0       # microns - MUST BE FLOAT
STRIPHEIGHT = 2000      # microns
GAP = 5.0               # microns - Strip Spacing - MUST BE FLOAT
PROXIMITY = 2.5         # 6.3 #12.6, microns - "c" in Charge Sharing Effect: 2nd Order, like drift field fringing
                        # set PROXIMITY=GAP/2 for no charge carrier loss between strips

# SENSOR PROPERTIES
propertiesDict = {'atomicNum':14, 'fano':0, 'energyPerHole':0, 'massDensity':0, 'massAttenuation':0, 'electronMobility':0, 'holeMobility':0,}
NIST_data.getElementProperties(propertiesDict, ENERGY, TEMP)
FANO = propertiesDict['fano']
ENERGYPERHOLE =  propertiesDict['energyPerHole']
MASSDENSITY = propertiesDict['massDensity']
MASSATTENUATION = propertiesDict['massAttenuation']
ELECTRONMOBILITY = propertiesDict['electronMobility']
HOLEMOBILITY = propertiesDict['holeMobility']

# EVENTS
numEvents = 401        # Number of Photons Generated
numCounts = 0           # Counter for the number of Detector Counts
plotEvent = False
printEvent = False
numLostPhotons = 0      # Counter for photons that are not absorbed
xPositionArray = []
zPositionArray = []
tStripArray = []
charge0array = []
charge1array = []
charge2array = []
charge3array = []
charge4array = []
c0size=0
c1size=0
c2size=0
c3size=0
c4size=0

zeta = []
event = 0

#-----------------------------------------------
# MAIN

for event in range(0,numEvents):

      if(printEvent == True):
            print '---------------------'
            print 'Event = ' + str(event)

      #-------------------------
      # Sample the xray hit position (x, y) or just x
      #x = numpy.random.normal(MEANBEAMX, SIGMABEAMX, size=1)
      #x = numpy.random.uniform(0, STRIPWIDTH/2+GAP/2)  	                 # microns
      x = (-STRIPWIDTH-GAP) + ((2*STRIPWIDTH)+(2*GAP))/(numEvents-1)*event   # scan from middle of strip1 to middle of strip3
      y = 0
      xPositionArray.append(x)
      if(printEvent == True): print 'x = ' + str(x)

      #-------------------------
      # Sample the absorption length using PDF (see notes) = 1/scale * exp(-x/scale)
      # NOTE: Photon depth == Hole Drift Distance!
      photonDepth = numpy.random.exponential(scale = 1/(MASSDENSITY * MASSATTENUATION))       # units: cm
      z = DRIFTMAX-photonDepth*10000                                                          # units: microns
      #z = math.floor(photonDepth*10000)                                                      # units: microns
      #z = 2500                                                                               # units: microns
      if(printEvent == True): print 'z = ' + str(z)
      if (z > DRIFTMAX) or (z <= 0):
              print 'Photon not absorbed by detector!'
              numLostPhotons += 1
              xPositionArray.pop()      # remove the current x position
              continue                  # NOTE: this cut pulls the average away from the expectation value (1/mu = 0.929mm)

      zPositionArray.append(z)

      #-------------------------
      # Sample the number of electron/hole pairs for Energy, from G(Energy, Fano) 
      numElectrons = math.floor(numpy.random.normal(ENERGY/ENERGYPERHOLE*1000, math.sqrt(FANO*ENERGYPERHOLE*ENERGY*1000)/ENERGYPERHOLE))
      if(printEvent == True): print 'numElectrons = ' + str(numElectrons)


      #-------------------------
      # Calculate the cloud width due to transverse diffusion
      sigmaCloud = math.sqrt(2 * (KBOLTZMAN * TEMP / ELEMCHARGE) * (BIASHEIGHT/BIASVOLTAGE) * z)


      #-------------------------
      # Sample a list of electron positions by sampling the Gauss(x, sigmaCloud), Gauss(y, sigmaCloud)
      positionArrayX = numpy.random.normal(x, sigmaCloud, size=numElectrons)
      positionArrayY = numpy.random.normal(y, sigmaCloud, size=numElectrons)


      #-------------------------
      # Calculate Drift Time = tDrift 
      tDrift = z / (HOLEMOBILITY * BIASVOLTAGE / BIASHEIGHT) * 10E-8	      # units: s
      tStripArray.append(tDrift)
      if(printEvent == True): print 'tDrift = ' + str(tDrift)


      #-------------------------
      # Sample a list of Equivalen Noise Charge (ENC) from Poisson Distribution with lamdba=SIGMAENC
      noiseArray = numpy.random.poisson(SIGMAENC, NPIXELS)

      #-------------------------
      # Hit Positions in each Strip Region
      # NOTE: I wrote it this way to allow for the more general case where charge carriers can be lost in a gap region between strips.
      #       The value of PROXIMITY allows you to control the drift field fringing into the gap region. 
      #       For example, if PROXIMITY=GAP/2, then no charge is lost.   
      stripArray0 = positionArrayX[ (positionArrayX > (-5*STRIPWIDTH/2-(2*GAP)))           & (positionArrayX < (-3*STRIPWIDTH/2-(2*GAP))) ]
      proxArray0  = positionArrayX[ (positionArrayX > (-3*STRIPWIDTH/2-(2*GAP)))           & (positionArrayX < (-3*STRIPWIDTH/2-(2*GAP)+PROXIMITY)) ]
      gapArray0   = positionArrayX[ (positionArrayX > (-3*STRIPWIDTH/2-(2*GAP)+PROXIMITY)) & (positionArrayX < (-3*STRIPWIDTH/2-GAP-PROXIMITY)) ]
      proxArray1  = positionArrayX[ (positionArrayX > (-3*STRIPWIDTH/2-GAP-PROXIMITY))     & (positionArrayX < (-3*STRIPWIDTH/2-GAP)) ]
      stripArray1 = positionArrayX[ (positionArrayX > (-3*STRIPWIDTH/2-GAP))               & (positionArrayX < (  -STRIPWIDTH/2-GAP)) ]
      proxArray2  = positionArrayX[ (positionArrayX > (  -STRIPWIDTH/2-GAP))               & (positionArrayX < (  -STRIPWIDTH/2-GAP+PROXIMITY)) ]
      gapArray1   = positionArrayX[ (positionArrayX > (  -STRIPWIDTH/2-GAP+PROXIMITY))     & (positionArrayX < (  -STRIPWIDTH/2-PROXIMITY)) ]
      proxArray3  = positionArrayX[ (positionArrayX > (  -STRIPWIDTH/2-PROXIMITY))         & (positionArrayX < (  -STRIPWIDTH/2)) ]
      stripArray2 = positionArrayX[ (positionArrayX > (  -STRIPWIDTH/2))                   & (positionArrayX < (   STRIPWIDTH/2)) ]
      proxArray4  = positionArrayX[ (positionArrayX > (   STRIPWIDTH/2))                   & (positionArrayX < (   STRIPWIDTH/2+PROXIMITY)) ]
      gapArray2   = positionArrayX[ (positionArrayX > (   STRIPWIDTH/2+PROXIMITY))         & (positionArrayX < (   STRIPWIDTH/2+GAP-PROXIMITY)) ]
      proxArray5  = positionArrayX[ (positionArrayX > (   STRIPWIDTH/2+GAP-PROXIMITY))     & (positionArrayX < (   STRIPWIDTH/2+GAP)) ]
      stripArray3 = positionArrayX[ (positionArrayX > (   STRIPWIDTH/2+GAP))               & (positionArrayX < ( 3*STRIPWIDTH/2+GAP)) ]
      proxArray6  = positionArrayX[ (positionArrayX > ( 3*STRIPWIDTH/2+GAP))               & (positionArrayX < ( 3*STRIPWIDTH/2+GAP+PROXIMITY)) ] 
      gapArray3   = positionArrayX[ (positionArrayX > ( 3*STRIPWIDTH/2+GAP+PROXIMITY))     & (positionArrayX < ( 3*STRIPWIDTH/2+(2*GAP)-PROXIMITY)) ]
      proxArray7  = positionArrayX[ (positionArrayX > ( 3*STRIPWIDTH/2+(2*GAP)-PROXIMITY)) & (positionArrayX < ( 3*STRIPWIDTH/2+(2*GAP))) ]
      stripArray4 = positionArrayX[ (positionArrayX > ( 3*STRIPWIDTH/2+(2*GAP)))           & (positionArrayX < ( 5*STRIPWIDTH/2+(2*GAP))) ]

      if(printEvent == True):    
            print 'hits on strip0 = ' + str(numpy.size(stripArray0))
            print 'hits in prox0  = ' + str(numpy.size(proxArray0))
            print 'misses in gap0 = ' + str(numpy.size(gapArray0))
            print 'hits in prox1  = ' + str(numpy.size(proxArray1))
            print 'hits on strip1 = ' + str(numpy.size(stripArray1))
            print 'hits in prox2  = ' + str(numpy.size(proxArray2))
            print 'misses in gap1 = ' + str(numpy.size(gapArray1))
            print 'hits in prox3  = ' + str(numpy.size(proxArray3))
            print 'hits on strip2 = ' + str(numpy.size(stripArray2))
            print 'hits in prox4  = ' + str(numpy.size(proxArray4))
            print 'misses in gap2 = ' + str(numpy.size(gapArray2))
            print 'hits in prox5  = ' + str(numpy.size(proxArray5))
            print 'hits on strip3 = ' + str(numpy.size(stripArray3))
            print 'hits in prox6  = ' + str(numpy.size(proxArray6))
            print 'misses in gap3 = ' + str(numpy.size(gapArray3))
            print 'hits in prox7  = ' + str(numpy.size(proxArray7))
            print 'hits on strip4 = ' + str(numpy.size(stripArray4))
            print 'total holes collected = ' + str(numpy.size(stripArray0) +numpy.size(proxArray0) +numpy.size(gapArray0) +numpy.size(proxArray1) 
                                            +numpy.size(stripArray1) +numpy.size(proxArray2) +numpy.size(gapArray1) +numpy.size(proxArray3) 
                                            +numpy.size(stripArray2) +numpy.size(proxArray4) +numpy.size(gapArray2) +numpy.size(proxArray5) 
                                            +numpy.size(stripArray3) +numpy.size(proxArray6) +numpy.size(gapArray3) +numpy.size(proxArray7) 
                                            +numpy.size(stripArray4))

      # Strip 0 collected holes + ENC 
      chargeCarriers0 = numpy.size(stripArray0) + numpy.size(proxArray0) + noiseArray[0]
      charge0array.append(chargeCarriers0)
      c0size=c0size+1

      # Strip 1 collected holes + ENC
      chargeCarriers1 = numpy.size(proxArray1) + numpy.size(stripArray1) + numpy.size(proxArray2) + noiseArray[1]
      charge1array.append(chargeCarriers1)
      c1size=c1size+1

      # Strip 2 collected holes + ENC
      chargeCarriers2 = numpy.size(proxArray3) + numpy.size(stripArray2) + numpy.size(proxArray4) + noiseArray[2]
      charge2array.append(chargeCarriers2)
      c2size=c2size+1

      # Strip 3 collected holes + ENC
      chargeCarriers3 = numpy.size(proxArray5) + numpy.size(stripArray3) + numpy.size(proxArray6) + noiseArray[3]
      charge3array.append(chargeCarriers3)
      c3size=c3size+1
    
      # Strip 4 collected holes + ENC
      chargeCarriers4 = numpy.size(proxArray7) + numpy.size(stripArray4) + noiseArray[4]
      charge4array.append(chargeCarriers4)
      c4size=c4size+1


      # Event Energy
      if(printEvent == True):
            print 'chargeCarriers0 = ' + str(chargeCarriers0)
            print 'chargeCarriers1 = ' + str(chargeCarriers1)
            print 'chargeCarriers2 = ' + str(chargeCarriers2)
            print 'chargeCarriers3 = ' + str(chargeCarriers3)
            print 'chargeCarriers4 = ' + str(chargeCarriers4)


      if(plotEvent == True):
              # Plot Electron hit/miss positions
              pylab.figure(1)
              pylab.scatter(positionArrayX, positionArrayY, s=5, facecolor='blue', alpha=0.15)  	      # s=size, alpha=Transparency
              pylab.title('Hole Cloud Projection on Readout Plane, ' + str(ENERGY) + ' keV')
              pylab.xlabel('x position (um)')
              pylab.ylabel('y position (um)')
              pylab.grid(True)
              
              # Draw strips, gaps, and proximity regions
              # NOTE: Rectangle((x,y), width, height, angle=0.0, **kwargs); where (x,y) is lower left corner
              currentAxis = pylab.gca()
              currentAxis.add_patch(Rectangle((-5*STRIPWIDTH/2-2*GAP, -2500), STRIPWIDTH, 5000, facecolor='blue', alpha=0.10))
              currentAxis.add_patch(Rectangle((-3*STRIPWIDTH/2-2*GAP, -2500), GAP, 5000, facecolor='yellow', alpha=0.10))
              currentAxis.add_patch(Rectangle((-3*STRIPWIDTH/2-GAP, -2500), STRIPWIDTH, 5000, facecolor='blue', alpha=0.10))
              currentAxis.add_patch(Rectangle((-STRIPWIDTH/2-GAP, -2500), GAP, 5000, facecolor='yellow', alpha=0.10))
              currentAxis.add_patch(Rectangle((-STRIPWIDTH/2, -2500), STRIPWIDTH, 5000, facecolor='blue', alpha=0.10))
              currentAxis.add_patch(Rectangle(( STRIPWIDTH/2, -2500), GAP, 5000, facecolor='yellow', alpha=0.10))
              currentAxis.add_patch(Rectangle(( STRIPWIDTH/2+GAP, -2500), STRIPWIDTH, 5000, facecolor='blue', alpha=0.10))
              currentAxis.add_patch(Rectangle(( 3*STRIPWIDTH/2+GAP, -2500), GAP, 5000, facecolor='yellow', alpha=0.10))
              currentAxis.add_patch(Rectangle(( 3*STRIPWIDTH/2+2*GAP, -2500), STRIPWIDTH, 5000, facecolor='blue', alpha=0.10))

              currentAxis.add_patch(Rectangle((-3*STRIPWIDTH/2-2*GAP, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))
              currentAxis.add_patch(Rectangle((-3*STRIPWIDTH/2-GAP-PROXIMITY, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))
              currentAxis.add_patch(Rectangle((-STRIPWIDTH/2-GAP, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))
              currentAxis.add_patch(Rectangle((-STRIPWIDTH/2-PROXIMITY, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))
              currentAxis.add_patch(Rectangle(( STRIPWIDTH/2, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))
              currentAxis.add_patch(Rectangle(( STRIPWIDTH/2+GAP-PROXIMITY, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))
              currentAxis.add_patch(Rectangle(( 3*STRIPWIDTH/2+GAP, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))
              currentAxis.add_patch(Rectangle(( 3*STRIPWIDTH/2+2*GAP-PROXIMITY, -2500), PROXIMITY, 5000, facecolor='red', alpha=0.20))

              # Set the Display Region
              pylab.xlim([-5*STRIPWIDTH/2-2*GAP, 5*STRIPWIDTH/2+2*GAP])
              pylab.ylim([-50, 50])
              pylab.show()


#--------------------------------------------
# Analysis

print '-------------------------------------'
print 'Simulation finished'
print 'Beam Energy = ' + str(ENERGY) + ' keV'
print 'Mass Desnsity = ' + str(MASSDENSITY)
print 'Mass Attenuation = ' + str(MASSATTENUATION)
print 'ENC rms = ' + str(SIGMAENC)
print 'Number of X-rays = ' + str(numEvents)
print 'X-rays not absorbed = ' + str(numLostPhotons)
print 'X-rays absorbed = ' + str(len(zPositionArray))
print 'Average z = ' + str(numpy.sum(zPositionArray)/len(zPositionArray))


# Histogram xPositionArray
#pylab.figure(4)
#counts, bins, ignored = pylab.hist(xPositionArray, 25, normed=False, facecolor='red', alpha=0.4)
#pylab.title('X-ray beam positions')
#pylab.xlabel('x (microns)')
#pylab.ylabel('Counts')
#pylab.grid(True)

# Histogram zPositionArray
#pylab.figure(5)
#binWidth = 50
#pylab.xlim([0-(binWidth/2), DRIFTMAX+(binWidth/2)])
#pylab.ylim(0,300)
#pylab.hist(zPositionArray, bins=numpy.arange((0-binWidth/2), (2501+binWidth/2), binWidth), normed=False, facecolor='blue', alpha=0.4)
#pylab.hist(zPositionArray, bins=numpy.arange((0), (2501), binWidth), normed=False, facecolor='blue', alpha=0.4, log=False, label='Energy = ' +str(ENERGY)+'keV')
#pylab.title('Hole Cluster Drift Distance, ' + str(ENERGY) + ' keV')
#pylab.xlabel('z (microns)')
#pylab.ylabel('Counts')
#pylab.grid(True)
#pylab.label = 'extra stuff'
##pylab.legend(loc='upper right', frameon=True, fancybox=True, shadow=False, mode='none', title='Stats' )
#pylab.show()


#--------------------------------------------
# Pete's Plots:

pylab.figure(2)
pylab.title('1D scan in x')
pylab.xlabel('x (microns)')
pylab.ylabel('holes collected')
print 'charge0array size = ' + str(len(charge0array))
print 'xPositionArray size = ' + str(len(xPositionArray))
pylab.plot(xPositionArray,charge0array, '-r', label='strip_0')
pylab.plot(xPositionArray,charge1array, '-b', label='strip_1')
pylab.plot(xPositionArray,charge2array, '-g', label='strip_2')
pylab.plot(xPositionArray,charge3array, '-c', label='strip_3')
pylab.plot(xPositionArray,charge4array, '-m', label='strip_4')

pylab.figure(3)
pylab.title('???')
pylab.plot(charge1array,charge2array)

pylab.figure(4)
pylab.title('zeta')
print "Plotting zeta"
print "c1size="
print c1size
print "charge1(1)" 
print charge1array[1] 
print "charge2(1)"
print charge2array[1]
#zeta=charge1array/(charge1array+charge2array)

for i in range(numEvents - numLostPhotons):
      p=float(charge1array[i])/float((charge1array[i]+charge2array[i]))
      #print p,charge1array[i],charge2array[i]
      zeta.append(float(charge1array[i])/float((charge1array[i]+charge2array[i])))

pylab.plot(xPositionArray,zeta)
pylab.show()
