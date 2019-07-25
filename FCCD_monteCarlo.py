#!/APSshare/anaconda/x86_64/bin/python2.7

#!/Users/dp/anaconda/bin/python2.7
#!/usr/bin/python


# Descrpition: A Monte Carlo of semiconducter detector physics
#              to investigate charge splitting.
#
# Authors: R.Woods
#    Date: 05/03/15
#          08/03/15
#	         08/27/15
#          05/15/18

import numpy
import time
import math
import pylab
from pylab import Rectangle
import sys
import NIST_data2 as NIST_data

#-----------------------------------------------
# PHYSICAL CONSTANTS
KBOLTZMAN = 1.38E-23        # J/K
ELEMCHARGE = 1.602E-19      # C

# BEAM CONSTANTS
RADSOURCE = 'Cd109'
ENERGY = 0.0            # keV
MEANBEAMX = 0			# microns
MEANBEAMY = 0			# microns
SIGMABEAMX = 1			# microns
SIGMABEAMY = 1			# microns

# DETECTOR CONSTANTS
TEMP = 270			    # Kelvin
BIASVOLTAGE = 30		# volts
BIASHEIGHT = 350		# microns
DRIFTMAX = 350			# microns
PIXELWIDTH = 30.0		# microns - MUST BE FLOAT
NPIXELS = 49            # number of pixels
NROWS = 7               # number of rows
NCOLUMNS = 7            # number of columns

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
numEvents = 10		    # Number of Photons Generated
numCounts = 0			# Counter for the number of Detector Counts
plotEvent = True        # Display each event on the pixel array
numLostPhotons = 0		# Counter for photons that are not absorbed
spectrumArray= []
xPositionArray = []
yPositionArray = []
zPositionArray = []
nClusterArray = []
tDriftArray = []
event = 0

#-----------------------------------------------
# MAIN
for event in range(0,numEvents):
      print '---------------------'
      print 'Event = ' + str(event)
    
      #-------------------------
      # Sample the Xray energy
      #ENERGY = NIST_data.getRadioactiveSourcePhoton(RADSOURCE)
      ENERGY = 8.0
      spectrumArray.append(ENERGY)
      #print 'xray energy = ' + str(ENERGY)

      #-------------------------
      # Sample the xray hit position (x, y) or just x
      #x = numpy.random.normal(MEANBEAMX, SIGMABEAMX, size=1)
      #x = (-PIXELWIDTH-GAP) + ((2*PIXELWIDTH)+(2*GAP))/(numEvents-1)*event          # scan from middle of strip1 to middle of strip3
      #x = int(numpy.random.uniform(-3*PIXELWIDTH/2, 3*PIXELWIDTH/2))  	             # microns
      #y = int(numpy.random.uniform(-3*PIXELWIDTH/2, 3*PIXELWIDTH/2))  	             # microns
      x = 2*PIXELWIDTH*event/(numEvents-1) - PIXELWIDTH 
      y = -2*PIXELWIDTH*event/(numEvents-1) + PIXELWIDTH
      xPositionArray.append(x)
      yPositionArray.append(y)
      print 'x = ' + str(x) + ', y = ' + str(y)

      #-------------------------
      # Sample the absorption length using PDF (see notes) = 1/scale * exp(-x/scale)
      # NOTE: Photon depth == Hole Drift Distance!
      NIST_data.getElementProperties(propertiesDict, ENERGY, TEMP)
      MASSDENSITY = propertiesDict['massDensity']
      MASSATTENUATION = propertiesDict['massAttenuation']
        
      photonDepth = numpy.random.exponential(scale = 1/(MASSDENSITY * MASSATTENUATION))       # units: cm
      z = DRIFTMAX-photonDepth*10000							      # units: microns
      #z = math.floor(photonDepth*10000)				                      # units: microns
      #z = 2500 								              # units: microns
      print 'z = ' + str(z)
      if (z > DRIFTMAX) or (z <= 0):
              print 'Photon not absorbed by sensor!'
              numLostPhotons += 1
              continue  				# NOTE: this cut pulls the average away from the expectation value (1/mu = 0.929mm)
      zPositionArray.append(z)

      #-------------------------
      # Sample the number of electron/hole pairs for Energy, from G(Energy, Fano) <------ ADD AMPLIFIER NOISE HERE! G(Energy, Fano^2+Amp^2)
      numElectrons = int(math.floor(numpy.random.normal(ENERGY/ENERGYPERHOLE*1000, math.sqrt(FANO*ENERGYPERHOLE*ENERGY*1000)/ENERGYPERHOLE)))
      print 'numElectrons = ' + str(numElectrons)
      nClusterArray.append(numElectrons)

      #-------------------------
      # Calculate the cloud width due to transverse diffusion
      sigmaCloud = math.sqrt(2 * (KBOLTZMAN * TEMP / ELEMCHARGE) * (BIASHEIGHT/BIASVOLTAGE) * z)

      #-------------------------
      # Sample a list of electron positions by sampling the Gauss(x, sigmaCloud), Gauss(y, sigmaCloud)
      electronArrayX = numpy.random.normal(x, sigmaCloud, size=numElectrons)
      electronArrayY = numpy.random.normal(y, sigmaCloud, size=numElectrons)

      #-------------------------
      # Calculate Drift Time = tDrift 
      tDrift = z / (HOLEMOBILITY * BIASVOLTAGE / BIASHEIGHT) * 10E-8	      # units: s
      tDriftArray.append(tDrift)
      print 'tDrift = ' + str(tDrift)

      #-------------------------
      # Sum the Hits on each Pixel  
      pixelCountArray = numpy.zeros(NPIXELS)
    
      # Loop over all elecrons in the distributions - brute force!
      for k in range(0, numElectrons):
          xPos = electronArrayX[k]
          yPos = electronArrayY[k]

          # Loop over all pixels
          for i in range(0, NCOLUMNS):
              for j in range(0, NROWS):
                  xEdge = -NROWS*PIXELWIDTH/2 + j*PIXELWIDTH
                  yEdge = -NCOLUMNS*PIXELWIDTH/2 + i*PIXELWIDTH
                
                  if( (xPos > xEdge) & (xPos < xEdge+PIXELWIDTH) & (yPos > yEdge) & (yPos < yEdge+PIXELWIDTH) ):
                      pixelCountArray[i*NCOLUMNS+j] += 1


      print 'Pixel Counts = ' + str(pixelCountArray)
      holesCollected = pixelCountArray.sum()       

      #-------------------------
      # Plot Pixels
      if (plotEvent == True):
              # Plot Electron hit/miss positions
              pylab.figure(1)
              pylab.scatter(electronArrayX, electronArrayY, s=5, facecolor='magenta', alpha=0.15)  # s=size, alpha=Transparency
              pylab.title('Electron Cloud Projection onto Readout Plane')
              pylab.xlabel('x position (um)')
              pylab.ylabel('y position (um)')
              pylab.grid(False)
              
              # Draw pixels
              # NOTE: Rectangle((x,y), width, height, angle=0.0, **kwargs); where (x,y) is lower left corner
              currentAxis = pylab.gca()

              for i in range(0, NCOLUMNS):
                  for j in range(0, NROWS):
                      xEdge = -NROWS*PIXELWIDTH/2 + j*PIXELWIDTH
                      yEdge = -NCOLUMNS*PIXELWIDTH/2 + i*PIXELWIDTH
                      currentAxis.add_patch(Rectangle((xEdge, yEdge), PIXELWIDTH, PIXELWIDTH, facecolor='blue', alpha=0.10))
                      currentAxis.annotate(str(int(pixelCountArray[i*NCOLUMNS+j])), (xEdge+PIXELWIDTH/2, yEdge+PIXELWIDTH/2), 
                                                      color='black', weight='normal', fontsize=12, ha='center', va='center')

              # Set the Display Region
              pylab.xlim([-NROWS*PIXELWIDTH/2, NROWS*PIXELWIDTH/2])
              pylab.ylim([-NCOLUMNS*PIXELWIDTH/2, NCOLUMNS*PIXELWIDTH/2])
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
print 'X-rays event energies = ' + str(spectrumArray)
