#!/APSshare/anaconda/x86_64/bin/python2.7

# Authors: R.Woods
#    Date: 08/27/15

import math

#propertiesDict = {'atomicNum':32; 'fano':0; 'energyPerHole':0; 'massDensity':0; 'massAttenuation':0; 'electronMobility':0; 'holeMobility':0; 'driftVelocity':0;}

def getElementProperties(propertiesDict, energy, temperature, biasfield):
	#----------
	# GERMANIUM
	if(propertiesDict['atomicNum'] == 32):
		print 'getting Germanium properties...'
		energy = energy * 0.001								# X-ray Energy (convert to MeV)
		propertiesDict['fano'] = 0.129							# Fano factor @ 77K
		propertiesDict['energyPerHole'] = 2.96  					# eV per e/hole pair
		propertiesDict['massDensity'] = 5.32						# g/cm^3 @ 300K
		propertiesDict['electronMobility'] = 4.9E7*math.pow(temperature, -1.66)		# cm^2/V*s 
		propertiesDict['holeMobility'] = 1.05E9*math.pow(temperature, -2.33) 		# 4.2E4 cm^2/V*s at 77K from G.F.Knoll

		# Mass Attenuation Data
		energyList = [1.00000E-03,1.10304E-03,1.21670E-03,1.21670E-03,1.23215E-03,1.24780E-03,1.24780E-03,1.32844E-03,1.41430E-03,1.41430E-03,1.50000E-03,2.00000E-03,3.00000E-03,4.00000E-03,5.00000E-03,6.00000E-03,8.00000E-03,1.00000E-02,1.11031E-02,1.11031E-02,1.50000E-02,2.00000E-02,3.00000E-02,4.00000E-02,5.00000E-02,6.00000E-02,8.00000E-02,1.00000E-01,1.50000E-01,2.00000E-01,3.00000E-01,4.00000E-01,5.00000E-01,6.00000E-01,8.00000E-01,1.00000E+00,1.25000E+00,1.50000E+00,2.00000E+00,3.00000E+00,4.00000E+00,5.00000E+00,6.00000E+00,8.00000E+00,1.00000E+01,1.50000E+01,2.00000E+01,]
		massAttenCoeffList = [1.89300E+03,1.50200E+03,1.19000E+03,4.38900E+03,4.73400E+03,4.97400E+03,6.69800E+03,6.34800E+03,5.55400E+03,6.28700E+03,5.47500E+03,2.71100E+03,9.61300E+02,4.49700E+02,2.47200E+02,1.50900E+02,6.89000E+01,3.74200E+01,2.81100E+01,1.98100E+02,9.15200E+01,4.22200E+01,1.38500E+01,6.20700E+00,3.33500E+00,2.02300E+00,9.50100E-01,5.55000E-01,2.49100E-01,1.66100E-01,1.13100E-01,9.32700E-02,8.21200E-02,7.45200E-02,6.42600E-02,5.72700E-02,5.10100E-02,4.65700E-02,4.08600E-02,3.52400E-02,3.27500E-02,3.15800E-02,3.10700E-02,3.10300E-02,3.15600E-02,3.34000E-02,3.52800E-02,]

		# Linear Interpolation between data points:
		for i, val in enumerate(energyList):
			if(energy < val):
				lowEnergy = energyList[i-1]
				highEnergy = energyList[i]
				lowMassAtten = massAttenCoeffList[i-1]
				highMassAtten = massAttenCoeffList[i]
				propertiesDict['massAttenuation'] = ( highMassAtten - lowMassAtten)/(highEnergy - lowEnergy) * (energy - lowEnergy) + lowMassAtten
				break

		# Drift Velocity Curve??? Can't find them yet.
		propertiesDict['driftVelocity'] = 5.0E6			# cm/s at T=77K, E=1000 V/cm

		return propertiesDict
		
	#----------
	# SILICON
	elif(propertiesDict['atomicNum'] == 14):
		print 'getting Silicon properties...'
		energy = energy * 0.001												# X-ray Energy (convert to MeV)
		propertiesDict['fano'] = 0.143											# Fano factor @ 77K
		propertiesDict['energyPerHole'] = 3.76  									# eV per e/hole pair
		propertiesDict['massDensity'] = 2.33										# g/cm^3 @ 300K
		propertiesDict['electronMobility'] = 2.1E4									# cm^2/V*s at 77K (1350 @ 300K)

		# Hole Mobility
		if(temperature == 77):
			propertiesDict['holeMobility'] = 1.1E4									# cm^2/V*s from G.F.Knoll
		elif(temperature == 300):
			propertiesDict['holeMobility'] = 480									# cm^2/V*s from G.F.Knoll

		# Mass Attenuation Data
		energyList = [1.0000E-03, 1.5000E-03, 1.8389E-03, 1.8389E-03, 2.0000E-03, 3.0000E-03, 4.0000E-03, 5.0000E-03, 6.0000E-03, 8.0000E-03, 1.0000E-02, 1.5000E-02, 2.0000E-02, 3.0000E-02, 4.0000E-02, 5.0000E-02, 6.0000E-02, 8.0000E-02, 1.0000E-01, 1.5000E-01, 2.0000E-01, 3.0000E-01, 4.0000E-01, 5.0000E-01, 6.0000E-01, 8.0000E-01, 1.0000E+00, 1.2500E+00, 1.5000E+00, 2.0000E+00, 3.0000E+00, 4.0000E+00, 5.0000E+00, 6.0000E+00, 8.0000E+00, 1.0000E+01, 1.5000E+01, 2.0000E+01,] 
		massAttenCoeffList = [1.57000E+03, 5.35500E+02, 3.09200E+02, 3.19200E+03, 2.77700E+03, 9.78400E+02, 4.52900E+02, 2.45000E+02, 1.47000E+02, 6.46800E+01, 3.38900E+01, 1.03400E+01, 4.46400E+00, 1.43600E+00, 7.01200E-01, 4.38500E-01, 3.20700E-01, 2.22800E-01, 1.83500E-01, 1.44800E-01, 1.27500E-01, 1.08200E-01, 9.61400E-02, 8.74800E-02, 8.07700E-02, 7.08200E-02, 6.36100E-02, 5.68800E-02, 5.18300E-02, 4.48000E-02, 3.67800E-02, 3.24000E-02, 2.96700E-02, 2.78800E-02, 2.57400E-02, 2.46200E-02, 2.35200E-02, 2.33800E-02,]

		# Linear Interpolation between data points:
		for i, val in enumerate(energyList):
			if(energy < val):
				lowEnergy = energyList[i-1]
				highEnergy = energyList[i]
				lowMassAtten = massAttenCoeffList[i-1]
				highMassAtten = massAttenCoeffList[i]
				propertiesDict['massAttenuation'] = ( highMassAtten - lowMassAtten)/(highEnergy - lowEnergy) * (energy - lowEnergy) + lowMassAtten
				break

		# Drift Velocity
		# Model from G.Ottavianim, C.Canali, A.Quaranta: "Electron and Hole Drift Velocity Measurement in Silicon and Their Empirical Relation to Electric Field and Tempeature", IEEE NS-22, 192(1975)
		propertiesDict['driftVelocity'] = 1.31E8 * math.pow(temperature, -2.2)* biasfield / math.pow( 1 + math.pow( biasfield/(1.24 * math.pow(temperature, 1.68)), 0.46*math.pow(temperature, 0.17)), 1/(0.46 * math.pow(temperature, 0.17)) )

		return propertiesDict
