
import numpy
import sys 
sys.path.append('Libs')
import profilesLib
import PRFG3Lib


# Setting the profiles 
Ubas = 5.6
scale = 400.0
zBottom = 0.0/scale
zTop = 1400.0/scale
nZ = 51

[z0,zMin,kt] = profilesLib.getECCoefficients("EC Category I")

parameters = dict()
parameters['seeds'] = numpy.array([45,101,212,333,490,700,1100])/scale
parameters['uProfile'] = profilesLib.UProfile(zBottom,zTop,nZ).ECProfile(Ubas,z0,kt,zMin,scale)
parameters['sProfile'] = profilesLib.SProfile(zBottom,zTop,nZ).ECProfile(Ubas,kt)
parameters['lProfile'] = profilesLib.LProfile(zBottom,zTop,nZ).ESDUProfile(z0,scale)
parameters['tProfile'] = profilesLib.TProfile().rampProfile(10.0)

# Setting the parameters for the seed
seedTurbParameters = dict()
seedTurbParameters['L'] = numpy.array([[1,0.6,0.4],[0.6,0.5,0.3],[0.35,0.25,0.2]])
seedTurbParameters['Covariance'] = numpy.diag(numpy.power([1,0.75,0.5],2))
seedTurbParameters['dEd'] = 5
seedTurbParameters['finalE'] = 0.95
        
mainSeed = PRFG3Lib.PRFG3Turb(seedTurbParameters)
mainSeed.PRFG()
mainSeed.plotInflowMetrics()

parameters['main seed'] = mainSeed

# Writing out
profilesLib.writeInflow(parameters,fileName='inflowGenerationDict')
print('Done!')
 



