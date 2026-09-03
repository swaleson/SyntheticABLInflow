import numpy

ofOpen = '''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  21.12                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
version 21.12;
format ascii;
class dictionary;
location "constant";
object random;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n'''

def fmt(val):
    """ Format  s number for output"""
    return "{:.3f}".format(val)

def getECCoefficients(category):
    """ Returns the eurocode coefficients """
    if category=="EC Category 0":
        z0 = 0.003
        zMin = 1.0
        kt = 0.16    #Not reported in EC Table 4.1 
    elif category=="EC Category I":
        z0 = 0.03
        zMin = 1.0
        kt = 0.18
    elif category=="EC Category II":
        z0 = 0.05
        zMin = 2.0
        kt = 0.19
    elif category=="EC Category III":
        z0 = 0.3
        zMin = 5.0
        kt = 0.22
    elif category=="EC Category IV":
        z0 = 0.7
        zMin = 7.86
        kt = 0.24
    return [z0,zMin,kt]

def getAIJCoefficients(category):
    """ Returns the AIJ coefficients """
    if category=="AIJ Category I":
        zb = 5.0
        zg = 250.0
        alpha = 0.1
    elif category=="AIJ Category II":
        zb = 5.0
        zg = 350.0
        alpha = 0.15
    elif category=="AIJ Category III":
        zb = 10.0
        zg = 450.0
        alpha = 0.2
    elif category=="AIJ Category IV":
        zb = 20.0
        zg = 550.0
        alpha = 0.27
    elif category=="AIJ Category V":
        zb = 30.0
        zg = 650.0
        alpha = 0.35
    return [zb,zg,alpha]


def writeInflow(parameters,fileName='inflowDict'):
    """ Writes the inflowDict """
    strOutList = [ofOpen]
    # Writing descriptions
    strOutList.append(parameters['tProfile'].description)
    strOutList.append(parameters['uProfile'].description)
    strOutList.append(parameters['sProfile'].description)
    strOutList.append(parameters['lProfile'].description)
    strOut = "\n// ".join(strOutList) + '\n'
    strOut = strOut + parameters["main seed"].description
    strOut = strOut + '\n\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n'
    # Writing profiles
    cutOffLevel = 'cutOff 0.5; \n\n'
    strOut = strOut  + cutOffLevel
    seedsList = list(parameters['seeds'])
    seedsStr = "\n".join([fmt(ss) for ss in seedsList])
    strOut = strOut + 'seeds (' + seedsStr + ');\n\n' 
    strOut = strOut + parameters['tProfile'].write() + '\n\n'
    strOut = strOut + parameters['uProfile'].write() + '\n\n'
    strOut = strOut + parameters['sProfile'].write() + '\n\n'
    strOut = strOut + parameters['lProfile'].write() + '\n\n'
    # Writing main seed
    kStr = "kMainSeed ( \n"
    kStr = kStr + str(parameters["main seed"].data["Turbulent field"]["k"].tolist()).replace(",","").replace("[","(").replace("]",")").replace(") (",")\n(")[1:-1]
    kStr = kStr + ");\n\n"
    pStr = "pMainSeed ( \n"
    pStr = pStr + str(parameters["main seed"].data["Turbulent field"]["p"].tolist()).replace(",","").replace("[","(").replace("]",")").replace(") (",")\n(")[1:-1]
    pStr = pStr + ");\n\n"
    qStr = "qMainSeed ( \n"
    qStr = qStr + str(parameters["main seed"].data["Turbulent field"]["q"].tolist()).replace(",","").replace("[","(").replace("]",")").replace(") (",")\n(")[1:-1]
    qStr = qStr + ");\n"
    strOut = strOut + kStr + pStr + qStr
    
    with open(fileName,'w') as ff:
        ff.write(strOut)

class AbstractProfile():
    """ Abstract profile class"""
    def __init__(self,zMin=0,zMax=300,nZ=101):
        """ Provides standard Lux profiles """
        self.parameters = dict()
        self.parameters['zCoo'] = numpy.linspace(zMin,zMax,nZ)
        self.values = None
        self.description = None
        
    def write(self):
        """ Writes the provided profile so that it can be added to the OF dict """
        if type(self.values) == type(None):
            return
        name = self.description.split('.')[0]
        strOutList = [name +' (']
        for iL in range(self.values.shape[0]):
            strOutList.append('(' + fmt(self.values[iL,0]) + ' ' + fmt(self.values[iL,1]) + ' ' + fmt(1.0) + ')')
        strOutList.append(');')
        return '\n'.join(strOutList)
    

class UProfile(AbstractProfile):
    """ Class to create mean(u) profiles """
    
    def uniformProfile(self,U):
        """ Constant profile """
        uu = numpy.ones(self.parameters['zCoo'].shape)*U
        prof = numpy.array([self.parameters['zCoo'],uu]).transpose()
        description = 'uProfile.uniformProfile(U='+str(U)+')'
        self.values = prof
        self.description = description
        return self
        
    def linearProfile(self,U0,dUdz):
        """ Constant profile """
        uu = self.parameters['zCoo']*dUdz + U0
        prof = numpy.array([self.parameters['zCoo'],uu]).transpose()
        description = 'uProfile.linearProfile(U0='+str(U0)+',dUdz='+str(dUdz)+')'
        self.values = prof
        self.description = description
        return self
    
    def powerLaw(self,U0,zRef,alpha,scale):
        """ Power-law profile"""
        uu = U0*numpy.power((self.parameters['zCoo']/zRef),alpha)
        uMin = uu.max()/20.0
        uu[uu<uMin] = uMin
        prof = numpy.array([self.parameters['zCoo'],uu]).transpose()
        description = 'uProfile.powerLaw(U0='+str(U0)+',zRef='+str(zRef)+',alpha='+str(alpha)+')'
        self.values = prof
        self.description = description
        return self
    
    def ECProfile(self,Ubas,z0,kt,zMin,scale):
        """ General EC profiles """
        uu = Ubas*kt*numpy.log(self.parameters['zCoo']*scale/z0)
        uZmin = Ubas*kt*numpy.log(zMin/z0)
        idMin = numpy.where((self.parameters['zCoo']*scale)<zMin)
        uu[idMin] = uZmin  
        prof = numpy.array([self.parameters['zCoo'],uu]).transpose()
        description = 'uProfile.ECProfile(Ubas='+str(Ubas)+',z0='+str(z0)+',kt='+str(kt)+',zMin='+str(zMin)+',scale='+str(scale)+')'
        self.values = prof
        self.description = description
        return self
    
    def AIJProfile(self,U,zb,zg,alpha,scale):
        """ General AIJ profiles """
        uu = U*1.7*numpy.power(self.parameters['zCoo']*scale/zg,alpha)
        uZmin = U*1.7*numpy.power(zb/zg,alpha)
        idMin = numpy.where((self.parameters['zCoo']*scale)<zb)
        uu[idMin] = uZmin
        prof = numpy.array([self.parameters['zCoo'],uu]).transpose()
        description = 'uProfile.AIJProfile(U='+str(U)+',zb='+str(zb)+',zg='+str(zg)+',alpha='+str(alpha)+',scale='+str(scale)+')'
        self.values = prof
        self.description = description
        return self
    
    
class SProfile(AbstractProfile):
    """ Class to create std(u) profiles """

    def uniformProfile(self,su):
        """ Constant profile """
        su = numpy.ones(self.parameters['zCoo'].shape)*su
        prof = numpy.array([self.parameters['zCoo'],su]).transpose()
        description = 'sProfile.uniformProfile(su='+str(su)+')'
        self.values = prof
        self.description = description
        return self
    
    def ECProfile(self,Ubas,kt):
        """ S profile for EC (constant) """
        su = numpy.ones(self.parameters['zCoo'].shape)*Ubas*kt
        prof = numpy.array([self.parameters['zCoo'],su]).transpose()
        description = 'sProfile.ECProfile(Ubas='+str(Ubas)+',kt='+str(kt)+')'
        self.values = prof
        self.description = description
        return self
    
class LProfile(AbstractProfile):
    """ Class to create Lux profiles """
    
    def uniformProfile(self,L):
        """ Constant profile """
        ll = numpy.ones(self.parameters['zCoo'].shape)*L
        prof = numpy.array([self.parameters['zCoo'],ll]).transpose()
        description = 'lProfile.uniformProfile(L='+str(L)+')'
        self.values = prof
        self.description = description
        return self
    
    def solari1993Profile(self,z0,scale):
        """ Lux profile after Solari 1993 """
        lu = 300.0*numpy.power((self.parameters['zCoo']*scale/300.0),(0.46+0.074*numpy.log(z0)))
        lMin = lu.max()/20.0
        lu[lu<lMin] = lMin
        prof = numpy.array([self.parameters['zCoo'],lu/scale]).transpose()
        description = 'lProfile.solari1993Profile(z0='+str(z0)+',scale='+str(scale)+')'
        self.values = prof
        self.description = description
        return self
    
    def ESDUProfile(self,z0,scale):
        """ Lux profile after Solari 1993 """
        lu = 25.0*numpy.power((self.parameters['zCoo']*scale),0.35)*z0**(-0.063)
        lMin = lu.max()/20.0
        lu[lu<lMin] = lMin
        prof = numpy.array([self.parameters['zCoo'],lu/scale]).transpose()
        description = 'lProfile.ESDUProfile(z0='+str(z0)+',scale='+str(scale)+')'
        self.values = prof
        self.description = description
        return self
    
class TProfile(AbstractProfile):
    """ Class to create T profiles """
    
    # Overwriting AbstractProfile
    def __init__(self):
        """ Provides standard Lux profiles """
        self.values = None
        self.description = None
    
    def constantProfile(self):
        """ Constant profile """
        prof = numpy.array([[0,1],[1.0E6,1]])
        description = 'tProfile.constantProfile()'
        self.values = prof
        self.description = description
        return self
        
    def rampProfile(self,tRamp):
        """ Constant profile """
        prof = numpy.array([[0,0],[tRamp,1]])
        description = 'tProfile.rampProfile(tRamp='+str(tRamp)+')'
        self.values = prof
        self.description = description
        return self
    
    
