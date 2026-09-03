import numpy
import copy
from functools import partial
import scipy.optimize
import scipy.integrate
from matplotlib import pyplot

def fmt(val):
    """ Format  s number for output"""
    return "{:.3f}".format(val)

class PRFG3Turb():
    
    def __init__(self,parameters):
        """ Initializing flow """

        self.parameters = copy.deepcopy(parameters)
        
        # Making sure they are unitary
        self.parameters['Covariance'] = self.parameters['Covariance']/self.parameters['Covariance'][0,0]
        self.parameters['L'] = self.parameters['L']/self.parameters['L'][0,0]

        # Setting marginal spectra
        Sux = partial(self.vonKarman,Ls=self.parameters['L'][0,0])
        Suy = partial(self.vonKarman,Ls=self.parameters['L'][0,1])
        Suz = partial(self.vonKarman,Ls=self.parameters['L'][0,2])
            
        Svx = partial(self.vonKarman,Ls=self.parameters['L'][1,0])
        Svy = partial(self.vonKarman,Ls=self.parameters['L'][1,1])
        Svz = partial(self.vonKarman,Ls=self.parameters['L'][1,2])
            
        Swx = partial(self.vonKarman,Ls=self.parameters['L'][2,0])
        Swy = partial(self.vonKarman,Ls=self.parameters['L'][2,1])
        Swz = partial(self.vonKarman,Ls=self.parameters['L'][2,2])
        
        self.spectra = [[Sux,Suy,Suz],[Svx,Svy,Svz],[Swx,Swy,Swz]]
        self.description = self.writeDescription()
    
    def writeDescription(self):
        """ Writes the parameters used for generation"""
        strOutList = []
        strOutList.append('// Covariance = [' + " ".join([fmt(vv) for vv in self.parameters['Covariance'][0,:]]) + ']')
        strOutList.append('//              [' + " ".join([fmt(vv) for vv in self.parameters['Covariance'][1,:]]) + ']')
        strOutList.append('//              [' + " ".join([fmt(vv) for vv in self.parameters['Covariance'][2,:]]) + ']')
        strOutList.append('// L          = [' + " ".join([fmt(vv) for vv in self.parameters['L'][0,:]]) + ']')
        strOutList.append('//              [' + " ".join([fmt(vv) for vv in self.parameters['L'][1,:]]) + ']')
        strOutList.append('//              [' + " ".join([fmt(vv) for vv in self.parameters['L'][2,:]]) + ']')
        strOutList.append('// dEd = ' + str(self.parameters['dEd']))
        strOutList.append('// finalE = ' + str(self.parameters['finalE']))
        return "\n".join(strOutList)
           
    # ------- Main PRFG3 routines 
    def PRFG(self):
        """ Samples the 3d spectra """
        kVectorX = self.sampleSpectra_1d(self.spectra[0])
        kVectorY = self.sampleSpectra_1d(self.spectra[1])
        kVectorZ = self.sampleSpectra_1d(self.spectra[2])
        
        ESampling,KSampling,varVS = self.sampleSpectra_3d(kVectorX,kVectorY,kVectorZ)
        
        KSamplingCorr,inside = self.correctSampledSpectra(ESampling,KSampling)
        
        allp, allq, allk = self.findPQK(ESampling,KSamplingCorr)
        
        self.data = dict()
        self.data["Turbulent field"] = dict()
        self.data["Turbulent field"]["p"] = allp
        self.data["Turbulent field"]["q"] = allq
        self.data["Turbulent field"]["k"] = allk
        
        self.rescaleFlow()
        self.inflowMetrics()
        
    def findPQK(self,ESampling,KSamplingCorr):
        """ Finds p, q and k """
        nK = KSamplingCorr.shape[0]
        
        allp = numpy.zeros(KSamplingCorr.shape)
        allq = numpy.zeros(KSamplingCorr.shape)
        allk = numpy.zeros(KSamplingCorr.shape)
        
        for iK in range(nK):
            
            #rSign = numpy.sign(numpy.random.randn(3))
            #kVector = numpy.multiply(KSamplingCorr[iK,:],rSign )
            kVector = KSamplingCorr[iK,:]
            eVector = ESampling[iK,:]
            
            #xi = numpy.random.randn(3)
            #xi = numpy.random.uniform(low=-1, high=1, size=[3])
            r1 = numpy.random.uniform(low=0, high=2.0*numpy.pi )
            r2 = numpy.random.uniform(low=-numpy.pi/2.0,high=numpy.pi/2.0)
            xi = numpy.array([numpy.cos(r2)*numpy.cos(r1),numpy.cos(r2)*numpy.sin(r1),numpy.sin(r2)])
            
            dir1 = numpy.cross(kVector,xi)/numpy.linalg.norm(numpy.cross(kVector,xi))
            dir2 = numpy.cross(kVector,dir1)/numpy.linalg.norm(numpy.cross(kVector,dir1))
             
            d1 = dir1[0]
            d2 = dir1[1]
            d3 = dir1[2]
            
            v1 = dir2[0]
            v2 = dir2[1]
            v3 = dir2[2]
            
            E1 = eVector[0]
            E2 = eVector[1]
            E3 = eVector[2]
            
            num = (2.0*d1**2*E3 - 2.0*d3**2*E1)*(-2.0*d2**2*d1*v1+2.0*d1**2*d2*v2) - (2*d1**2*E2-2.0*d2**2*E1)*(-2.0*d3**2*d1*v1+2.0*d1**2*d3*v3)
            den = (-d3**2*v1**2+d1**2*v3**2)*(-2.0*d2**2*d1*v1+2.0*d1**2*d2*v2) +(d2**2*v1**2-d1**2*v2**2)*(-2.0*d3**2*d1*v1+2.0*d1**2*d3*v3)
            bq = numpy.lib.scimath.sqrt(num/den)
            aq = (2.0*d1**2*E2-2.0*d2**2*E1+bq**2*(d2**2*v1**2-d1**2*v2**2))/(bq*(-2.0*d2**2*d1*v1+2.0*d1**2*d2*v2))
            ap = numpy.lib.scimath.sqrt((2.0*E1-(aq*d1+bq*v1)**2) /(d1**2))
            
            alphasSol = [ap,aq,bq]
            
            p = alphasSol[0]*dir1
            q = alphasSol[1]*dir1 + alphasSol[2]*dir2
            
            # Checking if imag part is present
            if numpy.linalg.norm(numpy.imag(numpy.array(alphasSol)))>0.0:    
                funToMin = lambda par: numpy.linalg.norm(numpy.divide((numpy.power(par[0]*dir1,2) + numpy.power(par[1]*dir1+par[2]*dir2,2) - 2.0*eVector),2.0*eVector))
                res = scipy.optimize.minimize(funToMin,[0,0,0],tol=1.0E-6,method = 'Nelder-Mead')
                alphasSol = res.x
                p = alphasSol[0]*dir1
                q = alphasSol[1]*dir1 + alphasSol[2]*dir2
            
            # Randomizing signs
            rSign = numpy.sign(numpy.random.randn(3))
            allk[iK,:] = numpy.multiply(kVector,rSign)
            p = numpy.multiply(p,rSign)
            q = numpy.multiply(q,rSign)
            
            if iK%2:
                allp[iK,:] = p
                allq[iK,:] = q
            else:
                allp[iK,:] = q
                allq[iK,:] = p 
                
        return allp,allq,allk
    
    # ------- Routines used to sample spectra 
    @staticmethod
    def vonKarman(k,Ls):
        """ von Karman spectrum """
        S = 4.0*Ls/(numpy.power(1.0+70.8*numpy.power(k*Ls,2),5.0/6.0))                
        return S
    
    def sampleSpectra_1d(self,Svec):
        """ Defines the spectral sampling """
        dEp = 1.0/self.parameters['dEd']
        finalE = self.parameters['finalE']
        nS = len(Svec)
        # Initializing sampling
        kVector = [0.0]
        Etot = numpy.zeros(nS)
        # Sampling and integrating energy for all U components at the same time
        while True:
            k = kVector[-1]
            dk = numpy.array([dEp/Sk(k) for Sk in Svec]).min()#*0.5*(1.0+numpy.random.uniform(low=0.1, high=1.1))
            kk = numpy.linspace(k,k+dk,100)
            
            for iS in range(nS):
                Etot[iS]  = Etot[iS] + numpy.trapz(Svec[iS](kk),kk);
            kVector.append(k+dk)
            if Etot.min()>=numpy.power(finalE,1.0/3.0):
                break   
        return numpy.array(kVector)
        
    def sampleSpectra_3d(self,kVectorX,kVectorY,kVectorZ):
        """ Samples 3d spectra """
        dkVectorX = kVectorX[1:]-kVectorX[:-1]
        dkVectorY = kVectorY[1:]-kVectorY[:-1]
        dkVectorZ = kVectorZ[1:]-kVectorZ[:-1]
        
        kVectorX = kVectorX[:-1] + 0.5*dkVectorX
        kVectorY = kVectorY[:-1] + 0.5*dkVectorY
        kVectorZ = kVectorZ[:-1] + 0.5*dkVectorZ
        
        nKx = kVectorX.shape[0]
        nKy = kVectorY.shape[0]
        nKz = kVectorZ.shape[0]
        
        varV = numpy.zeros(3)
        
        # Gauss integration points (using 3x3 points for every cell in k-space)
        xlg = numpy.array([-0.774596669241484,0.0,0.774596669241484])
        wlg = numpy.array([0.555555555555555,0.888888888888889,0.555555555555555])
        
        ESampling = list()
        KSampling = list()
        
        # Starting integration cycle
        for iKx in range(nKx):
            dKx = dkVectorX[iKx]
            Kmx = kVectorX[iKx]
            for iKy in range(nKy):
                dKy = dkVectorY[iKy]
                Kmy = kVectorY[iKy]
                for iKz in range(nKz):
                    dKz = dkVectorZ[iKz]
                    Kmz = kVectorZ[iKz]
                    
                    Xg = Kmx + xlg*dKx/2.0
                    Yg = Kmy + xlg*dKy/2.0
                    Zg = Kmz + xlg*dKz/2.0

                    [XX,YY,ZZ] = numpy.meshgrid(Xg,Yg,Zg)
                    [XXw,YYw,ZZw] = numpy.meshgrid(wlg*(dKx/2.0),wlg*(dKy/2.0),wlg*(dKz/2.0))
                    
                    dE = numpy.zeros(3)
                    for iU in range(3):
                        Sux = self.spectra[iU][0](XX)
                        Suy = self.spectra[iU][1](YY)
                        Suz = self.spectra[iU][2](ZZ)
                        Su = self.parameters['Covariance'][iU,iU]*numpy.multiply(numpy.multiply(Sux,Suy),Suz)
                        wwM = numpy.multiply(numpy.multiply(XXw,YYw),ZZw)
                        dE[iU] = numpy.multiply(Su,wwM).sum()
                    varV = varV + dE
                    ESampling.append(dE)
                    drKx = (numpy.random.rand()-0.5)*dKx
                    drKy = (numpy.random.rand()-0.5)*dKy
                    drKz = (numpy.random.rand()-0.5)*dKz
                    KSampling.append([max(Kmx+drKx,drKx/2.0),max(Kmy+drKy,drKy/2.0),max(Kmz+drKz,drKz/2.0)])

        ESampling = numpy.array(ESampling)
        KSampling = numpy.array(KSampling)
        return ESampling,KSampling,varV


    # ------- Routines used to ensure realizability
    def checkDomain(self,k,E):
        """ Checks if the point is inside the existance domain """
        u = E[0]/(k[1]*k[2])**2
        v = E[1]/(k[0]*k[2])**2
        w = E[2]/(k[0]*k[1])**2
        resu = u**2 + v**2 + w**2 - 2.0*(u*v+u*w+v*w)
        return resu
    
    def existanceGrad(self,y,metT,eVect):
        """ Trajectory to arrive to existance domain """
        funVal = self.checkDomain(y,eVect)
        epsToUse = 100.0*(numpy.finfo(float).eps)
        
        if funVal<(-epsToUse):
            dydt = y*0.0
        else:
            ddy = numpy.linalg.norm(y)*1.0E-6
            
            grad1 = (self.checkDomain(y+numpy.array([1.0,0.0,0.0])*ddy,eVect)-funVal)/ddy
            grad2 = (self.checkDomain(y+numpy.array([0.0,1.0,0.0])*ddy,eVect)-funVal)/ddy
            grad3 = (self.checkDomain(y+numpy.array([0.0,0.0,1.0])*ddy,eVect)-funVal)/ddy
            
            dydt = -numpy.array([grad1,grad2,grad3])
            
            # Making orthogonal to energy
            dydt = numpy.linalg.inv(metT).dot(dydt)
            
        return dydt
    
    def correctSampledSpectra(self,ESampling,KSampling):
        """ Corrects in order to verify existance condition """ 
        KSamplingCorr = numpy.zeros(KSampling.shape)
        nK = KSampling.shape[0]
        nInside = 0.0
        
        for iK in range(nK):
            
            if iK in numpy.arange(1,nK,int(nK/10)):
                print(*'*',end='')
            kVec = KSampling[iK,:]
            eVec = ESampling[iK,:]
            
            kCorr = kVec/numpy.linalg.norm(KSampling[iK,:])
            eCorr = eVec/numpy.linalg.norm(ESampling[iK,:])
            
            conE = self.checkDomain(kCorr,eCorr)
            # Point is already inside
            if conE<-100.0*(numpy.finfo(float).eps):
                KSamplingCorr[iK,:] = KSampling[iK,:]
                nInside = nInside + 1
            # Point must be correctd
            else:
                
                wVec = 1.0/kVec
                metT = numpy.diag(wVec)

                fun = lambda t,kCorr,self: self.existanceGrad(kCorr,metT,eCorr)
                
                rOde = scipy.integrate.ode(fun).set_integrator('dopri5')
                rOde.set_initial_value(kCorr)
                rOde.set_f_params(self)
                
                kvNorm = numpy.sqrt(metT.dot(kVec).transpose().dot(metT.dot(kVec)))
                kCorrNorm = numpy.sqrt(metT.dot(kCorr).transpose().dot(metT.dot(kCorr)))
                
                dt = 0.3*kCorrNorm/100.0
                epsToUse = 1000.0*(numpy.finfo(float).eps) 
                iC = 0    
                # Iteratively correcting up to realizability
                while (rOde.successful() and iC<100 and self.checkDomain(kCorr,eCorr)>(-10.0*epsToUse)): 
                    kCorr = rOde.integrate(rOde.t+dt)
                    iC = iC+1
                
                kCorrNorm = numpy.sqrt(metT.dot(kCorr).transpose().dot(metT.dot(kCorr)))
                
                kMult = kvNorm/kCorrNorm
                kCorr = kMult*kCorr
                
                conE = self.checkDomain(kCorr,eVec)
                KSamplingCorr[iK,:] = kCorr

                if conE<0.0:
                    nInside = nInside + 1
        
        print('\nDone!')
        return KSamplingCorr,nInside/nK
    
    
    # ------- Routines used for final rescaling
    def rescaleFlow(self):
        """ Rescale flow properties """
        # Rescaling variance    
        metricDict = self.inflowMetrics()
        allp = self.data["Turbulent field"]["p"]
        allq = self.data["Turbulent field"]["q"]
        
        # Rescaling covariance, it compensates only for very small sampling errors
        multV = numpy.average(numpy.sqrt(numpy.diag(metricDict['Covariance Target'])/numpy.diag(metricDict['Covariance'])))
        allp = allp*multV
        allq = allq*multV 
        self.data["Turbulent field"]["p"] = allp 
        self.data["Turbulent field"]["q"] = allq 
        
        # Rescaling integral scales, partially compensating for k-corrections
        metricDict = self.inflowMetrics()
        # Here it is possible to define different criteria
        LMult1 = numpy.average(metricDict['L']/metricDict['L Target'])
        LMult2 = metricDict['L'][0,0]/metricDict['L Target'][0,0]
        LMult = (LMult1+2.0*LMult2)/3.0
        # Applying rescaling
        self.data["Turbulent field"]["k"] = self.data["Turbulent field"]["k"]*LMult
        
    def getIntegralScale(self,nu,nd,Lext,plot=False):
        """ Calculates the integral scale of u_nu along nd"""
        nu = numpy.array(nu)/numpy.linalg.norm(nu)
        nd = numpy.array(nd)/numpy.linalg.norm(nd)
        eVec = 0.5*(numpy.power(self.data["Turbulent field"]["p"].dot(nu),2) + numpy.power(self.data["Turbulent field"]["q"].dot(nu),2))
        kVec = self.data["Turbulent field"]["k"].dot(nd)
        Ls = self.getIntegralScale1d(eVec,kVec,Lext,plot=plot)
        return Ls
           
    def getIntegralScale1d(self,eVec,kVec,Lext,plot=False):
        """ Calculates the integral scale """
        xVec = numpy.linspace(0.0,15.0*Lext,600)
        corrVec = xVec*0.0
        nK = kVec.shape[0]
        varTot =  eVec.sum()
        for iK in range(nK):
            corrVec = corrVec + eVec[iK]*numpy.cos(2.0*numpy.pi*kVec[iK]*xVec)
        corrVec = corrVec/varTot
        corrId = numpy.where(corrVec>=0.15)[0][-1]
        try:
            corrVecS = corrVec[:corrId]
            xVecS = xVec[:corrId]
        except:
            corrVecS = corrVec
            xVecS = xVec
        expDec = lambda LL: numpy.linalg.norm(numpy.exp(-xVecS/LL[0]) - corrVecS)
        LL = scipy.optimize.minimize(expDec,x0=[Lext,],tol=1.0E-3,method = 'Nelder-Mead').x[0]
        if plot:
            pyplot.figure()
            pyplot.plot(xVecS,corrVecS,'-b',linewidth=2)
            pyplot.plot(xVecS,numpy.exp(-xVecS/LL),'--r',linewidth=2)
        return LL
        
    def inflowMetrics(self):
        """ Measures inflow metrics"""
        allp = self.data["Turbulent field"]["p"]
        allq = self.data["Turbulent field"]["q"]
        allk = self.data["Turbulent field"]["k"]
        
        nx = numpy.array([1,0,0])
        ny = numpy.array([0,1,0])
        nz = numpy.array([0,0,1])
        Lux = self.getIntegralScale(nx,nx,self.parameters['L'][0,0])
        Luy = self.getIntegralScale(nx,ny,self.parameters['L'][0,1])
        Luz = self.getIntegralScale(nx,nz,self.parameters['L'][0,2])
        
        Lvx = self.getIntegralScale(ny,nx,self.parameters['L'][1,0])
        Lvy = self.getIntegralScale(ny,ny,self.parameters['L'][1,1])
        Lvz = self.getIntegralScale(ny,nz,self.parameters['L'][1,2])
        
        Lwx = self.getIntegralScale(nz,nx,self.parameters['L'][2,0])
        Lwy = self.getIntegralScale(nz,ny,self.parameters['L'][2,1])
        Lwz = self.getIntegralScale(nz,nz,self.parameters['L'][2,2])
        
        LMatr = numpy.array([[Lux,Luy,Luz],[Lvx,Lvy,Lvz],[Lwx,Lwy,Lwz]])
        LMatrT = self.parameters['L']
        

        covMatrix = numpy.zeros([3,3])
        for iW in range(allk.shape[0]):
            covMatrix = covMatrix + 0.5*(numpy.outer(allp[iW,:],allp[iW,:]) + numpy.outer(allq[iW,:],allq[iW,:]))
        
        inflowMetricsDict = dict()
        inflowMetricsDict['Covariance'] = covMatrix
        inflowMetricsDict['Covariance Target'] = copy.deepcopy(self.parameters['Covariance'])
        inflowMetricsDict['L'] = LMatr
        inflowMetricsDict['L Target'] = LMatrT
        
        self.data["Metrics"] = inflowMetricsDict
        
        return inflowMetricsDict

        
    def plotInflowMetrics(self):
        # Plotting inflow metrics
        inflowMetricsDict = self.data["Metrics"]
        print('------------- INFLOW METRICS --------------')
        print('Number of velocity waves: ' + str(self.data["Turbulent field"]["k"].shape[0]))
        print('-------------------------------------------')
        print('Target covariance:')
        print(inflowMetricsDict['Covariance Target'] )
        print('Obtained covariance:')
        print(inflowMetricsDict['Covariance'])
        print('Covariance % errors:')
        print(100.0*(inflowMetricsDict['Covariance Target']-inflowMetricsDict['Covariance'])/inflowMetricsDict['Covariance Target'])
        print('-------------------------------------------')        
        print('Target L:')
        print(inflowMetricsDict['L Target'])
        print('Obtained L:')
        print(inflowMetricsDict['L'])
        print('L % errors:')
        print(100.0*numpy.divide(inflowMetricsDict['L']-inflowMetricsDict['L Target'],inflowMetricsDict['L Target']))
        print('-------------------------------------------')        
 
      
# # Used for testing
# parameters = dict()
# parameters['Lu'] = numpy.array([1,0.6,0.5])
# parameters['Lv'] = numpy.array([0.6,0.5,0.5])
# parameters['Lw'] = numpy.array([0.5,0.5,0.3])
# parameters['Vars'] = numpy.array(numpy.power([1,0.75,0.5],2))
# parameters['dEd'] = 5
# parameters['finalE'] = 0.9
        
# thisFlow = PRFG3Turb(parameters)
# thisFlow.PRFG3()
# thisFlow.plotInflowMetrics()