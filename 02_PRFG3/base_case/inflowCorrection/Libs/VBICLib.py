""" Python code to create VBIC method data structure for OF """
import numpy
import copy
import os
import io
import scipy.sparse
import pandas

version = '4.0u'

ofDictstart = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  6.0                                   |
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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"""

def readOFMesh(folder):
    """ Reads the mesh from OF """
    # Reading patch 
    with open(folder+os.sep+'boundary','r') as ff:
        boundary = ff.read()
    # Reading all faces
    print("Importing faces...")
    with open(folder+os.sep+'faces','r') as ff:
        faces = ff.read()
    facesData = faces.split("(\n")[1].split("\n)\n")[0]
    facesData = facesData.replace("3(","3(0 ").replace("("," ").replace(")","")
    facesDataStr = io.StringIO(facesData)
    facesDataStr.seek(0,0)
    try:
        allFaces = pandas.io.parsers.read_csv(facesDataStr,sep=' ',skiprows=0,header=None,skipinitialspace=True,usecols=[0,1,2,3,4],dtype=numpy.float).values
    except:
        print('Using Python parsing to read faces...')
        allFaces = pandas.io.parsers.read_csv(facesDataStr,sep=' ',skiprows=0,header=None,skipinitialspace=True,usecols=[0,1,2,3,4],dtype=numpy.float,engine='python').values

    # Reading all points
    print("Importing points...")
    with open(folder+os.sep+'points','r') as ff:
        points = ff.read()
    pointsData = points.split("(\n")[1].split("\n)\n")[0]
    pointsData = pointsData.replace("("," ").replace(")","")
    pointsDataStr = io.StringIO(pointsData)
    facesDataStr.seek(0,0)
    allPoints = pandas.io.parsers.read_csv(pointsDataStr,sep=' ',skiprows=0,header=None,skipinitialspace=True,dtype=numpy.float).values

    print("OF mesh was imported...")
    
    dictOut = dict()
    dictOut["boundary"] = boundary
    dictOut["faces"] = allFaces.astype(numpy.int)
    dictOut["points"] = allPoints
    return dictOut

def getPatch(ofMesh,patchName):
    """ Gets patch and related info """
    # Getting patch
    bs = ofMesh["boundary"].split(patchName)[1]
    bs = bs.split("}")[0].strip()
    nFaces = int(bs.split("nFaces")[1].split(";")[0])
    startFace = int(bs.split("startFace")[1].split(";")[0])
    # Setting patch local reference system
    allFaces = ofMesh["faces"]
    allPoints = ofMesh["points"]
    patchFaces = allFaces[startFace:startFace+nFaces,:].astype(int)
    l1 = allPoints[patchFaces[:,2]] - allPoints[patchFaces[:,1]]
    l2 = allPoints[patchFaces[:,-1]] - allPoints[patchFaces[:,1]]
    normals = -numpy.cross(l1,l2)
    normals = normals/numpy.linalg.norm(normals,axis=1)[:,None]
    nn = numpy.average(normals,axis=0)
    t1 = l1[0,:]/numpy.linalg.norm(l1[0,:])
    t2 = numpy.cross(nn,t1)
    print("Adopted patch local reference system (t1,t2,n) is...")
    print(str(t1)+ "    " + str(t2)+ "    " + str(nn))
    # Point used as base, no need to be in the centre
    pOri =  numpy.average(allPoints[allFaces[startFace:startFace+nFaces,1].astype(int)],axis=0)
    
    print("Transforming points to local reference system")
    pointTrans = numpy.zeros([allPoints.shape[0],2])
    pointTrans[:,0] = (allPoints-pOri).dot(t1)
    pointTrans[:,1] = (allPoints-pOri).dot(t2)
    
    dictOut = dict()
    dictOut['local ref'] = dict()
    dictOut['faces'] = patchFaces
    dictOut['points'] = allPoints
    dictOut['local ref']['t1'] = t1
    dictOut['local ref']['t2'] = t2
    dictOut['local points'] = pointTrans
    
    return dictOut


def shapeFun(xi,eta):
    """ Shape functions """
    N = numpy.array([0.25*(1.-xi)*(1.-eta), 0.25*(1.+xi)*(1.-eta), 0.25*(1.+xi)*(1.+eta), 0.25*(1.-xi)*(1.+eta)])
    dN = numpy.array([[-0.25+0.25*eta, -0.25+0.25*xi], [0.25-0.25*eta, -0.25-0.25*xi], [0.25+0.25*eta, 0.25+0.25*xi], [-0.25-0.25*eta, 0.25-0.25*xi]])
    return [N,dN]

def gaussRule(npg):
    """ Gauss integration """
    if npg==1:
        xpg = [0.0]
        wpg = [2.0]
    if npg==2:
        xpg = [-0.5773502691896258,0.5773502691896258]
        wpg = [1.0,1.0]
    elif npg==3:
        xpg = [-0.7745966692414834,0.0,0.7745966692414834]
        wpg = [5.0/9.0,8.0/9.0,5.0/9.0]
    return [xpg,wpg]

def getLocalGradientTrans(coo):
    """ Gets matrices to calculate gradients from dofs"""
    coo = copy.deepcopy(coo)
    # Building ref
    oo = coo[0,:]
    localCoo = copy.deepcopy(coo)
    localCooX = localCoo[:,0] - oo[0]
    localCooY = localCoo[:,1] - oo[1]
    
    NdN = shapeFun(0.0,0.0)
    dN = NdN[1]
    
    dx = numpy.matmul(localCooX,dN)
    dy = numpy.matmul(localCooY,dN)
    
    Jm = numpy.array([dx,dy]) # dx/d(xi)
    Jm1 = numpy.linalg.inv(Jm) # dxi/d(x)
    gradU = Jm1.dot(dN.T)
    return gradU

def assembleGradMatrix(nodes,faces):
    """ Assembles matrices for gradient calculation """
    nFaces = faces.shape[0]
    nNodes = nodes.shape[0]
    gradV1 = scipy.sparse.dok_matrix((nFaces,nNodes))
    gradV2 = scipy.sparse.dok_matrix((nFaces,nNodes))
    for iF in range(nFaces):
        thisF = faces[iF,:]
        coo = nodes[thisF,:]
        locGrad = getLocalGradientTrans(coo)
        gradV1[iF,thisF] = locGrad[0,:]
        gradV2[iF,thisF] = locGrad[1,:]
    return [gradV1,gradV2]

def getLocalKLaplace(coo):
    """ Assembling local stiffness matrix for Laplace operator """
    coo = copy.deepcopy(coo)
    
    # Building ref ATTENTION MAYBE REMOVE
    oo = coo[0,:]
    
    localCoo = copy.deepcopy(coo)
    localCooX = localCoo[:,0] - oo[0]
    localCooY = localCoo[:,1] - oo[1]
    
    npg = 2
    gpi = gaussRule(npg)
    xpg = gpi[0]
    wpg = gpi[1]
    
    Ke = numpy.zeros([4,4])
    
    for iX in range(npg):
        for iY in range(npg):
            
            NdN = shapeFun(xpg[iX],xpg[iY])
            dN = NdN[1]
            
            dx = numpy.matmul(localCooX,dN)
            dy = numpy.matmul(localCooY,dN)
            
            Jm = numpy.array([dx,dy]) # dx/d(xi)
            # Openfoam ordering is clockwise
            Jq = -numpy.linalg.det(Jm)
            Jm1 = numpy.linalg.inv(Jm) # dxi/d(x)
            
            B = numpy.matmul(dN,Jm1).T
            Ke = Ke + wpg[iX]*wpg[iY]*Jq*numpy.matmul(B.T,B)
    
    return Ke

def assembleLaplace(nodes,faces):
    """ Assembles the laplace operator """
    nFaces = faces.shape[0]
    nNodes = nodes.shape[0]
    globK = scipy.sparse.dok_matrix((nNodes,nNodes))
    for iF in range(nFaces):
        thisF = faces[iF,:]
        coo = nodes[thisF,:]
        locK = getLocalKLaplace(coo)
        globK[numpy.ix_(thisF,thisF)] += locK
    return globK

def SMatrix(matrix,order=1.0):
    """ Converts to sparse matrix for writing"""
    sMatrixCoo = numpy.array(scipy.sparse.find(matrix)).T
    idxOrd = numpy.argsort(order*(sMatrixCoo[:,0]*1.0E9+sMatrixCoo[:,1]))
    sMatrix = sMatrixCoo[idxOrd,:]
    return sMatrix

def smatrix2String(name,matrix,order=1.0):
    """ Converts smatrix to string """
    sMatrix = SMatrix(matrix,order)
    strOut = [name + ' (']
    fmt = '{:5e}'
    for ss in sMatrix:
        strOut.append('('+fmt.format(ss[0])+' '+fmt.format(ss[1])+' '+fmt.format(ss[2])+')')
    strOut.append(');')
    strOut = '\n'.join(strOut)
    return strOut

def matrix2String(name,matrix):
    """ Converts matrix to string """
    nR = matrix.shape[0]
    strOut = [name + ' (']
    fmt = '{:5e}'
    for iR in range(nR):
        thisRD = matrix[iR,:]
        thisR = [fmt.format(dd) for dd in thisRD]
        thisLine = "(" + " ".join(thisR) + ")"
        strOut.append(thisLine)
    strOut.append(');')
    strOut = '\n'.join(strOut)
    return strOut

def matrix2List(name,matrix):
    """ Converts matrix to list """
    nR = matrix.shape[0]
    nC = matrix.shape[1]
    matrixL = matrix.reshape([1,nR*nC])[0]
    nRL = nR*nC
    strOut = [name + ' (']
    fmt = '{:5e}'
    for iR in range(nRL):
        thisv = matrixL[iR]
        thisR = fmt.format(thisv)
        strOut.append(thisR)
    strOut.append(');')
    strOut = '\n'.join(strOut)
    return strOut
    
def renumberOFPatch(ofPatch):
    """ Changes OF numering """
    allFaces = ofPatch['faces']
    allPointsOri = ofPatch['points']
    allPoints = ofPatch['local points']
    faces = allFaces[:,1::] 
    allFacesNodes = numpy.sort(numpy.unique(faces))
    nodesOri = allPointsOri[allFacesNodes.astype(numpy.int),:]
    nodes = allPoints[allFacesNodes.astype(numpy.int),:]
    nUniNodes = nodes.shape[0]
    # Renumbering connectivity
    for iN in range(nUniNodes):
        thisN = allFacesNodes[iN]
        faces[faces==thisN]=iN
        
    # Creating model connectivity for reordering
    nNodes = nodes.shape[0]
    Kmod = scipy.sparse.dok_matrix((nNodes,nNodes))
    for iF in range(faces.shape[0]):
        Kmod[numpy.ix_(faces[iF],faces[iF])] = 1.0
    ckOrder = scipy.sparse.csgraph.reverse_cuthill_mckee(scipy.sparse.csc_matrix(Kmod),symmetric_mode=True)
    # Reordering
    rNodesOri = nodesOri[ckOrder,:]
    rNodes = nodes[ckOrder,:]
    rFaces = numpy.zeros(faces.shape).astype(int)
    for iN in range(nNodes):
        rFaces[faces==ckOrder[iN]]=iN
        
    ofPatch['Renumbered'] = dict()
    ofPatch['Renumbered']['faces'] = rFaces.astype(numpy.int)
    ofPatch['Renumbered']['nodes ori'] = rNodesOri
    ofPatch['Renumbered']['nodes'] = rNodes
    return ofPatch
    
def getEdges(faces):
    """ Gets the edge list """
    edges = numpy.zeros([0,4])
    nF = faces.shape[0]
    for iF in range(nF):
        thisface = faces[iF]
        e1 = [thisface[0],thisface[1]]
        e2 = [thisface[1],thisface[2]]
        e3 = [thisface[2],thisface[3]]
        e4 = [thisface[3],thisface[0]]
        theseEdges = [e1,e2,e3,e4]
        for iE in range(len(theseEdges)):
            thisEd = theseEdges[iE]
            id1 = edges[:,0]==thisEd[1]
            id2 = edges[:,1]==thisEd[0]
            id21 = numpy.logical_and(id1,id2).any()
            if id21:
                idEdge = numpy.where(numpy.multiply(id1,id2))[0][0]
                edges[idEdge,3] = iF 
            else:
                thisEdAdd = numpy.array([thisEd[0],thisEd[1],iF,-1])
                edges = numpy.vstack([edges,thisEdAdd])
            
    return edges.astype(int)    

def bcCalculationStructure(nodes,edges,faces):
    """ Prepare structure needed to calculate bc fluxes """
    nE = edges.shape[0]
    bcEdges = []
    bcNorm = []
    for iE in range(nE):
        thisEdge = edges[iE,:]
        if thisEdge[3]==-1:
            edgeFace = faces[thisEdge[2]]
            ll = nodes[thisEdge[1],:] - nodes[thisEdge[0],:] 
            eCentre = numpy.average(nodes[thisEdge[[0,1]],:],axis=0)
            fCentre = numpy.average(nodes[edgeFace,:],axis=0)
            nn = eCentre-fCentre
            nn = nn - nn.dot(ll)/numpy.linalg.norm(ll)*nn
            nn = nn/numpy.linalg.norm(nn)
            leng = numpy.linalg.norm(ll)
            bcEdges.append([thisEdge[0],thisEdge[1],leng])
            bcNorm.append([nn[0],nn[1],0.0])
    bcStruDict = dict()
    bcStruDict['bcEdges'] = numpy.array(bcEdges)
    bcStruDict['bcNorm'] = numpy.array(bcNorm)
    return bcStruDict

def plotEdges(nodes,edges):
    """ Plots the edges """
    import matplotlib.pyplot as plt
    fig = plt.gcf()
    ax = fig.gca(projection='3d')
    nE = edges.shape[0]
    for iE in range(nE):
        thisE = edges[iE,:].astype(int)
        xx = [nodes[thisE[0],0],nodes[thisE[1],0]]
        yy = [nodes[thisE[0],1],nodes[thisE[1],1]]
        zz = [nodes[thisE[0],2],nodes[thisE[1],2]]
        if thisE[3]==-1:
            ax.plot(xx,yy,zz,color='r')
        else:
            ax.plot(xx,yy,zz,color='k')
            
def writeOFDict(fname,inList):
    """ Writes all data to OF """
    strOut = [ofDictstart]
    corrLevel =   'corr 1.0;            // 0.0 - no correction, 1.0 - correction'
    strOut.append(corrLevel)

    for ss in inList:
            strOut.append(ss)
    with open(fname,'w') as ff:
        ff.write("\n".join(strOut))