""" Python code to create VBIC method data structure for OF """
import sys 
sys. path. append('Libs')
import os
import numpy
import scipy.sparse
import scipy.sparse.linalg
import VBICLib

""" INPUT HERE """
folder = os.getcwd() + os.sep + 'polyMesh'   # Path to folder containing the openFoam mesh 
patchName = 'inlet'                          # Name of the inflow patch
dictName = 'inflowCorrectionDict'               # Path of output file (VBICDict)
showMesh = False                          # Plot mesh?




"""-----------------------------------------------------------------------------------"""
""" ---------------------------- END OF USER EDITABLE PART ---------------------------"""
ofMesh = VBICLib.readOFMesh(folder)
ofPatch = VBICLib.getPatch(ofMesh,patchName)
ofPatch = VBICLib.renumberOFPatch(ofPatch)

nodesOri = ofPatch['Renumbered']['nodes ori']
nodes = ofPatch['Renumbered']['nodes']
faces = ofPatch['Renumbered']['faces']

# Laplace operator
print("Assembling laplace...")
laplaceK = VBICLib.assembleLaplace(nodes,faces)
print("Assembling gradients...")
gradMat = VBICLib.assembleGradMatrix(nodes,faces)
print("Getting edges...")
edges = VBICLib.getEdges(faces)

# Enforcing symmery (redundant)
K = 0.5*(laplaceK + laplaceK.T)
# Restraining K (Just setting level)
K[0,:] = 0.0
K[:,0] = 0.0
K[0,0] = 1.0

print("LU decomposition...")
LU = scipy.sparse.linalg.splu(scipy.sparse.csc_matrix(K))

Lo = LU.L - scipy.sparse.diags(LU.L.diagonal())
Ld = LU.L.diagonal()

Uo = LU.U - scipy.sparse.diags(LU.U.diagonal())
Ud = LU.U.diagonal()

print("Assembling BCs...")
bcStru = VBICLib.bcCalculationStructure(nodes,edges,faces)

# Calculating face centres
faceCentres = 0.25*(nodesOri[faces[:,0]] + nodesOri[faces[:,1]] + nodesOri[faces[:,2]] + nodesOri[faces[:,3]])

if showMesh:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    plt.figure()
    fig = plt.gcf()
    VBICLib.plotEdges(nodesOri,edges)
    ax = fig.gca(projection='3d')
    bcCoo = 0.5*(nodesOri[bcStru['bcEdges'][:,0].astype(int)] + nodesOri[bcStru['bcEdges'][:,1].astype(int)])
    bcNorm = numpy.outer(bcStru['bcNorm'][:,0],ofPatch['local ref']['t1']) + numpy.outer(bcStru['bcNorm'][:,1],ofPatch['local ref']['t2'])
    ax.scatter(bcCoo[:,0], bcCoo[:,1], bcCoo[:,2])
    ax.quiver(bcCoo[:,0], bcCoo[:,1], bcCoo[:,2], bcNorm[:,0], bcNorm[:,1], bcNorm[:,2], length=0.1, normalize=True)

print("Preparing to export in OF...")
dataToOF = []
fmt = '{:5e}'
v2Str = "t1 (" + fmt.format(ofPatch['local ref']['t1'][0]) + "  " + fmt.format(ofPatch['local ref']['t1'][1]) + "  " + fmt.format(ofPatch['local ref']['t1'][2]) + ");"
dataToOF.append(v2Str)
v3Str = "t2 (" + fmt.format(ofPatch['local ref']['t2'][0]) + "  " + fmt.format(ofPatch['local ref']['t2'][1]) + "  " + fmt.format(ofPatch['local ref']['t2'][2]) + ");"
dataToOF.append(v3Str)

dataToOF.append(VBICLib.smatrix2String('LoOp',Lo))
dataToOF.append(VBICLib.smatrix2String('UoOp',Uo,order=-1.0))
dataToOF.append(VBICLib.smatrix2String('grad1Op',gradMat[0]))
dataToOF.append(VBICLib.smatrix2String('grad2Op',gradMat[1]))

dataToOF.append(VBICLib.matrix2List('LdOp',numpy.array([Ld])))
dataToOF.append(VBICLib.matrix2List('UdOp',numpy.array([Ud])))

dataToOF.append(VBICLib.matrix2List('LUpc',numpy.array([LU.perm_c])))
dataToOF.append(VBICLib.matrix2List('LUpr',numpy.array([LU.perm_r])))

dataToOF.append(VBICLib.matrix2List('faces',faces))
dataToOF.append(VBICLib.matrix2String('nodes',nodesOri))
dataToOF.append(VBICLib.matrix2String('faceCentres',faceCentres))
dataToOF.append(VBICLib.matrix2String('bcEdges',bcStru['bcEdges']))
dataToOF.append(VBICLib.matrix2String('bcNorm',bcStru['bcNorm']))

print("Writing dict...")
VBICLib.writeOFDict(dictName,dataToOF)
print("Done!")