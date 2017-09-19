import os, sys
import scipy as sp
import scipy.linalg as lin
from utils_seqomm import *

border = '#'+'-'*30
if (len(sys.argv) != 2):
    sys.exit("Missing arguments. There should be one arguments: caseName.")
else:
    print '\n'+border+'-'*30, '\n#',' '*20, 'SEQOMM\n', border+'-'*30
caseName = sys.argv[1]
dirPath = './data/'+caseName
savePath = './outputs/'+caseName
os.system('mkdir -pv '+savePath)
sys.path.append(dirPath)
import userDefinedParameters as param


N = param.numSimulationSamples
G, X = loadSimulations(dirPath+'/'+param.simDir, N)
dim = X.shape[1]
numMom = param.numMoments
measPath = param.measDir+'/noise_'+str(param.noiseLevel)
vol = param.stochVolume
#--- Generate OMM.in input file
generateOMMInput(dim, N, numMom, param.simDir, measPath, param.maxIterOMM, vol, param.alphaTol, param.pseudoInverseTol, dirPath)

# maxIterGlob = 5
# maxNx = 10
# caseName = sys.argv[1]
# dirPath = '../../../Data/'+caseName
# savePath = '/home/ROCQ/reo/etixier/local/Data'
# # TODO: get this from the input file DE.in
# dim = 2
# numMom = 1
# vol = 100.
# N = 512
# noiseLevel = 60
# askAtEveryIteration = False
# #-----------------------------------------

# M = loadMeasuredMoments(dirPath, numMom, noiseLevel)
# G, X = loadSimulations(dirPath, N)
# Gp = computePowersOfG(G, numMom)

# finalRho = []; finalResidual = []; S = []
# for iterGlob in range(maxIterGlob):
#     answer = promptForContinuation(dirPath, iterGlob, askAtEveryIteration)
#     if (not answer): break
#     rhoList = []; resNormList = []
#     for nx in range(1,maxNx+1):
#         os.system('python sgm/extractSelection.py '+dirPath+' '+str(iterGlob)+' '+str(nx))
#         os.system('./runOMM.sh '+dirPath)
#         try:
#             pdf = sp.loadtxt(savePath+'/pdf.txt')
#         except IOError:
#             sys.exit("pdf.txt file not found. ./omm/computePDF has probably failed...")
#             print "The program will now stop."
#             break
#         rho = pdf[:,dim]
#         R = computeGlobalResidual(Gp, M, rho, vol, N)
#         print 'iter =',iterGlob+1,'| nx =', nx, ' -- ||R|| =', format(lin.norm(R),'2.3e')
#         resNormList.append(lin.norm(R))
#         rhoList.append(rho)

#     iBest = sp.argmin(sp.array(resNormList))
#     print 'OMM sub-iterations are over. Saving the best PDF for next global iteration'
#     print '- Summary:'
#     print '-   best ever: nx = ',iBest+1,'||R|| =',resNormList[iBest]
#     S.append(iBest)
#     bestRho = rhoList[iBest]
#     savePDF(dirPath, iterGlob, bestRho, X)
#     deleteSelection(dirPath, iterGlob)

#     finalRho.append(bestRho)
#     finalResidual.append(resNormList[iBest])

# print finalResidual
# S = sp.hstack((sp.array(S),sp.array(finalResidual)))


# sp.savetxt('summary.out',S)
