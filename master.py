import os, sys
import scipy as sp
import scipy.linalg as lin
from utils_seqomm import *

border = '#'+'-'*30
if (len(sys.argv) != 2):
    sys.exit("Missing arguments. There should be one argument: caseName.")
else:
    print '\n'+border+'-'*30, '\n#',' '*20, 'SEQOMM\n', border+'-'*30
caseName = sys.argv[1]
dirPath = './data/'+caseName
savePath = './outputs/'+caseName+'/PDFs'
os.system('mkdir -pv '+savePath)
sys.path.append(dirPath)
import userDefinedParameters as param


N = param.numSimulationSamples
G, X = loadSimulations(dirPath+'/'+param.simDir, N, param.normalizeData)
dim = X.shape[1]
numMom = param.numMoments
noiseLevel = param.noiseLevel
measPath = param.measDir+'/noise_'+str(noiseLevel)
vol = param.stochVolume
maxIterGlob = param.maxIterGlob
maxNx = param.maxClusters
#--- Generate OMM.in input file
generateOMMInput(dim, N, numMom, param.simDir+param.normalizeData*'/normalized', measPath, param.maxIterOMM, vol, param.alphaTol, param.pseudoInverseTol, dirPath)

askAtEveryIteration = False
#-----------------------------------------

M = loadMeasuredMoments(dirPath, numMom, noiseLevel)
Gp = computePowersOfG(G, numMom)

finalRho = []; finalResidual = []; S = []
for iterGlob in range(maxIterGlob):
    answer = promptForContinuation(caseName, iterGlob, askAtEveryIteration)
    if (not answer): break
    rhoList = []; resNormList = []
    for nx in range(1,maxNx+1):
        os.system('python ./src/pds/extractSelection.py '+caseName+' '+str(iterGlob)+' '+str(nx))
        os.system('./runOMMAlone.sh ' + caseName)
        try:
            pdf = sp.loadtxt('./pdf.txt')
        except IOError:
            sys.exit("pdf.txt file not found. ./omm/computePDF has probably failed...")
            print "The program will now stop."
            break
        rho = pdf[:,dim]
        R = computeGlobalResidual(Gp, M, rho, vol, N)
        print ' iter =',iterGlob+1,'| nx =', nx, ' -- ||R|| =', format(lin.norm(R),'2.2e')
        resNormList.append(lin.norm(R))
        rhoList.append(rho)

    iBest = sp.argmin(sp.array(resNormList))
    print 'OMM sub-iterations are over. Saving the best PDF for next global iteration'
    print '- Summary:'
    print '-   best ever: nx = ',iBest+1,'||R|| =', format(resNormList[iBest],'2.2e')
    S.append(iBest)
    bestRho = rhoList[iBest]
    savePDF(savePath, iterGlob, bestRho, X)
    deleteSelection(caseName, iterGlob)

    finalRho.append(bestRho)
    finalResidual.append(resNormList[iBest])

#print finalResidual
S = sp.hstack((sp.array(S),sp.array(finalResidual)))


sp.savetxt('summary.out',S)
