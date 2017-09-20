import os, sys
import scipy as sp
sys.path.append('./src/utils')
from ioBin import *
import inspect

def printDict(dict, pad):
    s = ''
    for tup in dict:
        s += ' - '+tup[0].ljust(pad)+': '+str(tup[1])+'\n'
    return s

def generateOMMInput(dim, N, m, simDir, measDir, maxIter, vol, alpha, tol, dirPath):
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    #print 'function name "%s"' % inspect.getframeinfo(frame)[2]
    contents = []
    for i in args:
        print "    %s = %s" % (i, values[i])
        contents.append( str(values[i])+'\n')
    fileName = dirPath+'/OMM.in'
    with open(fileName, 'w') as file:
        file.writelines(contents[:-1])
    return 0

def promptForContinuation(dirPath, iterGlob, ask):
    sgm = checkSGM(dirPath, iterGlob)
    if (sgm): return True
    if (not ask):
        os.system('python ./sgm/activeClustering.py '+dirPath+' '+str(iterGlob))
        return True
    is_valid=0
    while not is_valid:
        boo = raw_input("Start a new global iteration? The SGM matrix and physical points clusters need to be recomputed. (yes/no): ")
        if (boo=="yes" or boo=="no"): is_valid=1
        else: print " Please answer by yes or no."
    if (boo=="no"):
        print "This was the last global iteration. The program will now stop."
        return False
    os.system('python ./sgm/activeClustering.py '+dirPath+' '+str(iterGlob))
    return True

def checkSGM(dirPath, iterGlob):
    fileName = dirPath+'/sgromm/selList_'+str(iterGlob)+'.txt'
    return os.path.isfile(fileName)

def loadMeasuredMoments(dirPath, numMom, noiseLevel):
    M = []
    for a in range(numMom):
        fileName = dirPath + '/measurements/noise_'+str(noiseLevel)+'/moment' + str(a+1) + '.txt'
        M.append(sp.loadtxt(fileName))
    return sp.vstack((M))

def loadSimulations(dirPath, N, normalized):
    X = readLargeBin(dirPath+'/collocation.bin')[:N,:]
    if (normalized): dirPath += '/normalized'
    G = readLargeBin(dirPath+'/data.bin')[:N,:]
    return G, X

def computePowersOfG(G, numMom):
    H = G.copy()
    Gp = [H.copy()]
    for m in range(1,numMom):
        H *= G
        Gp.append(H.copy())
    return Gp

def computeGlobalResidual(Gp, M, rho, vol, N):
    numMom = M.shape[0]
    #z = (vol/N)*sp.sum(rho) # PDF norm (should be 1)
    #appMom = sp.vstack([(vol/N)*sp.dot((G**(a+1)).T, rho) for a in range(numMom)])
    appMom = sp.vstack([(vol/N)*sp.dot((Gp[a]).T, rho) for a in range(numMom)])
    R = appMom - M
    return R

def savePDF(dirPath, iter_glob, rho, X):
    pdf = sp.vstack((X.T,rho)).T
    print pdf.shape
    fileName = dirPath+'/sgromm/pdf_'+str(iter_glob+1)+'.txt'
    sp.savetxt(fileName, pdf)

def deleteSelection(dirPath, iter_glob):
    if (iter_glob>0):
        fileName = dirPath+'/sgromm/selList_'+str(iter_glob)+'.txt'
        os.system('rm -fv '+fileName)
