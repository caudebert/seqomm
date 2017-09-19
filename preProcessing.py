import os, sys
sys.path.append('./src/utils')
import scipy as sp
from ioBin import *


def saveMoment(Gm, m, sd, var_G, dir, correction):
    z = sp.mean(Gm,0)
    # Noise correction assuming it is a zero-mean additive noise with known variance
    if (correction):
        if (m==0):   # moment order 1: do nothing
            pass
        elif (m==1): # moment order 2: remove E(eps^2)
            z -= sd**2
        elif (m==2): # moment order 3: remove 3.*E(g)*E(eps^2)
            z -= 3.*z*(sd**2)
        elif (m==3): # moment order 4: remove 6*E(g^2)*E(eps^2) + 3*E(sd^4)
            z -= (6.*var_G*(sd**2) + 3.*(sd**4))
    sp.savetxt(dir+'/moment'+str(m+1)+'.txt',z)
    return 0

#-----------------------------
# LOAD USER-DEFINED PARAMETERS
#-----------------------------
border = '#'+'-'*30
if (len(sys.argv) != 3):
    sys.exit("Missing arguments. There should be two arguments: caseName and mode.")
else:
    print '\n'+border+'-'*30, '\n#',' '*20, 'PRE-PROCESSING\n', border+'-'*30
caseName = sys.argv[1]
mode = int(sys.argv[2]) # mode=0: simulations only, mode=1: measurements only, mode=2:everything
dirName = './data/'+caseName
sys.path.append(dirName)
import userDefinedParameters as param
maxMom = param.maxMoments
corr = param.correction
noiseLevels = sp.array(param.SNRList)
#---
noiseStds = (8.*sp.exp(0.1*sp.log(10)*noiseLevels))**(-.5)
noiseNames = ['nn']+[str(int(n)) for n in noiseLevels]


# Pre-process simulations
if (mode==0 or mode==2):
    print border, '\n# Pre-process simulations \n', border
    os.system('mkdir -pv '+dirName+'/simulations/normalized')
    originalFile = dirName+'/simulations/data.bin'
    normalizedFile = dirName+'/simulations/normalized/data.bin'
    G = readLargeBin(originalFile)
    a, b = G.min(), G.max()
    G = (G-a)/(b-a)
    writeLargeBin(normalizedFile, G)
    os.system('cd '+ dirName+'/simulations/normalized && ln -sf ../collocation.bin . && cd - > /dev/null')
    sp.savetxt(dirName+'/simulations/normalized/normalization.txt', sp.array([a,b]))
    print '  Done.'
# Pre-process measurements
if (mode>=1):
    print border, '\n# Pre-process measurements \n', border
    if (mode==1):
        normFile =dirName+'/simulations/normalized/normalization.txt'
        if (not os.path.isfile(normFile)):
            sys.exit('Error. Missing normalization file. Make sure you normalized the simulations.')
        a, b = sp.loadtxt(normFile)
    
    measurementsFile = dirName+'/measurements/data.bin'
    G = readLargeBin(measurementsFile)
    N = G.shape[0]
    G = (G-a)/(b-a)

    # Zero-noise (should be for debugging only)
    dir = dirName+'/measurements/noise_'+noiseNames[0]
    os.system('mkdir -pv '+ dir)
    Gm = G.copy()
    for m in range(0,maxMom):
        if (m>0): Gm*= G
        sp.savetxt(dir+'/moment'+str(m+1)+'.txt',sp.mean(Gm,0))

    # Noisy measurements
    for k in range(1, len(noiseNames)):
        dir = dirName+'/measurements/noise_'+noiseNames[k]
        os.system('mkdir -pv '+ dir)
        sd = noiseStds[k-1]
        Y = G + sp.random.normal(0.,sd, G.shape)
        print "  SNR =", int(noiseLevels[k-1]),'dB: adding noise with std =',"{:.2e}".format(sd)
        Gm = Y.copy()
        var_G = 0.*Gm
        for m in range(0,maxMom):
            if (m>0): Gm*= Y
            if (m==1): var_G = sp.mean(Gm,0) # variance of G for kurtosis correction
            saveMoment(Gm, m, sd, var_G, dir, corr)
    print '  Done.'


print border+'-'*30+'\n'
