import os,sys
import scipy as sp
import scipy.linalg as lin
sys.path.append('./')
from utils_seqomm import *
from sklearn.cluster import AgglomerativeClustering

def readInputFile(inputFile):
  with open(inputFile,'r') as f: lines = f.readlines()
  caseName = lines[0][:-1]; print "|  caseName: ", caseName,
  N = int(lines[2]); print " | N = ", N,
  maxK = int(lines[4]); print " | maxK = ", maxK, " | ",
  ev_threshold = float(lines[5]); print "ev_threshold ", ev_threshold, " | ",
  vol = float(lines[6]); print "vol ", vol
  #glob_iter = int(lines[7]); print " | iter = ", glob_iter
  return caseName, N, maxK, ev_threshold, vol#, glob_iter

def readPDF(caseName, glob_iter, N, vol):
  saveDir = './outputs/'+caseName+'/PDFs'
  if (glob_iter==0):
    th = readLargeBin('./data/'+caseName+'/simulations/collocation.bin')
    dim = th.shape[1]
    rho = th[:,0]*0. + (1./vol)
    os.system('mkdir -pv '+saveDir)
    sp.savetxt(saveDir+'/pdf_0.txt',rho)
    rho = rho[:N]
  else:
    pdfName = saveDir+'/pdf_'+str(glob_iter)+'.txt'
    pdf = sp.loadtxt(pdfName)[:N,:]
    dim = pdf.shape[1]-1
    rho = pdf[:,dim]
  return rho, dim

def checkDerivatives(caseName):
  fileName = './outputs/'+caseName+'/derivatives/dgdx1.bin'
  der = os.path.isfile(fileName)
  if (der): return True
  is_valid=0
  while not is_valid:
    boo = raw_input("It appears the derivatives have not yet been computed. Should we compute them ? This is donce once and for all. (yes/no): ")
    if (boo=="yes" or boo=="no"): is_valid=1
    else: print " Please answer by yes or no."
  if (boo=="no"):
    sys.exit('Error. We cannot continue without computing the derivatives. The program stops.')
    return False
  return False


border = '#'+'-'*30
if (len(sys.argv) != 3):
    sys.exit("Missing arguments. There should be 2 arguments: caseName and global iteration.")
else:
    print '\n'+border+'-'*30, '\n#',' '*20, 'ACTIVE CLUSTERING\n', border+'-'*30
caseName = sys.argv[1]
glob_iter = int(sys.argv[2])
dirName = './data/'+caseName
sys.path.append(dirName)
import userDefinedParameters as param

N = param.numSimulationSamples
maxK = param.maxClusters
ev_threshold = param.traceThreshold
vol = param.stochVolume

paramList = [('case name', caseName), ('glob. iter.', glob_iter), ('max clusters', maxK), ('trace threshold', ev_threshold)]
print printDict(paramList,20)

# PDF at previous iteration
rho, dim = readPDF(caseName, glob_iter, N, vol)


# Load derivatives
der = checkDerivatives(caseName)
if (not der): os.system('python ./src/pds/computeDerivatives.py '+caseName)
print (' loading the derivatives... ='),;sys.stdout.flush()
fileName =  './outputs/'+caseName+'/derivatives/dgdx1.bin'
dgdx = readLargeBin(fileName)[:,:N]
nT = dgdx.shape[0]
DGDX = sp.zeros((nT,N,dim))
DGDX[:,:,0] = dgdx
for d in xrange(dim-1):
  print ('='),;sys.stdout.flush()
  fileName = './outputs/'+caseName+'/derivatives/dgdx' + str(d+2) + '.bin'
  dgdx = readLargeBin(fileName)[:,:N]
  DGDX[:,:,d+1] = dgdx
print ' done.'


trace = sp.zeros((nT,))
E = []
print (' looping on physical points... '), ;sys.stdout.flush()
for i in xrange(nT):
  if (i%100==0): print ('='),;sys.stdout.flush()
  dfdx = DGDX[i,:,:]
  C = (1./N)*sp.dot(rho*dfdx.T,dfdx)
  lam,ee = lin.eigh(C)
  trace[i] = sp.sum(lam)
  E.append(ee[:,sp.argmax(lam)])
print ' done.'

idx = sp.nonzero(trace >= (trace.max())*ev_threshold)[0]
numAS = len(idx)
print '  '+str(numAS) + '/' + str(nT) + ' physical DOFs selected after thresholding.'
maxK = min(maxK,numAS)
X = sp.array(E)[idx,:]
D = 1.-sp.absolute(sp.dot(X,X.T))

print ' agglomerative clustering:'
print '  - initializing clustering tree...',
est = AgglomerativeClustering(linkage="complete",affinity='precomputed',compute_full_tree=True,memory='./outputs/'+caseName+'/cache')
print ' done.'

S = []

print '  - agglomerative clustering for K = 1 to '+str(maxK)+'...',
for nK in range(1,maxK+1):
  #print nK
  est.set_params(n_clusters = nK)
  est.fit(D)
  labels = est.labels_

  selec = []
  for k in range(nK):
    idk = sp.nonzero(labels==k)[0]
    best_k = idx[idk[sp.argmax(trace[idx[idk]])]]
    selec.append(best_k)
  S.append(selec)
print ' done.'

selDir = './outputs/'+caseName+'/DOFSelection'
os.system('mkdir -pv '+selDir)
print ' write list of clusters...',
f = open(selDir+'/selList_'+str(glob_iter)+'.txt', 'w')
for s in S:
    f.writelines(["%u " % item for item in s])
    f.write("\n")
f.close()
print ' done.\n'
print border+'-'*30+'\n'
