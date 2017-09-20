import os,sys
import itertools as ito
import scipy as sp
import scipy.linalg as lin
from sklearn.neighbors import KDTree
sys.path.append('./')
from utils_seqomm import *
import multiprocessing

def collMatrix(u, u0):
  a = 2./(sp.amax(u,0)-sp.amin(u,0))
  b = -(sp.amax(u,0)+sp.amin(u,0))/(sp.amax(u,0)-sp.amin(u,0))
  x = a*u+b; x0 = a*u0 + b # rescaling so that -1 <= x,x0 <= 1
  nx = x.shape[0]
  x = sp.hstack(( sp.ones((nx,1)),x )); x0 = sp.hstack((1.,x0))
  P = sp.zeros((nx, numCol))
  D = sp.zeros((numCol,dim))
  for col, t in enumerate(colList):
    P[:,col] = x[:,t[0]]*x[:,t[1]]
    for d in range(dim):
      if   (t.count(d+1)==0): D[col,d] = 0.
      elif (t.count(d+1)==1): D[col,d] = x0[t[1-t.index(d+1)]]
      elif (t.count(d+1)==2): D[col,d] = 2.*x0[d+1]
  D *= a # Because of the x = a*u+b rescaling
  return P,D

def mp_deriv(it):
  dfdx = sp.zeros((nC,dim))
  n2Err = sp.zeros((nC,))
  for i in range(nC):
    #print th[i,:].shape, dim
    x0 = sp.reshape(th[i,:],(1,dim))
    #print x0.shape
    dist,idx = tree.query(x0,k=K)
    x0 = sp.reshape(x0,(dim,))
    #print idx.shape, K
    idx = sp.reshape(idx,(K,))
    x = th[idx,:]
    y = G[idx,it]
    P,D = collMatrix(x,x0)
    q = lin.lstsq(P,y)[0]
    n2Err[i] = lin.norm(sp.dot(P,q) - y)/lin.norm(y)
    dfdx[i,:] = sp.dot(q,D)
  return (dfdx,n2Err)

border = '#'+'-'*30
if (len(sys.argv) != 2):
    sys.exit("Missing arguments. There should be one argument: caseName.")
else:
    print '\n'+border+'-'*30, '\n#',' '*10, 'COMPUTE DERIVATIVES (w.r.t. parameters)\n', border+'-'*30
caseName = sys.argv[1]
dirName = './data/'+caseName
sys.path.append(dirName)
import userDefinedParameters as param
nProcs = param.numProcs
nC = param.numSimulationSamples
neighbour_regularization = param.neighborReg
G, th = loadSimulations(dirName+'/'+param.simDir, nC, param.normalizeData)
dim = th.shape[1]
tree = KDTree(th)
nT = G.shape[1]

paramList = [('case name',caseName), ('stoch. dim.', dim), ('nb samples', nC), ('nb phys. DOFs', nT), ('nb procs', nProcs)]
print printDict(paramList,20)



colList = list(ito.combinations_with_replacement(range(dim+1),2))
numCol = len(colList)
K = numCol + neighbour_regularization # neighbour_regularization should be >=0

if __name__ == "__main__": 
  p = multiprocessing.Pool(nProcs)

print '\n  Computing the derivatives for '+str(nT)+' physical DOFs...',
sys.stdout.flush()
P = p.map(mp_deriv,range(nT))
print ' done.'

output = sp.asarray([p[0] for p in P])
n2Err = sp.asarray([p[1] for p in P])

saveDir = './outputs/' + caseName+'/derivatives'
os.system('mkdir -pv ' + saveDir)
for d in range(dim):
  fileName = saveDir+'/dgdx'+str(d+1)+'.bin'
  print '  Writing to ',fileName
  i = writeLargeBin(fileName,output[:,:,d])


sp.savetxt(saveDir+'/mean_fitting_error.txt',sp.mean(n2Err,1))
print border+'-'*30+'\n'
