import os, sys
import scipy as sp

caseName = sys.argv[1]
selDir = './outputs/'+caseName+'/DOFSelection'
j = int(sys.argv[2])
nk = int(sys.argv[3])
S = sp.genfromtxt(selDir+'/selList_'+str(j)+'.txt',skip_header=nk-1,max_rows=1).astype('int')
if (len(S.shape)==0): S = sp.array([S])
#print str(len(S))+' pysical points were selected.'
sp.savetxt('./data/'+caseName+'/selectedDOFs.txt',sp.reshape(S,(len(S),1)),fmt='%u')
