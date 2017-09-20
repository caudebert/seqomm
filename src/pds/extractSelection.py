import os, sys
import scipy as sp

caseName = sys.argv[1]
j = int(sys.argv[2])
nk = int(sys.argv[3])
S = sp.genfromtxt(caseName+'/sgromm/selList_'+str(j)+'.txt',skip_header=nk-1,max_rows=1).astype('int')
if (len(S.shape)==0): S = sp.array([S])
#print str(len(S))+' pysical points were selected.'
sp.savetxt(caseName+'/selectedTimeSteps.txt',sp.reshape(S,(len(S),1)),fmt='%u')
