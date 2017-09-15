####### USER DEFINED PARAMETERS #########
# Test case: demo

# General settings
numSimulationSamples = 512
maxMoments = 3
numMoments = 1
stochVolume = 100. # stochastic volume, i.e. volume of the parameter domain: \int_\Theta 1*d\theta. @TODO: put this in user-defined quadrature weights for non-MC quadrature rules support. 
maxIterGlob = 3

# preProcessing settings
SNRList = [10., 20., 30., 40., 50., 60.] # signal-to-noise ratio list for synthetic measurements (in dB)

# OMM (observable moment matching) settings
maxIterOMM = 100 # Max iterations of Newton method
alphaTol = 1.e-3
pseudoInverseTol = 1.e-8

# PDS (physical DOF selection) settings
maxClusters = 62 # Maximum number of clusters. Must be <= total number of DOFs
traceThreshold = 1.e-2 # Threshold on the SGM (sensitivity Gram matrix) trace

# Options
correction = True # Activate measurements noise correction in the empirical moments. Works only if the noise is additive, zero-mean with known variance
normalizeData = True
