### Welcome to our seqomm project !

#### Disclaimer: work in progress, some parts of the code are missing...

<p align="center"> 
<img src="./images/seqomm.png" title="what is this?" alt="Front page
illustration" width=auto height="300">
</p>

### What is this?
"seqomm" stands for SEQuential Observable Moment Matching (OMM).
Given a computational model, the OMM method provides an estimation of its
uncertain parameters probability density function (PDF) from the knowledge of the
statistical moments of some observable quantities.
By default, the moment matching is applied to every observable degree of
freedom (DOF) of the model but it may become computationally challenging when
complex models are considered (say, 3-D PDE model with millions of DOFs).
To alleviate the computational cost of OMM, we propose an algorithm, referred to as Physical DOF Selection (PDS), that
selects the model DOFs where the moments are matched.

The repository is organized as follows:

 * ./data
    
 * ./src
 
   * ./omm
   * ./pds

### What do I need?
You will need the following tools for the present project:
    * C++ compile
    * Python version >= 2.7

You will also need the following external libraries/packages:
    * Eigen 3: C++ (header-only) library (Eigen 2 might work but hasn't been tested)
    * GSL/GSLCBLAS: GSL implementation (in C) of BLAS routines (included with most OS X/Linux distributions)
    * scikit-learn: machine-learning library implemented in Python

### What do I need to do?
There is not much to do if you want to run the demo. You only need to compile
the OMM code which is written in C++. We provide a template BASH script for the
compilation in `./compileOMM.sh`


### Try the demo
