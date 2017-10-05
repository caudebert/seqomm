# seqomm
### Welcome to our seqomm project !

<p align="center"> 
<img src="./images/seqomm.png" title="what is this?" alt="Front page
illustration" width=auto height="270">
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
All the methods and algorithms are presented in [this preprint](https://hal.archives-ouvertes.fr/hal-01391254).

The repository is organized as follows:

 * `./data`: This where the user data are stored. A demo test case is already
   there. TO use your own data, read the `README.md` file.
    
 * `./src`: This is where the different sources are. Some are written in C++,
   others in Python and the ensemble is driven by the Python script
   `master.py` in the project root directory.
 
   * `./src/omm`: OMM sources, more info in `README.md` file
   * `./src/pds`: PDS sources, more info in `README.md` file
   * `./src/utils`: some utilities (I/O, formatting, ...)

### What do I need?
You will need the following tools for the present project:

   * C++ compiler
   * Python version >= 2.7

You will also need the following external libraries/packages:

   * [Eigen 3](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download): C++ (header-only) library (Eigen 2 might work but hasn't been tested)
   * GSL/GSLCBLAS: GSL implementation (in C) of BLAS routines (included with most OS X/Linux distributions)
   * [scikit-learn](http://scikit-learn.org/stable/install.html): machine-learning library implemented in Python

### What do I need to do?
There is not much to do if you want to run the demo. You only need to compile
the OMM code which is written in C++. We provide a template BASH script for the
compilation in `./compileOMM.sh`.

### Try the demo
Make sure you have compiled what needs to be compiled (see the above paragraph).
Then, execute the following commands in your console:

```console
python preProcessing.py demo 2
python master.py demo
```
