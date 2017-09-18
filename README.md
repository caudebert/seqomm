### Welcome to our seqomm project !

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
