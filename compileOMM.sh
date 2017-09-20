#!/bin/bash

cd src/omm
g++ -O3 -I /path/to/gsl/include -I /path/to/eigen/include main.cpp DE.cpp -o computePDF -L /path/to/gsl/lib -lgsl -lgslcblas -lm
cd -
