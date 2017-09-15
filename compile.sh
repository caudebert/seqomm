#!/bin/bash
homeDir=/home/ROCQ/reo/etixier
cd src/omm
g++ -O3 -I $homeDir/.local/gsl-2.4/ -I $homeDir/.local/include/eigen3 main.cpp DE.cpp -o computePDF -lgsl -lgslcblas -lm
cd -
