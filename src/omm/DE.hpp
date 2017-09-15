
/*=============================================================================
 This file is part of the code DEMMA
 Density Estimation using the Moment Matching Approach
 Copyright (C) 2014, INRIA

 DEMMA is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by the Free
 Software Foundation; either version 2.1 of the License, or (at your option)
 any later version.

 DEMMA is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 more details.

 You should have received a copy of the GNU Lesser General Public License
 along with DEMMA. If not, see http://www.gnu.org/licenses/.
 =============================================================================*/

/*!
 \file DE.hpp
 \authors J.-F. Gerbeau, D. Lombardi, E. Tixier
 \date 05/04/2017
 \brief
 */


#ifndef DE_HPP
#define DE_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <unistd.h>
#include <locale.h>
#include <fstream>
#include <stdio.h>
#include <new>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <map>
#include "Eigen/Dense"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>


using namespace std;
using namespace Eigen;

// Mother class

class DE {

protected:
  unsigned int m_numStochPts;
  unsigned int m_numMoments;
  unsigned int m_numTimeSnap;
  unsigned int m_physDim;
  unsigned int m_stochDim;
  int m_numStochRec;
  int m_maxIter;
  string m_modelDir, m_modelName;
  string m_momentDir, m_momentName;
  string m_timeStepsDir;
  string m_extension;
  bool m_isInit;
  double m_tolInv;
  double m_prevRes;
  bool m_user_lambda;
  double m_volume;
  double m_tolEta;
  double m_alpha;

  double m_minVal, m_maxVal;
  vector<vector<double>*>* m_g;
  vector<vector<double>*> m_expMoments; 
  vector<vector<double>*> m_lagrangeMulti; 
  vector<vector<double>*> m_theta;    
  vector<double> m_rho;
  vector<unsigned int> m_selectedTimeSteps;
  VectorXd m_lambda;
  double m_mu;
  string m_fileFormat; 

public:
  
  // Constructor & Destructor:
  DE();
  ~DE();
  
  
  // Access functions:
//  inline unsigned int numStochPts(){return m_numStochPts;}
//  inline unsigned int numMoments(){return m_numMoment;}
//  inline unsigned int numTimeSnap(){return m_numTimeSnap;}

  // Functions
  bool isConverged(double);
  bool isStagnating(double, double);
  // Methods:
  void initialize(const string&);
  void readFromFile(string&, vector<double>&);
  void readFromText(string&, vector<double>&);
  void readData();
  void dataBinaryRead();
  void collocationBinaryRead();
  void readModelOutputs();
  void findSaddlePoint();
  void computePseudoInverse(MatrixXd, double, int, MatrixXd&, double&, int&);
  void computeSVD(MatrixXd, MatrixXd&, MatrixXd&, MatrixXd&);
  void saveResult2();
  void computeStatistics();
  void initializeLagrangeMultipliers(VectorXd&, double&);
  // Auxiliary
  void fillMatrices(MatrixXd&, VectorXd&);
  void fillMatrices2(MatrixXd&, VectorXd&);
  string i2s(int);
  int s2i(string);
  double s2d(string);
};


#endif

