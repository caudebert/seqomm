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
 \file DE.cpp
 \authors J.-F. Gerbeau, D. Lombardi, E. Tixier
 \date 05/04/2017
 \brief
 */


#include "DE.hpp"
using namespace std;
using namespace Eigen;

// 1 - Implementation of the mother class UQpp

// 1.0 Constructor & Destructor
DE::DE(){
  m_g       = NULL;
  m_isInit  = false;
  m_tolInv  = 0.;
  m_prevRes = 1.0e5;
  m_tolEta = 1.e-3;
  m_alpha = 0.;
}

DE::~DE(){

  if(m_g != NULL){
    for(size_t i=0; i<m_g->size(); i++){
      delete (*m_g)[i];
    }
    delete m_g;
  }

  for(size_t i=0; i<m_theta.size(); i++){
    if( m_theta[i] != NULL){ 
      delete ( m_theta[i] );
    }
  }

  for(size_t i=0; i<m_expMoments.size(); i++){
    if( m_expMoments[i] != NULL){ 
      delete ( m_expMoments[i] );
    }
  }

  for(size_t i=0; i<m_lagrangeMulti.size(); i++){
    if( m_lagrangeMulti[i] != NULL){ 
      delete ( m_lagrangeMulti[i] );
    }
  }

}



// 1.1 fInitialize 
void DE::initialize(const string& logDir){

  string tmpString;
  vector<string> theWholeFile;
  
  string fname = logDir + "/OMM.in";
  cout << endl;
  cout << "Opening the input file: " << fname <<endl;
  
  ifstream logFile;
  logFile.open(fname.c_str());  
  if (!logFile) {puts("File not found"); exit(1);}
 
  while(!logFile.eof()){
    getline(logFile,tmpString);
    theWholeFile.push_back(tmpString); 
  }
  logFile.close();
  
  cout << "File has been read." << endl;
  cout << endl;
  
  int tmpInt;

  tmpInt = s2i(theWholeFile[0]);
  m_stochDim = static_cast<unsigned int>(tmpInt); // Stochastic dimension = number of uncertain parameters

  tmpInt = s2i(theWholeFile[1]);
  m_numStochPts = static_cast<unsigned int>(tmpInt); // Number of collocation points (= quadrature points = parameter samples). Must be <= number of rows in data.bin
  
  tmpInt = s2i(theWholeFile[2]);
  m_numMoments = static_cast<unsigned int>(tmpInt); // Number of moments orders to be matched.
  
  m_modelDir      = logDir + "/" + theWholeFile[3]; // Directory of data.bin and collocation.bin
  m_modelName     = "g_"; // Deprecated
  m_momentDir     = logDir + "/" + theWholeFile[4]; // Directory of measurments
  m_momentName    = "moment"; // File prefix for measured moments (e.g. "moment")
  m_timeStepsDir  = logDir; // Where to get selectedTimeSteps.txt
  m_extension     = ".txt"; // Extension of data files in ascii format
  m_maxIter       = s2i(theWholeFile[5]); // Maximum number of iterations in the Newton method
  m_fileFormat    = "binary"; // Only binary for this release
  m_volume        = s2d(theWholeFile[6]); // Stochastic volume. Necessary if you want your PDF to integrate to 1.
  m_tolEta        = s2d(theWholeFile[7]); // Tolerance on the representation error
  m_alpha         = s2d(theWholeFile[8]); // Tolerance on the pseudo-inverse calculation
  if (m_fileFormat != "ascii" && m_fileFormat != "binary"){
    puts("Error: Only ascii and binary files are currently supported. Check your OMM.in file.");
  }
  // Initialization completed
  m_isInit = true;
}

// 1.2 Read experimental moments (from real or synthetized data)
void DE::readData(){

  puts( "Warning : Make sure outputs from the model are correctly normalized. Values should be ranging from 0 to 1.");

  cout << "Starting reading the experimental moments...";

  m_expMoments.resize(m_numMoments);
  for (unsigned int i=0; i<m_numMoments; i++){
    m_expMoments[i] = new vector<double>;
    string tmpFile = m_momentDir + "/" + m_momentName + i2s(i+1) + m_extension;
    readFromFile ( tmpFile, *(m_expMoments[i]) );
  }
  m_numTimeSnap = m_expMoments[0]->size();

  cout << " Done." << endl;
  //cout << " Note : Number of time steps = " << m_numTimeSnap << endl;


  cout << "Reading the stochastic collocation points...";
  // Read the parameters samples
  m_theta.resize(m_stochDim);
  m_rho.resize(m_numStochPts);  

  if (m_fileFormat=="ascii"){
    for (unsigned int i=0; i<m_stochDim; i++){
      m_theta[i] = new vector<double>;
      string tmpFile = m_modelDir + "/theta_" + i2s(i+1) + ".txt";
      readFromFile ( tmpFile, *(m_theta[i]) );
    }
  }
  else if(m_fileFormat=="binary"){
    collocationBinaryRead();
  }
  cout << " Done." << endl;
  cout << "Reading the outputs of the model...";
  // Read the corresponding outputs of the model
  m_g = new vector<vector<double>* >;
  m_g->resize(m_numStochPts);
 
  if (m_fileFormat=="ascii"){
    for(unsigned int j=0; j<m_numStochPts; j++){
      (*m_g)[j] = new vector<double>;
      string idSim = i2s(j);  
      string tmpFile = m_modelDir + "/" + m_modelName + idSim + m_extension;
      readFromFile ( tmpFile, *((*m_g)[j]) );
    }
  }
  else if(m_fileFormat=="binary"){
    dataBinaryRead();
  }

  cout << " Done." << endl;
  unsigned int Nt = (*m_g)[0]->size();

  //cout << " Note : Number of stochastic points = " << m_numStochPts << endl;
  
  if (Nt != m_numTimeSnap){
    puts("Error: Experimental moments and outputs of the model do not have the same length (time-wise). Please make sure the outputs of the model you generated are consistent with the experimental data.");
    exit(1);
  }

  string tmpFile = m_timeStepsDir + "/selectedTimeSteps.txt";
  ifstream timeFile(tmpFile.c_str());
  if (!timeFile){
    cout << "No file " << tmpFile << " was found. The whole signal will be taken into account.";
  }

//  vector<unsigned int> selectedTimeSteps;
  string line;
  if (timeFile){
    while(timeFile >> line){
      unsigned int t;
      istringstream parseStr(line);
      parseStr >> t;
      m_selectedTimeSteps.push_back(t);
    }
    m_numTimeSnap  = m_selectedTimeSteps.size();
    if (m_numTimeSnap>Nt){
      puts("Error: Number of selected time steps cannot exceed size of raw data. Program stops.");
      exit(1);
    }
  }
  else{
    for (int t=0; t<m_numTimeSnap; t++){
      m_selectedTimeSteps.push_back(t);
    }
  }

  cout << m_numTimeSnap << " time steps were selected." << endl;

  for (unsigned int m=0; m<m_numMoments; m++){
    vector<double> tmpVec = (*(m_expMoments[m]));
    for (int i=0; i<m_numTimeSnap; i++){
      int index = m_selectedTimeSteps[i];
      (*(m_expMoments[m]))[i] = tmpVec[index];
    }
    m_expMoments[m]->resize(m_numTimeSnap);
  }



  for (unsigned int j=0; j<m_numStochPts; j++){
    vector<double> tmpVec =  (*(*m_g)[j]);
    for (int i=0; i<m_numTimeSnap; i++){
      int index = m_selectedTimeSteps[i];
      (*(*m_g)[j])[i] = tmpVec[index];
    }
    (*m_g)[j]->resize(m_numTimeSnap);
  }

  m_lagrangeMulti.resize(m_numMoments);
  for (unsigned int i=0; i<m_numMoments; i++){
    m_lagrangeMulti[i] = new vector<double>;
    m_lagrangeMulti[i] -> resize(m_numTimeSnap);
  }

  cout << endl << "------------------------------" << endl;
  cout << "DATA SUMMARY:" << endl << "  Stoch. Dim. = " << m_stochDim << " | Nc = " << m_numStochPts << " | " << "Nx = " << m_numTimeSnap << " | Nm = " << m_numMoments <<  " | Hessian size = " << m_numMoments*m_numTimeSnap+1 << endl;
  cout << "------------------------------" << endl << endl;

}

//------------------------------------------------------------------------------------------------

// 1.3 Solve the inverse problem
void DE::findSaddlePoint(){

  double ti = clock();

  const int nT = m_numTimeSnap;
  const int nM = m_numMoments;
  const int nS = m_numStochPts;
  const double K = (1./nS); // Quadrature coefficient
  const int maxIter = m_maxIter; 
  cout << "Solving the inverse problem: " << endl;  

  VectorXd rho(nS), logrho(nS);
  int nH = nT*nM + 1;
  MatrixXd G (nS, nH-1);
  VectorXd expMom(nH-1);
  VectorXd lambda = VectorXd::Constant(nH-1, 0.); 
  double mu = 0.;

  MatrixXd Ufull= MatrixXd::Constant(nS,nH-1,0.); MatrixXd Vfull= MatrixXd::Constant(nH-1,nH-1,0.); 
  MatrixXd Sfull = MatrixXd::Constant(nH-1,nH-1,0.);
  MatrixXd U = MatrixXd::Constant(nS,nH-1,0.); MatrixXd V = MatrixXd::Constant(nH-1,nH-1,0.); 
  MatrixXd S = MatrixXd::Constant(nH-1,nH-1,0.);
  MatrixXd diagRhoU = MatrixXd::Constant(nS,nH-1,0.);

  fillMatrices(G,expMom);


  if (G.rows()>= G.cols()){
    computeSVD(G,Ufull,Vfull,Sfull);
  }
  else {
    MatrixXd Utmp= MatrixXd::Constant(nH-1,nS,0.); MatrixXd Vtmp= MatrixXd::Constant(nS,nS,0.); 
    MatrixXd Stmp = MatrixXd::Constant(nS,nS,0.);
    computeSVD(G.transpose(),Utmp,Vtmp,Stmp);
    cout << "OK" << endl;
    cout << Utmp.rows() << " " << Utmp.cols() << endl;
    cout << Vtmp.rows() << " " << Vtmp.cols() << endl;
    cout << Stmp.rows() << " " << Stmp.cols() << endl;
    Ufull.leftCols(nS) = Vtmp;
    Sfull.topLeftCorner(nS,nS) = Stmp;
    Vfull.leftCols(nS) = Utmp;
  }

/*
  cout << Ufull.block(0,0,10,10) << endl<< endl;
  cout << Vfull.block(0,0,10,10) << endl<< endl;
  cout << Sfull.block(0,0,10,10) << endl<< endl;
*/
  // Truncation procedure
  cout << scientific; cout.precision(3);
  // Truncation procedure
  int trunc = nH-1;
  for (int tr=1; tr<nH; tr++){ 
    U.leftCols(tr) = Ufull.leftCols(tr);
    V.leftCols(tr) = Vfull.leftCols(tr);
    S.topLeftCorner(tr,tr) = Sfull.topLeftCorner(tr,tr);
    double eta = ( expMom.transpose() - expMom.transpose() * V * V.transpose()).norm() / expMom.norm();
    //double eta = ( MatrixXd::Identity(nH-1,nH-1) - V * V.transpose() ).norm();
    //cout << tr << " " << eta << endl;
    if (eta < m_tolEta){
      trunc = tr;
      break;
    }
  }
  cout << "Truncation after " << trunc << " singular values." << endl;

  U.conservativeResize(nS,trunc); 
  V.conservativeResize(nH-1,trunc); 
  S.conservativeResize(trunc,trunc);

  MatrixXd H(trunc+1,trunc+1); MatrixXd H1(trunc,trunc); MatrixXd H2(trunc,1); MatrixXd H3(1,1);
  VectorXd x(trunc+1);
  MatrixXd P(trunc+1,trunc+1);
  VectorXd r(trunc); double s = 0.; VectorXd R(trunc+1);
  VectorXd beta = VectorXd::Constant(trunc,0.);


  //initializeLagrangeMultipliers(lambda, mu);

  //--------------------

  double prevRes = 1.e10;
  cout << scientific; cout.precision(3);
  for (int iter=0; iter<maxIter; iter++){ 

    logrho = (U*S*beta).array() + mu;
    for (int i=0; i<nS; i++){
      rho(i) = exp(logrho(i));

    }

    double entropy = -m_volume*K*(rho.transpose())*logrho;
    r = expMom.transpose()*V - K*rho.transpose()*U*S;
    //cout << r.rows() << " " << r.cols() << endl;
    s = 1. - K*rho.sum();
    R.head(trunc) = r; R(trunc) = s;
    //cout << s << " " << r.norm() << endl;
    cout << "Iter " << iter << ": ||R|| = " << R.norm();
    cout << " | Entropy = " << entropy;
    if (R.norm()<1.e-12){
      cout << " ||R||<tol. Stop." << endl;
      break;
    }
    else if ( isStagnating(prevRes,R.norm()) ){
      cout << " Residual is stagnating. Stop." << endl;
      break;
    }
    prevRes = R.norm();
    //H1 = -K * S * (U.transpose()) * diagrho * U * S;
    diagRhoU = (rho*MatrixXd::Constant(1,trunc, 1.));
    diagRhoU = (diagRhoU.array() * U.array());
    H1 = -K * S * (U.transpose()) * diagRhoU * S;
    //cout << H1.rows() << " " << H1.cols() << endl;

    H2 = -K * S * (U.transpose()) * rho;
    //cout << H2.rows() << " " << H2.cols() << endl;
    
    H3(0,0) = -K * rho.sum();
    H.topLeftCorner(trunc,trunc) = H1;
    H.topRightCorner(trunc,1) = H2;
    H.bottomLeftCorner(1,trunc) = H2.transpose();
    H.bottomRightCorner(1,1) = H3; 

    //cout << H.rows() << " " << H.cols() << endl;

    // Levenberg-Marquardt
    //MatrixXd D = MatrixXd::Identity(trunc+1,trunc+1); D(nH-1,nH-1) = 0.;
    //H += m_alpha*D;
    double cond=0.; int kount = 0;
    computePseudoInverse(H, m_alpha, trunc+1, P, cond, kount);
    //cout << P << endl;
    
    x = P*R;

    cout << " | |x|= " << x.norm() << " | cond= " << cond<< " | sigma= " << kount <<"/"<<trunc+1 << endl;
    //cout << "nH = " << nH << endl;
    beta -= x.head(trunc);
    mu   -= (x.tail(1))(0);

  }

  rho /= (K*rho.sum());
  for (int i=0; i<nS; i++){
    m_rho[i] = rho(i);
  }

  lambda = V * beta;
  m_lambda = lambda;
  m_mu = mu;

}

// 1.4.2 read Text Files
void DE::computePseudoInverse(MatrixXd H, double tolinv, int NH, MatrixXd& P, double& cond, int& kount){

//  ofstream debug("debug3.txt");

  gsl_matrix * gA_t = gsl_matrix_alloc (NH, NH);
  for (int i=0; i<NH; i++){
    for (int j=0; j<NH; j++){
      double val = H(j,i);
      gsl_matrix_set (gA_t, i, j, val); // gA_t is the transpose of H.
    }
//    debug << endl;
  }
	
  gsl_matrix * U = gsl_matrix_alloc (NH, NH);
  gsl_matrix * V= gsl_matrix_alloc (NH, NH);
  gsl_vector * S = gsl_vector_alloc (NH);


  // Computing the SVD of the transpose of A
  // The matrix 'gA_t' will contain 'U' after the function is called
  gsl_vector * work = gsl_vector_alloc (NH);
  gsl_linalg_SV_decomp (gA_t, V, S, work);
  gsl_vector_free(work);
  gsl_matrix_memcpy (U, gA_t);

  // Cut the lowest singular values
  gsl_matrix * SI = gsl_matrix_calloc (NH, NH);
  double s0 = gsl_vector_get (S, 0);
  double smin =   gsl_vector_get (S, NH-1);
  cond = s0 / smin;
  for (int i = 0; i < NH; i++) {
    double si = gsl_vector_get (S, i);
    if ((si/s0)>=tolinv){
      //cout << i << " " << si << endl;
      gsl_matrix_set (SI, i, i, 1.0/si);
      kount++;
    }
  }	
  //cout << " | Hessian cond. = \E[31;1m" << gsl_vector_get (S, 0)/gsl_vector_get (S, NH-1) << "\E[m";
	
  gsl_matrix * VT = gsl_matrix_alloc (NH,NH);
  gsl_matrix_transpose_memcpy (VT, V);					// Tranpose of V
		
		
  //THE PSEUDOINVERSE//
  //----------------------------------------------------------
  //Computation of the pseudoinverse of trans(A) as pinv(A) = U·inv(S).trans(V)   with trans(A) = U.S.trans(V)
  //----------------------------------------------------------
  gsl_matrix * SIpVT = gsl_matrix_alloc (NH, NH);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, SI, VT, 0.0, SIpVT);		// Calculating  inv(S).trans(V)
                 	
  gsl_matrix * pinv = gsl_matrix_alloc (NH,NH);	// Calculating  U·inv(S).trans(V)
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                	1.0, U, SIpVT,
                	0.0, pinv);

  for (int i=0; i<NH; i++){
    for (int j=0; j<NH; j++){
      double pij = gsl_matrix_get (pinv, i, j);
      P(i,j) = pij;
    }
  }

  gsl_matrix_free(pinv);
  gsl_matrix_free(VT);
  gsl_matrix_free(SI);
  gsl_matrix_free(SIpVT);
  gsl_matrix_free(gA_t);
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(S);

//  debug.close();

}
//----------------------------------------------------------
//------------------------------------------------------------------------------------------------

//----------------------------------------------------------
// 1.4.2 read Text Files
void DE::computeSVD(MatrixXd A, MatrixXd& Umat, MatrixXd& Vmat, MatrixXd& Smat){

//  WARNING: M>=N
  int M = A.rows();
  int N = A.cols();
  gsl_matrix * gA = gsl_matrix_alloc (M,N);
  for (int i=0; i<M; i++){
    for (int j=0; j<N; j++){
      double val = A(i,j);
      gsl_matrix_set (gA, i, j, val); //
    }
  }
	
  gsl_matrix * U = gsl_matrix_alloc (M,N);
  gsl_matrix * V= gsl_matrix_alloc (N,N);
  gsl_vector * S = gsl_vector_alloc (N);


  // Computing the SVD of the transpose of A
  // The matrix 'gA' will contain 'U' after the function is called
  gsl_vector * work = gsl_vector_alloc (N);
  gsl_linalg_SV_decomp (gA, V, S, work);
  gsl_vector_free(work);
  gsl_matrix_memcpy (U, gA);

  double val = 0.;
  for (int j=0; j<N; j++){
    Smat(j,j) = gsl_vector_get (S, j);
    for (int i=0; i<M; i++){
      Umat(i,j) = gsl_matrix_get (U, i, j);
    }
    for (int i=0; i<N; i++){
      Vmat(i,j) = gsl_matrix_get (V, i, j);
    }
  }



  gsl_matrix_free(gA);
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(S);

//  debug.close();

}

// 1.5 Save results in text files
void DE::saveResult2(){

  ofstream outPDF("./pdf.txt");

  double stochVolume = m_volume;

  outPDF.precision(15);
  outPDF << scientific;
  for (int s=0; s<m_numStochPts; s++){
    for (unsigned int dim=0; dim<m_stochDim; dim++){
      outPDF << (*(m_theta[dim]))[s] << " ";
    }
    outPDF << m_rho[s]/stochVolume << endl;
  }
  outPDF.close();

  // ofstream outLambda("lambda.txt");
  // outLambda.precision(15);
  // outLambda << scientific << m_lambda << endl << m_mu << endl;
  // outLambda.close();


}

// 1.6 Compute the parameters statistics (mean, sd, etc.)
void DE::computeStatistics(){

  const int nS = m_numStochPts;
  const int nD = m_stochDim;
  VectorXd rho(nS);
  MatrixXd theta(nS,nD);
  const int nM = max<int>(2,m_numMoments);
  MatrixXd moment(nM,nD);

  for (int s=0; s<nS; s++){
    for (unsigned int d=0; d<nD; d++){
      theta(s,d) = (*(m_theta[d]))[s];
    }
    rho(s) = m_rho[s];
  }

  MatrixXd multi = theta;
  for (int m=0; m<nM; m++){
    moment.row(m) = (1./nS)*(rho.transpose())*multi;
    multi = multi.array() * theta.array();
  }

  //ofstream outFile("statistics.txt"); outFile << scientific ; outFile.precision(12);
  cout << endl;
  cout << fixed << "ESTIMATED PARAMETER STATISTICS :" << endl;
  cout.precision(6);
  cout << "m1(mean) ";
  for (int m=1; m<nM; m++){cout << "m" << m+1 << "       ";}
  cout << "std.";
  cout << endl;
  for (int d=0; d<nD; d++){
    //double mean = moment(0,d);
    double sd   = sqrt( moment(1,d) - pow(moment(0,d),2) );
    //cout << "TH" << d+1 << ": mean = " << mean << ", sd = " << sd;
    for (int m=0; m<nM; m++){
      cout << moment(m,d) << " ";
    }
    cout << sd;
    cout << endl;
//    outFile << mean << " " << sd << endl;
  }
  cout << endl;
  //outFile << moment << endl;
  //outFile.close();
}

// Test if the residual is converged
bool DE::isConverged(double res){

  if (res != res) { 
    puts("Error: Residual has a NaN norm. Program stops. Check the tolerance parameters.");
    exit(1);
  }

  const double tol = 1.0e-4;
  //bool yes = ( abs((res-m_prevRes)/m_prevRes) < tol );//|| (res<1.0e-8);
  bool yes = ( abs((res-m_prevRes)/m_prevRes) < tol || (res<1.0e-8) );
  //bool yes = res<1.e-2;
  if (yes) {cout << endl << endl << "\E[32;1m Residual has converged or it is stagnating.\E[m" << endl;}

  return yes;
}

// Test if the residual is converged
bool DE::isStagnating(double old_res, double new_res){

  if (new_res != new_res) { 
    puts("Error: Residual has a NaN norm. Program stops. Check the tolerance parameters.");
    exit(1);
  }

  const double tol = 1.0e-4;
  //bool yes = ( abs((res-m_prevRes)/m_prevRes) < tol );//|| (res<1.0e-8);
  bool yes = ( abs((old_res-new_res)/new_res) < tol );
  //bool yes = res<1.e-2;
  if (yes) {cout << endl << endl << "\E[32;1m Residual is stagnating.\E[m" << endl;}

  return yes;
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// AUXILIARY

// Fill Eigen matrices and vectors
void DE::fillMatrices(MatrixXd& G, VectorXd& M){
  
  for (int s=0; s<m_numStochPts; s++){
    for (int m=0; m<m_numMoments; m++){
      for (int t=0; t<m_numTimeSnap; t++){
        int j = m*m_numTimeSnap + t;
        double gVal = (*(*m_g)[s])[t];
        G(s,j) = pow(gVal,m+1);
      }
    }
  }

  for (int m=0; m<m_numMoments; m++){
    for (int t=0; t<m_numTimeSnap; t++){
      int j = m*m_numTimeSnap + t;
      M(j) = (*(m_expMoments[m]))[t];
    }
  }

}

// Fill Eigen matrices and vectors
void DE::fillMatrices2(MatrixXd& G, VectorXd& M){
  
  for (int s=0; s<m_numStochPts; s++){
    for (int m=0; m<m_numMoments; m++){
      for (int t=0; t<m_numTimeSnap; t++){
        int j = m*m_numTimeSnap + t;
        double gVal = (*(*m_g)[s])[t];
        G(s,j) = pow(gVal,m+1);
      }
    }
    G(s,m_numMoments*m_numTimeSnap) = 1.;
  }

  for (int m=0; m<m_numMoments; m++){
    for (int t=0; t<m_numTimeSnap; t++){
      int j = m*m_numTimeSnap + t;
      M(j) = (*(m_expMoments[m]))[t];
    }
  }
  M(m_numMoments*m_numTimeSnap) = 1.;

}
//--------------------------------------------------------------------------------------------------
// Initialize Lagrange Multipliers
void DE::initializeLagrangeMultipliers(VectorXd& lambda, double& mu){
  int nx = m_numMoments*m_numTimeSnap + 1;

  vector<double> user_lambda;
  string userFile = "./user_lambda.txt";
  readFromFile ( userFile, user_lambda );
  if (user_lambda.size() != nx){
    cout << endl << "Error: user_lambda.txt has a wrong size." << endl;
    exit(1);
  }

  for (int i=0; i<nx-1; i++){
    lambda(i) = user_lambda[i];
  }
  mu = user_lambda[nx-1];

}
//--------------------------------------------------------------------------------------------------
// 1.4.2 read Text Files
void DE::readFromFile(string& fName, vector<double>& field){
  
  readFromText(fName, field);

}

// 1.4.2 read Text Files
void DE::dataBinaryRead(){
  string fileName = m_modelDir + "/data.bin";
  ifstream binaryFile (fileName.c_str(), ios::binary); 
  if (!binaryFile){puts("Error: could not open collocation file."); exit(1);}
  int dummy;
  binaryFile.read((char*)&dummy,sizeof(int));
  cout << endl<<"  checking binary file header: Nc = " << m_numStochPts << "/" << dummy;
  if (dummy < m_numStochPts){puts("Error in model outputs file: number of asked samples greater than number of available ones."); exit(1);}
  binaryFile.read((char*)&dummy,sizeof(int));
  cout << ", nT = " << dummy << ". ";
  if (dummy != m_numTimeSnap){puts("Error in model outputs file: number of physical points not consistent with experimental data."); exit(1);}
  //m_g = new vector<vector<double>* >;
  //m_g->resize(m_numSamples);
  for(unsigned int i=0; i<m_numStochPts; i++){
    (*m_g)[i] = new vector<double>;
    (*m_g)[i]->resize(m_numTimeSnap);
    binaryFile.read( (char*) &(*((*m_g)[i]))[0] , ((*m_g)[i]->size())*sizeof(double) );
  }
  binaryFile.close();
}

// 1.4.2 read Text Files
void DE::collocationBinaryRead(){
  string fileName = m_modelDir + "/collocation.bin";
  ifstream binaryFile (fileName.c_str(), ios::binary);
  if (!binaryFile){puts("Error: could not open collocation file."); exit(1);}
  int dummy;
  binaryFile.read((char*)&dummy,sizeof(int));
  cout << endl << "  checking binary file header : max_nC = " << dummy;
  if (dummy < m_numStochPts){puts("Error in collocation file: number of asked samples greater than number of available ones."); exit(1);}
  binaryFile.read((char*)&dummy,sizeof(int));
  cout << ", sD = " << dummy << ". ";
  if (dummy != m_stochDim){puts("Error in collocation file: number of columns not consistent with stochastic dimension."); exit(1);}
  for(unsigned int d=0; d<m_stochDim; d++){
    m_theta[d] = new vector<double>;
    m_theta[d]->resize(m_numStochPts);
  }
  double val;
  for(unsigned int i=0; i<m_numStochPts; i++){
    for(unsigned int d=0; d<m_stochDim; d++){
      binaryFile.read( (char*) &val , sizeof(double) );
      (*(m_theta[d]))[i] = val;
    }
  }
  binaryFile.close();
}


// 1.4.2 read Text Files
void DE::readFromText(string& fName, vector<double>& field){

  string tmpString;
  field.clear();
  
  ifstream textFile;
  textFile.open(fName.c_str());  
  if(!textFile){
    puts("Error in file name!");
    cout << fName << endl;
    puts("The above file was not found or could not be opened. The program will now stop.");
    cout << endl;
    exit(1);
  }
  
  int cc = 0;
  while(textFile >> tmpString){
    double tmpVal;
    istringstream parseStr(tmpString);
    parseStr >> tmpVal;
    field.push_back(tmpVal);
    cc = cc + 1;
  }
  textFile.close();
}

// int<=>string by using streams
string DE::i2s(int number){
  stringstream tmpS;
  tmpS << number;
  return tmpS.str();
}

int DE::s2i(string text){
  int value;
  istringstream (text) >> value;
  return value;
}

// string=>double by using streams
double DE::s2d(string text){
  double value;
  stringstream (text) >> value;
  return value;
}
