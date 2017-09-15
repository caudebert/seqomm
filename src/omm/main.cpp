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
 \file main.cpp
 \authors J.-F. Gerbeau, D. Lombardi, E. Tixier
 \date 05/04/2017
 \brief
 */

#include "DE.hpp"

using namespace std;

int main(int argc, char ** argv) {


  if (argc != 2){puts("Error : Wrong number of arguments. Specify your case name."); return 1;}
  const string caseName = argv[1];  

  cout << "Estimating the PDF from your data using the OMM algorithm." << endl;
  const string dirOfLog = "./" + caseName;
  
  DE estimation;
  estimation.initialize(dirOfLog);
  estimation.readData();
  estimation.findSaddlePoint();
  estimation.saveResult2();
  estimation.computeStatistics();

  return 0;
}
