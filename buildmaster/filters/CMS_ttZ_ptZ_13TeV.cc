/*
Reference:
!!!!!THIS HEADER MUST BE WRITTEN!!!!!!!!!!
*/

#include "CMS_ttZ_ptZ_13TeV.h"

void CMS_ttZ_ptZ_13TeVFilter::ReadData()
{
  fstream f1;
  
  //Central values and total uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_ttZ_ptZ_13TeV/CMS-TOP-18-009-ptZ.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values and total uncertainties
  string line;

  for(int i=0; i<2; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData; i++)
    {
      char comma;
      double sys_up, sys_do;
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin3[i] = 13000.;
      lstream >> fKin2[i] >> comma
	      >> fData[i] >> comma
	      >> fStat[i] >> comma
	      >> sys_up >> comma
	      >> sys_do;

      double sys = (sys_up+sys_do)/2.;
      fData[i] = fData[i] + (sys_up-sys_do)/2.;
      
      fSys[i][0].add  = sqrt(pow(sys,2.)-pow(fStat[i],2.));
      fSys[i][0].mult = fSys[i][0].add/fData[i] * 100.;
      fSys[i][0].type = ADD;
      fSys[i][0].name = "UNCORR";
    }

  f1.close();

}
