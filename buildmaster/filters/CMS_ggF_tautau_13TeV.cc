/*
Reference:
   [CMS PAS HIG-18-032]
   A measurement of the inclusive cross section for the production of a Higgs 
   boson decaying into a pair of tau leptons, H→ττ. The measurement is based on
   pp collision data collected by the CMS experiment in 2016 and 2017 
   corresponding to an integrated luminosity of 77.4 fb−1 at a center-of-mass 
   energy of 13 TeV.

There are 5 simplified template cross sections for gg->H.

The implementation is based on Figures 10-11 and on Table 9.
*/

#include "CMS_ggF_tautau_13TeV.h"

void CMS_ggF_tautau_13TeVFilter::ReadData()
{
  fstream f1;
  fstream f2;
  fstream f3;

  //Central values and total uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_ggF_tautau_13TeV/data.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Correlations between the 6 data points
  stringstream datafile_corr("");
  datafile_corr << dataPath()
		<< "rawdata/CMS_ggF_tautau_13TeV/corr.txt";
  f2.open(datafile_corr.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_corr.str() << endl;
      exit(-1);
    }

  //SM predictions
  stringstream datafile_SM("");
  datafile_SM << dataPath()
		<< "rawdata/CMS_ggF_tautau_13TeV/SMpredictions.txt";
  f3.open(datafile_SM.str().c_str(), ios::in);

  if (f3.fail())
    {
      cerr << "Error opening data file " << datafile_SM.str() << endl;
      exit(-1);
    }

  //Read central values and total uncertainties
  string line;

  double* Sys = new double[fNData];
  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];
  double* SM_cov = new double[fNData];
  
  for(int i=0; i<fNData; i++)
    {
   
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 0.;
      double sysp, sysm;
      lstream >> fData[i] >> sysp >> sysm;

      corrmat[i] = new double[fNData];
      syscor[i]  = new double[fNData];
      getline(f2,line);
      istringstream kstream(line);

      for(int j=0; j<fNData; j++)
	{
	  kstream >> corrmat[i][j];
	}
      
      getline(f3,line);
      istringstream jstream(line);
      double SS_cv, SS_er;
      jstream >> SS_cv >> SS_er;
      SM_cov[i] = SS_er;

      fData[i] = fData[i] * SS_cv;
      fStat[i] = 0.;
      Sys[i]   = (sysp-sysm)/2.* fData[i];  
    }

  //Generate covariance matrix from correlation matrix
  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  corrmat[i][j] = corrmat[i][j]*Sys[i]*Sys[j];
	}
    }

  //Generate artificial systematics from covariance matrix
  if(!genArtSys(fNData,corrmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNSys; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
    }
  
  //Generate SM covariance matrix
  for(int i=0; i< fNData; i++)
    {
      for (int j=0; j< fNData; j++)
	{
	  cout << SM_cov[i]*SM_cov[j] << "   ";
	}
      cout << endl;
    }
  
  f1.close();

}
