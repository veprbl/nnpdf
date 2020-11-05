/*
Reference:
   [18....]

MUST REWRITE THE DESCRIPTION OF THIS DATA SET


   Combined measurements of Higgs boson production and decay using up to
   80 fb-1 of proton–proton collision data at √s = 13 TeV collected with 
   the ATLAS experiment
   Phys. Rev. D101 012002

There are 16 signal strengths in the following order
ggF gamma gamma
    Z Z
    W W
    tau tau
VBF gamma gamma
    Z Z
    W W
    tau tau
    b bbar
VH  gamma gamma
    Z Z
    b bar
tth gamma gamma
    V V
    tau tau
    b b

The implementation is based on Fig.5 (including uncertainties) and on Fig.6 
(the correlation matrix for the total uncertainty).
*/

#include "ATLAS_ggF_13TeV_2015.h"

//1) Distribution differential in pTH
void ATLAS_ggF_13TeV_2015_pTHFilter::ReadData()
{
  fstream f1;

  //Central values and total uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_ggF_13TeV_2015/HEPData-ins1674946-v1-Table_1.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values and total uncertainties
  string line;

  for(int i=0; i<12; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData; i++)
    {
   
      getline(f1,line);
      istringstream lstream(line);
      fKin2[i] = 0.;
      fKin3[i] = 13000.;
      double bin1, bin2, statp, statm, systp, systm;
      char comma;
      lstream >> bin1 >> comma
	      >> bin2 >> comma
	      >> fData[i] >> comma
	      >> statp >> comma
	      >> statm >> comma
	      >> systp >> comma
	      >> systm;

      fKin1[i] = (bin2-bin1)/2.;
      fData[i] = fData[i] + (statp+statm)/2. + (systp+systm)/2.;
      fStat[i] = (statp-statm)/2.;
      fSys[i][0].add  = (systp-systm)/2.;
      fSys[i][0].mult = fSys[i][0].add*1e2/fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
    }
  
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }

  double* SM_cov = new double[fNData];
  
  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      double ddum, SM_cv, SM_erp, SM_erm;
      char comma;
      lstream >> ddum >> comma
	      >> ddum >> comma
	      >> SM_cv >> comma
	      >> SM_erp >> comma
	      >> SM_erm;

      SM_cv = SM_cv + (SM_erp+SM_erm)/2.;
      SM_cov[i] = (SM_erp-SM_erm)/2.;
      cout << SM_cv << endl;
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

//2) Distribution differential in yH
void ATLAS_ggF_13TeV_2015_yHFilter::ReadData()
{
  fstream f1;

  //Central values and total uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_ggF_13TeV_2015/HEPData-ins1674946-v1-Table_2.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values and total uncertainties
  string line;

  for(int i=0; i<12; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData; i++)
    {
   
      getline(f1,line);
      istringstream lstream(line);
      fKin2[i] = 0.;
      fKin3[i] = 13000.;
      double bin1, bin2, statp, statm, systp, systm;
      char comma;
      lstream >> bin1 >> comma
	      >> bin2 >> comma
	      >> fData[i] >> comma
	      >> statp >> comma
	      >> statm >> comma
	      >> systp >> comma
	      >> systm;

      fKin1[i] = (bin2-bin1)/2.;
      fData[i] = fData[i] + (statp+statm)/2. + (systp+systm)/2.;
      fStat[i] = (statp-statm)/2.;
      fSys[i][0].add  = (systp-systm)/2.;
      fSys[i][0].mult = fSys[i][0].add*1e2/fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
    }
  
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }

  double* SM_cov = new double[fNData];
  
  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      double ddum, SM_cv, SM_erp, SM_erm;
      char comma;
      lstream >> ddum >> comma
	      >> ddum >> comma
	      >> SM_cv >> comma
	      >> SM_erp >> comma
	      >> SM_erm;

      SM_cv = SM_cv + (SM_erp+SM_erm)/2.;
      SM_cov[i] = (SM_erp-SM_erm)/2.;
      cout << SM_cv << endl;
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
