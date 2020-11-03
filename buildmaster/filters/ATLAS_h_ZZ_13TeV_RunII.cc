/*
Reference:
   [2004.03969]
   Inclusive and differential fiducial cross sections of the Higgs boson are 
   measured in the H→ZZ∗→4ℓH \rightarrow ZZ^{*} \rightarrow 4\ell H→ZZ∗→4ℓ 
   (ℓ=e,μ\ell = e,\mu ℓ=e,μ) decay channel. The results are based on 
   proton−proton collision data produced at the Large Hadron Collider at a 
   centre-of-mass energy of 13 TeV and recorded by the ATLAS detector from 2015 
   to 2018, equivalent to an integrated luminosity of 139 \hbox {fb}^{-1}. 
   Eur.Phys.J.C 80 (2020) 10, 942

The implementation is based on Tabs 19a (central values and ucnertainties) and 
19b (correlation matrix) of the relevant hepdata entry, see
https://www.hepdata.net/record/ins1790439
*/

#include "ATLAS_h_ZZ_13TeV_RunII.h"

void ATLAS_h_ZZ_13TeV_RunIIFilter::ReadData()
{
  fstream f1;
  fstream f2;

  //Central values and total uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_h_ZZ_13TeV_RunII/HEPData-ins1790439-v1-Figure_19a.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Correlations between the 10 data points
  stringstream datafile_corr("");
  datafile_corr << dataPath()
		<< "rawdata/ATLAS_h_ZZ_13TeV_RunII/HEPData-ins1790439-v1-Figure_19b.csv";
  f2.open(datafile_corr.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_corr.str() << endl;
      exit(-1);
    }

  //Read central values and total uncertainties
  string line;

  double* Sys = new double[fNData];
  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];

  for(int i=0; i<9; i++)
    getline(f1,line);

  for(int i=0; i<8; i++)
    getline(f2,line); 
  
  for(int i=0; i<fNData; i++)
    {
      double ptmin, ptmax;
      char comma;
      double statp, statm, systp, systm;
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 130000;
      lstream >> fKin1[i] >> comma
	      >> ptmin >> comma
	      >> ptmax >> comma
	      >> fData[i] >> comma
	      >> statp >> comma
	      >> statm >> comma
	      >> systp >> comma
	      >> systm;

      fStat[i] = (statp-statm)/2.;
      Sys[i]   = (systp-systm)/2.;
      fData[i] = fData[i] +  (statp+statm)/2. + (systp+systm)/2.;
	     
      corrmat[i] = new double[fNData];
      syscor[i]  = new double[fNData];
      getline(f2,line);
      getline(f2,line);

      for(int j=0; j<fNData; j++)
	{
	  string sdum;
	  getline(f2,line);
	  istringstream kstream(line);
	  kstream >> sdum >> comma >> sdum >> corrmat[i][j];
	}
      
      getline(f2,line);
      getline(f2,line);
      getline(f2,line);
      getline(f2,line);
      getline(f2,line);
    }

  getline(f1,line);
  getline(f1,line);

  double* SM_cov = new double[fNData];
  
  for(int i=0; i<fNData; i++)
    {
      double ddum;
      char comma;
      getline(f1,line);
      istringstream lstream(line);
      double SM_cv, SM_erp, SM_erm;
      lstream >> ddum >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> SM_cv >> comma
	      >> SM_erp >> comma
	      >> SM_erm;
      
      SM_cov[i] = (SM_erp - SM_erm)/2.;
      SM_cv = SM_cv + (SM_erp + SM_erm)/2.;
      cout << SM_cv << endl;
      
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
  f2.close();

}
