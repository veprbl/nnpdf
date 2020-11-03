/*
Reference:
   [CMS PAS HIG-18-029]
   Measurement of Higgs boson production via gluon fusion and vector boson 
   fusion in the diphoton decay channel at sqrt(s)=13 TeV

There are 6 simplified template cross sections for gg->H.

The implementation is based on Table 6 (including uncertainties).
*/

#include "CMS_ggF_aa_13TeV.h"

void CMS_ggF_aa_13TeVFilter::ReadData()
{
  fstream f1;

  //Central values and total uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_ggF_aa_13TeV/data.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values and total uncertainties
  string line;

  double* SM_cov = new double[fNData];
  
  for(int i=0; i<fNData; i++)
    {
   
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 0.;
      double statp, statm;
      double sysp, sysm;
      double theop, theom, SM_cv, ratio;
      lstream >> SM_cv >> fData[i] >> ratio
	      >> statp >> statm
	      >> sysp >> sysm
	      >> theop >> theom;

      statp = statp * SM_cv;
      statm = statm * SM_cv;
      sysp  = sysp  * SM_cv;
      sysm  = sysm  * SM_cv;
      theop = theop/ratio * SM_cv;
      theom = theom/ratio * SM_cv;

      fData[i] = fData[i] + (statp+statm)/2. + (sysp+sysm)/2.;
      fStat[i] = (statp-statm)/2.;
      fSys[i][0].add  = (sysp-sysm)/2.;
      fSys[i][0].mult = fSys[i][0].add*1e2/fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";

      SM_cov[i] = (theop-theom)/2.;
      SM_cv = SM_cv + (theop+theom)/2.;

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
