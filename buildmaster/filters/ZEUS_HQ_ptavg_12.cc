/*
Reference: 1010.6167
Hepdata: https://www.hepdata.net/record/ins875006
Published in Eur.Phys.J.C 70 (2010) 965-982
DESY-HERA.  Double-differential dijet cross sections in neutral current deep 
inelastic ep scattering have been measured with the ZEUS detector using an 
integrated luminosity of 374 pb^-1. The measurement was performed at large 
values of the photon virtuality, Q^2, between 125 and 20000 GeV^2. The jets 
were reconstructed with the k_T cluster algorithm in the Breit reference frame 
and selected by requiring their transverse energies in the Breit frame, 
E_T,B^jet, to be larger than 8 GeV. In addition, the invariant mass of the 
dijet system, M_jj, was required to be greater than 20 GeV.

The data is taken from Hepdata, specifically from Tabs. 13-18. The cross
section is differential in pT in bins of ET and Q2.
*/

#include "ZEUS_HQ_ptavg_12.h"

void ZEUS_HQ_ptavg_12Filter::ReadData()
{
  const int ntab=6;

  const double Q2[6]={187.5, 375.0, 750.0, 1500.0, 3500.0, 10500.0};

  int count = 0;
  int k;
  
  for(int i=0; i<ntab; i++)
    {
      int index = 13+i;
      fstream f1;
      stringstream datafile("");
      datafile << dataPath() << "rawdata/" << fSetName
	       << "/HEPData-ins875006-v1-Table_" << index << ".csv";
      f1.open(datafile.str().c_str(), ios::in);

      if (f1.fail())
	{
	  cerr << "Error opening data file " << datafile.str() << endl;
	  exit(-1);
	}
      
      string line;

      for(int j=0; j<18; j++)
	{
	  getline(f1,line);
	}

      int p;
      
      if(i<4)
	p=4;
      else
	p=3;
      
      for(int j=0; j<p; j++)
	{
	  k = count;
	  getline(f1,line);
	  istringstream lstream(line);
	  char comma;
	  double Et, Etmin, Etmax, ddum;
	  double sys1p, sys1m, sys2p, sys2m;
	  double shift, delta;
	  lstream >> Et >> comma
		  >> Etmin >> comma
		  >> Etmax >> comma
		  >> fData[k] >> comma
		  >> fStat[k] >> comma
		  >> ddum  >> comma
		  >> sys1p >> comma
		  >> sys1m >> comma
		  >> sys2p >> comma
		  >> sys2m;

	  fKin1[k] = Et;
	  fKin2[k] = Q2[i]; //Q2
	  fKin3[k] = 318.;  //GeV

	  symmetriseErrors(sys1p, -sys1m, &shift, &delta);
	  fSys[k][0].add = delta;
	  fData[k] += shift;
	  fSys[k][0].type = ADD;
	  fSys[k][0].name = "UNCORR";

	  symmetriseErrors(sys2p, -sys2m, &shift, &delta);
	  fSys[k][1].add = delta;
	  fData[k] += shift;
	  fSys[k][1].type = MULT;
	  fSys[k][1].name = "CORR";

	  fSys[k][0].mult = fSys[k][0].add/fData[k] * 100;
	  fSys[k][1].mult = fSys[k][1].add/fData[k] * 100;

	  count++;
	}
      
      f1.close();

    }

}
