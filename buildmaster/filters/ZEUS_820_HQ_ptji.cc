/*
Reference: hep-ex/0208037
Hepdata: https://www.hepdata.net/record/ins593409
Published in Phys.Lett.B 547 (2002) 164-180
DESY-HERA. Measurement of the inclusive jet differential cross sections in 
neutral current deep inelastic scattering of 27.5 GeV positrons from 820 GeV 
protons with boson virtualities Q**2 > 125 GeV. The data were taken during the 
1996-97 HERA running period and have an integrated luminsoity of 38.6 pb-1. 
Jets areidentified in the Breit frame using the longitudinal kT cluster 
algorithm.

The data is taken from Hepdata, specifically from Tabs. 4-5-6. The cross
section is differential in pT in bins of ET and Q2.
*/

#include "ZEUS_820_HQ_ptji.h"

void ZEUS_820_HQ_ptjiFilter::ReadData()
{
  const int ntab=3;

  const double Q2[6]={187.5, 375.0, 750.0, 1500.0, 3500.0, 7500.0};
  
  for(int i=0; i<ntab; i++)
    {
      int index = 4+i;
      fstream f1;
      stringstream datafile("");
      datafile << dataPath() << "rawdata/" << fSetName
	       << "/HEPData-ins593409-v1-Table_" << index << ".csv";
      f1.open(datafile.str().c_str(), ios::in);

      if (f1.fail())
	{
	  cerr << "Error opening data file " << datafile.str() << endl;
	  exit(-1);
	}
      
      string line;

      for(int j=0; j<17; j++)
	{
	  getline(f1,line);
	}

      for(int j=0; j<5; j++)
	{
	  int k = 10*i+j;
	  cout << i << "  " << j << "   " << k << endl;
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
	  fKin2[k] = Q2[2*i]; //Q2
	  fKin3[k] = 300.;    //GeV

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
	  
	}

      const int count=5;
      
      for(int j=0; j<9; j++)
	{
	  getline(f1,line);
	}

      for(int j=0; j<5; j++)
	{
	  int k = 10*i+j+count;
	  cout << i << "  " << j << "   " << k << endl;
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
	  fKin2[k] = Q2[2*i+1]; //Q2 
	  fKin3[k] = 300.;      //GeV

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
	  
	}
      
      f1.close();

    }

}
