/*
!!!!!!Data set description needed here !!!!!!!! 
*/

#include "ZEUS_DISJETS.h"

void ZEUS_DISJETSFilter::ReadData()
{
  const int ntab = 6;

  for(int i=0; i<ntab; i++)
    {
      int index = 5+i;
      fstream f1;
      stringstream datafile("");
      datafile << dataPath() << "rawdata/" << fSetName
	       << "/Table" << index << ".csv";
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
	  int k = 5*i+j;
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

	  fKin1[k] = 0.;
	  fKin2[k] = 0.;
	  fKin3[k] = 0.;

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
