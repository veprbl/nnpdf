/*
Reference: 0911.5678
Hepdata: https://www.hepdata.net/record/ins838435
Published in Eur.Phys.J.C 67 (2010) 1-24
DESY-HERA.  The production of jets is studied in deep-inelastic e+p scattering 
at low negative four momentum transfer squared 5<Q2<100 GeV and at inelasticity 
0.2<y<0.7 using data recorded by the H1 detector at HERA in the years 1999 and 
2000, corresponding to an integrated luminosity of 43.5 pb-1. Inclusive jet and
2-jet cross sections are measured as a function of Q2 and jet transverse 
momentum. The 2-jet cross section is also measured as a function of the proton momentum fraction.

The data is taken from Hepdata, specifically from Tabs. 7-8. The cross
section is differential in pT in bins of Q2.
*/

#include "H1_HeraI_LQ_ptij.h"

//1JET
void H1_HeraI_LQ_ptij_1JETFilter::ReadData()
{
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/H1_HeraI_LQ_ptji/HEPData-ins838435-v1-Table_7.csv";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  string line;
  
  for(int j=0; j<14; j++)
    {
      getline(f1,line);
    }
  
  for(int k=0; k<fNData; k++)
    {
      getline(f1,line);
      istringstream lstream(line);
      char comma;
      double ddum;
      lstream >> fKin1[k] >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> fKin2[k] >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> fData[k] >> comma
	      >> fStat[k] >> comma >> comma
	      >> ddum >> comma >> comma
	      >> fSys[k][0].mult >> comma >> comma
	      >> ddum >> comma >> comma
	      >> fSys[k][1].mult >> comma >> comma
	      >> ddum >> comma >> comma	
	      >> fSys[k][2].mult >> comma >> comma
	      >> ddum >> comma >> comma
	      >> fSys[k][3].mult >> comma >> comma
	      >> ddum >> comma >> comma
	      >> fSys[k][4].mult >> comma >> comma
	      >> ddum >> comma;
      
      fStat[k] = fStat[k]/100.*fData[k];
      fKin3[k] = 319.;    //GeV

      for (int j=0; j<fNSys; j++)
	{
	  fSys[k][j].add = fSys[k][j].mult/100.*fData[k];
	  fSys[k][j].type = ADD;
	}

      fSys[k][5].mult = 1.5;
      fSys[k][5].add = fSys[k][5].mult*fData[k]/100.;
      fSys[k][5].type = MULT;
      
      fSys[k][0].name = "UNCORR";
      fSys[k][1].name = "MODEL_H1_HeraI_LQ";
      fSys[k][2].name = "ElEnSC_H1_HeraI_LQ";
      fSys[k][3].name = "ElPolAngle_H1_HeraI_LQ";
      fSys[k][4].name = "HadEnSC_H1_HeraI_LQ";
      fSys[k][5].name = "HERA1_LUMI";

    }

  f1.close();
  
}

//2JET
void H1_HeraI_LQ_ptij_2JETFilter::ReadData()
{
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/H1_HeraI_LQ_ptji/HEPData-ins838435-v1-Table_8.csv";
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
  
  for(int k=0; k<fNData; k++)
    {
      getline(f1,line);
      istringstream lstream(line);
      char comma;
      double ddum;
      lstream >> fKin1[k] >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> fKin2[k] >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> fData[k] >> comma
	      >> fStat[k] >> comma >> comma
	      >> ddum >> comma >> comma
	      >> fSys[k][0].mult >> comma >> comma
	      >> ddum >> comma >> comma
	      >> fSys[k][1].mult >> comma >> comma
	      >> ddum >> comma >> comma	
	      >> fSys[k][2].mult >> comma >> comma
	      >> ddum >> comma >> comma
	      >> fSys[k][3].mult >> comma >> comma
	      >> ddum >> comma >> comma
	      >> fSys[k][4].mult >> comma >> comma
	      >> ddum >> comma;

      fStat[k] = fStat[k]/100.*fData[k];
      fKin3[k] = 319.;    //GeV

      for (int j=0; j<fNSys; j++)
	{
	  fSys[k][j].add = fSys[k][j].mult/100.*fData[k];
	  fSys[k][j].type = ADD;
	}

      fSys[k][5].mult = 1.5;
      fSys[k][5].add = fSys[k][5].mult*fData[k]/100.;
      fSys[k][5].type = MULT;
      
      fSys[k][0].name = "UNCORR";
      fSys[k][1].name = "MODEL_H1_HeraI_LQ";
      fSys[k][2].name = "ElEnSC_H1_HeraI_LQ";
      fSys[k][3].name = "ElPolAngle_H1_HeraI_LQ";
      fSys[k][4].name = "HadEnSC_H1_HeraI_LQ";
      fSys[k][5].name = "HERA1_LUMI";

    }

  f1.close();
  
}
