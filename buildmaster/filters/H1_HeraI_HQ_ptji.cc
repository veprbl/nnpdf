/*
Reference: 0706.3722
Hepdata: https://www.hepdata.net/record/ins753951
Published in Phys.Lett.B 653 (2007) 134
DESY-HERA. Measurement of inclusive jet production in neutral current DIS of 
27.5 GeV positrons with 920 GeV protons equivalent to a centre of mass energy 
of 319 GeV. The data with an integrated luminosity of 65.4 pb-1 were analysed 
using the inclusive Kt jet finding algorithm in the Breit frame of reference, 
with single and double differential distributions being presented as functions 
of Q**2 (& gt; 150 Gev**2) and ET.

The data is taken from the paper, specifically from Tab.1. The cross
section is differential in pT in bins of Q2.
*/

#include "H1_HeraI_HQ_ptij.h"

//1JET
void H1_HeraI_HQ_ptij_1JETFilter::ReadData()
{

  double Q2[fNData]={175.,175.,175.,175.,
		     235.,235.,235.,235.,
		     335.,335.,335.,335.,
		     550.,550.,550.,550.,
		     2850.,2850.,2850.,2850.,
		     10000.,10000.,10000.,10000};
  
  double pT[fNData]={9.,15.5,24.,40.,
		     9.,15.5,24.,40.,
		     9.,15.5,24.,40.,
		     9.,15.5,24.,40.,
		     9.,15.5,24.,40.,
		     9.,15.5,24.,40.};
    
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/H1_HeraI_HQ_ptji/data.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  string line;
  
  for(int k=0; k<fNData; k++)
    {
      getline(f1,line);
      istringstream lstream(line);
      double ddum;
      lstream >> fData[k] 
	      >> fStat[k]
	      >> ddum
	      >> fSys[k][0].mult
	      >> ddum
	      >> fSys[k][1].mult
	      >> fSys[k][2].mult
	      >> fSys[k][3].mult
	      >> fSys[k][4].mult;
      
      fStat[k] = fStat[k]/100.*fData[k];
      fKin3[k] = 319.;    //GeV
      fKin1[k] = Q2[k];
      fKin2[k] = pT[k];
      
      for (int j=0; j<fNSys; j++)
	{
	  fSys[k][j].add = fSys[k][j].mult/100.*fData[k];
	  fSys[k][j].type = ADD;
	}

      fSys[k][5].mult = 1.5;
      fSys[k][5].add = fSys[k][5].mult*fData[k]/100.;
      fSys[k][5].type = MULT;
      
      fSys[k][0].name = "UNCORR";
      fSys[k][1].name = "MODEL_H1_HeraI_HQ";
      fSys[k][2].name = "ElEnSC_H1_HeraI_HQ";
      fSys[k][3].name = "ElPolAngle_H1_HeraI_HQ";
      fSys[k][4].name = "HadEnSC_H1_HeraI_HQ";
      fSys[k][5].name = "HERA1_LUMI";

    }

  f1.close();
  
}

