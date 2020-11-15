/*
Experiment: CMS
Eprint: 1907.08330
Published in Eur.Phys.J.C 80 (2020) 370
Hepdata: https://www.hepdata.net/record/ins1744604
Description: Measurement of the differential cross sections for t-channel 
single top quark and antiquark production in proton-proton collisions at a 
centre-of-mass energy of 13 TeV by the CMS experiment at the LHC. From a data 
set corresponding to an integrated luminosity of 35.9 fb−1, events containing 
one muon or electron and two or three jets are analysed. The cross section is 
measured as a function of the top quark rapidity. See table 3 in Hepdata.
*/

#include "CMS_SINGLETOP_TCH_DIFF_13TEV.h"

void CMS_SINGLETOP_TCH_DIFF_13TEV_TTBAR_RAPFilter::ReadData()
{
  //Opening files
  fstream f1;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
           << "rawdata/CMS_SINGLETOP_TCH_DIFF_13TEV/HEPData-ins1744604-v1-Table_3.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Read central values
  string line;
  for (int i = 0; i < 11; i++)
  {
    getline(f1, line);
  }

  for (int i = 0; i < fNData; i++)
  {
    double pt_top, dud;
    char comma;

    getline(f1, line);
    istringstream lstream(line);
    lstream >> pt_top >> comma
            >> dud >> comma
            >> dud >> comma
            >> fData[i] >> comma
            >> fStat[i] >> comma
	    >> dud;

    fKin1[i] = pt_top; //mtt [GeV]
    fKin2[i] = Mt * Mt;
    fKin3[i] = 13000; //sqrt(s) [GeV]

    for(int j=0; j<fNSys; j++)
      {
	lstream >> comma >> fSys[i][j].add
		>> comma >> dud;
	
	fSys[i][j].mult = fSys[i][j].add/fData[i]*100;
	fSys[i][j].type = MULT;
	fSys[i][j].name = "CORR";
	if(j==fNSys-1)
	  fSys[i][j].name="CMSLUMI15";
      }
  }
  
  f1.close();
}


