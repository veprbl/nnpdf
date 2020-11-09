/*
Experiment: CMS
Eprint: arXiv:1908.07305
Published in: Eur. Phys. J. C79 (2019) 1028
Hepdata: https://www.hepdata.net/record/95758?version=1
Description: Absolute differential cross-section as a function of 
$m^{t\bar{t}}$ at parton level in the resolved topology. The values 
used here are obtained from Table 837 of the Hepdata entry. The kinematic 
binning is equivalent to the that of the CMS dilepton measurement 
[arXiv:1811.06625]. There are 208 ssytematic ucnertainties.
*/

#include "ATLAS_TTB_DIFF_13TEV_LJ.h"

void ATLAS_TTB_DIFF_13TEV_LJ_TTMFilter::ReadData()
{
  //Opening files
  fstream f1;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
           << "rawdata/ATLAS_TTB_DIFF_13TEV_LJ/HEPData-ins1750330-v1-Table_837.csv";
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
            >> fStat[i] >> comma >> comma >> comma
	    >> dud >> comma;

    fKin1[i] = pt_top; //mtt [GeV]
    fKin2[i] = Mt * Mt;
    fKin3[i] = 13000; //sqrt(s) [GeV]
    fStat[i] = fStat[i]/100.*fData[i];

    for(int j=0; j<fNSys; j++)
      {
	lstream >> comma >> fSys[i][j].mult >> comma;
	
	fSys[i][j].mult = fSys[i][j].mult/sqrt(2.);
	fSys[i][j].add  = fSys[i][j].mult/100.*fData[i];
	fSys[i][j].type = MULT;
	fSys[i][j].name = "CORR";
      }
  }
  
  f1.close();
}


