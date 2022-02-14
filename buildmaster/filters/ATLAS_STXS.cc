#include "ATLAS_STXS.h"

//1) 2020 Higgs data

void ATLAS_STXS_2020Filter::ReadData()
{

  // Opening files
  fstream rData, rCorr;

  // rapidity distribution
  stringstream DataFile("");
  DataFile << dataPath() << "rawdata/" << fSetName 
	    << "/data_2020.csv";
  rData.open(DataFile.str().c_str(), ios::in);

  if (rData.fail()) {
    cerr << "Error opening data file " << DataFile.str() << endl;
    exit(-1);
  }

  // correlation matrix
  stringstream DataFileCorr("");
  DataFileCorr << dataPath() << "rawdata/" << fSetName 
	       << "/cov_2020.csv";
  rCorr.open(DataFileCorr.str().c_str(), ios::in);

  if (rCorr.fail()) {
    cerr << "Error opening data file " << DataFileCorr.str() << endl;
    exit(-1);
  }

  // Starting filter
  double s = 13000;

  string line;

  for (int i = 0; i < fNData; i++)
  {
    getline(rData,line);                  
    istringstream lstream(line); 

    fKin1[i] = 1;          // dummy variable
    fKin2[i] = 1;          // dummy variable
    fKin3[i] = s;          // sqrt(s)

    lstream >> fData[i];        
    fStat[i] = 0;          // only total covariance matrix is provided,accounting for both stat and sys
  }
  
  // Defining covariance matrix
  double** covmat = new double*[fNData];
  for (int i = 0; i < fNData; i++) 
    covmat[i] = new double[fNData];
 
  // Reading Covariance Matrix
  for (int i = 0; i < fNData; i++){
    for (int j = 0; j < fNData; j++){
	    rCorr >> covmat[i][j];
	}
  }
  
  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];
  
  if(!genArtSys(fNData,covmat,syscor))
    {
      cerr << " in " << fSetName 
	   << " : cannot generate artificial systematics" << endl;
      exit(-1);
    }
  
  // Copy the artificial systematics in the fSys matrix
  for (int i = 0; i < fNData; i++)
    for (int l = 0; l < fNSys; l++)
    {
      fSys[i][l].add  = syscor[i][l];
      fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
      fSys[i][l].type = ADD;  
      fSys[i][l].name = "CORR";
    }
  
  rData.close();
  rCorr.close();
  
  for(int i = 0; i < fNData; i++) 
    delete[] covmat[i];
  delete[] covmat;
  
  for(int i = 0; i < fNData; i++) 
    delete[] syscor[i];
  delete[] syscor;
  
}