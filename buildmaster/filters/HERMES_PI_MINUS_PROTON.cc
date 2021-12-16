#include "HERMES.h"

//1) PROTON, PI-

void HERMES_PI_MINUS_PROTONFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/HERMES/hermes.proton.zQ2-3D.zQ2-proj.vmsub.mults_piminus.list";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/HERMES/hermes.proton.zQ2-3D.zQ2-proj.vmsub.covmat_mults.list";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << covfile.str() << endl;
      exit(-1);
    }

  //Read central values
  string line;
  for (int i=0; i<36; i++)
    {
      getline(f1,line);
      getline(f2,line);
    }

  for(int i=0; i<fNData; i++)
    {
      int idum;
      double ddum;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum >> fData[i] >> ddum >> fStat[i];
      fKin1[i]=0.0;
      fKin2[i]=0.0;
      fKin3[i]=0.0;	
    }

  //Read covariance matrix
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  double row, col, ddum;
	  getline(f2,line);
	  istringstream lstream(line);
	  lstream >> row >> col >> ddum >> covmat[i][j];
	  fSys[i][j].name = "CORR";
	  fSys[i][j].type = ADD;
	}
    }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  //Print results on screen
  cout << "  values:" << endl;
  
  for(int i=0; i<fNData; i++)
    {
      cout << "  - errors:" << endl;
      cout << "    - {label: unc, value: " << fStat[i] << "}" << endl;
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult  = fSys[i][j].add/fData[i];
	  cout << "    - {label: unc, value: " << fSys[i][j].mult << "}" << endl;
	}
      cout << "    value: " << fData[i] << endl;
    }
  
  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];
  
  delete[] syscor;
  
  for(int i=0; i<fNData; i++)
    delete[] covmat[i];
  
  delete[] covmat;
  
  f1.close();
  f2.close();
  
}
