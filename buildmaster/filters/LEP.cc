/*


Add here a description of the data set




*/

#include "LEP.h"

void LEP_eeWW_182GeV::ReadData()
{
  fstream f1;

  //Central value and uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/LEP_WW/eeww_182GeV.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values and uncertainties
  string line;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin1[i] = 182.0;     //GeV
      fKin1[i] = 0.;
      lstream >> fData[i] >> fStat[i] >> fSys[i][0].add;
      fSys[i][0].mult = fSys[i][0].add*1e2/fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
    }

  f1.close();
  
}

void LEP_eeWW_189GeV::ReadData()
{
  fstream f1;

  //Central value and uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/LEP_WW/eeww_189GeV.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values and uncertainties
  string line;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin1[i] = 189.0;     //GeV
      fKin1[i] = 0.;
      lstream >> fData[i] >> fStat[i] >> fSys[i][0].add;
      fSys[i][0].mult = fSys[i][0].add*1e2/fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
    }

  f1.close();
  
}

void LEP_eeWW_198GeV::ReadData()
{
  fstream f1;

  //Central value and uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/LEP_WW/eeww_198GeV.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values and uncertainties
  string line;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin1[i] = 198.0;     //GeV
      fKin1[i] = 0.;
      lstream >> fData[i] >> fStat[i] >> fSys[i][0].add;
      fSys[i][0].mult = fSys[i][0].add*1e2/fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
    }

  f1.close();
  
}

void LEP_eeWW_206GeV::ReadData()
{
  fstream f1;

  //Central value and uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/LEP_WW/eeww_206GeV.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values and uncertainties
  string line;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin1[i] = 206.0;     //GeV
      fKin1[i] = 0.;
      lstream >> fData[i] >> fStat[i] >> fSys[i][0].add;
      fSys[i][0].mult = fSys[i][0].add*1e2/fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
    }

  f1.close();
  
}

void LEP_eeWW_all::ReadData()
{
  fstream f1;

  //Central value and uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/LEP_WW/eeww_all.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values and uncertainties
  string line;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin1[i] = 0.;     //GeV
      fKin1[i] = 0.;
      lstream >> fData[i] >> fStat[i]
	      >> fSys[i][0].add
	      >> fSys[i][1].add
	      >> fSys[i][2].add;
      
      fSys[i][0].mult = fSys[i][0].add*1e2/fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";

      fSys[i][1].mult = fSys[i][1].add*1e2/fData[i];
      fSys[i][1].type = ADD;
      fSys[i][1].name = "UNCORR";

      fSys[i][2].mult = fSys[i][2].add*1e2/fData[i];
      fSys[i][2].type = ADD;
      fSys[i][2].name = "UNCORR";
      
    }

  f1.close();
  
}
