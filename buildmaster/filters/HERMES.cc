//Add some information here in the header file




#include "HERMES.h"

void HERMES_PIminus_proton_decFilter::ReadData()
{
  fstream f1;

  //Central values and ttal uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/HERMES/hermes.proton.zQ2-3D.zQ2-proj.vmsub.mults_piminus.list";
  f1.open(datafile.str().c_str(),ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  string line;
  
  //Read centra values and total ucnertainties
  for(int i=0; i<36; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 0.;

      int idum;
      
      lstream >> idum >> fData[i] >> fStat[i] >> fSys[i][0].add;
      
    }

  f1.close();

  //Write relevant results on file
  ofstream f2;
  f2.open ("HERMES_PIminus_proton_dec.txt");

  for(int i=0; i<fNData; i++)
    {
      f2 << "  - errors:" << endl;
      f2 << "    - {label: unc, value: " << sqrt(fStat[i]*fStat[i]+fSys[0][i].add*fSys[0][i].add)
	 << "}" << endl;
      f2 << "    value: " << fData[i] << endl;
    }
    
  f2.close();
}

void HERMES_PIplus_proton_decFilter::ReadData()
{
  fstream f1;

  //Central values and ttal uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/HERMES/hermes.proton.zQ2-3D.zQ2-proj.vmsub.mults_piplus.list";
  f1.open(datafile.str().c_str(),ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  string line;
  
  //Read centra values and total ucnertainties
  for(int i=0; i<36; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 0.;

      int idum;
      
      lstream >> idum >> fData[i] >> fStat[i] >> fSys[i][0].add;
      
    }

  f1.close();

  //Write relevant results on file
  ofstream f2;
  f2.open ("HERMES_PIplus_proton_dec.txt");

  for(int i=0; i<fNData; i++)
    {
      f2 << "  - errors:" << endl;
      f2 << "    - {label: unc, value: " << sqrt(fStat[i]*fStat[i]+fSys[0][i].add*fSys[0][i].add)
	 << "}" << endl;
      f2 << "    value: " << fData[i] << endl;
    }
    
  f2.close();
}

void HERMES_PIminus_deuteron_decFilter::ReadData()
{
  fstream f1;

  //Central values and ttal uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/HERMES/hermes.deuteron.zQ2-3D.zQ2-proj.vmsub.mults_piminus.list";
  f1.open(datafile.str().c_str(),ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  string line;
  
  //Read centra values and total ucnertainties
  for(int i=0; i<36; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 0.;

      int idum;
      
      lstream >> idum >> fData[i] >> fStat[i] >> fSys[i][0].add;
      
    }

  f1.close();

  //Write relevant results on file
  ofstream f2;
  f2.open ("HERMES_PIminus_deuteron_dec.txt");

  for(int i=0; i<fNData; i++)
    {
      f2 << "  - errors:" << endl;
      f2 << "    - {label: unc, value: " << sqrt(fStat[i]*fStat[i]+fSys[0][i].add*fSys[0][i].add)
	 << "}" << endl;
      f2 << "    value: " << fData[i] << endl;
    }
    
  f2.close();
}

void HERMES_PIplus_deuteron_decFilter::ReadData()
{
  fstream f1;

  //Central values and ttal uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/HERMES/hermes.deuteron.zQ2-3D.zQ2-proj.vmsub.mults_piplus.list";
  f1.open(datafile.str().c_str(),ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  string line;
  
  //Read centra values and total ucnertainties
  for(int i=0; i<36; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 0.;

      int idum;
      
      lstream >> idum >> fData[i] >> fStat[i] >> fSys[i][0].add;
      
    }

  f1.close();

  //Write relevant results on file
  ofstream f2;
  f2.open ("HERMES_PIplus_deuteron_dec.txt");

  for(int i=0; i<fNData; i++)
    {
      f2 << "  - errors:" << endl;
      f2 << "    - {label: unc, value: " << sqrt(fStat[i]*fStat[i]+fSys[0][i].add*fSys[0][i].add)
	 << "}" << endl;
      f2 << "    value: " << fData[i] << endl;
    }
    
  f2.close();
}
