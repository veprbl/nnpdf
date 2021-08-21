/*
Reference: 1611.03421
Hepdata: https://www.hepdata.net/record/ins1496981
Published in Eur.Phys.J.C 77 (2017) 215
A precision measurement of jet cross sections in neutral current deep-inelastic 
scattering for photon virtualities and inelasticities is presented, using data 
taken with the H1 detector at HERA, corresponding to an integrated luminosity 
of 290 pb-1. Double-differential inclusive jet and dijet cross sections are 
measured simultaneously and are presented as a function of jet transverse 
momentum observables and as a function of Q2. Jet cross sections normalised to 
the inclusive neutral current DIS cross section in the respective Q2-interval 
are also determined. 

The data is taken from Hepdata
*/

#include "H1_LQ_ptji.h"

//1JET (absolute)
void H1_LQ_ptji_1JETFilter::ReadData()
{
  const int ntab=8;
  string nametab[ntab];
  nametab[0]="InclusivejetsforQ^2=5.5-8.0GeV^2.csv";
  nametab[1]="InclusivejetsforQ^2=8.0-11.0GeV^2.csv";
  nametab[2]="InclusivejetsforQ^2=11.0-16.0GeV^2.csv";
  nametab[3]="InclusivejetsforQ^2=16.0-22.0GeV^2.csv";
  nametab[4]="InclusivejetsforQ^2=22.0-30.0GeV^2.csv";
  nametab[5]="InclusivejetsforQ^2=30.0-42.0GeV^2.csv";
  nametab[6]="InclusivejetsforQ^2=42.0-60.0GeV^2.csv";
  nametab[7]="InclusivejetsforQ^2=60.0-80.0GeV^2.csv";

  double Q2[ntab]={6.75,9.5,13.5,19.5,26.,36.,51.,70.};
  int k=0;
  const int nptbin=6;
  
  for(int i=0; i<ntab; i++)
    {

      fstream f1;
      stringstream datafile("");
      datafile << dataPath() << "rawdata/H1_LQ_ptji/"
	       << nametab[i];
      f1.open(datafile.str().c_str(), ios::in);
      
      if (f1.fail())
	{
	  cerr << "Error opening data file " << datafile.str() << endl;
	  exit(-1);
	}
      
      string line;
      
      for(int j=0; j<16; j++)
	{
	  getline(f1,line);
	}
      
      for(int j=0; j<nptbin; j++)
	{
	  double sysp[fNSys], sysm[fNSys];
	  getline(f1,line);
	  istringstream lstream(line);
	  char comma;
	  double ddum;
	  lstream >> fKin2[k] >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> fData[k] >> comma
		  >> fStat[k] >> comma >> comma
		  >> ddum >> comma >> comma;
	  
	  for(int l=0; l<fNSys-1; l++)
	    {
	      lstream >> sysp[l] >> comma >> comma
		      >> sysm[l] >> comma >> comma;
	    }
	  
	  lstream >> sysp[fNSys-1] >> comma >> comma
		  >> sysm[fNSys-1] >> comma;	    
	  
	  fKin1[k] = Q2[i];
	  fStat[k] = fStat[k]/100.*fData[k];
	  fKin3[k] = 319.;    //GeV
	  
	  //Symmetrise errors
	  for(int l=0; l<fNSys; l++)
	    {
	      double shift, delta;
	      sysp[l] = sysp[l]/100.*fData[k];
	      sysm[l] = sysm[l]/100.*fData[k];
	      symmetriseErrors(sysp[l], sysm[l], &delta, &shift);
	      
	      fData[k] += shift;
	      fSys[k][l].add = delta;
	      fSys[k][l].mult = fSys[k][l].add/fData[k]*100.;
	      fSys[k][l].type = MULT;
	    }
	  
	  fSys[k][0].name = "UNCORR";
	  fSys[k][1].name = "Model_H1_LQ";
	  fSys[k][2].name = "ModelRW_H1_LQ";
	  fSys[k][3].name = "JES_H1_LQ";
	  fSys[k][4].name = "RCES_H1_LQ";
	  fSys[k][5].name = "ElEn_H1_LQ";
	  fSys[k][6].name = "ElTh_H1_LQ";
	  fSys[k][7].name = "Lumi_H1_LQ";
	  fSys[k][8].name = "LArn_H1_LQ";
	  fSys[k][9].name = "UNCORR";
	  fSys[k][10].name = "RadErr_H1_LQ";
	  
	  k++;
	}
      
      f1.close();
      
    }
}

//2JET (absolute)
void H1_LQ_ptji_2JETFilter::ReadData()
{
  const int ntab=8;
  string nametab[ntab];
  nametab[0]="DijetsforQ^2=5.5-8.0GeV^2.csv";
  nametab[1]="DijetsforQ^2=8.0-11.0GeV^2.csv";
  nametab[2]="DijetsforQ^2=11.0-16.0GeV^2.csv";
  nametab[3]="DijetsforQ^2=16.0-22.0GeV^2.csv";
  nametab[4]="DijetsforQ^2=22.0-30.0GeV^2.csv";
  nametab[5]="DijetsforQ^2=30.0-42.0GeV^2.csv";
  nametab[6]="DijetsforQ^2=42.0-60.0GeV^2.csv";
  nametab[7]="DijetsforQ^2=60.0-80.0GeV^2.csv";

  double Q2[ntab]={6.75,9.5,13.5,19.5,26.,36.,51.,70.};
  int k=0;
  const int nptbin=6;
  
  for(int i=0; i<ntab; i++)
    {

      fstream f1;
      stringstream datafile("");
      datafile << dataPath() << "rawdata/H1_LQ_ptji/"
	       << nametab[i];
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
      
      for(int j=0; j<nptbin; j++)
	{
	  double sysp[fNSys], sysm[fNSys];
	  getline(f1,line);
	  istringstream lstream(line);
	  char comma;
	  double ddum;
	  lstream >> fKin2[k] >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> fData[k] >> comma
		  >> fStat[k] >> comma >> comma
		  >> ddum >> comma >> comma;
	  
	  for(int l=0; l<fNSys-1; l++)
	    {
	      lstream >> sysp[l] >> comma >> comma
		      >> sysm[l] >> comma >> comma;
	    }
	  
	  lstream >> sysp[fNSys-1] >> comma >> comma
		  >> sysm[fNSys-1] >> comma;	    
	  
	  fKin1[k] = Q2[i];
	  fStat[k] = fStat[k]/100.*fData[k];
	  fKin3[k] = 319.;    //GeV

	  //Symmetrise errors
	  for(int l=0; l<fNSys; l++)
	    {
	      double shift, delta;
	      sysp[l] = sysp[l]/100.*fData[k];
	      sysm[l] = sysm[l]/100.*fData[k];
	      symmetriseErrors(sysp[l], sysm[l], &delta, &shift);
	      
	      fData[k] += shift;
	      fSys[k][l].add = delta;
	      fSys[k][l].mult = fSys[k][l].add/fData[k]*100.;
	      fSys[k][l].type = MULT;
	    }
      
	  fSys[k][0].name = "UNCORR";
	  fSys[k][1].name = "Model_H1_LQ";
	  fSys[k][2].name = "ModelRW_H1_LQ";
	  fSys[k][3].name = "JES_H1_LQ";
	  fSys[k][4].name = "RCES_H1_LQ";
	  fSys[k][5].name = "ElEn_H1_LQ";
	  fSys[k][6].name = "ElTh_H1_LQ";
	  fSys[k][7].name = "Lumi_H1_LQ";
	  fSys[k][8].name = "LArn_H1_LQ";
	  fSys[k][9].name = "UNCORR";
	  fSys[k][10].name = "RadErr_H1_LQ";
	  
	  k++;
	}
      
      f1.close();

    }
}

//1JET (normalised)
void H1_LQ_ptji_1JET_NORMFilter::ReadData()
{
  const int ntab=8;
  string nametab[ntab];
  nametab[0]="NormalisedinclusivejetsforQ^2=5.5-8.0GeV^2.csv";
  nametab[1]="NormalisedinclusivejetsforQ^2=8.0-11.0GeV^2.csv";
  nametab[2]="NormalisedinclusivejetsforQ^2=11.0-16.0GeV^2.csv";
  nametab[3]="NormalisedinclusivejetsforQ^2=16.0-22.0GeV^2.csv";
  nametab[4]="NormalisedinclusivejetsforQ^2=22.0-30.0GeV^2.csv";
  nametab[5]="NormalisedinclusivejetsforQ^2=30.0-42.0GeV^2.csv";
  nametab[6]="NormalisedinclusivejetsforQ^2=42.0-60.0GeV^2.csv";
  nametab[7]="NormalisedinclusivejetsforQ^2=60.0-80.0GeV^2.csv";

  double Q2[ntab]={6.75,9.5,13.5,19.5,26.,36.,51.,70.};
  int k=0;
  const int nptbin=6;
  
  for(int i=0; i<ntab; i++)
    {

      fstream f1;
      stringstream datafile("");
      datafile << dataPath() << "rawdata/H1_LQ_ptji/"
	       << nametab[i];
      f1.open(datafile.str().c_str(), ios::in);
      
      if (f1.fail())
	{
	  cerr << "Error opening data file " << datafile.str() << endl;
	  exit(-1);
	}
      
      string line;
      
      for(int j=0; j<16; j++)
	{
	  getline(f1,line);
	}
      
      for(int j=0; j<nptbin; j++)
	{
	  double sysp[fNSys], sysm[fNSys];
	  getline(f1,line);
	  istringstream lstream(line);
	  char comma;
	  double ddum;
	  lstream >> fKin2[k] >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> fData[k] >> comma
		  >> fStat[k] >> comma >> comma
		  >> ddum >> comma >> comma;
	  
	  for(int l=0; l<fNSys-1; l++)
	    {
	      lstream >> sysp[l] >> comma >> comma
		      >> sysm[l] >> comma >> comma;
	    }
	  
	  lstream >> sysp[fNSys-1] >> comma >> comma
		  >> sysm[fNSys-1] >> comma;	    
	  
	  fKin1[k] = Q2[i];
	  fStat[k] = fStat[k]/100.*fData[k];
	  fKin3[k] = 319.;    //GeV
	  
	  //Symmetrise errors
	  for(int l=0; l<fNSys; l++)
	    {
	      double shift, delta;
	      sysp[l] = sysp[l]/100.*fData[k];
	      sysm[l] = sysm[l]/100.*fData[k];
	      symmetriseErrors(sysp[l], sysm[l], &delta, &shift);
	      
	      fData[k] += shift;
	      fSys[k][l].add = delta;
	      fSys[k][l].mult = fSys[k][l].add/fData[k]*100.;
	      fSys[k][l].type = MULT;
	    }
	  
	  fSys[k][0].name = "UNCORR";
	  fSys[k][1].name = "Model_H1_LQ";
	  fSys[k][2].name = "ModelRW_H1_LQ";
	  fSys[k][3].name = "JES_H1_LQ";
	  fSys[k][4].name = "RCES_H1_LQ";
	  fSys[k][5].name = "ElEn_H1_LQ";
	  fSys[k][6].name = "ElTh_H1_LQ";
	  fSys[k][7].name = "Lumi_H1_LQ";
	  fSys[k][8].name = "LArn_H1_LQ";
	  fSys[k][9].name = "UNCORR";
	  fSys[k][10].name = "RadErr_H1_LQ";
	  
	  k++;
	}
      
      f1.close();
      
    }
}

//2JET (normalised)
void H1_LQ_ptji_2JET_NORMFilter::ReadData()
{
  const int ntab=8;
  string nametab[ntab];
  nametab[0]="NormaliseddijetsforQ^2=5.5-8.0GeV^2.csv";
  nametab[1]="NormaliseddijetsforQ^2=8.0-11.0GeV^2.csv";
  nametab[2]="NormaliseddijetsforQ^2=11.0-16.0GeV^2.csv";
  nametab[3]="NormaliseddijetsforQ^2=16.0-22.0GeV^2.csv";
  nametab[4]="NormaliseddijetsforQ^2=22.0-30.0GeV^2.csv";
  nametab[5]="NormaliseddijetsforQ^2=30.0-42.0GeV^2.csv";
  nametab[6]="NormaliseddijetsforQ^2=42.0-60.0GeV^2.csv";
  nametab[7]="NormaliseddijetsforQ^2=60.0-80.0GeV^2.csv";

  double Q2[ntab]={6.75,9.5,13.5,19.5,26.,36.,51.,70.};
  int k=0;
  const int nptbin=6;
  
  for(int i=0; i<ntab; i++)
    {

      fstream f1;
      stringstream datafile("");
      datafile << dataPath() << "rawdata/H1_LQ_ptji/"
	       << nametab[i];
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
      
      for(int j=0; j<nptbin; j++)
	{
	  double sysp[fNSys], sysm[fNSys];
	  getline(f1,line);
	  istringstream lstream(line);
	  char comma;
	  double ddum;
	  lstream >> fKin2[k] >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> fData[k] >> comma
		  >> fStat[k] >> comma >> comma
		  >> ddum >> comma >> comma;
	  
	  for(int l=0; l<fNSys-1; l++)
	    {
	      lstream >> sysp[l] >> comma >> comma
		      >> sysm[l] >> comma >> comma;
	    }
	  
	  lstream >> sysp[fNSys-1] >> comma >> comma
		  >> sysm[fNSys-1] >> comma;	    
	  
	  fKin1[k] = Q2[i];
	  fStat[k] = fStat[k]/100.*fData[k];
	  fKin3[k] = 319.;    //GeV

	  //Symmetrise errors
	  for(int l=0; l<fNSys; l++)
	    {
	      double shift, delta;
	      sysp[l] = sysp[l]/100.*fData[k];
	      sysm[l] = sysm[l]/100.*fData[k];
	      symmetriseErrors(sysp[l], sysm[l], &delta, &shift);
	      
	      fData[k] += shift;
	      fSys[k][l].add = delta;
	      fSys[k][l].mult = fSys[k][l].add/fData[k]*100.;
	      fSys[k][l].type = MULT;
	    }
      
	  fSys[k][0].name = "UNCORR";
	  fSys[k][1].name = "Model_H1_LQ";
	  fSys[k][2].name = "ModelRW_H1_LQ";
	  fSys[k][3].name = "JES_H1_LQ";
	  fSys[k][4].name = "RCES_H1_LQ";
	  fSys[k][5].name = "ElEn_H1_LQ";
	  fSys[k][6].name = "ElTh_H1_LQ";
	  fSys[k][7].name = "Lumi_H1_LQ";
	  fSys[k][8].name = "LArn_H1_LQ";
	  fSys[k][9].name = "UNCORR";
	  fSys[k][10].name = "RadErr_H1_LQ";
	  
	  k++;
	}
      
      f1.close();

    }
}
