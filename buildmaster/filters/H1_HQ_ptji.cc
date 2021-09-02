/*
Reference: 1406.4709
Hepdata: https://www.hepdata.net/record/ins1301218
Published in Eur.Phys.J.C 75 (2015) 65, 2015.
Double differential jet cross sections are obtained using a regularised 
unfolding procedure. They are presented as a function of Q**2 and the 
transverse momentum of the jet, PT(JET), and as a function of Q**2 and the 
proton longitudinal momentum fraction, XI, carried by the parton participating 
in the hard interaction.

The data is taken from the paper (Tables 8-9)
*/

#include "H1_HQ_ptji.h"

//1JET (absolute)
void H1_HQ_ptji_1JETFilter::ReadData()
{
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/H1_HQ_ptji/Inclusive_Jets.dat";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  string line;

  double Q2[fNData]={175.,175., 175., 175.,
		     235.,235.,235.,235.,
		     345.,345.,345.,345.,
		     550.,550.,550.,550.,
		     2850.,2850.,2850.,2850.,
		     10000.,10000.,10000.,10000};
  double pT[fNData]={9., 14.5, 24., 40.,
		     9., 14.5, 24., 40.,
		     9., 14.5, 24., 40.,
		     9., 14.5, 24., 40.,
		     9., 14.5, 24., 40.,
		     9., 14.5, 24., 40.};
  
  for(int i=0; i<fNData; i++)
    {
      const int nsysasy=5;
      double ddum;
      double sysp[nsysasy], sysm[nsysasy];
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fData[i]
	      >> fStat[i]
	      >> ddum
	      >> fSys[i][0].mult;

      fSys[i][6].mult = 0.6; //Lar noise
      fSys[i][7].mult = 2.9; //normalisation

      fKin1[i] = Q2[i];
      fKin2[i] = pT[i];
      fKin3[i] = 319.;
      
      for(int j=0; j<nsysasy; j++)
	{
	  double delta, shift;
	  lstream >> sysp[j] >> sysm[j];
	  sysp[j] = sysp[j]*fData[i]/100.;
	  sysm[j] = sysm[j]*fData[i]/100.;
	  symmetriseErrors(sysp[j], sysm[j], &delta, &shift);
	  fData[i] += shift;
	  fSys[i][j+1].add = delta;
	  fSys[i][j+1].mult = fSys[i][j+1].add/fData[i] * 100.;
	}

      fStat[i] = fStat[i] * fData[i]/100.;
      fSys[i][0].add = fSys[i][0].mult * fData[i]/100.;
      fSys[i][6].add = fSys[i][6].mult * fData[i]/100.;
      fSys[i][7].add= fSys[i][7].mult * fData[i]/100.;
      
      fSys[i][0].name = "Model_H1_HQ";
      fSys[i][1].name = "JES_H1_HQ";
      fSys[i][2].name = "RCES_H1_HQ";
      fSys[i][3].name = "El_H1_HQ";
      fSys[i][4].name = "Theta_H1_HQ";
      fSys[i][5].name = "ID_H1_HQ";
      fSys[i][6].name = "LarNoise_H1_HQ";
      fSys[i][7].name = "Lumi_H1_HQ";

      for(int j=0; j<fNSys; j++)
	{
	  fSys[i][j].type = MULT;
	}
    }

  f1.close();

}

//2JET (absolute)
void H1_HQ_ptji_2JETFilter::ReadData()
{
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/H1_HQ_ptji/Di_Jets.dat";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  string line;

  double Q2[fNData]={175.,175., 175., 175.,
		     235.,235.,235.,235.,
		     345.,345.,345.,345.,
		     550.,550.,550.,550.,
		     2850.,2850.,2850.,2850.,
		     10000.,10000.,10000.,10000};
  double pT[fNData]={9., 14.5, 24., 40.,
		     9., 14.5, 24., 40.,
		     9., 14.5, 24., 40.,
		     9., 14.5, 24., 40.,
		     9., 14.5, 24., 40.,
		     9., 14.5, 24., 40.};
  
  for(int i=0; i<fNData; i++)
    {
      const int nsysasy=5;
      double ddum;
      double sysp[nsysasy], sysm[nsysasy];
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fData[i]
	      >> fStat[i]
	      >> ddum
	      >> fSys[i][0].mult;

      fSys[i][6].mult = 0.6; //Lar noise
      fSys[i][7].mult = 2.9; //normalisation

      fKin1[i] = Q2[i];
      fKin2[i] = pT[i];
      fKin3[i] = 319.;
      
      for(int j=0; j<nsysasy; j++)
	{
	  double delta, shift;
	  lstream >> sysp[j] >> sysm[j];
	  sysp[j] = sysp[j]*fData[i]/100.;
	  sysm[j] = sysm[j]*fData[i]/100.;
	  symmetriseErrors(sysp[j], sysm[j], &delta, &shift);
	  fData[i] += shift;
	  fSys[i][j+1].add = delta;
	  fSys[i][j+1].mult = fSys[i][j+1].add/fData[i] * 100.;
	}

      fStat[i] = fStat[i] * fData[i]/100.;
      fSys[i][0].add = fSys[i][0].mult * fData[i]/100.;
      fSys[i][6].add = fSys[i][6].mult * fData[i]/100.;
      fSys[i][7].add= fSys[i][7].mult * fData[i]/100.;
      
      fSys[i][0].name = "Model_H1_HQ";
      fSys[i][1].name = "JES_H1_HQ";
      fSys[i][2].name = "RCES_H1_HQ";
      fSys[i][3].name = "El_H1_HQ";
      fSys[i][4].name = "Theta_H1_HQ";
      fSys[i][5].name = "ID_H1_HQ";
      fSys[i][6].name = "LarNoise_H1_HQ";
      fSys[i][7].name = "Lumi_H1_HQ";

      for(int j=0; j<fNSys; j++)
	{
	  fSys[i][j].type = MULT;
	}
    }

  f1.close();

}
