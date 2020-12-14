// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class NTVNUDMNFilter: public CommonData
{
public: NTVNUDMNFilter():
  CommonData("NTVNUDMN") { ReadData(); }

private:
  void ReadData();
};

class NTVNBDMNFilter: public CommonData
{
public: NTVNBDMNFilter():
  CommonData("NTVNBDMN") { ReadData(); }

private:
  void ReadData();
};

class NTVNUDMN_PDF4LHCFilter: public CommonData
{
public: NTVNUDMN_PDF4LHCFilter():
  CommonData("NTVNUDMN_PDF4LHC") { ReadData(); }

private:
  void ReadData();
};

class NTVNBDMN_PDF4LHCFilter: public CommonData
{
public: NTVNBDMN_PDF4LHCFilter():
  CommonData("NTVNBDMN_PDF4LHC") { ReadData(); }

private:
  void ReadData();
};
