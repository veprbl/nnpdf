// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLAS_hxsec_RunII
 *  \brief  CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class LEP_eeWW_182GeVFilter: public CommonData
{
public: LEP_eeWW_182GeVFilter():
  CommonData("LEP_eeWW_182GeV") { ReadData(); }

private:
  void ReadData();
};

class LEP_eeWW_189GeVFilter: public CommonData
{
public: LEP_eeWW_189GeVFilter():
  CommonData("LEP_eeWW_189GeV") { ReadData(); }

private:
  void ReadData();
};

class LEP_eeWW_198GeVFilter: public CommonData
{
public: LEP_eeWW_198GeVFilter():
  CommonData("LEP_eeWW_198GeV") { ReadData(); }

private:
  void ReadData();
};

class LEP_eeWW_206GeVFilter: public CommonData
{
public: LEP_eeWW_206GeVFilter():
  CommonData("LEP_eeWW_206GeV") { ReadData(); }

private:
  void ReadData();
};

class LEP_eeWW_allFilter: public CommonData
{
public: LEP_eeWW_allFilter():
  CommonData("LEP_eeWW_all") { ReadData(); }

private:
  void ReadData();
};
