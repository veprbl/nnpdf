// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class HERMES
 *  \brief HERMES CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class HERMES_PIplus_protonFilter: public CommonData
{
public: HERMES_PIplus_protonFilter():
  CommonData("HERMES_PIplus_proton") { ReadData(); }

private:
  void ReadData();
};

class HERMES_PIminus_protonFilter: public CommonData
{
public: HERMES_PIminus_protonFilter():
  CommonData("HERMES_PIminus_proton") { ReadData(); }

private:
  void ReadData();
};

class HERMES_PIplus_deuteronFilter: public CommonData
{
public: HERMES_PIplus_deuteronFilter():
  CommonData("HERMES_PIplus_deuteron") { ReadData(); }

private:
  void ReadData();
};

class HERMES_PIminus_deuteronFilter: public CommonData
{
public: HERMES_PIminus_deuteronFilter():
  CommonData("HERMES_PIminus_deuteron") { ReadData(); }

private:
  void ReadData();
};
