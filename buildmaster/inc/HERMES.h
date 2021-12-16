// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class HERMES_PI_PLUS_PROTONFilter: public CommonData
{
public: HERMES_PI_PLUS_PROTONFilter():
  CommonData("HERMES_PI_PLUS_PROTON") { ReadData(); }

private:
  void ReadData();
};

class HERMES_PI_MINUS_PROTONFilter: public CommonData
{
public: HERMES_PI_MINUS_PROTONFilter():
  CommonData("HERMES_PI_MINUS_PROTON") { ReadData(); }

private:
  void ReadData();
};

class HERMES_PI_PLUS_DEUTERONFilter: public CommonData
{
public: HERMES_PI_PLUS_DEUTERONFilter():
  CommonData("HERMES_PI_PLUS_DEUTERON") { ReadData(); }

private:
  void ReadData();
};

class HERMES_PI_MINUS_DEUTERONFilter: public CommonData
{
public: HERMES_PI_MINUS_DEUTERONFilter():
  CommonData("HERMES_PI_MINUS_DEUTERON") { ReadData(); }

private:
  void ReadData();
};

class HERMES_KA_PLUS_PROTONFilter: public CommonData
{
public: HERMES_KA_PLUS_PROTONFilter():
  CommonData("HERMES_KA_PLUS_PROTON") { ReadData(); }

private:
  void ReadData();
};

class HERMES_KA_MINUS_PROTONFilter: public CommonData
{
public: HERMES_KA_MINUS_PROTONFilter():
  CommonData("HERMES_KA_MINUS_PROTON") { ReadData(); }

private:
  void ReadData();
};

class HERMES_KA_PLUS_DEUTERONFilter: public CommonData
{
public: HERMES_KA_PLUS_DEUTERONFilter():
  CommonData("HERMES_KA_PLUS_DEUTERON") { ReadData(); }

private:
  void ReadData();
};

class HERMES_KA_MINUS_DEUTERONFilter: public CommonData
{
public: HERMES_KA_MINUS_DEUTERONFilter():
  CommonData("HERMES_KA_MINUS_DEUTERON") { ReadData(); }

private:
  void ReadData();
};

