// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLAS_ggF_13TeV
 *  \brief  CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_ggF_13TeV_2015_pTHFilter: public CommonData
{
public: ATLAS_ggF_13TeV_2015_pTHFilter():
  CommonData("ATLAS_ggF_13TeV_2015_pTH") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_ggF_13TeV_2015_yHFilter: public CommonData
{
public: ATLAS_ggF_13TeV_2015_yHFilter():
  CommonData("ATLAS_ggF_13TeV_2015_yH") { ReadData(); }

private:
  void ReadData();
};

