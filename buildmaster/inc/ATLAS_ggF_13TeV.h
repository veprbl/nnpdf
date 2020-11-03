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

class ATLAS_ggF_13TeVFilter: public CommonData
{
public: ATLAS_ggF_13TeVFilter():
  CommonData("ATLAS_ggF_13TeV") { ReadData(); }

private:
  void ReadData();
};
