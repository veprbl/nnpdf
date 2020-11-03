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

class CMS_ggF_aa_13TeVFilter: public CommonData
{
public: CMS_ggF_aa_13TeVFilter():
  CommonData("CMS_ggF_aa_13TeV") { ReadData(); }

private:
  void ReadData();
};
