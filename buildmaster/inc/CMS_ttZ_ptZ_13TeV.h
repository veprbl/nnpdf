// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class CMS_ttZ_ptZ_13TeV
 *  \brief  CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMS_ttZ_ptZ_13TeVFilter: public CommonData
{
public: CMS_ttZ_ptZ_13TeVFilter():
  CommonData("CMS_ttZ_ptZ_13TeV") { ReadData(); }

private:
  void ReadData();
};
