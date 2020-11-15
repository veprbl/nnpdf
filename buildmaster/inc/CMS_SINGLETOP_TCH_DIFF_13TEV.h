// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

//Normalised distributions

class CMS_SINGLETOP_TCH_DIFF_13TEV_TTBAR_RAPFilter: public CommonData
{
 public: CMS_SINGLETOP_TCH_DIFF_13TEV_TTBAR_RAPFilter():
  CommonData("CMS_SINGLETOP_TCH_DIFF_13TEV_TTBAR_RAP") { ReadData(); }

 private:
  void ReadData();
};
