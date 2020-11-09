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

class ATLAS_TTB_DIFF_13TEV_LJ_TTMFilter: public CommonData
{
 public: ATLAS_TTB_DIFF_13TEV_LJ_TTMFilter():
  CommonData("ATLAS_TTB_DIFF_13TEV_LJ_TTM") { ReadData(); }

 private:
  void ReadData();
};
