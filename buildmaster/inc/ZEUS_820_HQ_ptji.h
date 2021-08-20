// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
//          Luca Rottoli,     luca.rottoli@physics.oc.ac.uk

#pragma once

#include "buildmaster_utils.h"

class ZEUS_820_HQ_ptjiFilter: public CommonData
{
public: ZEUS_820_HQ_ptjiFilter():
  CommonData("ZEUS_820_HQ_ptji") { ReadData(); }

private:
  void ReadData();
};
