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

class ZEUS_HQ_ptavg_12Filter: public CommonData
{
public: ZEUS_HQ_ptavg_12Filter():
  CommonData("ZEUS_HQ_ptavg_12") { ReadData(); }

private:
  void ReadData();
};
