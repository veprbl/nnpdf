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

class ZEUS_DISJETSFilter: public CommonData
{
public: ZEUS_DISJETSFilter():
  CommonData("ZEUS_DISJETS") { ReadData(); }

private:
  void ReadData();
};
