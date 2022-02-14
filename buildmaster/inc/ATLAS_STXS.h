#pragma once

#include "buildmaster_utils.h"

// ********* Filters **************

class ATLAS_STXS_2020Filter: public CommonData
{
public: ATLAS_STXS_2020Filter():
  CommonData("ATLAS_STXS_2020") { ReadData(); }

private:
  void ReadData();
};