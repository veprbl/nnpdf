#pragma once

#include "buildmaster_utils.h"

// ********* Filters **************

class LHCB_DMESON_R_13_5Filter: public CommonData
{ public: LHCB_DMESON_R_13_5Filter():
  CommonData("LHCB_DMESON_R_13_5") { ReadData(); }

private:
  void ReadData();
}; 

class LHCB_DMESON_N7Filter: public CommonData
{ public: LHCB_DMESON_N7Filter():
  CommonData("LHCB_DMESON_N7") { ReadData(); }

private:
  void ReadData();
}; 

class LHCB_DMESON_N5Filter: public CommonData
{ public: LHCB_DMESON_N5Filter():
  CommonData("LHCB_DMESON_N5") { ReadData(); }

private:
  void ReadData();
}; 

class LHCB_DMESON_N13Filter: public CommonData
{ public: LHCB_DMESON_N13Filter():
  CommonData("LHCB_DMESON_N13") { ReadData(); }

private:
  void ReadData();
};