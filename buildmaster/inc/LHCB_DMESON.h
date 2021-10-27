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

class LHCB_DMESON_R_nuclear_forwardFilter: public CommonData
{ public: LHCB_DMESON_R_nuclear_forwardFilter():
  CommonData("LHCB_DMESON_R_nuclear_forward") { ReadData(); }

private:
  void ReadData();
}; 

class LHCB_DMESON_R_nuclear_backwardFilter: public CommonData
{ public: LHCB_DMESON_R_nuclear_backwardFilter():
  CommonData("LHCB_DMESON_R_nuclear_backward") { ReadData(); }

private:
  void ReadData();
}; 