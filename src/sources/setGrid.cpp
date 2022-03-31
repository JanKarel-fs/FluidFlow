#include "setGrid.hpp"

void setGrid(Grid& g, const Settings& setting) {
  switch (setting.grid_type) {
  case 1:
    g = Grid_gamm(setting.mCells, setting.nCells, setting.ghostCells, setting.nodeWeightType);
    break;
  case 2:
    g = Grid_xy(setting.name1, setting.name2, setting.ghostCells, setting.nodeWeightType);
    break;
  default:
    cout << "No such grid type is supported!" << endl;
    exit(0);
  }
}
