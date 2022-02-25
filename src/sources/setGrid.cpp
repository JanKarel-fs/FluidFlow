#include "setGrid.hpp"

void setGrid(Grid& g, const Settings& setting) {
  switch (setting.grid_type) {
  case 1:
    g = Grid_gamm(setting.mCells, setting.nCells, setting.ghostCells);
    break;
  default:
    cout << "No such grid type is supported!" << endl;
    exit(0);
  }
}
