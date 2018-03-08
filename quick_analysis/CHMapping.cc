#include <iostream>
#include <fstream>
#include <string>
#include "CHMapping.h"

using namespace std;

void CHMapping::read_elec_map(){
  ifstream file("./quick_analysis/src_txtfile/map_CERN_Hexaboard_OneLayers_May2017.txt");
  string line;
  int chip,ch,iu,iv,type,dum,formatCH;
  while(true){
    getline(file,line);
    if( file.eof() ) break;
    file >> chip >> ch >> dum >> dum >> dum >> iu >> iv >> type;
    formatCH = (chip-1)*32+ch/2;
    if( formatCH == 94 ) type = -1; // there is no chip 2 , ch 60
    CH_type[formatCH] = type;
    CH_iu  [formatCH] = iu;
    CH_iv  [formatCH] = iv;
  }
  file.close();
}
void CHMapping::read_CH_x_y(){
  ifstream file("./quick_analysis/src_txtfile/CH_map.txt");
  string line;
  int chip,ch,type,formatCH;
  double posx,posy;
  while(true){
    getline(file,line);
    if( file.eof() ) break;
    file >> chip >> ch >> posx >> posy >> type;
    formatCH = chip*32+ch/2;
    CH_x[formatCH] = posx;
    CH_y[formatCH] = posy;}
  file.close();
}

