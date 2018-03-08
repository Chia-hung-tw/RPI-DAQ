#ifndef CHMapping_h
#define CHMapping_h

class CHMapping{
  
 public:

  double CH_x[128];
  double CH_y[128];
  int    CH_type[128];

  CHMapping(){
    read_elec_map();
    read_CH_x_y(); };
  
 private:
  int CH_iu[128];
  int CH_iv[128];
  void read_elec_map();
  void read_CH_x_y();
  
};

#endif
