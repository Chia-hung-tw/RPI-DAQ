#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "CHMapping.h"

//================================= Options ===================================
// Masking time sample 9 ~ 12, since they are effect by trigger signal
const bool MASKTS = true;     // Mask trigger signal
const int  EXPECTED_TS  = 11;

const bool DrawConnectCH = true; // false to draw Uncnnected CH
const int  HGorLGtodraw = 2; //0 for HG, 1 for LG, 2 for both
const int  SCAtodraw = 0;
//=============================================================================

#define RAWSIZE 30787
#define totSCA 15  // 13SCAs + TOT + TOA
#define nSCA 13    // Real SCAs
using namespace std;

unsigned char raw[RAWSIZE];
unsigned int ev[4][1924];
int dati[4][128][totSCA];
int dati_sum[4][128][nSCA];
double dati_sumsq[4][128][nSCA];
int mem_counter[4][128][nSCA];
int evt_counter;

float avg_HG  [4][64][nSCA];
float sigma_HG[4][64][nSCA];
float avg_LG  [4][64][nSCA];
float sigma_LG[4][64][nSCA];

int decode_raw();
int format_channels();
int roll_position(); //Output is the location of TS 0
int init(); // initialize the sum and sumsq array
int channel_sum(int rollpos);
void ped_and_noise();
void dump_PandN(string fileName);
void dump_frame();   // hexagonal shape and chip,ch info
void dump_HGxynoise_avg(); // for plotting
void dump_LGxynoise_avg(); // for plotting
void dump_ped_check();     // dump HG pedestal for plotting
void dump_gnuplotconfig(string fileN); // tell gnuplot which file is the input
void Inj_ana();

int main(){

  evt_counter = 0;
  
  string logcontent;
  ifstream logfile("./runinfo.log");
  logfile >> logcontent;
  char fileName[100];
  sprintf(fileName,"%s",logcontent.c_str());
  ifstream file(fileName);
  cout << "input file: " << fileName << endl;
  if (file.is_open()){
    init();
    while(true){
      if(evt_counter % 100 == 0)
	cout << "processing event " << evt_counter << "..." << endl;
      unsigned char testeof;
      file >> testeof;
      if( file.eof() ) break;
      else{
	testeof = raw[0];
	for(int i = 1 ; i < RAWSIZE; ++i)
	  file >> raw[i];
	decode_raw();
	format_channels();
	int rollpos = roll_position();
	//Sum and Sum square calculator (use to calculate RMS and sigma)
	channel_sum(rollpos);
	evt_counter++;
      }
    }
  }
  file.close();
  cout << "Evt = " << evt_counter << endl;
  ped_and_noise();
  dump_PandN(fileName);

  //plotting stuff
  //dump_frame();
  if(HGorLGtodraw == 0)
    dump_HGxynoise_avg();
  else if(HGorLGtodraw == 1)
    dump_LGxynoise_avg();
  else{
    dump_HGxynoise_avg();
    dump_LGxynoise_avg();}
  dump_gnuplotconfig(fileName);
  dump_ped_check();
  Inj_ana();
}

int decode_raw(){
    int i, j, k, m;
    unsigned char x;
	unsigned int t;
	unsigned int bith, bit11, bit10, bit9, bit8, bit7, bit6, bit5, bit4, bit3, bit2, bit1, bit0;
	for( i = 0; i < 1924; i = i+1){
		for (k = 0; k < 4; k = k + 1){
			ev[k][i] = 0;
		}
   	}
    for( i = 0; i < 1924; i = i+1){
        for (j=0; j < 16; j = j+1){
            x = raw[1 + i*16 + j];
            x = x & 15;
            for (k = 0; k < 4; k = k + 1){
            	ev[k][i] = ev[k][i] | (unsigned int) (((x >> (3 - k) ) & 1) << (15 - j));
            }
        }
    }
    
    /*****************************************************/
    /*    Gray to binary conversion                      */
    /*****************************************************/
   	for(k = 0; k < 4 ; k = k +1 ){
   		for(i = 0; i < 128*totSCA; i = i + 1){
   			bith = ev[k][i] & 0x8000;
   			t = ev[k][i] & 0x7fff;
   			bit11 = (t >> 11) & 1;
        	bit10 = bit11 ^ ((t >>10) &1);
        	bit9 = bit10 ^ ((t >>9) &1);
        	bit8 = bit9 ^ ((t >>8) &1);
        	bit7 = bit8 ^ ((t >>7) &1);
        	bit6 = bit7 ^ ((t >>6) &1);
        	bit5 = bit6 ^ ((t >>5) &1);
        	bit4 = bit5 ^ ((t >>4) &1);
        	bit3 = bit4 ^ ((t >>3) &1);
        	bit2 = bit3 ^ ((t >>2) &1);
        	bit1 = bit2 ^ ((t >>1) &1);
        	bit0 = bit1 ^ ((t >>0) &1);
        	ev[k][i] =  bith | ((bit11 << 11) + (bit10 << 10) + (bit9 << 9) + (bit8 << 8) + (bit7 << 7) + (bit6 << 6) + (bit5 << 5) + (bit4 << 4) + (bit3  << 3) + (bit2 << 2) + (bit1  << 1) + bit0);
        }
    }
}

int format_channels(){
    int chip, hit, ch;
	
    for (chip =0; chip < 4; chip = chip +1 ){
        for (ch = 0; ch < 128; ch = ch +1 ){
            for (hit = 0 ; hit <totSCA ; hit = hit +1){
                dati[chip][ch][hit] = ev[chip][hit*128+ch] & 0x0FFF;
            }
        }
    }
    return(0);
}

int roll_position(){
  unsigned int roll_check;  //Just to check if 4 chip has same rollmask
  int chip,first,sec,rollpos;
  bool skip_evt;
  for (chip =0; chip < 4; chip = chip +1 ){
    unsigned int roll;
    roll = ev[chip][1920] & 0x1FFF;
    //printf("roll pos = %x\n",roll);
    if(chip == 0) roll_check = roll;
    else if( roll_check != roll ){
      cout << "Problematic event!( No. " << evt_counter <<") Chip " << chip
	   << "has different rollMask! Skip event!" << endl;
      skip_evt = true;
    }

    unsigned char bits[13];
    first = -1; // first is actually second XD
    sec   = -1;
    for(int bit = 0; bit < 13 ; ++bit){
      bits[bit] = (roll >> bit) & 1;}
    for(int bit = 0 ; bit < 13; ++bit) {if((int)bits[bit] == 1) first = bit;}
    for(int bit = 0 ; bit < 13; ++bit) {
      if((int)bits[bit] == 1 && first != bit) sec = bit;}
    if(first == 12 && sec == 11) rollpos = 0;
    else if(first == 12 && sec == 0)  rollpos = 1;
    else rollpos = first+1;
    
    //for(int bit = 0; bit < 13 ; ++bit)      cout << (int)bits[bit] << " ";
    //cout << first << " , " << sec << ", rollpos = " << rollpos << endl;
    //getchar();
   }
  if(skip_evt){ return (-1); }
  return rollpos;
}
int init(){
  int chip, hit, ch;
  for (chip =0; chip < 4; chip = chip +1 ){
    for (ch = 0; ch < 128; ch = ch +1 ){
      for (hit = 0 ; hit < nSCA ; hit = hit +1){
	dati_sum[chip][ch][hit]    = 0;
	dati_sumsq[chip][ch][hit]  = 0;
	mem_counter[chip][ch][hit] = 0;
      }
    }
  }
  return(0);

}
int channel_sum(int rollpos){
  int chip, hit, ch;
  for (chip =0; chip < 4; chip = chip +1 ){
    for (ch = 0; ch < 128; ch = ch +1 ){
      int check = 0;
      for (hit = 0 ; hit < nSCA ; hit = hit +1){
	if(MASKTS){
	  int ts;
	  if(hit >= rollpos) ts = hit - rollpos;
	  else ts = 12+hit-rollpos+1;
	  if (ts >= EXPECTED_TS)  continue;}
	
	check++;
	mem_counter[chip][ch][hit]++;
	dati_sum[chip][ch][hit] += dati[chip][ch][hit];
	dati_sumsq[chip][ch][hit] += dati[chip][ch][hit]*dati[chip][ch][hit];
      }
      if(MASKTS && check != EXPECTED_TS)
	cout << "rollposition is weird! "<< check << " TS was saved!" << endl;
    }
  }
  return(0);
}

void ped_and_noise(){
  int chip, hit, ch;
  for (chip =0; chip < 4; chip = chip +1 ){
    for (ch = 0; ch < 128; ch = ch +1 ){
      for (hit = 0 ; hit < nSCA ; hit = hit +1){
	//LG
	if( ch < 64 ) {
	  int correct_CH = 63-ch;
	  //	  cout << "ch: " << correct_CH << ",SCA " << hit << ", mem = " << mem_counter[chip][ch][hit] << endl;
	  avg_LG  [chip][correct_CH][hit] = (float)dati_sum[chip][ch][hit] / mem_counter[chip][ch][hit];
	  sigma_LG[chip][correct_CH][hit] = (float)dati_sumsq[chip][ch][hit] / mem_counter[chip][ch][hit] - avg_LG  [chip][correct_CH][hit]*avg_LG  [chip][correct_CH][hit];
	  sigma_LG[chip][correct_CH][hit] = sqrt(sigma_LG[chip][correct_CH][hit]);
	  //cout << (float)dati_sumsq[chip][ch][hit] / evt_counter << " , " << avg_LG  [chip][correct_CH][hit] << " , N = " << evt_counter << endl;
	  //getchar();
	}
	//HG
	else{
	  int correct_CH = 127-ch;
	  avg_HG  [chip][correct_CH][hit] = (float)dati_sum[chip][ch][hit] / mem_counter[chip][ch][hit];
	  sigma_HG[chip][correct_CH][hit] = (float)dati_sumsq[chip][ch][hit] / mem_counter[chip][ch][hit] - avg_HG  [chip][correct_CH][hit]*avg_HG  [chip][correct_CH][hit];
	  sigma_HG[chip][correct_CH][hit] = sqrt(sigma_HG[chip][correct_CH][hit]);
	}
      }
    }
  }  
}

void dump_PandN(string fileName){
  int start = fileName.find_last_of("/");
  int end   = fileName.find(".raw");
  string Name = fileName.substr(start+1,end-start-1);

  char outtitle[100];
  sprintf(outtitle,"./ANA_out/%s_HG.txt",Name.c_str());
  ofstream fileHG(outtitle);
  cout << "Analysis output for HG at: " << outtitle << endl;
  sprintf(outtitle,"./ANA_out/%s_LG.txt",Name.c_str());
  ofstream fileLG(outtitle);
  cout << "Analysis output for LG at: " << outtitle << endl;
  fileHG << "CHIP\tCH\t";
  fileLG << "CHIP\tCH\t";
  for(int i = 0; i < nSCA; ++i){
    fileHG << "SCA " << i << " ";
    fileLG << "SCA " << i << " ";}
  fileHG << "\n";
  fileLG << "\n";
  int chip, hit, ch;
  for (chip =0; chip < 4; chip = chip +1 ){
    for (ch = 0; ch < 64; ch = ch +1 ){
      fileHG << chip << "\t" << ch << "\t";
      fileLG << chip << "\t" << ch << "\t";
      for (hit = 0 ; hit < nSCA ; hit = hit +1){
	fileHG << fixed << setprecision(2) << avg_HG[chip][ch][hit] << " ";
	fileLG << fixed << setprecision(2) << avg_LG[chip][ch][hit] << " ";}
      fileHG << "\n";
      fileLG << "\n";
      fileHG << chip << "\t" << ch << "\t";
      fileLG << chip << "\t" << ch << "\t";
      for (hit = 0 ; hit < nSCA ; hit = hit +1){
	fileHG << fixed << setprecision(2) << sigma_HG[chip][ch][hit] << " ";
	fileLG << fixed << setprecision(2) << sigma_LG[chip][ch][hit] << " ";}
      fileHG << "\n";
      fileLG << "\n";      
    }
  }
  fileHG.close();
  fileLG.close();
}

void dump_frame(){
  ifstream file("./quick_analysis/src_txtfile/poly_frame.txt");
  string line;
  int iu,iv,mem;
  int MAX_MEM = 6;
  double posx[MAX_MEM],posy[MAX_MEM];
  for(int header = 0; header < 4; ++header )     getline(file,line);
  ofstream outfile("./gnuplot/gnu_frame.dat");
  outfile << "# X\tY" << endl;
  
  while(true){
    getline(file,line);
    if( file.eof() ) break;
    file >> iu >> iv >> mem;
    
    for(int i = 0; i < mem ; ++i){
      getline(file,line);
      file >> posx[i] >> posy[i];
      outfile << posx[i] << "\t" << posy[i] << endl;}
    outfile << posx[0] << "\t" << posy[0] << endl;
    outfile << endl;
  }
  file.close();
  outfile.close();

}

void dump_HGxynoise_avg(){
  CHMapping chmap;
  char title[100];
  sprintf(title,"./gnuplot/HGSCA%d.dat",SCAtodraw);
  ofstream outfile(title);
  
  outfile << "# X\tY\tavg\tsigma\tchip\tch" << endl;
  double norm_fac = -1;
  
  int chip, ch, formatCH;
  float max = -1;
  
  for (chip =0; chip < 4; chip = chip +1 ){
    for (ch = 0; ch < 64; ch = ch +1 ){
      formatCH = chip*32+ch/2;
      // remove unconnected CH and the inexist CH (chip 2 ch 60)
      if( DrawConnectCH && (ch %2 != 0 || formatCH == 94) ) continue;
      if( !DrawConnectCH && (ch %2 == 0 || formatCH == 94) ) continue;
      outfile << setprecision(6) << chmap.CH_x[formatCH]
	      << "\t" << chmap.CH_y[formatCH] << "\t"
	      << setprecision(3) << avg_HG[chip][ch][SCAtodraw] << "\t"
	      << sigma_HG[chip][ch][SCAtodraw] << "\t"
	      << chip << "\t" << ch << endl;

    }
  }
      outfile.close();  
}

void dump_LGxynoise_avg(){
  CHMapping chmap;
  char title[100];
  sprintf(title,"./gnuplot/LGSCA%d.dat",SCAtodraw);
  ofstream outfile(title);
  
  outfile << "# X\tY\tavg\tsigma\tchip\tch" << endl;
  double norm_fac = -1;
  
  int chip, ch, formatCH;
  float max = -1;
  
  for (chip =0; chip < 4; chip = chip +1 ){
    for (ch = 0; ch < 64; ch = ch +1 ){
      formatCH = chip*32+ch/2;
      // remove unconnected CH and the inexist CH (chip 2 ch 60)
      if( ch %2 != 0 || formatCH == 94) continue;
      outfile << setprecision(6) << chmap.CH_x[formatCH]
	      << "\t" << chmap.CH_y[formatCH] << "\t"
	      << setprecision(3) << avg_LG[chip][ch][SCAtodraw] << "\t"
	      << sigma_LG[chip][ch][SCAtodraw] << "\t"
	      << chip << "\t" << ch << endl;

    }
  }
      outfile.close();  
}

void dump_ped_check(){

  ofstream file("./gnuplot/ped_check.dat");
  file << "#Draw SCA " << SCAtodraw << " to check..."<< endl;
  file << "#CH num\tHGADC\tLGADC" << endl;
  int chip,ch;
  for(chip = 0; chip < 4;++chip){
    for(ch = 0; ch < 64; ++ch){      
      file << chip*64+ch << "\t" << avg_HG[chip][ch][SCAtodraw]
	   << "\t" << avg_LG[chip][ch][SCAtodraw] << endl;
    }
  }
  file.close();
}

void dump_gnuplotconfig(string fileN){
  int start = fileN.find_last_of("/");
  int end   = fileN.find(".raw");
  string Name = fileN.substr(start+1,end-start-1);
  
  ofstream file("./gnuplot/plotconfig.txt");
  file << Name << endl;
  if(DrawConnectCH) file << "CC" << endl;
  else file << "UC" << endl;
  if(HGorLGtodraw == 0)      file << "HG" << endl;
  else if(HGorLGtodraw == 1) file << "LG" << endl;
  else file << "HL" << endl;
  file << SCAtodraw << endl;
  file.close();
}

void Inj_ana(){
  int N_inj  = 1;
  int tar_CH = 36;
  int tar_TS = 4;
  
  string logcontent;
  ifstream logfile("./runinfo.log");
  logfile >> logcontent;
  char fileName[100];
  sprintf(fileName,"%s",logcontent.c_str());
  ifstream file(fileName);

  //Find inj.txt
  int end = logcontent.find(".raw");
  string fn = logcontent.substr(0,end);
  sprintf(fileName,"%s_Inj.txt",fn.c_str());
  ifstream file_inj(fileName);
  if(!file_inj.is_open()){
    cout << "Can't find file " << fileName << ", exit Injection analysis!"
	 << "\n This run may not be an injection run!" << endl;
    return;}
  string dum_line;
  for(int i = 0 ; i < 5; ++i) {
    if(i == 1) file_inj >> tar_CH;
    getline(file_inj,dum_line);    
  }
  ofstream fout("gnuplot/dac_adc.dat");
  fout << "#evt\tdac\t";
  for(int i = 0 ; i < 4 ; ++i){
    fout << "c" << i << ",C" << tar_CH << "HG\t";
    fout << "c" << i << ",C" << tar_CH << "LG\t";}
  fout << "\n";

  cout << "dumping CH " << tar_CH << " TS " << tar_TS << " for all 4 chips to "
       << "gnuplot/dac_adc.dat ...." << endl;
  int chip,ch,hit;
  
  if (file.is_open()){
    while(true){
      unsigned char testeof;
      file >> testeof;
      if( file.eof() ) break;
      else{
	testeof = raw[0];
	for(int i = 1 ; i < RAWSIZE; ++i)
	  file >> raw[i];
	decode_raw();
	format_channels();
	int rollpos = roll_position();
	float sig_evt[4][N_inj*2];

	for (chip =0; chip < 4; chip = chip +1 ){
	  int fake_CH_HG = 127 - tar_CH;
	  int fake_CH_LG = 63 - tar_CH;
	  for (hit = 0 ; hit < nSCA ; hit = hit +1){
	    int ts;
	    if(hit >= rollpos) ts = hit - rollpos;
	    else ts = 12+hit-rollpos+1;
	    if (ts == tar_TS){
	      sig_evt[chip][0] = dati[chip][fake_CH_HG][hit];
	      //sig_evt[chip][0] =
	      //dati[chip][fake_CH_HG][hit] - avg_HG [chip][tar_CH][hit];
	      sig_evt[chip][1] = dati[chip][fake_CH_LG][hit];
	      //sig_evt[chip][1] =
	      //dati[chip][fake_CH_LG][hit] - avg_LG [chip][tar_CH][hit];
	    }}}

	int nevt,dac;
	file_inj >> nevt >> dac;
	fout << nevt << "\t" << dac << "\t";
	for(int i = 0 ; i < 4 ; ++i){
	  fout << setprecision(4) << sig_evt[i][0] << "\t"
	       << sig_evt[i][1] << "\t";
	}
	//TOTS
	for(int i = 0 ; i < 4 ; ++i){
	  int fake_CH_HG = 127 - tar_CH;
	  fout << setprecision(4) << dati[i][fake_CH_HG][14] << "\t";}
	fout << endl;
      }
    }
  }
  file.close();
  file_inj.close();
  fout.close();
  
}
