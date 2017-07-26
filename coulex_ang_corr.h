// coulex_ang_dist.h

//ROOT stuff
#include "TRandom3.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TFile.h"
#include "TApplication.h"
#include "TH1.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;

#define S32K    32768
#define NSPECT  10
#define RAD2DEG 57.2957795131
#define DEG2RAD 0.01745329251
#define AMU     931.4940954
#define MXS     -82.4393359375 // mass excess in MeV for 84Kr
#define EMASS   0.5109989461   // electron mass in MeV
#define PI      3.14159265359
#define PI2     1.57079632679
#define CONST   1.1283791671   // 4/sqrt(4*pi)

// functions
double calcWeight(int);
double getXi(double,double); // from the Coulex Geant4 code

// calculate VJM(theta,xi) for given input
// J,M,theta,xi using bilinear interpolation
double calcVJMthxi(int,int,double,double); 

// components of weight calculation
double calcUJM(int,int,double,double);
double calcWJM(int,int,double,double);

double getCoeff(int);
double wignerSmallD(int,int,int);
double associatedLegendre(int,int,double);
void readBranchData();
int factorial(int);

// for output of spectra with new weights
void addTreeDataToOutHist(int,double);
void readGroupMap();
double FWHM_response(double);


// ROOT magic for loading libraries
TApplication *theApp;

// file I/O
FILE *config;
char cfgstr[256],str1[256],str2[256];
char inp_filename[256],out_filename[256],tree_name[256],sort_path[256];
char prIn_name[256],prOut_name[256],rrOut_name[256],prDec_name[256];
char gx_name[256],gy_name[256],gz_name[256];
double Ap,Zp,Ar,Zr; // projectile and recoil A, Z
double E0; // projectile excitation energy in MeV
double A,DEprime;
char str[256];
char group_file[256],pos_path[256],col_path[256],csi_path[256]; // group sorting stuff
bool fwhmResponse=false;  // whether to do energy convolution
double fwhmF,fwhmG,fwhmH; // energy convolution parameters
double sort_scaling=1.0;  // gate scaling factor

// trees to extract data from
TTree *tree;
// leaves
TLeaf *lprIn_E,*lprIn_px,*lprIn_py,*lprIn_pz,*lprIn_x,*lprIn_y,*lprIn_z,*lprIn_b; // projectile before reaction
TLeaf *lprOut_E,*lprOut_px,*lprOut_py,*lprOut_pz,*lprOut_x,*lprOut_y,*lprOut_z,*lprOut_b,*lprOut_w; // proj after reaction
TLeaf *lprDec_E,*lprDec_px,*lprDec_py,*lprDec_pz,*lprDec_x,*lprDec_y,*lprDec_z,*lprDec_b,*lprDec_df; // proj decay
TLeaf *lgx,*lgy,*lgz; // gamma first interaction position (x,y,z)

// pre-computed tables
TRandom3 *randGen;

// for output of group spectra
TLeaf *sortLeaf,*gateLeaf;
TBranch *sortBranch,*gateBranch;
double newWeight;
int  group_map[17][5][25];
TLeaf *posLeaf,*colLeaf,*csiLeaf;
TBranch *posBranch,*colBranch,*csiBranch;
float dOutHist[NSPECT][S32K];
