#ifndef analysisClass_h
#define analysisClass_h

#include "baseClass.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <memory>

using namespace std;

class analysisClass : public baseClass {
public :
  analysisClass(string * inputList, string * cutFile, string * treeName,  string *outputFileName=0, string * cutEfficFile=0);
  virtual ~analysisClass();
  void Loop();

  // Define Functions
  //

  TFile* HBMmapfile;
  TFile* HEMmapfile;
  TFile* HBPmapfile;
  TFile* HEPmapfile;
  //
  TTree* HBMmaptree;
  TTree* HEMmaptree;
  TTree* HBPmaptree;
  TTree* HEPmaptree;
  //
  int    iEta_HBM, iPhi_HBM, RBX_HBM, RM_HBM;
  int    iEta_HEM, iPhi_HEM, RBX_HEM, RM_HEM;
  int    iEta_HBP, iPhi_HBP, RBX_HBP, RM_HBP;
  int    iEta_HEP, iPhi_HEP, RBX_HEP, RM_HEP;
  //
  void   ReadHCALmaps( );
  void   CloseHCALmaps( );
  //
  int HBM_EtaPhitoRBXrm( int, int );
  int HEM_EtaPhitoRBXrm( int, int );
  int HBP_EtaPhitoRBXrm( int, int );
  int HEP_EtaPhitoRBXrm( int, int );
  //
  int EtaPhitoRBXrm( int, int );

};

#endif

#ifdef analysisClass_cxx

#endif // #ifdef analysisClass_cxx
