#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
using namespace std;

void analysisClass::printBits(unsigned int num){
  unsigned int size = sizeof(num);
  unsigned int maxPow = 1<<(size*8-1);
  //printf("MAX POW : %u\n",maxPow);
  int i=0;
  for(;i<size*8;++i){
    // print last bit and shift left.
    printf("%u ",num&maxPow ? 1 : 0);
    num = num<<1;
  }
}

int analysisClass::getBit( unsigned int num,  int N ){
  unsigned int size = sizeof(num);
  unsigned int maxPow = 1<<(size*8-1);
  //
  int i=0;
  for(;i<size*8;++i){
    if( i==size*8-N && ( num&maxPow ? 1 : 0  ) ) return 1;
    num = num<<1;
  }
  return 0;
}

/*
int analysisClass::RBXnoiseCtr( TH1D*&histo ){
  // this function checks for RBX noise in a 1D histo stored using the iEta/iPhi/iDepth-to-RBXrm mapping.
  // ex: HBMphiHitsPerEvent_Histo
  //
  int badHPDperRBXctr=0;
  //
  for( unsigned int irbx=1; irbx<=18; irbx++){
    if( histo.GetBinContent( int( (irbx-1)*4+1 ) )>filterCut ) badHPDperRBXctr++;
    if( histo.GetBinContent( int( (irbx-1)*4+2 ) )>filterCut ) badHPDperRBXctr++;
    if( histo.GetBinContent( int( (irbx-1)*4+3 ) )>filterCut ) badHPDperRBXctr++;
    if( histo.GetBinContent( int( (irbx-1)*4+4 ) )>filterCut ) badHPDperRBXctr++;
  }

  return badHPDperRBXctr;
}
*/
