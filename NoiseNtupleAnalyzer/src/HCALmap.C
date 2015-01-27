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

void analysisClass::ReadHCALmaps( ){

  //input file
  HBMmapfile=new TFile("/afs/cern.ch/user/h/hsaka/public/HCALmaps/HBM_EtaPhiRBXRM.root","READ");
  HEMmapfile=new TFile("/afs/cern.ch/user/h/hsaka/public/HCALmaps/HEM_EtaPhiRBXRM.root","READ");
  HBPmapfile=new TFile("/afs/cern.ch/user/h/hsaka/public/HCALmaps/HBP_EtaPhiRBXRM.root","READ");
  HEPmapfile=new TFile("/afs/cern.ch/user/h/hsaka/public/HCALmaps/HEP_EtaPhiRBXRM.root","READ");
  //
  HBMmaptree=(TTree*)(HBMmapfile->Get("tree"));
  HEMmaptree=(TTree*)(HEMmapfile->Get("tree"));
  HBPmaptree=(TTree*)(HBPmapfile->Get("tree"));
  HEPmaptree=(TTree*)(HEPmapfile->Get("tree"));
  //
  HBMmaptree->SetBranchAddress("iEta",  &iEta_HBM  );
  HBMmaptree->SetBranchAddress("iPhi",  &iPhi_HBM  );
  HBMmaptree->SetBranchAddress("iDepth",&iDepth_HBM);
  HBMmaptree->SetBranchAddress("RBX" ,  &RBX_HBM   );
  HBMmaptree->SetBranchAddress("RM"  ,  &RM_HBM    );
  //
  HEMmaptree->SetBranchAddress("iEta",  &iEta_HEM  );
  HEMmaptree->SetBranchAddress("iPhi",  &iPhi_HEM  );
  HEMmaptree->SetBranchAddress("iDepth",&iDepth_HEM);
  HEMmaptree->SetBranchAddress("RBX" ,  &RBX_HEM   );
  HEMmaptree->SetBranchAddress("RM"  ,  &RM_HEM    );
  //
  HBPmaptree->SetBranchAddress("iEta",  &iEta_HBP  );
  HBPmaptree->SetBranchAddress("iPhi",  &iPhi_HBP  );
  HBPmaptree->SetBranchAddress("iDepth",&iDepth_HBP);
  HBPmaptree->SetBranchAddress("RBX" ,  &RBX_HBP   );
  HBPmaptree->SetBranchAddress("RM"  ,  &RM_HBP    );
  //
  HEPmaptree->SetBranchAddress("iEta",  &iEta_HEP  );
  HEPmaptree->SetBranchAddress("iPhi",  &iPhi_HEP  );
  HEPmaptree->SetBranchAddress("iDepth",&iDepth_HEP);
  HEPmaptree->SetBranchAddress("RBX" ,  &RBX_HEP   );
  HEPmaptree->SetBranchAddress("RM"  ,  &RM_HEP    );
  //
}


int analysisClass::HBM_EtaPhitoRBXrm( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HBMmaptree->GetEntries(); ientry++){
    HBMmaptree->GetEntry(ientry);
    if( iEta_HBM==abs(ieta) && iPhi_HBM==iphi && iDepth_HBM==idepth ){  output=int((RBX_HBM-int(1))*int(4)+RM_HBM); break; }
  }
  //
  return output;
}
int analysisClass::HBM_EtaPhitoRBX( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HBMmaptree->GetEntries(); ientry++){
    HBMmaptree->GetEntry(ientry);
    if( iEta_HBM==abs(ieta) && iPhi_HBM==iphi && iDepth_HBM==idepth ){  output=int(RBX_HBM); break; }
  }
  //
  return output;
}

int analysisClass::HEM_EtaPhitoRBXrm( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HEMmaptree->GetEntries(); ientry++){
    HEMmaptree->GetEntry(ientry);
    if( iEta_HEM==abs(ieta) && iPhi_HEM==iphi && iDepth_HEM==idepth ){  output=int((RBX_HEM-int(1))*int(4)+RM_HEM); break; }
  }
  //
  return output;
}
int analysisClass::HEM_EtaPhitoRBX( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HEMmaptree->GetEntries(); ientry++){
    HEMmaptree->GetEntry(ientry);
    if( iEta_HEM==abs(ieta) && iPhi_HEM==iphi && iDepth_HEM==idepth ){  output=int(RBX_HEM); break; }
  }
  //
  return output;
}

int analysisClass::HBP_EtaPhitoRBXrm( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HBPmaptree->GetEntries(); ientry++){
    HBPmaptree->GetEntry(ientry);
    if( iEta_HBP==ieta && iPhi_HBP==iphi && iDepth_HBP==idepth ){  output=int((RBX_HBP-int(1))*int(4)+RM_HBP); break; }
  }
  //
  return output;
}
int analysisClass::HBP_EtaPhitoRBX( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HBPmaptree->GetEntries(); ientry++){
    HBPmaptree->GetEntry(ientry);
    if( iEta_HBP==ieta && iPhi_HBP==iphi && iDepth_HBP==idepth ){  output=int(RBX_HBP); break; }
  }
  //
  return output;
}


int analysisClass::HEP_EtaPhitoRBXrm( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HEPmaptree->GetEntries(); ientry++){
    HEPmaptree->GetEntry(ientry);
    if( iEta_HEP==ieta && iPhi_HEP==iphi && iDepth_HEP==idepth ){  output=int((RBX_HEP-int(1))*int(4)+RM_HEP); break; }
  }
  //
  return output;
}
int analysisClass::HEP_EtaPhitoRBX( int ieta, int iphi, int idepth ){
  //
  int output=999;
  //
  for(unsigned int ientry=0; ientry<HEPmaptree->GetEntries(); ientry++){
    HEPmaptree->GetEntry(ientry);
    if( iEta_HEP==ieta && iPhi_HEP==iphi && iDepth_HEP==idepth ){  output=int(RBX_HEP); break; }
  }
  //
  return output;
}

int analysisClass::EtaPhitoRBXrm( int ieta, int iphi, int idepth ){
  //
  int output=888;
  //
  if(ieta>  0 && ieta< 17 && idepth<3 ) output=HBP_EtaPhitoRBXrm(ieta,iphi,idepth); //HBP
  if(ieta<  0 && ieta>-17 && idepth<3 ) output=HBM_EtaPhitoRBXrm(ieta,iphi,idepth); //HBM
  if(ieta== 16 && idepth==3 ) output=HEP_EtaPhitoRBXrm(ieta,iphi,idepth); //HEP
  if(ieta==-16 && idepth==3 ) output=HEM_EtaPhitoRBXrm(ieta,iphi,idepth); //HEM
  if(ieta> 16 && ieta< 30 ) output=HEP_EtaPhitoRBXrm(ieta,iphi,idepth); //HEP
  if(ieta<-16 && ieta>-30 ) output=HEM_EtaPhitoRBXrm(ieta,iphi,idepth); //HEM
  //
  return output;
}


void analysisClass::CloseHCALmaps( ){

  HBMmapfile->Close();
  HEMmapfile->Close();
  HBPmapfile->Close();
  HEPmapfile->Close();
  //

}


