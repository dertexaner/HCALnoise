#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
//
#include "HCALmap.C"

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

  std::cout << "analysisClass::analysisClass(): ends " << std::endl;
}

analysisClass::~analysisClass()
{
  std::cout << "analysisClass::~analysisClass(): begins " << std::endl;

  std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   
    
   if (fChain == 0) return;

   //////////read HCALmaps
   ReadHCALmaps( );


   //////////book histos here
   //
   //HB phiHits: each phi segment has 16+2 etas feeding into it (similarly for plus and minus sides)
   //If at a given phi, an eta channel has nonzero energy, it is considered 
   //
   TH1D* HBMphiHitsPerEvent_Histo = new TH1D("HBMphiHitsPerEvent_Histo","HBMphiHitsPerEvent_Histo",72, 0.5, 72.5);
   HBMphiHitsPerEvent_Histo->Sumw2();
   TH1D* HBPphiHitsPerEvent_Histo = new TH1D("HBPphiHitsPerEvent_Histo","HBPphiHitsPerEvent_Histo",72, 0.5, 72.5);
   HBPphiHitsPerEvent_Histo->Sumw2();
   TH1D* HBMphiHits_Histo = new TH1D("HBMphiHits_Histo","HBMphiHits_Histo",72, 0.5, 72.5); //minus side
   HBMphiHits_Histo->Sumw2();
   TH1D* HBPphiHits_Histo = new TH1D("HBPphiHits_Histo","HBPphiHits_Histo",72, 0.5, 72.5); //plus side
   HBPphiHits_Histo->Sumw2();

   TH1D* HEMphiHitsPerEvent_Histo = new TH1D("HEMphiHitsPerEvent_Histo","HEMphiHitsPerEvent_Histo",72, 0.5, 72.5);
   HEMphiHitsPerEvent_Histo->Sumw2();
   TH1D* HEPphiHitsPerEvent_Histo = new TH1D("HEPphiHitsPerEvent_Histo","HEPphiHitsPerEvent_Histo",72, 0.5, 72.5);
   HEPphiHitsPerEvent_Histo->Sumw2();
   TH1D* HEMphiHits_Histo = new TH1D("HEMphiHits_Histo","HEMphiHits_Histo",72, 0.5, 72.5); //minus side
   HEMphiHits_Histo->Sumw2();
   TH1D* HEPphiHits_Histo = new TH1D("HEPphiHits_Histo","HEPphiHits_Histo",72, 0.5, 72.5); //plus side
   HEPphiHits_Histo->Sumw2();

   TH1D* BarrelHPDoccupancy_Histo = new TH1D("BarrelHPDoccupancy_Histo","BarrelHPDoccupancy_Histo",19,-0.5,18.5);
   BarrelHPDoccupancy_Histo->Sumw2();
   TH1D* EndcapHPDoccupancy_Histo = new TH1D("EndcapHPDoccupancy_Histo","EndcapHPDoccupancy_Histo",19,-0.5,18.5);
   EndcapHPDoccupancy_Histo->Sumw2();

   TH1D* HBMHPDoccupancy_Histo = new TH1D("HBMHPDoccupancy_Histo","HBMHPDoccupancy_Histo",19,-0.5,18.5);
   HBMHPDoccupancy_Histo->Sumw2();
   TH1D* HBPHPDoccupancy_Histo = new TH1D("HBPHPDoccupancy_Histo","HBPHPDoccupancy_Histo",19,-0.5,18.5);
   HBPHPDoccupancy_Histo->Sumw2();
   TH1D* HEMHPDoccupancy_Histo = new TH1D("HEMHPDoccupancy_Histo","HEMHPDoccupancy_Histo",19,-0.5,18.5);
   HEMHPDoccupancy_Histo->Sumw2();
   TH1D* HEPHPDoccupancy_Histo = new TH1D("HEPHPDoccupancy_Histo","HEPHPDoccupancy_Histo",19,-0.5,18.5);
   HEPHPDoccupancy_Histo->Sumw2();

   TH1D* badHPDperPileup_Histo  =  new TH1D("badHPDperPileup_Histo","badHPDperPileup_Histo",10,0.5,50.5);
   badHPDperPileup_Histo->Sumw2();
   TH1D* NormalizedbadHPDperPileup_Histo  =  new TH1D("NormalizedbadHPDperPileup_Histo","NormalizedbadHPDperPileup_Histo",10,0.5,50.5);
   NormalizedbadHPDperPileup_Histo->Sumw2();
   TH1D* Pileup_Histo  =  new TH1D("Pileup_Histo","Pileup_Histo",10,0.5,50.5);
   Pileup_Histo->Sumw2();

   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////
     ///Stuff to be done every event

     //for( int ch=0; ch<PulseCount; ch++){
     //std::cout<<"input2: "<<IEta[ch]<<"/"<<IPhi[ch]<<" :: "<<EtaPhitoRBXrm(IEta[ch],IPhi[ch])<<std::endl;
     //}

     HBMphiHitsPerEvent_Histo->Reset("MICES");
     HBPphiHitsPerEvent_Histo->Reset("MICES");
     HEMphiHitsPerEvent_Histo->Reset("MICES");
     HEPphiHitsPerEvent_Histo->Reset("MICES");
     //
     for( int ch=0; ch<PulseCount; ch++){
       if(IEta[ch]>  0 && IEta[ch]< 17 ) HBPphiHitsPerEvent_Histo->Fill(HBP_EtaPhitoRBXrm(IEta[ch],IPhi[ch]),1); //HBP
       if(IEta[ch]<  0 && IEta[ch]>-17 ) HBMphiHitsPerEvent_Histo->Fill(HBM_EtaPhitoRBXrm(IEta[ch],IPhi[ch]),1); //HBM
       if(IEta[ch]> 16 && IEta[ch]< 30 ) HEPphiHitsPerEvent_Histo->Fill(HEP_EtaPhitoRBXrm(IEta[ch],IPhi[ch]),1); //HEP
       if(IEta[ch]<-16 && IEta[ch]>-30 ) HEMphiHitsPerEvent_Histo->Fill(HEM_EtaPhitoRBXrm(IEta[ch],IPhi[ch]),1); //HEM
     }


     int badHPDctr=0;
     const int filtercut=17;
     for( unsigned int ihpd=1; ihpd<73; ihpd++){
       EndcapHPDoccupancy_Histo->Fill(HEPphiHitsPerEvent_Histo->GetBinContent(ihpd));
       EndcapHPDoccupancy_Histo->Fill(HEMphiHitsPerEvent_Histo->GetBinContent(ihpd));
       HEPHPDoccupancy_Histo->Fill(HEPphiHitsPerEvent_Histo->GetBinContent(ihpd));
       HEMHPDoccupancy_Histo->Fill(HEMphiHitsPerEvent_Histo->GetBinContent(ihpd));
       BarrelHPDoccupancy_Histo->Fill(HBMphiHitsPerEvent_Histo->GetBinContent(ihpd));
       BarrelHPDoccupancy_Histo->Fill(HBPphiHitsPerEvent_Histo->GetBinContent(ihpd));
       HBMHPDoccupancy_Histo->Fill(HBMphiHitsPerEvent_Histo->GetBinContent(ihpd));
       HBPHPDoccupancy_Histo->Fill(HBPphiHitsPerEvent_Histo->GetBinContent(ihpd));

       if( HEPphiHitsPerEvent_Histo->GetBinContent(ihpd)>filtercut ){
	 HEPphiHits_Histo->Fill(ihpd,1);
	 //print "endcap+  hdp ",HEPphiHitsPerEvent_Histo->GetBinContent(ihpd)," hits";
	 badHPDctr+=1;
       }
       if( HEMphiHitsPerEvent_Histo->GetBinContent(ihpd)>filtercut ){
	 HEMphiHits_Histo->Fill(ihpd,1);
	 //print "endcap-  hdp ",HEMphiHitsPerEvent_Histo->GetBinContent(ihpd)," hits";
	 badHPDctr+=1;
       }
       if( HBPphiHitsPerEvent_Histo->GetBinContent(ihpd)>filtercut ){
	 HBPphiHits_Histo->Fill(ihpd,1);
	 //print "barrel+  hdp ",HBPphiHitsPerEvent_Histo->GetBinContent(ihpd)," hits";
	 badHPDctr+=1;
       }
       if( HBMphiHitsPerEvent_Histo->GetBinContent(ihpd)>filtercut ){
	 HBMphiHits_Histo->Fill(ihpd,1);
	 //print "barrel-  hdp ",HBMphiHitsPerEvent_Histo->GetBinContent(ihpd)," hits";
	 badHPDctr+=1;
       }
     }
      
     badHPDperPileup_Histo->Fill(NumberOfGoodPrimaryVertices,badHPDctr);
     Pileup_Histo->Fill(NumberOfGoodPrimaryVertices,1);

     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   NormalizedbadHPDperPileup_Histo->Divide(badHPDperPileup_Histo, Pileup_Histo,1,1);


   //TFile* newoutfile = new TFile("/afs/cern.ch/user/h/hsaka/HCALwork/CMSSW_5_3_7_patch6/src/HcalNoise/HcalNoiseAnalyzer/Ccode/rootNtupleAnalyzerV2/data/output/newrootFile.root","RECREATE");

   //////////write histos 
   output_root_->cd();
   //newoutfile->cd();

   std::cout<<"close2B"<<std::endl;

   HBMphiHits_Histo->Write();
   HBPphiHits_Histo->Write();
   
   HEMphiHits_Histo->Write();
   HEPphiHits_Histo->Write();
   
   BarrelHPDoccupancy_Histo->Write();
   EndcapHPDoccupancy_Histo->Write();

   HBPHPDoccupancy_Histo->Write();
   HBMHPDoccupancy_Histo->Write();
   HEPHPDoccupancy_Histo->Write();
   HEMHPDoccupancy_Histo->Write();

   badHPDperPileup_Histo->Write();
   NormalizedbadHPDperPileup_Histo->Write();
   Pileup_Histo->Write();
 
   CloseHCALmaps( );

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
