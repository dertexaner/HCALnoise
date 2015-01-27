#define analysisClass_cxx
#include "analysisClass.h"
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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

   //////////debug flag
   bool debug_  = true;
   bool debug2_ = false;
   bool debug3_ = false;

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

   TH1D* barrelbadHPDperPileup_Histo  =  new TH1D("barrelbadHPDperPileup_Histo","barrelbadHPDperPileup_Histo",20,0.5,100.5);
   barrelbadHPDperPileup_Histo->Sumw2();
   TH1D* barrelallHPDperPileup_Histo  =  new TH1D("barrelallHPDperPileup_Histo","barrelallHPDperPileup_Histo",20,0.5,100.5);
   barrelallHPDperPileup_Histo->Sumw2();
   TH1D* endcapbadHPDperPileup_Histo  =  new TH1D("endcapbadHPDperPileup_Histo","endcapbadHPDperPileup_Histo",20,0.5,100.5);
   endcapbadHPDperPileup_Histo->Sumw2();
   TH1D* endcapallHPDperPileup_Histo  =  new TH1D("endcapallHPDperPileup_Histo","endcapallHPDperPileup_Histo",20,0.5,100.5);
   endcapallHPDperPileup_Histo->Sumw2();

   TH1D* badHPDperPileup_Histo  =  new TH1D("badHPDperPileup_Histo","badHPDperPileup_Histo",20,0.5,100.5);
   badHPDperPileup_Histo->Sumw2();
   TH1D* allHPDperPileup_Histo  =  new TH1D("allHPDperPileup_Histo","allHPDperPileup_Histo",20,0.5,100.5);
   allHPDperPileup_Histo->Sumw2();

   TH1D* barrelNormbadHPDperPileup_Histo  =  new TH1D("barrelNormbadHPDperPileup_Histo","barrelNormbadHPDperPileup_Histo",20,0.5,100.5);
   barrelNormbadHPDperPileup_Histo->Sumw2();
   TH1D* endcapNormbadHPDperPileup_Histo  =  new TH1D("endcapNormbadHPDperPileup_Histo","endcapNormbadHPDperPileup_Histo",20,0.5,100.5);
   endcapNormbadHPDperPileup_Histo->Sumw2();

   TH1D* NormbadHPDperPileup_Histo  =  new TH1D("NormbadHPDperPileup_Histo","NormbadHPDperPileup_Histo",20,0.5,100.5);
   NormbadHPDperPileup_Histo->Sumw2();

   TH1D* Pileup_Histo  =  new TH1D("Pileup_Histo","Pileup_Histo",20,0.5,100.5);
   Pileup_Histo->Sumw2();

   TH1D* verifyHPDhits_Histo  =  new TH1D("verifyHPDhits_Histo","verifyHPDhits_Histo",4,0.5,4.5);
   verifyHPDhits_Histo->Sumw2();

   TH3D* plusEtaPhiEnergyMap_Histo  = new TH3D( "plusEtaPhiEnergyMap_Histo",  "plusEtaPhiEnergyMap_Histo", 29,0.5,29.5, 72,0.5,72.5, 3,0.5,3.5 );
   TH3D* minusEtaPhiEnergyMap_Histo = new TH3D("minusEtaPhiEnergyMap_Histo", "minusEtaPhiEnergyMap_Histo", 29,0.5,29.5, 72,0.5,72.5, 3,0.5,3.5 );


   /////////initialize variables
   int pileup;

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=6000; jentry<nentries-1; jentry++) {
     //
     //if( jentry%2==0 ) continue; // skip even events
     //
     // load event "jentry"
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     //
     if(jentry < 10 || (jentry-1)%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     if( (jentry+1)%15000 == 0 ) break;

     ////////////////////// User's code starts here ///////////////////////
     ///Stuff to be done every event

     if( debug_ ){
       if( HPDHits>17 ) std::cout<<"This event has HPDHits: "<<HPDHits<<std::endl;
     }

     if( debug3_ ){
       for( int ch=0; ch<PulseCount; ch++){
	 std::cout<<"Eta/Phi: "<<IEta[ch]<<"/"<<IPhi[ch]<<" :: RBXrm: "<<EtaPhitoRBXrm(IEta[ch],IPhi[ch])<<std::endl;
       }
     }

     plusEtaPhiEnergyMap_Histo->Reset("MICES");
     minusEtaPhiEnergyMap_Histo->Reset("MICES");

     for( int ch=0; ch<PulseCount; ch++){
       //
       if( debug2_ ){ if( ch==0 ){ std::cout<<"Ch=0 for event  : "<<IEta[ch]<<" "<<IPhi[ch]<<" "<<Energy[ch]<<std::endl; } }
       //
       if( IEta[ch]>0 ) plusEtaPhiEnergyMap_Histo->Fill(IEta[ch],IPhi[ch],Depth[ch],Energy[ch]);
       if( IEta[ch]<0 ) minusEtaPhiEnergyMap_Histo->Fill(abs(IEta[ch]),IPhi[ch],Depth[ch],Energy[ch]);
       //
     }

     pileup=0;
     pileup+=NumberOfGoodPrimaryVertices; // add pileup from jentry


     /*
     // load event "jentry+1"
     ientry = LoadTree(jentry+1);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry+1);   nbytes += nb;
     //
     for( int ch=0; ch<PulseCount; ch++){
       //
       if( debug2_ ){ if( ch==0 ){ std::cout<<"Ch=0 for event+1: "<<IEta[ch]<<" "<<IPhi[ch]<<" "<<Energy[ch]<<std::endl; std::cout<<std::endl; } }
       //
       if( IEta[ch]>0 ) plusEtaPhiEnergyMap_Histo->Fill(IEta[ch],IPhi[ch],Depth[ch],Energy[ch]);
       if( IEta[ch]<0 ) minusEtaPhiEnergyMap_Histo->Fill(abs(IEta[ch]),IPhi[ch],Depth[ch],Energy[ch]);
     }

     pileup+=NumberOfGoodPrimaryVertices; // add pileup from jentry+1
     */

     //map from 3D histo (for even mix) back to 4 histos.
     HBMphiHitsPerEvent_Histo->Reset("MICES");
     HBPphiHitsPerEvent_Histo->Reset("MICES");
     HEMphiHitsPerEvent_Histo->Reset("MICES");
     HEPphiHitsPerEvent_Histo->Reset("MICES");
     // ------------
     for( int ie=1; ie<=16; ie++){
       for( int ip=1; ip<=72; ip++){
	 for( int id=1; id<=2; id++){
	   //
	   //--- plus side
	   if( plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>1.5 ){
	     HBPphiHitsPerEvent_Histo->Fill(HBP_EtaPhitoRBXrm(ie,ip),1); //HBP
	   }
	   //
	   //--- minus side
	   if( minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>1.5 ){
	     HBMphiHitsPerEvent_Histo->Fill(HBM_EtaPhitoRBXrm(ie,ip),1); //HBM
	   }
	   //
	 }
       }
     }
     // ------------
     for( int ie=16; ie<=29; ie++){
       for( int ip=1; ip<=72; ip++){
	 for( int id=1; id<=3; id++){
	   //
	   if( ie==16 && id<3 ) continue; //eta 16 has first 2 depths in HB
	   //
	   //--- plus side
	   if( plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>1.5 ){
	     HEPphiHitsPerEvent_Histo->Fill(HEP_EtaPhitoRBXrm(ie,ip),1); //HEP
	   }
	   //
	   //--- minus side
	   if( minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>1.5 ){
	     HEMphiHitsPerEvent_Histo->Fill(HEM_EtaPhitoRBXrm(ie,ip),1); //HEM
	   }
	   //
	 }
       }
     }


     int badHPDctr=0;
     int barrelbadHPDctr=0;
     int endcapbadHPDctr=0;
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
	 if( debug_ ) std::cout<<"HEP:  HPD with "<<HEPphiHitsPerEvent_Histo->GetBinContent(ihpd)<<" hits"<<std::endl;
	 badHPDctr++;
	 endcapbadHPDctr++;
       }
       if( HEMphiHitsPerEvent_Histo->GetBinContent(ihpd)>filtercut ){
	 HEMphiHits_Histo->Fill(ihpd,1);
	 if( debug_ ) std::cout<<"HEM:  HPD with "<<HEMphiHitsPerEvent_Histo->GetBinContent(ihpd)<<" hits"<<std::endl;
	 badHPDctr++;
	 endcapbadHPDctr++;
       }
       if( HBPphiHitsPerEvent_Histo->GetBinContent(ihpd)>filtercut ){
	 HBPphiHits_Histo->Fill(ihpd,1);
	 if( debug_ ) std::cout<<"HBP:  HPD with "<<HBPphiHitsPerEvent_Histo->GetBinContent(ihpd)<<" hits"<<std::endl;
         badHPDctr++;
         barrelbadHPDctr++;
       }
       if( HBMphiHitsPerEvent_Histo->GetBinContent(ihpd)>filtercut ){
	 HBMphiHits_Histo->Fill(ihpd,1);
	 if( debug_ ) std::cout<<"HBM:  HPD with "<<HBMphiHitsPerEvent_Histo->GetBinContent(ihpd)<<" hits"<<std::endl;
	 badHPDctr++;
	 barrelbadHPDctr++;
       }
     }
      
     barrelbadHPDperPileup_Histo->Fill(pileup,barrelbadHPDctr);
     barrelallHPDperPileup_Histo->Fill(pileup,144);//36*4  36 RBXs, 4 RMs

     endcapbadHPDperPileup_Histo->Fill(pileup,endcapbadHPDctr);
     endcapallHPDperPileup_Histo->Fill(pileup,144);//36*4  36 RBXs, 4 RMs

     badHPDperPileup_Histo->Fill(pileup,badHPDctr);
     allHPDperPileup_Histo->Fill(pileup,288);//72*4  72 RBXs, 4 RMs

     Pileup_Histo->Fill(pileup,1);

     if( HPDHits> 17 && badHPDctr> 0 ) verifyHPDhits_Histo->Fill(1);//ok
     if( HPDHits<=17 && badHPDctr> 0 ) verifyHPDhits_Histo->Fill(2);//bad
     if( HPDHits> 17 && badHPDctr==0 ) verifyHPDhits_Histo->Fill(3);//bad
     if( HPDHits<=17 && badHPDctr==0 ) verifyHPDhits_Histo->Fill(4);//ok



     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   barrelNormbadHPDperPileup_Histo->Divide(barrelbadHPDperPileup_Histo, barrelallHPDperPileup_Histo,1,1,"cl=0.683 b(1,1) mode");
   endcapNormbadHPDperPileup_Histo->Divide(endcapbadHPDperPileup_Histo, endcapallHPDperPileup_Histo,1,1,"cl=0.683 b(1,1) mode");
   NormbadHPDperPileup_Histo->Divide(badHPDperPileup_Histo, allHPDperPileup_Histo,1,1,"cl=0.683 b(1,1) mode");


   //////////write histos 
   output_root_->cd();

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

   barrelbadHPDperPileup_Histo->Write();
   barrelallHPDperPileup_Histo->Write();
   endcapbadHPDperPileup_Histo->Write();
   endcapallHPDperPileup_Histo->Write();
   badHPDperPileup_Histo->Write();
   allHPDperPileup_Histo->Write();

   barrelNormbadHPDperPileup_Histo->Write();
   endcapNormbadHPDperPileup_Histo->Write();
   NormbadHPDperPileup_Histo->Write();

   Pileup_Histo->Write();

   verifyHPDhits_Histo->Write();
 
   CloseHCALmaps( );

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
