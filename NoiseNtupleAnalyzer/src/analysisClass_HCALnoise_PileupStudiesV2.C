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
   bool debug_  = false;
   bool debug2_ = false;
   bool debug3_ = false;
   bool debug4_ = false;
   bool debug5_ = false;

   //////////initialize variables
   filterCut = 10;        // No. of pixels with allowed activity in an HPD.
   pixelEnergyCut = 1.5;  // Energy treshold (GeV) of a pixel to be considered as active.
   int pileup;

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

   TH1D* barrelbadHPDperPileup_Histo  =  new TH1D("barrelbadHPDperPileup_Histo","barrelbadHPDperPileup_Histo",60,0.5,300.5);
   barrelbadHPDperPileup_Histo->Sumw2();
   TH1D* barrelallHPDperPileup_Histo  =  new TH1D("barrelallHPDperPileup_Histo","barrelallHPDperPileup_Histo",60,0.5,300.5);
   barrelallHPDperPileup_Histo->Sumw2();
   TH1D* endcapbadHPDperPileup_Histo  =  new TH1D("endcapbadHPDperPileup_Histo","endcapbadHPDperPileup_Histo",60,0.5,300.5);
   endcapbadHPDperPileup_Histo->Sumw2();
   TH1D* endcapallHPDperPileup_Histo  =  new TH1D("endcapallHPDperPileup_Histo","endcapallHPDperPileup_Histo",60,0.5,300.5);
   endcapallHPDperPileup_Histo->Sumw2();

   TH1D* badHPDperPileup_Histo  =  new TH1D("badHPDperPileup_Histo","badHPDperPileup_Histo",60,0.5,300.5);
   badHPDperPileup_Histo->Sumw2();
   TH1D* allHPDperPileup_Histo  =  new TH1D("allHPDperPileup_Histo","allHPDperPileup_Histo",60,0.5,300.5);
   allHPDperPileup_Histo->Sumw2();

   TH1D* barrelNormbadHPDperPileup_Histo  =  new TH1D("barrelNormbadHPDperPileup_Histo","barrelNormbadHPDperPileup_Histo",60,0.5,300.5);
   barrelNormbadHPDperPileup_Histo->Sumw2();
   TH1D* endcapNormbadHPDperPileup_Histo  =  new TH1D("endcapNormbadHPDperPileup_Histo","endcapNormbadHPDperPileup_Histo",60,0.5,300.5);
   endcapNormbadHPDperPileup_Histo->Sumw2();

   TH1D* NormbadHPDperPileup_Histo  =  new TH1D("NormbadHPDperPileup_Histo","NormbadHPDperPileup_Histo",60,0.5,300.5);
   NormbadHPDperPileup_Histo->Sumw2();

   TH1D* Pileup_Histo  =  new TH1D("Pileup_Histo","Pileup_Histo",60,0.5,300.5);
   Pileup_Histo->Sumw2();

   TH1D* OfficialDecision_Histo  =  new TH1D("OfficialDecision_Histo","OfficialDecision_Histo",2,-0.5,1.5);
   OfficialDecision_Histo->Sumw2();
   TH1D* PulseCount_Histo  =  new TH1D("PulseCount_Histo","PulseCount_Histo",6001,-0.5,6000.5);
   PulseCount_Histo->Sumw2();
   TH1D* ChannelEnergy_Histo  =  new TH1D("ChannelEnergy_Histo","ChannelEnergy_Histo",2010,-5,1000);
   ChannelEnergy_Histo->Sumw2();
   TH1D* ChannelEnergyBarrel_Histo  =  new TH1D("ChannelEnergyBarrel_Histo","ChannelEnergyBarrel_Histo",2010,-5,1000);
   ChannelEnergyBarrel_Histo->Sumw2();
   TH1D* ChannelEnergyEndcap_Histo  =  new TH1D("ChannelEnergyEndcap_Histo","ChannelEnergyEndcap_Histo",2010,-5,1000);
   ChannelEnergyEndcap_Histo->Sumw2();

   TH1D* ChargeTS4TS5_Histo  =  new TH1D("ChargeTS4TS5_Histo","ChargeTS4TS5_Histo",9000,-500,8500);
   ChargeTS4TS5_Histo->Sumw2();
   TH1D* ChargeTS0TS9_Histo  =  new TH1D("ChargeTS0TS9_Histo","ChargeTS0TS9_Histo",9000,-500,8500);
   ChargeTS0TS9_Histo->Sumw2();
   TH1D* Energy_Histo  =  new TH1D("Energy_Histo","Energy_Histo",8200,-100,4000);
   Energy_Histo->Sumw2();

   TH1D* HPDHitsCalculated_Histo = new TH1D("HPDHitsCalculated_Histo","HPDHitsCalculated_Histo",21,-0.5,20.5);
   HPDHitsCalculated_Histo->Sumw2();
   TH1D* HPDHitsOfficial_Histo   = new TH1D("HPDHitsOfficial_Histo","HPDHitsOfficial_Histo",21,-0.5,20.5);
   HPDHitsOfficial_Histo->Sumw2();

   TH1D* verifyHPDhits_Histo  =  new TH1D("verifyHPDhits_Histo","verifyHPDhits_Histo",4,0.5,4.5);
   verifyHPDhits_Histo->Sumw2();

   TH1D* RBXnoiseCounter_Histo  =  new TH1D("RBXnoiseCounter_Histo","RBXnoiseCounter_Histo",2,0.5,2.5);
   RBXnoiseCounter_Histo->Sumw2();

   TH3D* plusEtaPhiEnergyMap_Histo  = new TH3D( "plusEtaPhiEnergyMap_Histo",  "plusEtaPhiEnergyMap_Histo", 29,0.5,29.5, 72,0.5,72.5, 3,0.5,3.5 );
   TH3D* minusEtaPhiEnergyMap_Histo = new TH3D("minusEtaPhiEnergyMap_Histo", "minusEtaPhiEnergyMap_Histo", 29,0.5,29.5, 72,0.5,72.5, 3,0.5,3.5 );

   TH2D* plusEtaPhiEnergyMapD1A_Histo  = new TH2D( "plusEtaPhiEnergyMapD1A_Histo",  "plusEtaPhiEnergyMapD1A_Histo", 29,0.5,29.5, 72,0.5,72.5);
   TH2D* plusEtaPhiEnergyMapD1B_Histo  = new TH2D( "plusEtaPhiEnergyMapD1B_Histo",  "plusEtaPhiEnergyMapD1B_Histo", 29,0.5,29.5, 72,0.5,72.5);
   TH2D* plusEtaPhiEnergyMapD1C_Histo  = new TH2D( "plusEtaPhiEnergyMapD1C_Histo",  "plusEtaPhiEnergyMapD1C_Histo", 29,0.5,29.5, 72,0.5,72.5);
   TH2D* plusEtaPhiEnergyMapD1D_Histo  = new TH2D( "plusEtaPhiEnergyMapD1D_Histo",  "plusEtaPhiEnergyMapD1D_Histo", 29,0.5,29.5, 72,0.5,72.5);

   TH1D* chBarrelOccupancyBefore_Histo  =  new TH1D("chBarrelOccupancyBefore_Histo","chBarrelOccupancyBefore_Histo",1501,-0.5,1500.5);
   chBarrelOccupancyBefore_Histo->Sumw2();
   TH1D* chEndcapOccupancyBefore_Histo  =  new TH1D("chEndcapOccupancyBefore_Histo","chEndcapOccupancyBefore_Histo",1501,-0.5,1500.5);
   chEndcapOccupancyBefore_Histo->Sumw2();
   TH1D* chBarrelOccupancyAfter_Histo  =  new TH1D("chBarrelOccupancyAfter_Histo","chBarrelOccupancyAfter_Histo",1501,-0.5,1500.5);
   chBarrelOccupancyAfter_Histo->Sumw2();
   TH1D* chEndcapOccupancyAfter_Histo  =  new TH1D("chEndcapOccupancyAfter_Histo","chEndcapOccupancyAfter_Histo",1501,-0.5,1500.5);
   chEndcapOccupancyAfter_Histo->Sumw2();

   /////////Mixing scenario
   // Set ONLY one of these to TRUE ... or all of them to FALSE to get "NoMixing"
   bool mixPlus1_=false; 
   bool mixPlus2_=false;
   bool mixPlus3_=false;

   /////////initialize variables

   //Long64_t nentries = fChain->GetEntriesFast();
   long nentries = fChain->GetEntries();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   //Long64_t startEvent=1;
   //Long64_t nbytes = 0, nb = 0;
   //for (Long64_t jentry=startEvent; jentry<nentries-3;jentry++) {
   long startEvent=1;
   long nbytes = 0, nb = 0;
   for ( long jentry=startEvent; jentry<nentries-5; jentry++ ) {
     //
     if( jentry%4!=1 )  continue; // [1] 2 3 4, [5] 6 7 8, [9] 10 11 12,... // Keep the total number of events constant thru out the mixing scenarios
     //
     //if(  mixPlus1_ && !mixPlus2_ && !mixPlus3_ &&   jentry%2==0 )  continue; // [1] 2, [3] 4, [5] 6, [7] 8,...
     //if( !mixPlus1_ &&  mixPlus2_ && !mixPlus3_ &&   jentry%3!=1 )  continue; // [1] 2 3, [4] 5 6, [7] 8 9,...
     //if( !mixPlus1_ && !mixPlus2_ &&  mixPlus3_ &&   jentry%4!=1 )  continue; // [1] 2 3 4, [5] 6 7 8, [9] 10 11 12,...
     //
     // Old setup
     //if( jentry%2==0 ) continue; // skip even events
     //if( jentry%3!=0 ) continue; // [0] 1 2, [3] 4 5, [6] 7 8,...
     //if( jentry%4!=0 ) continue; // [0] 1 2 3, [4] 5 6 7, [8] 9 10 11,...
     //
     //Long64_t ientry = LoadTree(jentry);
     long ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     //std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////
     ///Stuff to be done every event

     // Event Filter:
     if( !OfficialDecision || NumberOfGoodPrimaryVertices<1 ) continue;

     if( debug5_ &&  Energy[1]>5){
       std::cout << "GeV & Charge"<< Energy[1]<<" & "<<(Charge[1][4]+Charge[1][5])<<std::endl;       
       std::cout << "Charge/Energy"<< (Charge[1][4]+Charge[1][5])/Energy[1]<<"  <<"<<std::endl;
     }

     //std::cout << "analysisClass::Loop(): jentry = BASE:: " << jentry <<"   PU:"<< NumberOfGoodPrimaryVertices<<" pulseCount: "<<PulseCount<<std::endl;
     //std::cout << "analysisClass::Loop(): jentry = BASE:: " << jentry <<"   PU:"<< NumberOfGoodPrimaryVertices<< std::endl;
	 
     // Some basic plots before any mixing..
     if(  OfficialDecision ) OfficialDecision_Histo->Fill(1);
     if( !OfficialDecision ) OfficialDecision_Histo->Fill(0);
     PulseCount_Histo->Fill(PulseCount);
     for( int ch=0; ch<PulseCount; ch++){
       ChannelEnergy_Histo->Fill(Energy[ch]);
       if( debug2_ ){ if( Energy[ch]<0.5 ) std::cout<<"Low Energy per Channel: ("<<jentry<<")  "<<Energy[ch]<<std::endl; }
     }

     double totalCharge=0;
     for( int ch=0; ch<PulseCount; ch++){
       ChargeTS4TS5_Histo->Fill( (Charge[ch][4]+Charge[ch][5]) );
       totalCharge=0;
       for( int iTS=0; iTS<10; iTS++ ){
	 totalCharge+=Charge[ch][iTS];
       }
       ChargeTS0TS9_Histo->Fill(totalCharge);
       Energy_Histo->Fill(Energy[ch]);
     }

     //There is also these stored without any mixing (splitting ChannelEnergy_Histo into two):
     // -- ChannelEnergyBarrel_Histo
     // -- ChannelEnergyEndcap_Histo
     // see below..

     int HPDHits_jentry=HPDHits;

     if( debug_ ){
       if( HPDHits_jentry>17 ) std::cout<<"V2 This event has HPDHits: "<<HPDHits_jentry<<" :: jentry "<<jentry<<std::endl;
     }

     if( debug2_ ){
       for( int ch=0; ch<PulseCount; ch++){
	 std::cout<<"Eta/Phi/Depth: "<<IEta[ch]<<"/"<<IPhi[ch]<<"/"<<Depth[ch]<<" :: RBXrm: "<<EtaPhitoRBXrm(IEta[ch],IPhi[ch],Depth[ch])<<std::endl;
       }
     }

   
     plusEtaPhiEnergyMap_Histo->Reset("MICES");
     minusEtaPhiEnergyMap_Histo->Reset("MICES");
     //
     plusEtaPhiEnergyMapD1A_Histo->Reset("MICES");
     plusEtaPhiEnergyMapD1B_Histo->Reset("MICES");
     plusEtaPhiEnergyMapD1C_Histo->Reset("MICES");
     plusEtaPhiEnergyMapD1D_Histo->Reset("MICES");
     

     for( int ch=0; ch<PulseCount; ch++){
       //
       if( debug3_ ){ if( ch==0 ){ std::cout<<"Ch=0 for event  : "<<IEta[ch]<<" "<<IPhi[ch]<<" "<<Energy[ch]<<std::endl; } }
       //
       if( IEta[ch]>0 ) plusEtaPhiEnergyMap_Histo->Fill(IEta[ch],IPhi[ch],Depth[ch],Energy[ch]);
       if( IEta[ch]<0 ) minusEtaPhiEnergyMap_Histo->Fill(abs(IEta[ch]),IPhi[ch],Depth[ch],Energy[ch]);
       //
       //Adding energy per channel in the first layer
       if( IEta[ch]>0 && Depth[ch]==1 ) plusEtaPhiEnergyMapD1A_Histo->Fill(IEta[ch],IPhi[ch],Energy[ch]);
       if( IEta[ch]>0 && Depth[ch]==1 ) plusEtaPhiEnergyMapD1B_Histo->Fill(IEta[ch],IPhi[ch],Energy[ch]);
       if( IEta[ch]>0 && Depth[ch]==1 ) plusEtaPhiEnergyMapD1C_Histo->Fill(IEta[ch],IPhi[ch],Energy[ch]);
       if( IEta[ch]>0 && Depth[ch]==1 ) plusEtaPhiEnergyMapD1D_Histo->Fill(IEta[ch],IPhi[ch],Energy[ch]);
     }

     //This piece had a typo in early version - First HCAL Noiw WG presentation (which was noticed and presented with the correct label in plots)
     //int occupiedChBarrelCtr=0;
     //int occupiedChEndcapCtr=0;
     //for( int ie=1; ie<=16; ie++){
     //for( int ip=1; ip<=72; ip++){
     //for( int id=1; id<=3; id++){
     //if( plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>1.5 ) occupiedChBarrelCtr++;
     //if( minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>1.5 ) occupiedChEndcapCtr++;
     //}
     //}
     //}
     //
     int occupiedChBarrelCtr=0;
     int occupiedChEndcapCtr=0;
     for( int ie=1; ie<=16; ie++){//loop over ieta in Barrel, over all iphis and idepths
       for( int ip=1; ip<=72; ip++){
	 for( int id=1; id<=3; id++){
	   if( ie==16 && id==3 ) continue; //eta 16 has first 2 depths in HB
	   if( plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut  ) occupiedChBarrelCtr++;
	   if( minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut ) occupiedChBarrelCtr++;
	   ChannelEnergyBarrel_Histo->Fill(plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id));
	   ChannelEnergyBarrel_Histo->Fill(minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id));
	 }
       }
     }
     for( int ie=16; ie<=29; ie++){//loop over ieta in Endcap, over all iphis and idepths
       for( int ip=1; ip<=72; ip++){
	 for( int id=1; id<=3; id++){
	   if( ie==16 && id<3 ) continue; //eta 16 has first 2 depths in HB
	   if( plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut  ) occupiedChEndcapCtr++;
	   if( minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut ) occupiedChEndcapCtr++;
           ChannelEnergyEndcap_Histo->Fill(plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id));
	   ChannelEnergyEndcap_Histo->Fill(minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id));
	 }
       }
     }

     if( debug4_ ) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;
     if( debug4_ ) std::cout<<"total no of Barrel channels with E>"<<pixelEnergyCut<<" GeV: "<<occupiedChBarrelCtr<<std::endl;
     if( debug4_ ) std::cout<<"total no of Endcap channels with E>"<<pixelEnergyCut<<" GeV: "<<occupiedChEndcapCtr<<std::endl;

     //Occupancies (hit channels) in Barrel and Endcap - before mixing
     chBarrelOccupancyBefore_Histo->Fill(occupiedChBarrelCtr);
     chEndcapOccupancyBefore_Histo->Fill(occupiedChEndcapCtr);

     if( NumberOfGoodPrimaryVertices<0 ) std::cout<<"NumberOfGoodPrimaryVerticesCheck: "<<NumberOfGoodPrimaryVertices<<std::endl;

     pileup=0;
     pileup+=NumberOfGoodPrimaryVertices; // add pileup from jentry

     if( mixPlus1_ || mixPlus2_ || mixPlus3_ ){
       // load event "jentry+1"
       ientry = LoadTree(jentry+1);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry+1);   nbytes += nb;
       //-- add event only if official decision is TRUE and pileup>0
       if( OfficialDecision && NumberOfGoodPrimaryVertices>0 ){
	 //
	 for( int ch=0; ch<PulseCount; ch++){
	   //
	   if( debug3_ ){ if( ch==0 ){ std::cout<<"Ch=0 for event+1: "<<IEta[ch]<<" "<<IPhi[ch]<<" "<<Energy[ch]<<std::endl; } }
	   //
	   if( IEta[ch]>0 ) plusEtaPhiEnergyMap_Histo->Fill(IEta[ch],IPhi[ch],Depth[ch],Energy[ch]);
	   if( IEta[ch]<0 ) minusEtaPhiEnergyMap_Histo->Fill(abs(IEta[ch]),IPhi[ch],Depth[ch],Energy[ch]);
	   //
	   if( IEta[ch]>0 && Depth[ch]==1 ) plusEtaPhiEnergyMapD1B_Histo->Fill(IEta[ch],IPhi[ch],Energy[ch]);
	   if( IEta[ch]>0 && Depth[ch]==1 ) plusEtaPhiEnergyMapD1C_Histo->Fill(IEta[ch],IPhi[ch],Energy[ch]);
	   if( IEta[ch]>0 && Depth[ch]==1 ) plusEtaPhiEnergyMapD1D_Histo->Fill(IEta[ch],IPhi[ch],Energy[ch]);
	 }
	 
	 if( NumberOfGoodPrimaryVertices<0 ) std::cout<<"NumberOfGoodPrimaryVerticesCheck: "<<NumberOfGoodPrimaryVertices<<std::endl;
	 pileup+=NumberOfGoodPrimaryVertices; // add pileup from jentry+1
	 std::cout << "analysisClass::Loop(): jentry =  ADD:: " << jentry+1 <<"   PU:"<< NumberOfGoodPrimaryVertices<< std::endl;
       }
     }


     if( mixPlus2_ || mixPlus3_ ){
       // load event "jentry+2"
       ientry = LoadTree(jentry+2);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry+2);   nbytes += nb;
       //-- add event only if official decision is TRUE and pileup>0
       if( OfficialDecision && NumberOfGoodPrimaryVertices>0 ){
	 for( int ch=0; ch<PulseCount; ch++){
	   //
	   if( debug3_ ){ if( ch==0 ){ std::cout<<"Ch=0 for event+2: "<<IEta[ch]<<" "<<IPhi[ch]<<" "<<Energy[ch]<<std::endl; } }
	   //
	   if( IEta[ch]>0 ) plusEtaPhiEnergyMap_Histo->Fill(IEta[ch],IPhi[ch],Depth[ch],Energy[ch]);
	   if( IEta[ch]<0 ) minusEtaPhiEnergyMap_Histo->Fill(abs(IEta[ch]),IPhi[ch],Depth[ch],Energy[ch]);
	   //
	   if( IEta[ch]>0 && Depth[ch]==1 ) plusEtaPhiEnergyMapD1C_Histo->Fill(IEta[ch],IPhi[ch],Energy[ch]);
	   if( IEta[ch]>0 && Depth[ch]==1 ) plusEtaPhiEnergyMapD1D_Histo->Fill(IEta[ch],IPhi[ch],Energy[ch]);
	 }
	 
	 pileup+=NumberOfGoodPrimaryVertices; // add pileup from jentry+2
	 std::cout << "analysisClass::Loop(): jentry =  ADD:: " << jentry+2 <<"   PU:"<< NumberOfGoodPrimaryVertices<< std::endl;
       }
     }

     if( mixPlus3_ ){
       // load event "jentry+3"
       ientry = LoadTree(jentry+3);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry+3);   nbytes += nb;
       //-- add event only if official decision is TRUE and pileup>0
       if( OfficialDecision && NumberOfGoodPrimaryVertices>0 ){
	 for( int ch=0; ch<PulseCount; ch++){
	   //
	   if( debug3_ ){ if( ch==0 ){ std::cout<<"Ch=0 for event+3: "<<IEta[ch]<<" "<<IPhi[ch]<<" "<<Energy[ch]<<std::endl; std::cout<<std::endl; } }
	   //
	   if( IEta[ch]>0 ) plusEtaPhiEnergyMap_Histo->Fill(IEta[ch],IPhi[ch],Depth[ch],Energy[ch]);
	   if( IEta[ch]<0 ) minusEtaPhiEnergyMap_Histo->Fill(abs(IEta[ch]),IPhi[ch],Depth[ch],Energy[ch]);
	   //
	   if( IEta[ch]>0 && Depth[ch]==1 ) plusEtaPhiEnergyMapD1D_Histo->Fill(IEta[ch],IPhi[ch],Energy[ch]);
	 }
	 
	 pileup+=NumberOfGoodPrimaryVertices; // add pileup from jentry+3
	 std::cout << "analysisClass::Loop(): jentry =  ADD:: " << jentry+3 <<"   PU:"<< NumberOfGoodPrimaryVertices<< std::endl;
       }
     }

     occupiedChBarrelCtr=0;
     occupiedChEndcapCtr=0;
     //This piece had a typo in early version - First HCAL Noiw WG presentation (which was noticed and presented with the correct label in plots)
     //for( int ie=1; ie<=16; ie++){
     //for( int ip=1; ip<=72; ip++){
     //for( int id=1; id<=3; id++){
     //if( plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>1.5 ) occupiedChBarrelCtr++;
     //if( minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>1.5 ) occupiedChEndcapCtr++;
     //}
     //}
     //}
     //
     for( int ie=1; ie<=16; ie++){//loop over ieta in Barrel, over all iphis and idepths
       for( int ip=1; ip<=72; ip++){
	 for( int id=1; id<=3; id++){
	   if( plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut  ) occupiedChBarrelCtr++;
	   if( minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut ) occupiedChBarrelCtr++;
	 }
       }
     }
     for( int ie=16; ie<=29; ie++){//loop over ieta in Endcap, over all iphis and idepths
       for( int ip=1; ip<=72; ip++){
	 for( int id=1; id<=3; id++){
	   if( ie==16 && id<3 ) continue; //eta 16 has first 2 depths in HB
	   if( plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut  ) occupiedChEndcapCtr++;
	   if( minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut ) occupiedChEndcapCtr++;
	 }
       }
     }
     if( debug4_ ) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;
     if( debug4_ ) std::cout<<"total no of Barrel channels with E>"<<pixelEnergyCut<<" GeV: "<<occupiedChBarrelCtr<<std::endl;
     if( debug4_ ) std::cout<<"total no of Endcap channels with E>"<<pixelEnergyCut<<" GeV: "<<occupiedChEndcapCtr<<std::endl;
     if( debug4_ ) std::cout<<std::endl;

     chBarrelOccupancyAfter_Histo->Fill(occupiedChBarrelCtr);
     chEndcapOccupancyAfter_Histo->Fill(occupiedChEndcapCtr);


     // Each RBX has four RMs, one HPD per RM.
     //
     // For each event: 
     // -- We map pulses (E>1.5 GeV) in iEta/iPhi/iDepth to RBX/RM to a single number "RBXrm", for each HBP, HBM, HEP, HEM.
     // -- RBXrm histograms are made for each of HBP, HBM, HEP, HEM.
     // -- These histograms are checked for 18 hits (full pixel occupancy).
     // -- These histograms are reset for every event.
     //
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
           if( plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut ){
             HBPphiHitsPerEvent_Histo->Fill(HBP_EtaPhitoRBXrm(ie,ip,id),1); //HBP
           }
           //
           //--- minus side
           if( minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut ){
             HBMphiHitsPerEvent_Histo->Fill(HBM_EtaPhitoRBXrm(ie,ip,id),1); //HBM
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
           if( plusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut ){
             HEPphiHitsPerEvent_Histo->Fill(HEP_EtaPhitoRBXrm(ie,ip,id),1); //HEP
           }
           //
           //--- minus side
           if( minusEtaPhiEnergyMap_Histo->GetBinContent(ie,ip,id)>pixelEnergyCut ){
             HEMphiHitsPerEvent_Histo->Fill(HEM_EtaPhitoRBXrm(ie,ip,id),1); //HEM
           }
           //
         }
       }
     }
     
     // dead channels
     // pending..

     // also check for 4-HPD RBX noie
     bool eventHasRBXnoise_=false;
     int  HBP_badHPDperRBXctr=0;
     int  HBM_badHPDperRBXctr=0;
     int  HEP_badHPDperRBXctr=0;
     int  HEM_badHPDperRBXctr=0;
     unsigned int  RBXwithNoise=0;
     //
     //if( badHPDctr>=4 ){//Do this check if event has at least four bad HPDs. //cant call this now, defined/used below..
     //
     for( unsigned int irbx=1; irbx<=18; irbx++){
       RBXwithNoise=irbx;
       HBP_badHPDperRBXctr=0; HBM_badHPDperRBXctr=0;
       HEP_badHPDperRBXctr=0; HEM_badHPDperRBXctr=0;
       for( unsigned int irm=1; irm<=4; irm++ ){//loop over RMs of the given RBX
	 if( HBPphiHitsPerEvent_Histo->GetBinContent( int( (irbx-1)*4+irm ) )>filterCut ) HBP_badHPDperRBXctr++;
	 if( HBMphiHitsPerEvent_Histo->GetBinContent( int( (irbx-1)*4+irm ) )>filterCut ) HBM_badHPDperRBXctr++;
	 if( HEPphiHitsPerEvent_Histo->GetBinContent( int( (irbx-1)*4+irm ) )>filterCut ) HEP_badHPDperRBXctr++;
	 if( HEMphiHitsPerEvent_Histo->GetBinContent( int( (irbx-1)*4+irm ) )>filterCut ) HEM_badHPDperRBXctr++;
       }
       if( HBP_badHPDperRBXctr==4 || HBM_badHPDperRBXctr==4 || HEP_badHPDperRBXctr==4 || HEM_badHPDperRBXctr==4 ){ eventHasRBXnoise_=true; break; }
     }
     //
     //}

     if(  eventHasRBXnoise_ ) std::cout<<"Event has RBX noise! RBX#"<<RBXwithNoise<<" :: jentry "<<jentry<<std::endl;
     if(  eventHasRBXnoise_ ) RBXnoiseCounter_Histo->Fill(2);// event has RBX noise
     if( !eventHasRBXnoise_ ) RBXnoiseCounter_Histo->Fill(1);// event has NO RBX noise


     // Overflow, underflow checks in iEta/iPhi to RBXrm assignments
     if( HBPphiHitsPerEvent_Histo->GetBinContent(0)>0 ) std::cout<<"HBP has undeflow: "<<HBPphiHitsPerEvent_Histo->GetBinContent(0)<<std::endl;
     if( HBMphiHitsPerEvent_Histo->GetBinContent(0)>0 ) std::cout<<"HBM has undeflow: "<<HBMphiHitsPerEvent_Histo->GetBinContent(0)<<std::endl;
     if( HEPphiHitsPerEvent_Histo->GetBinContent(0)>0 ) std::cout<<"HEP has undeflow: "<<HEPphiHitsPerEvent_Histo->GetBinContent(0)<<std::endl;
     if( HEMphiHitsPerEvent_Histo->GetBinContent(0)>0 ) std::cout<<"HEM has undeflow: "<<HEMphiHitsPerEvent_Histo->GetBinContent(0)<<std::endl;
     //
     if( HBPphiHitsPerEvent_Histo->GetBinContent(73)>0 ) std::cout<<"HBP has overflow: "<<HBPphiHitsPerEvent_Histo->GetBinContent(73)<<std::endl;
     if( HBMphiHitsPerEvent_Histo->GetBinContent(73)>0 ) std::cout<<"HBM has overflow: "<<HBMphiHitsPerEvent_Histo->GetBinContent(73)<<std::endl;
     if( HEPphiHitsPerEvent_Histo->GetBinContent(73)>0 ) std::cout<<"HEP has overflow: "<<HEPphiHitsPerEvent_Histo->GetBinContent(73)<<std::endl;
     if( HEMphiHitsPerEvent_Histo->GetBinContent(73)>0 ) std::cout<<"HEM has overflow: "<<HEMphiHitsPerEvent_Histo->GetBinContent(73)<<std::endl;


     // Event filter BEGIN
     //if( pileup<1 ) continue;
     //if( eventHasRBXnoise_ ) continue;
     // Event filter END


     int badHPDctr=0;
     int barrelbadHPDctr=0;
     int endcapbadHPDctr=0;
     //
     for( unsigned int ihpd=1; ihpd<73; ihpd++){
       EndcapHPDoccupancy_Histo->Fill(HEPphiHitsPerEvent_Histo->GetBinContent(ihpd));
       EndcapHPDoccupancy_Histo->Fill(HEMphiHitsPerEvent_Histo->GetBinContent(ihpd));
       HEPHPDoccupancy_Histo->Fill(HEPphiHitsPerEvent_Histo->GetBinContent(ihpd));
       HEMHPDoccupancy_Histo->Fill(HEMphiHitsPerEvent_Histo->GetBinContent(ihpd));
       BarrelHPDoccupancy_Histo->Fill(HBMphiHitsPerEvent_Histo->GetBinContent(ihpd));
       BarrelHPDoccupancy_Histo->Fill(HBPphiHitsPerEvent_Histo->GetBinContent(ihpd));
       HBMHPDoccupancy_Histo->Fill(HBMphiHitsPerEvent_Histo->GetBinContent(ihpd));
       HBPHPDoccupancy_Histo->Fill(HBPphiHitsPerEvent_Histo->GetBinContent(ihpd));

       //
       if( HEPphiHitsPerEvent_Histo->GetBinContent(ihpd)>filterCut ){
	 HEPphiHits_Histo->Fill(ihpd,1);
	 if( debug_ ) std::cout<<"HEP:  HPD with "<<HEPphiHitsPerEvent_Histo->GetBinContent(ihpd)<<" hits :: jentry "<<jentry<<std::endl;
	 badHPDctr++;
	 endcapbadHPDctr++;
       }
       if( HEMphiHitsPerEvent_Histo->GetBinContent(ihpd)>filterCut ){
	 HEMphiHits_Histo->Fill(ihpd,1);
	 if( debug_ ) std::cout<<"HEM:  HPD with "<<HEMphiHitsPerEvent_Histo->GetBinContent(ihpd)<<" hits :: jentry "<<jentry<<std::endl;
	 badHPDctr++;
	 endcapbadHPDctr++;
       }
       if( HBPphiHitsPerEvent_Histo->GetBinContent(ihpd)>filterCut ){
	 HBPphiHits_Histo->Fill(ihpd,1);
	 if( debug_ ) std::cout<<"HBP:  HPD with "<<HBPphiHitsPerEvent_Histo->GetBinContent(ihpd)<<" hits :: jentry "<<jentry<<std::endl;
         badHPDctr++;
         barrelbadHPDctr++;
       }
       if( HBMphiHitsPerEvent_Histo->GetBinContent(ihpd)>filterCut ){
	 HBMphiHits_Histo->Fill(ihpd,1);
	 if( debug_ ) std::cout<<"HBM:  HPD with "<<HBMphiHitsPerEvent_Histo->GetBinContent(ihpd)<<" hits :: jentry "<<jentry<<std::endl;
	 badHPDctr++;
	 barrelbadHPDctr++;
       }
     }

     //debugging
     //for( unsigned int ihpd=1; ihpd<73; ihpd++){
     //  std::cout<<"HBM ihp("<<ihpd<<"): "<<HBMphiHitsPerEvent_Histo->GetBinContent(ihpd)<<std::endl;
     // }

     if( debug_ && HPDHits_jentry>17 && ( HEPphiHitsPerEvent_Histo->GetMaximum()!=18 && HEMphiHitsPerEvent_Histo->GetMaximum()!=18 &&
			 HBPphiHitsPerEvent_Histo->GetMaximum()!=18 && HBMphiHitsPerEvent_Histo->GetMaximum()!=18 ) ){
       // This is a discrepancy between already stored "HPDHits_jentry" variable, and our check here..
       std::cout<<"Max HPD pixel occupancy in HEP: "<<HEPphiHitsPerEvent_Histo->GetMaximum()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HEM: "<<HEMphiHitsPerEvent_Histo->GetMaximum()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HBP: "<<HBPphiHitsPerEvent_Histo->GetMaximum()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HBM: "<<HBMphiHitsPerEvent_Histo->GetMaximum()<<std::endl;
       //
       std::cout<<"Max HPD pixel occupancy in HEP, RBXrm no: "<<HEPphiHitsPerEvent_Histo->GetMaximumBin()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HEM, RBXrm no: "<<HEMphiHitsPerEvent_Histo->GetMaximumBin()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HBP, RBXrm no: "<<HBPphiHitsPerEvent_Histo->GetMaximumBin()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HBM, RBXrm no: "<<HBMphiHitsPerEvent_Histo->GetMaximumBin()<<std::endl;
     }

     //
     //Calculated Maximum:
     int calcMaxHPDHits=0;
     if( calcMaxHPDHits< HEPphiHitsPerEvent_Histo->GetMaximum() ) calcMaxHPDHits=HEPphiHitsPerEvent_Histo->GetMaximum();
     if( calcMaxHPDHits< HEMphiHitsPerEvent_Histo->GetMaximum() ) calcMaxHPDHits=HEMphiHitsPerEvent_Histo->GetMaximum();
     if( calcMaxHPDHits< HBPphiHitsPerEvent_Histo->GetMaximum() ) calcMaxHPDHits=HBPphiHitsPerEvent_Histo->GetMaximum();
     if( calcMaxHPDHits< HBMphiHitsPerEvent_Histo->GetMaximum() ) calcMaxHPDHits=HBMphiHitsPerEvent_Histo->GetMaximum();
     //
     if( mixPlus1_==false && mixPlus2_==false && mixPlus3_==false ){ //if no mixing..
       HPDHitsOfficial_Histo->Fill( HPDHits_jentry );
     }
     HPDHitsCalculated_Histo->Fill( calcMaxHPDHits );
     //
     if( HPDHits_jentry != calcMaxHPDHits ){
       std::cout<<std::endl;
       std::cout<<"jentry: "<< jentry<<std::endl;
       std::cout<<"NominalMET: "<<NominalMET[0]<<std::endl;
       std::cout<<"HBHE SumET: "<<HBSumET+HESumET<<std::endl;
       std::cout<<"HPDNoOtherHits: "<<HPDNoOtherHits<<std::endl;
       std::cout<<"Max HPD pixel occupancy Official: "<<HPDHits_jentry<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HEP: "<<HEPphiHitsPerEvent_Histo->GetMaximum()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HEM: "<<HEMphiHitsPerEvent_Histo->GetMaximum()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HBP: "<<HBPphiHitsPerEvent_Histo->GetMaximum()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HBM: "<<HBMphiHitsPerEvent_Histo->GetMaximum()<<std::endl;
       for( int iphi=1; iphi<=72; iphi++){
	 if( HBPphiHitsPerEvent_Histo->GetBinContent(iphi)> 5 ) std::cout<<"HBP iphi: "<<iphi<<" -> "<<HBPphiHitsPerEvent_Histo->GetBinContent(iphi)<<std::endl;
	 if( HBMphiHitsPerEvent_Histo->GetBinContent(iphi)> 5 ) std::cout<<"HBM iphi: "<<iphi<<" -> "<<HBMphiHitsPerEvent_Histo->GetBinContent(iphi)<<std::endl;
	 if( HEPphiHitsPerEvent_Histo->GetBinContent(iphi)> 5 ) std::cout<<"HEP iphi: "<<iphi<<" -> "<<HEPphiHitsPerEvent_Histo->GetBinContent(iphi)<<std::endl;
	 if( HEMphiHitsPerEvent_Histo->GetBinContent(iphi)> 5 ) std::cout<<"HEM iphi: "<<iphi<<" -> "<<HEMphiHitsPerEvent_Histo->GetBinContent(iphi)<<std::endl;
       }
       //
       std::cout<<"Max HPD pixel occupancy in HEP, RBXrm no: "<<HEPphiHitsPerEvent_Histo->GetMaximumBin()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HEM, RBXrm no: "<<HEMphiHitsPerEvent_Histo->GetMaximumBin()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HBP, RBXrm no: "<<HBPphiHitsPerEvent_Histo->GetMaximumBin()<<std::endl;
       std::cout<<"Max HPD pixel occupancy in HBM, RBXrm no: "<<HBMphiHitsPerEvent_Histo->GetMaximumBin()<<std::endl;
     }

     //
     barrelbadHPDperPileup_Histo->Fill(pileup,barrelbadHPDctr);
     barrelallHPDperPileup_Histo->Fill(pileup,144);//36*4  36 RBXs, 4 RMs

     endcapbadHPDperPileup_Histo->Fill(pileup,endcapbadHPDctr);
     endcapallHPDperPileup_Histo->Fill(pileup,144);//36*4  36 RBXs, 4 RMs

     badHPDperPileup_Histo->Fill(pileup,badHPDctr);
     allHPDperPileup_Histo->Fill(pileup,288);//72*4  72 RBXs, 4 RMs

     if( pileup<0 ) std::cout<<"pileupCheck: "<<pileup<<std::endl;
     Pileup_Histo->Fill(pileup,1);

     if( HPDHits_jentry> 17 && badHPDctr> 0 ) verifyHPDhits_Histo->Fill(1);//ok
     if( HPDHits_jentry<=17 && badHPDctr> 0 ) verifyHPDhits_Histo->Fill(2);//bad -- overestimate
     if( HPDHits_jentry> 17 && badHPDctr==0 ) verifyHPDhits_Histo->Fill(3);//bad -- underestimate
     if( HPDHits_jentry<=17 && badHPDctr==0 ) verifyHPDhits_Histo->Fill(4);//ok

     ////////////////////// User's code ends here ///////////////////////
     //break;
   } // End loop over events

   barrelNormbadHPDperPileup_Histo->Divide(barrelbadHPDperPileup_Histo, barrelallHPDperPileup_Histo,1,1,"cl=0.683 b(1,1) mode");
   endcapNormbadHPDperPileup_Histo->Divide(endcapbadHPDperPileup_Histo, endcapallHPDperPileup_Histo,1,1,"cl=0.683 b(1,1) mode");
   NormbadHPDperPileup_Histo->Divide(badHPDperPileup_Histo, allHPDperPileup_Histo,1,1,"cl=0.683 b(1,1) mode");

   //////////write histos 
   output_root_->cd();

   OfficialDecision_Histo->Write();
   PulseCount_Histo->Write();
   ChannelEnergy_Histo->Write();
   ChannelEnergyBarrel_Histo->Write();
   ChannelEnergyEndcap_Histo->Write();

   HPDHitsCalculated_Histo->Write();
   HPDHitsOfficial_Histo->Write();

   ChargeTS4TS5_Histo->Write();
   ChargeTS0TS9_Histo->Write();
   Energy_Histo->Write();

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

   RBXnoiseCounter_Histo->Write();

   chBarrelOccupancyBefore_Histo->Write();
   chEndcapOccupancyBefore_Histo->Write();
   chBarrelOccupancyAfter_Histo->Write();
   chEndcapOccupancyAfter_Histo->Write();

   plusEtaPhiEnergyMapD1A_Histo->Write();
   plusEtaPhiEnergyMapD1B_Histo->Write();
   plusEtaPhiEnergyMapD1C_Histo->Write();
   plusEtaPhiEnergyMapD1D_Histo->Write();


   CloseHCALmaps( );

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
