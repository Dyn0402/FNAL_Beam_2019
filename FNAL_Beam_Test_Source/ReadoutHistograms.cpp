#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TObject.h"
#include "TKey.h"
#include "TTree.h"
#include "TObject.h"
#include "TIterator.h"
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <string.h>


using namespace std;

const int HCalNChannel = 16;
const int ECalNChannel = 16;
const int HodNChannel = 16;


typedef struct {
	
	int eventno;
	float HCalRawADC[HCalNChannel];
	float ECalRawADC[ECalNChannel];
	float Sc1RawADC;
	float Sc2RawADC;
	float MonRawADC;
	float Ce1RawADC;
	float Ce2RawADC;
	float PbGRawADC;
	float Veto1RawADC;
	float Veto2RawADC;
	float HodRawADC[HodNChannel];
    float junk[12];
	
	bool hhigh, ok, CalOk;
	unsigned char xmult, ymult, iSx, iSy, iEx, iEy, iEl;
	float eCalSum, hCalSum, eCalX, eCalY, eCalXloc, eCalYloc;
	
} Event;

void ReadoutHistograms(string fname, string outfilename) {
    
    //MAP HCAL= first 16 channels
    // ECAL = second 16 channels
    //Hod=third 16 channells
    //Sc1=48, Sc2=49, Ce1=50, Ce2=51
	
	float SingleEventRawADC[64];
	Event evn;
	
	int eventno;
	float HCalRawADC[HCalNChannel];
	for(int i=0; i<HCalNChannel; i++) HCalRawADC[i] = 0;
	float ECalRawADC[ECalNChannel];
	for(int i=0; i<ECalNChannel; i++) ECalRawADC[i] = 0;
    float junk[12]={0.};
	float Sc1RawADC = 0;
	float Sc2RawADC = 0;
	float Ce1RawADC = 0;
	float Ce2RawADC = 0;
	float HodRawADC[HodNChannel];
	for(int i=0; i<HodNChannel; i++) HodRawADC[i] = 0;
	
	//Channel lists modify here for configuration changes
	int HCalChannelList[HCalNChannel];
	for(int i=0; i<HCalNChannel; i++) HCalChannelList[i] = i;
	
	int ECalChannelList[ECalNChannel];
	for(int i=0; i<ECalNChannel; i++) ECalChannelList[i] = i+16;
	
	int HodChannelList[HodNChannel];
	for(int i=0; i<HodNChannel; i++) HodChannelList[i] = i+32;//i+48
	
	int Sc1ChannelList = 48;//32
	int Sc2ChannelList = 49;//40
	int Ce1ChannelList = 50; //Inner 35
	int Ce2ChannelList = 51; //Outer 37

	int index = 0;
	long int number;
	int thisevent;

	TFile *f = new TFile(outfilename.data(), "RECREATE");
	TTree *T = new TTree("T", "event data");
	
	//branches
	T->Branch("eventno", &evn.eventno, "eventno/I");
	T->Branch("HCalRawADC", evn.HCalRawADC, "HCalRawADC[16]/F");
	T->Branch("ECalRawADC", evn.ECalRawADC, "ECalRawADC[16]/F");
    T->Branch("HodRawADC", evn.HodRawADC, "HodRawADC[16]/F");
	T->Branch("Sc1RawADC", &evn.Sc1RawADC, "Sc1RawADC/F");
	T->Branch("Sc2RawADC", &evn.Sc2RawADC, "Sc2RawADC/F");

	T->Branch("Ce1RawADC", &evn.Ce1RawADC, "Ce1RawADC/F");
	T->Branch("Ce2RawADC", &evn.Ce2RawADC, "Ce2RawADC/F");
    T->Branch("junk",&evn.junk,"junk[12]/F");


	bool endofevent = false;
	ifstream indata;
	indata.open(fname);
	indata.clear();
	
	while(!indata.eof()) {
		indata >> number;

		
		if(index%65==0) {
			thisevent = number;
			if((index/65)%50000==0) cout << "Processing event: " << index/65 << endl;
		}
		else {
			SingleEventRawADC[index%65-1] = number;
			if(index%65==64) endofevent = true;
		}
		
		//end of single event in loop add to appropriate arrays
		if(endofevent) {
			eventno = thisevent;

			//for(int i=0; i<48; i++) cout << SingleEventRawADC[i] << endl;
						
			for(int i=0; i<HCalNChannel; i++) {
				HCalRawADC[i] = SingleEventRawADC[HCalChannelList[i]];
			}
			
			for(int i=0; i<ECalNChannel; i++) {
				ECalRawADC[i] = SingleEventRawADC[ECalChannelList[i]];
			}
			
			Sc1RawADC = SingleEventRawADC[Sc1ChannelList];
			Sc2RawADC = SingleEventRawADC[Sc2ChannelList];

			Ce1RawADC = SingleEventRawADC[Ce1ChannelList];
			Ce2RawADC = SingleEventRawADC[Ce2ChannelList];

			
			for(int i=0; i<HodNChannel; i++) {
				HodRawADC[i] = SingleEventRawADC[HodChannelList[i]];
			}
            
            for(int i=0;i<12;i++){
                junk[i]=SingleEventRawADC[i+52];
               // cout<<junk[i]<<endl;
            }
			
			//Write to event struct
			evn.eventno = eventno;
			for(int i=0; i<HCalNChannel; i++) {
				evn.HCalRawADC[i] = HCalRawADC[i];
			}
			for(int i=0; i<ECalNChannel; i++) {
				evn.ECalRawADC[i] = ECalRawADC[i];
			}
			evn.Sc1RawADC = Sc1RawADC;
			evn.Sc2RawADC = Sc2RawADC;

			evn.Ce1RawADC = Ce1RawADC;
			evn.Ce2RawADC = Ce2RawADC;

			for(int i=0; i<HodNChannel; i++) {
				evn.HodRawADC[i] = HodRawADC[i];
			}
            for(int i=0;i<12;i++){
                evn.junk[i]=junk[i];
               // cout<<evn.junk[i]<<endl;
                
            }
			
			T->Fill();
			
			for(int i=0; i<64; i++) SingleEventRawADC[i] = 0;
			endofevent = false;
		}
		
		index++;

	}
	indata.close();	
	
	f->Write();
	f->Close();
}
