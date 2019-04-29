///*
// * write_64.cpp
// *
// *  Created on: Apr 15, 2019
// *      Author: dylan
// */
//
//
//
//
//#include <vector>
//
//#include "TFile.h"
//#include "TTree.h"
//#include "TCanvas.h"
//#include "TH1.h"
//#include "TF1.h"
//#include "TSpectrum.h"
//
//#include "config.h"
//
//using namespace std;
//
//
//void write_64(TTree *tree, vector<int> events, vector<int> events2) {
//	//Set branches from tree.
//	float hCal[config::hCal_channels];
//	float eCal[config::eCal_channels];
//	float hod[config::hod_channels];
//	float sc1, sc2, ce1, ce2;
//	float junk[config::junk_channels];
//	TBranch *b_event_num = tree->GetBranch("eventno");
//	TBranch *b_hCal = tree->GetBranch("HCalRawADC");
//	TBranch *b_eCal = tree->GetBranch("ECalRawADC");
//	TBranch *b_hod = tree->GetBranch("HodRawADC");
//	TBranch *b_sc1 = tree->GetBranch("Sc1RawADC");
//	TBranch *b_sc2 = tree->GetBranch("Sc2RawADC");
//	TBranch *b_ce1 = tree->GetBranch("Ce1RawADC");
//	TBranch *b_ce2 = tree->GetBranch("Ce2RawADC");
//	TBranch *b_junk = tree->GetBranch("junk");
//	b_hCal->SetAddress(hCal);
//	b_eCal->SetAddress(eCal);
//	b_hod->SetAddress(hod);
//	b_sc1->SetAddress(&sc1);
//	b_sc2->SetAddress(&sc2);
//	b_ce1->SetAddress(&ce1);
//	b_ce2->SetAddress(&ce2);
//	b_junk->SetAddress(junk);
//
//	//Set channel histograms.
//	TH1D *channel_hists[64];
//	TH1D *channel_hists2[64];
//	string name;
//	for(int i=0; i<64; i++) {
//		name = "Channel " + to_string(i);
//		channel_hists[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
//		name = "Channel " + to_string(i) + " Noise";
//		channel_hists2[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
//		channel_hists2[i]->SetLineColor(kRed);
//	}
//
//	for(auto event_index:events) {
//		b_hCal->GetEntry(event_index);
//		b_eCal->GetEntry(event_index);
//		b_hod->GetEntry(event_index);
//		b_sc1->GetEntry(event_index);
//		b_sc2->GetEntry(event_index);
//		b_ce1->GetEntry(event_index);
//		b_ce2->GetEntry(event_index);
//		b_junk->GetEntry(event_index);
//
//		//Fill channel histograms
//		for(int i = 0; i<config::hCal_channels; i++) { //hCal
//			channel_hists[i]->Fill(hCal[i]);
//		}
//		for(int i = 0; i<config::eCal_channels; i++) { //eCal
//			channel_hists[i+16]->Fill(eCal[i]);
//		}
//		for(int i = 0; i<config::hod_channels; i++) { //hod
//			channel_hists[i+32]->Fill(hod[i]);
//		}
//		channel_hists[config::sc1_channel]->Fill(sc1);
//		channel_hists[config::sc2_channel]->Fill(sc2);
//		channel_hists[config::ce1_channel]->Fill(ce1);
//		channel_hists[config::ce2_channel]->Fill(ce2);
//		for(int i = 0; i<config::junk_channels; i++) { //junk
//			channel_hists[i+config::junk_start_channel]->Fill(junk[i]);
//		}
//
//	}
//
//
//	for(auto event_index:events2) {
//		b_hCal->GetEntry(event_index);
//		b_eCal->GetEntry(event_index);
//		b_hod->GetEntry(event_index);
//		b_sc1->GetEntry(event_index);
//		b_sc2->GetEntry(event_index);
//		b_ce1->GetEntry(event_index);
//		b_ce2->GetEntry(event_index);
//		b_junk->GetEntry(event_index);
//
//		//Fill channel histograms
//		for(int i = 0; i<config::hCal_channels; i++) { //hCal
//			channel_hists2[i]->Fill(hCal[i]);
//		}
//		for(int i = 0; i<config::eCal_channels; i++) { //eCal
//			channel_hists2[i+16]->Fill(eCal[i]);
//		}
//		for(int i = 0; i<config::hod_channels; i++) { //hod
//			channel_hists2[i+32]->Fill(hod[i]);
//		}
//		channel_hists2[config::sc1_channel]->Fill(sc1);
//		channel_hists2[config::sc2_channel]->Fill(sc2);
//		channel_hists2[config::ce1_channel]->Fill(ce1);
//		channel_hists2[config::ce2_channel]->Fill(ce2);
//		for(int i = 0; i<config::junk_channels; i++) { //junk
//			channel_hists2[i+config::junk_start_channel]->Fill(junk[i]);
//		}
//
//	}
//
//	TCanvas *ADC1_canvas = new TCanvas("ADC1" , "ADC1", config::canvas_x, config::canvas_y);
//	ADC1_canvas->cd();
//	ADC1_canvas->Divide(4,4);
//
//	for(int i=0; i<config::hCal_channels; i++) {
//		ADC1_canvas->cd(i+1);
//		channel_hists[i+config::hCal_start_channel]->Draw();
//		channel_hists2[i+config::hCal_start_channel]->Draw();
//	}
//
//
//
//}
//
//
//
//void write_64_2(TTree *tree, vector<vector<int>> events, vector<int> events2) {
//	//Set branches from tree.
//	float hCal[config::hCal_channels];
//	float eCal[config::eCal_channels];
//	float hod[config::hod_channels];
//	float sc1, sc2, ce1, ce2;
//	float junk[config::junk_channels];
//	TBranch *b_event_num = tree->GetBranch("eventno");
//	TBranch *b_hCal = tree->GetBranch("HCalRawADC");
//	TBranch *b_eCal = tree->GetBranch("ECalRawADC");
//	TBranch *b_hod = tree->GetBranch("HodRawADC");
//	TBranch *b_sc1 = tree->GetBranch("Sc1RawADC");
//	TBranch *b_sc2 = tree->GetBranch("Sc2RawADC");
//	TBranch *b_ce1 = tree->GetBranch("Ce1RawADC");
//	TBranch *b_ce2 = tree->GetBranch("Ce2RawADC");
//	TBranch *b_junk = tree->GetBranch("junk");
//	b_hCal->SetAddress(hCal);
//	b_eCal->SetAddress(eCal);
//	b_hod->SetAddress(hod);
//	b_sc1->SetAddress(&sc1);
//	b_sc2->SetAddress(&sc2);
//	b_ce1->SetAddress(&ce1);
//	b_ce2->SetAddress(&ce2);
//	b_junk->SetAddress(junk);
//
//	//Set channel histograms.
//	TH1D *channel_hists[64];
//	TH1D *channel_hists2[64];
//	string name;
//	for(int i=0; i<64; i++) {
//		name = "Channel " + to_string(i);
//		channel_hists[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
//		name = "Channel " + to_string(i) + " Noise";
//		channel_hists2[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
//		channel_hists2[i]->SetLineColor(kRed);
//	}
//
//	for(auto event_index:events) {
//		b_hCal->GetEntry(event_index);
//		b_eCal->GetEntry(event_index);
//		b_hod->GetEntry(event_index);
//		b_sc1->GetEntry(event_index);
//		b_sc2->GetEntry(event_index);
//		b_ce1->GetEntry(event_index);
//		b_ce2->GetEntry(event_index);
//		b_junk->GetEntry(event_index);
//
//		//Fill channel histograms
//		for(int i = 0; i<config::hCal_channels; i++) { //hCal
//			channel_hists[i]->Fill(hCal[i]);
//		}
//		for(int i = 0; i<config::eCal_channels; i++) { //eCal
//			channel_hists[i+16]->Fill(eCal[i]);
//		}
//		for(int i = 0; i<config::hod_channels; i++) { //hod
//			channel_hists[i+32]->Fill(hod[i]);
//		}
//		channel_hists[config::sc1_channel]->Fill(sc1);
//		channel_hists[config::sc2_channel]->Fill(sc2);
//		channel_hists[config::ce1_channel]->Fill(ce1);
//		channel_hists[config::ce2_channel]->Fill(ce2);
//		for(int i = 0; i<config::junk_channels; i++) { //junk
//			channel_hists[i+config::junk_start_channel]->Fill(junk[i]);
//		}
//
//	}
//
//
//	for(auto event_index:events2) {
//		b_hCal->GetEntry(event_index);
//		b_eCal->GetEntry(event_index);
//		b_hod->GetEntry(event_index);
//		b_sc1->GetEntry(event_index);
//		b_sc2->GetEntry(event_index);
//		b_ce1->GetEntry(event_index);
//		b_ce2->GetEntry(event_index);
//		b_junk->GetEntry(event_index);
//
//		//Fill channel histograms
//		for(int i = 0; i<config::hCal_channels; i++) { //hCal
//			channel_hists2[i]->Fill(hCal[i]);
//		}
//		for(int i = 0; i<config::eCal_channels; i++) { //eCal
//			channel_hists2[i+16]->Fill(eCal[i]);
//		}
//		for(int i = 0; i<config::hod_channels; i++) { //hod
//			channel_hists2[i+32]->Fill(hod[i]);
//		}
//		channel_hists2[config::sc1_channel]->Fill(sc1);
//		channel_hists2[config::sc2_channel]->Fill(sc2);
//		channel_hists2[config::ce1_channel]->Fill(ce1);
//		channel_hists2[config::ce2_channel]->Fill(ce2);
//		for(int i = 0; i<config::junk_channels; i++) { //junk
//			channel_hists2[i+config::junk_start_channel]->Fill(junk[i]);
//		}
//
//	}
//
//	TCanvas *ADC1_canvas = new TCanvas("ADC1" , "ADC1", config::canvas_x, config::canvas_y);
//	ADC1_canvas->cd();
//	ADC1_canvas->Divide(4,4);
//
//	for(int i=0; i<config::hCal_channels; i++) {
//		ADC1_canvas->cd(i+1);
//		channel_hists[i+config::hCal_start_channel]->Draw();
//		channel_hists2[i+config::hCal_start_channel]->Draw();
//	}
//
//
//
//}
