/*
 * pedestals.cpp
 *
 *  Created on: Apr 8, 2019
 *      Author: dylan
 */

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"

#include "../FNAL_Beam_Test_Source/config.h"

using namespace std;


vector<vector<float>> get_pedestals(TTree *tree, vector<int> ped_events) {
	//Set branches from tree.
	float hCal[config::hCal_channels];
	float eCal[config::eCal_channels];
	float hod[config::hod_channels];
	float sc1, sc2, ce1, ce2;
	float junk[config::junk_channels];
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_hCal = tree->GetBranch("HCalRawADC");
	TBranch *b_eCal = tree->GetBranch("ECalRawADC");
	TBranch *b_hod = tree->GetBranch("HodRawADC");
	TBranch *b_sc1 = tree->GetBranch("Sc1RawADC");
	TBranch *b_sc2 = tree->GetBranch("Sc2RawADC");
	TBranch *b_ce1 = tree->GetBranch("Ce1RawADC");
	TBranch *b_ce2 = tree->GetBranch("Ce2RawADC");
	TBranch *b_junk = tree->GetBranch("junk");
	b_hCal->SetAddress(hCal);
	b_eCal->SetAddress(eCal);
	b_hod->SetAddress(hod);
	b_sc1->SetAddress(&sc1);
	b_sc2->SetAddress(&sc2);
	b_ce1->SetAddress(&ce1);
	b_ce2->SetAddress(&ce2);
	b_junk->SetAddress(junk);

	//Set channel histograms.
	TH1D *channel_hists[64];
	string name;
	for(int i=0; i<64; i++) {
		name = "Channel " + to_string(i);
		channel_hists[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
	}

	//Fill channel histograms.
	int event_index = 0;
	int ped_event_index = 0;
	int next_ped_event = ped_events[ped_event_index];
	while(b_event_num->GetEntry(event_index)) {
		//Check filters
		if(event_index == next_ped_event) {
			next_ped_event = ped_events[++ped_event_index];

			//Get branch entries
			b_hCal->GetEntry(event_index);
			b_eCal->GetEntry(event_index);
			b_hod->GetEntry(event_index);
			b_sc1->GetEntry(event_index);
			b_sc2->GetEntry(event_index);
			b_ce1->GetEntry(event_index);
			b_ce2->GetEntry(event_index);
			b_junk->GetEntry(event_index);

			//Fill channel histograms
			for(int i = 0; i<config::hCal_channels; i++) { //hCal
				channel_hists[i]->Fill(hCal[i]);
			}
			for(int i = 0; i<config::eCal_channels; i++) { //eCal
				channel_hists[i+16]->Fill(eCal[i]);
			}
			for(int i = 0; i<config::hod_channels; i++) { //hod
				channel_hists[i+32]->Fill(hod[i]);
			}
			channel_hists[config::sc1_channel]->Fill(sc1);
			channel_hists[config::sc2_channel]->Fill(sc2);
			channel_hists[config::ce1_channel]->Fill(ce1);
			channel_hists[config::ce2_channel]->Fill(ce2);
			for(int i = 0; i<config::junk_channels; i++) { //junk
				channel_hists[i+config::junk_start_channel]->Fill(junk[i]);
			}
		}
		event_index++;
	}

	//Fit channel pedestals
	vector<vector<float>> ped_fits = {{},{}};
	Double_t *x_peaks, *y_peaks;
	int num_peaks;
	double x_peak, y_peak;
	TF1 *gauss = new TF1("Ped_Gauss", "gaus");
	for(int i=0; i<config::hCal_channels; i++) { //hCal
		TSpectrum *spec = new TSpectrum(4); //Hard-coded, fix
		num_peaks = spec->Search(channel_hists[i+config::hCal_start_channel], 10.0, "nobackground", config::ped_spec_thresh);
		x_peaks = spec->GetPositionX();
		y_peaks = spec->GetPositionY();
		x_peak = x_peaks[0];
		y_peak = y_peaks[0];
		for(int i=1; i<num_peaks; i++) {
			if(x_peaks[i] < x_peak) {
				x_peak = x_peaks[i];
				y_peak = y_peaks[i];
			}
		}
		gauss->SetParameters(y_peak, x_peak, 10.0);
		gauss->SetRange(x_peak - 20.0, x_peak + 20.0);
		channel_hists[i+config::hCal_start_channel]->Fit("Ped_Gauss", "RQ");
		ped_fits[0].push_back(gauss->GetParameter(1));
		ped_fits[1].push_back(gauss->GetParameter(2));
		channel_hists[i+config::hCal_start_channel]->Write();
		delete spec;
	}
	for(int i=0; i<config::eCal_channels; i++) { //eCal
		TSpectrum *spec = new TSpectrum(4); //Hard-coded, fix
		num_peaks = spec->Search(channel_hists[i+config::eCal_start_channel], 3.0, "nobackground", config::ped_spec_thresh);
		x_peaks = spec->GetPositionX();
		y_peaks = spec->GetPositionY();
		x_peak = x_peaks[0];
		y_peak = y_peaks[0];
		for(int i=1; i<num_peaks; i++) {
			if(x_peaks[i] < x_peak) {
				x_peak = x_peaks[i];
				y_peak = y_peaks[i];
			}
		}
		gauss->SetParameters(y_peak, x_peak, 3.0);
		gauss->SetRange(x_peak - 6.0, x_peak + 6.0);
		channel_hists[i+config::eCal_start_channel]->Fit("Ped_Gauss", "RQ");
		ped_fits[0].push_back(gauss->GetParameter(1));
		ped_fits[1].push_back(gauss->GetParameter(2));
		channel_hists[i+config::eCal_start_channel]->Write();
		delete spec;
	}
	for(int i=0; i<config::hod_channels; i++) { //hod
		TSpectrum *spec = new TSpectrum(4); //Hard-coded, fix
		num_peaks = spec->Search(channel_hists[i+config::hod_start_channel], 10.0, "nobackground", config::ped_spec_thresh);
		x_peaks = spec->GetPositionX();
		y_peaks = spec->GetPositionY();
		x_peak = x_peaks[0];
		y_peak = y_peaks[0];
		for(int i=1; i<num_peaks; i++) {
			if(x_peaks[i] < x_peak) {
				x_peak = x_peaks[i];
				y_peak = y_peaks[i];
			}
		}
		gauss->SetParameters(y_peak, x_peak, 10.0);
		gauss->SetRange(x_peak - 20.0, x_peak + 20.0);
		channel_hists[i+config::hod_start_channel]->Fit("Ped_Gauss", "RQ");
		ped_fits[0].push_back(gauss->GetParameter(1));
		ped_fits[1].push_back(gauss->GetParameter(2));
		channel_hists[i+config::hod_start_channel]->Write();
		delete spec;
	}
	for(int i=0; i<16; i++) { //ADC 4
		TSpectrum *spec = new TSpectrum(4); //Hard-coded, fix
		num_peaks = spec->Search(channel_hists[i+config::sc1_channel], 5.0, "nobackground", config::ped_spec_thresh);
		x_peaks = spec->GetPositionX();
		y_peaks = spec->GetPositionY();
		x_peak = x_peaks[0];
		y_peak = y_peaks[0];
		for(int i=1; i<num_peaks; i++) {
			if(x_peaks[i] < x_peak) {
				x_peak = x_peaks[i];
				y_peak = y_peaks[i];
			}
		}
		gauss->SetParameters(y_peak, x_peak, 5.0);
		gauss->SetRange(x_peak - 10.0, x_peak + 10.0);
		channel_hists[i+config::sc1_channel]->Fit("Ped_Gauss", "RQ");
		ped_fits[0].push_back(gauss->GetParameter(1));
		ped_fits[1].push_back(gauss->GetParameter(2));
		channel_hists[i+config::sc1_channel]->Write();
		delete spec;
	}

	return(ped_fits);
}
