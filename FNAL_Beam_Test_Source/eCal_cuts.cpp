/*
 * eCal_cuts.cpp
 *
 *  Created on: Apr 7, 2019
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


vector<int> get_led_events(TTree *tree, vector<int> events) {
	//Set branches from tree.
	float eCal_ADC[config::eCal_channels];
	float sc1;
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_eCal_ADC = tree->GetBranch("ECalRawADC");
	TBranch *b_sc1 = tree->GetBranch("Sc1RawADC");
	b_eCal_ADC->SetAddress(eCal_ADC);
	b_sc1->SetAddress(&sc1);

	//Set eCal histograms.
	TH1D *eCal_hists[config::eCal_channels];
	TH1D *eCal_upper_hists[config::eCal_channels];
	string name, name_upper;
	for(int i=0; i<config::eCal_channels; i++) {
		name = config::eCal_hist_name + to_string(i);
		name_upper = name + " LED range";
		eCal_hists[i] = new TH1D(name.data(), name.data(), config::bins, config::low_bin, config::high_bin);
		eCal_upper_hists[i] = new TH1D(name_upper.data(), name_upper.data(), config::bins, config::low_bin, config::high_bin);
	}

	//Fill eCal histograms.
	int event_index = 0;
	int next_event_index = 0;
	int next_event = events[next_event_index];
	while(b_event_num->GetEntry(event_index)) {
		if(event_index == next_event) {
			next_event = events[++next_event_index];
			b_eCal_ADC->GetEntry(event_index);
			for(int i=0; i<config::eCal_channels; i++) {
				eCal_hists[i]->Fill(eCal_ADC[i]);
				if(eCal_ADC[i] > config::eCal_led_low_ADC) {
					eCal_upper_hists[i]->Fill(eCal_ADC[i]);
				}
			}
		}
		event_index++;
	}

	//Fit peak and find led cuts.
	vector<vector<double>> led_cuts;
	vector<double> cut = {0.0, 0.0};
	TF1 *gauss = new TF1(config::eCal_fit_name.data(), "gaus");
	double mean, rms;
	Double_t *x_peak, *y_peak;
	for(int i:config::eCal_cut_channels) {
		TSpectrum *spec = new TSpectrum(config::eCal_spec_peaks);
		spec->Search(eCal_upper_hists[i], config::eCal_spec_sigma, "nobackground", config::eCal_spec_thresh);
		x_peak = spec->GetPositionX();
		y_peak = spec->GetPositionY();
		gauss->SetParameters(y_peak[0], x_peak[0], config::eCal_spec_sigma);
		gauss->SetRange(x_peak[0]-config::eCal_spec_sigma*config::eCal_fit_range, x_peak[0]+config::eCal_spec_sigma*config::eCal_fit_range);
		eCal_hists[i]->Fit(config::eCal_fit_name.data(), "RQ");
		mean = gauss->GetParameter(1);
		rms = gauss->GetParameter(2);
		cut[0] = mean - rms * config::eCal_rms_range;
		cut[1] = mean + rms * config::eCal_rms_range;
		led_cuts.push_back(cut);
		eCal_hists[i]->Write();
		delete spec;
	}


//	//Get number of noise coincidence for each event and determine noise events.
	vector<int> led_events;
	int coincidence = 0;
	int j = 0;
	event_index = 0;
	//Should I only check events or all events here?
	while(b_event_num->GetEntry(event_index)) {
		b_eCal_ADC->GetEntry(event_index);
		coincidence = 0;
		j = 0;
		for(int i:config::eCal_cut_channels) {
			if(eCal_ADC[i] > led_cuts[j][0] && eCal_ADC[i] < led_cuts[j][1]) {
				coincidence++;
			}
			j++;
		}
		if(coincidence >= config::eCal_coinc_thresh) {
			led_events.push_back(event_index);
		}
		event_index++;
	}

	//Delete objects from heap.
	for(int i=0; i<config::eCal_channels; i++) {
		delete eCal_hists[i];
		delete eCal_upper_hists[i];
	}
	delete gauss;


	return(led_events);
}





vector<int> get_led_events2(TTree *tree, vector<int> events) {
	//Set branches from tree.
	float eCal_ADC[config::eCal_channels];
	float sc1;
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_eCal_ADC = tree->GetBranch("ECalRawADC");
	TBranch *b_sc1 = tree->GetBranch("Sc1RawADC");
	b_eCal_ADC->SetAddress(eCal_ADC);
	b_sc1->SetAddress(&sc1);

	//Set eCal histograms.
	TH1D *eCal_hists[config::eCal_channels];
	TH1D *eCal_upper_hists[config::eCal_channels];
	string name, name_upper;
	for(int i=0; i<config::eCal_channels; i++) {
		name = config::eCal_hist_name + to_string(i);
		name_upper = name + " LED range";
		eCal_hists[i] = new TH1D(name.data(), name.data(), config::bins/4, config::low_bin, config::high_bin);
		eCal_upper_hists[i] = new TH1D(name_upper.data(), name_upper.data(), config::bins/4, config::low_bin, config::high_bin);
	}

	//Fill eCal histograms.
	int event_index = 0;
	int next_event_index = 0;
	int next_event = events[next_event_index];
	while(b_event_num->GetEntry(event_index)) {
		if(event_index == next_event) {
			next_event = events[++next_event_index];
			b_eCal_ADC->GetEntry(event_index);
			for(int i=0; i<config::eCal_channels; i++) {
				eCal_hists[i]->Fill(eCal_ADC[i]);
				if(eCal_ADC[i] > config::eCal_led_low_ADC) {
					eCal_upper_hists[i]->Fill(eCal_ADC[i]);
				}
			}
		}
		event_index++;
	}

	//Fit peak and find led cuts.
	vector<vector<double>> led_cuts;
	vector<double> cut = {0.0, 0.0};
	TF1 *gauss = new TF1(config::eCal_fit_name.data(), "gaus");
	TSpectrum *specs[config::eCal_channels];
	double mean, rms;
	int num_peaks;
	Double_t *x_peak, *y_peak;
	for(int i:config::eCal_cut_channels) {
		specs[i] = new TSpectrum(config::eCal_spec_peaks);
		num_peaks = specs[i]->Search(eCal_upper_hists[i], config::eCal_spec_sigma, "nobackground", config::eCal_spec_thresh);
		x_peak = specs[i]->GetPositionX();
		cout << i << "\t" << num_peaks << "\t" << x_peak[0] << endl;
		y_peak = specs[i]->GetPositionY();
		gauss->SetParameters(y_peak[0], x_peak[0], config::eCal_spec_sigma);
		gauss->SetRange(x_peak[0]-config::eCal_spec_sigma*config::eCal_fit_range, x_peak[0]+config::eCal_spec_sigma*config::eCal_fit_range);
		eCal_hists[i]->Fit(gauss, "RQ");
		mean = gauss->GetParameter(1);
		rms = gauss->GetParameter(2);
		cut[0] = mean - rms * config::eCal_rms_range;
		cut[1] = mean + rms * config::eCal_rms_range;
		led_cuts.push_back(cut);
		eCal_hists[i]->Write();
	}


//	//Get number of noise coincidence for each event and determine noise events.
	vector<int> led_events;
	int coincidence = 0;
	int j = 0;
	event_index = 0;
	//Should I only check events or all events here?
	while(b_event_num->GetEntry(event_index)) {
		b_eCal_ADC->GetEntry(event_index);
		coincidence = 0;
		j = 0;
		for(int i:config::eCal_cut_channels) {
			if(eCal_ADC[i] > led_cuts[j][0] && eCal_ADC[i] < led_cuts[j][1]) {
				coincidence++;
			}
			j++;
		}
		if(coincidence >= config::eCal_coinc_thresh) {
			led_events.push_back(event_index);
		}
		event_index++;
	}

	//Delete objects from heap.
	for(int i=0; i<config::eCal_channels; i++) {
		delete eCal_hists[i];
		delete eCal_upper_hists[i];
	}
	for(int i:config::eCal_cut_channels) {
		delete specs[i];
	}
	delete gauss;


	return(led_events);
}
