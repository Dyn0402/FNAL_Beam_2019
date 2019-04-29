/*
 * ce1_cuts.cpp
 *
 *  Created on: Apr 11, 2019
 *      Author: dylan
 */


#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TLine.h>

#include "../FNAL_Beam_Test_Source/config.h"

using namespace std;


vector<int> get_electron_events(TTree *tree, vector<int> beam_events, vector<vector<float>>ped_fits) {
	//Set branches from tree.
	float ce1;
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_ce1 = tree->GetBranch("Ce1RawADC");
	b_ce1->SetAddress(&ce1);

	//Set and fill sc1 histogram.
	TH1D *ce1_hist = new TH1D(config::ce1_electron_hist_name.data(), config::ce1_electron_hist_name.data(), config::bins, config::low_bin, config::high_bin);

	int event_index = 0;
	int beam_event_index = 0;
	int next_beam_event = beam_events[beam_event_index];
	while(b_event_num->GetEntry(event_index)) {
		if(event_index == next_beam_event) {
			next_beam_event = beam_events[++beam_event_index];
			b_ce1->GetEntry(event_index);
			ce1_hist->Fill(ce1);
		}
		event_index++;
	}

	ce1_hist->Write();

	//Determine background cut based on pedestal.
	double cut = ped_fits[0][config::ce1_channel] + ped_fits[1][config::ce1_channel]*config::ce1_electron_sigmas;

	//Set and save canvas showing fits and cuts.
	TCanvas *ce1_canvas = new TCanvas(config::ce1_electron_hist_name.data(), config::ce1_electron_hist_name.data(), config::canvas_x, config::canvas_y);

	TLine *cut_line = new TLine(cut, 0, cut, ce1_hist->GetMaximum());
	cut_line->SetLineColor(kRed);
	ce1_hist->Draw();
	cut_line->Draw();
	ce1_canvas->Write();

	//Determine electron events
	vector<int> electron_events;
	event_index = 0;
	beam_event_index = 0;
	next_beam_event = beam_events[beam_event_index];
	while(b_event_num->GetEntry(event_index)) {
		if(event_index == next_beam_event) {
			next_beam_event = beam_events[++beam_event_index];
			b_ce1->GetEntry(event_index);
			if(ce1 > cut) {
				electron_events.push_back(event_index);
			}
		}
		event_index++;
	}


	//Delete objects from heap
	delete ce1_canvas;
	delete cut_line;
	delete ce1_hist;


	return(electron_events);
}



vector<int> get_hadron_events(TTree *tree, vector<int> beam_events, vector<vector<float>>ped_fits) {
	//Set branches from tree.
	float ce1;
	TBranch *b_event_num = tree->GetBranch("eventno");
	TBranch *b_ce1 = tree->GetBranch("Ce1RawADC");
	b_ce1->SetAddress(&ce1);

	//Set and fill sc1 histogram.
	TH1D *ce1_hist = new TH1D(config::ce1_hadron_hist_name.data(), config::ce1_hadron_hist_name.data(), config::bins, config::low_bin, config::high_bin);

	int event_index = 0;
	int beam_event_index = 0;
	int next_beam_event = beam_events[beam_event_index];
	while(b_event_num->GetEntry(event_index)) {
		if(event_index == next_beam_event) {
			next_beam_event = beam_events[++beam_event_index];
			b_ce1->GetEntry(event_index);
			ce1_hist->Fill(ce1);
		}
		event_index++;
	}

	ce1_hist->Write();

	//Determine background cut based on pedestal.
	double cut = ped_fits[0][config::ce1_channel] + ped_fits[1][config::ce1_channel]*config::ce1_hadron_sigmas;

	//Set and save canvas showing fits and cuts.
	TCanvas *ce1_canvas = new TCanvas(config::ce1_hadron_hist_name.data(), config::ce1_hadron_hist_name.data(), config::canvas_x, config::canvas_y);

	TLine *cut_line = new TLine(cut, 0, cut, ce1_hist->GetMaximum());
	cut_line->SetLineColor(kRed);
	ce1_hist->Draw();
	cut_line->Draw();
	ce1_canvas->Write();

	//Determine electron events
	vector<int> hadron_events;
	event_index = 0;
	beam_event_index = 0;
	next_beam_event = beam_events[beam_event_index];
	while(b_event_num->GetEntry(event_index)) {
		if(event_index == next_beam_event) {
			next_beam_event = beam_events[++beam_event_index];
			b_ce1->GetEntry(event_index);
			if(ce1 < cut) {
				hadron_events.push_back(event_index);
			}
		}
		event_index++;
	}


	//Delete objects from heap
	delete ce1_canvas;
	delete cut_line;
	delete ce1_hist;


	return(hadron_events);
}
