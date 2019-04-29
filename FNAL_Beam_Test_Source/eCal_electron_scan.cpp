/*
 * eCal_electron_scan.cpp
 *
 *  Created on: Apr 14, 2019
 *      Author: dylan
 */

#include <vector>

#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TLine.h>

#include "config.h"

using namespace std;


const float ecal_weight_masha[16]={.974462,.954643,.946785,1.01851,1.00689,.920467,.90186,1.0301,1.11047,1.03755,1.03422,1.10631,1.07239,.898038,1.0314,.955892};
const float hcal_weight_masha[16]={0.81323318,    0.685412787,    0.949554328,    0.722023522,    0.725346692,    1.434315419,    1.178238512,    1.309768637,    1.074029024,    1.070514826,    1.031548337,    1.203621835,    0.932082456,    0.836555741,    0.859508861,    1.174245842};

void get_eCal_sum(TTree *tree, vector<int> events, vector<vector<float>> ped_fits) {
	int event_num;
	float eCal[config::eCal_channels];
	TBranch *b_eCal = tree->GetBranch("ECalRawADC");
	b_eCal->SetAddress(eCal);

	//Set eCal Sum histogram.
	TH1D *eCal_sum = new TH1D("eCal_Sum", "Sum of eCal Electron Energy", config::bins/2, config::low_bin, config::high_bin);

	double ele_energy;
	for(auto event_index:events) {
		b_eCal->GetEntry(event_index);

		ele_energy = 0.0;
		for(int i=0; i<config::eCal_channels; i++) {
			if(eCal[i] - ped_fits[0][i+config::eCal_start_channel] > 3 * ped_fits[1][i+config::eCal_start_channel]) {
				ele_energy += (eCal[i] - ped_fits[0][i+config::eCal_start_channel]) / ecal_weight_masha[i];
			}
		}
		eCal_sum->Fill(ele_energy);
	}

	//Set gauss TF1 and fit election energy sum peak.
	TF1 *ele_energy_fit = new TF1("Electron_energy", "gaus");
	eCal_sum->Fit(ele_energy_fit);

	//Write sum of eCal electron energy histogram to file.
	eCal_sum->Write();

	//Delete objects from heap.
	delete eCal_sum;

}


vector<int> get_hCal_zero(TTree* tree, vector<int> events, vector<vector<float>> ped_fits) {
	float hCal[config::hCal_channels];
	TBranch *b_hCal = tree->GetBranch("HCalRawADC");
	b_hCal->SetAddress(hCal);

	vector<int> hCal_zero;
	float hCal_energy;
	for(auto event_index:events) {
		b_hCal->GetEntry(event_index);

		hCal_energy = 0.0;
		for(int i=0; i<config::hCal_channels; i++) {
			if(hCal[i] - ped_fits[0][i+config::hCal_start_channel] > 3 * ped_fits[1][i+config::hCal_start_channel]) {
				hCal_energy += (hCal[i] - ped_fits[0][i+config::hCal_start_channel]) / hcal_weight_masha[i];
			}
		}
		if(hCal_energy == 0.0) {
			hCal_zero.push_back(event_index);
		}
	}

	return(hCal_zero);
}


vector<int> get_cut(TTree* tree, vector<int> events, string branch, double lower, double upper) {
	float val_ADC;
	TBranch *b_branch = tree->GetBranch(branch.data());;
//	TBranch *b_branch = tree->GetBranch("Sc1RawADC");
	b_branch->SetAddress(&val_ADC);

	//Set hist to display cut.
	TH1D *cut_hist = new TH1D((branch + "_cut_hist").data(), (branch + "_cut").data(), config::bins, config::low_bin, config::high_bin);

	vector<int> filtered_events;
	for(auto event_index:events) {
		b_branch->GetEntry(event_index);
		cut_hist->Fill(val_ADC);
		if(val_ADC > lower && val_ADC < upper) {
			filtered_events.push_back(event_index);
		}
	}

	//Set canvas to display cut.
	TCanvas *cut_can = new TCanvas(((branch + "_cut_can").data()), ((branch + "_cut_can").data()), config::canvas_x, config::canvas_y);

	TLine *cut_line_lower = new TLine(lower, 0, lower, cut_hist->GetMaximum());
	cut_line_lower->SetLineColor(kRed);
	TLine *cut_line_upper = new TLine(upper, 0, upper, cut_hist->GetMaximum());
	cut_line_upper->SetLineColor(kRed);
	cut_hist->Draw();
	cut_line_lower->Draw();
	cut_line_upper->Draw();
	cut_can->Write();

	//Delete heap objects.
	delete cut_can;
	delete cut_line_lower;
	delete cut_line_upper;
	delete cut_hist;

	return(filtered_events);
}

vector<int> get_cut_cal(TTree* tree, vector<int> events, string branch, int channel, double lower, double upper) {
	float val_ADC[config::eCal_channels];
	TBranch *b_branch = tree->GetBranch(branch.data());;
//	TBranch *b_branch = tree->GetBranch("Sc1RawADC");
	b_branch->SetAddress(val_ADC);

	//Set hist to display cut.
	TH1D *cut_hist = new TH1D((branch + to_string(channel) + "_cut_hist").data(), (branch  + to_string(channel) + "_cut").data(), config::bins, config::low_bin, config::high_bin);

	vector<int> filtered_events;
	for(auto event_index:events) {
		b_branch->GetEntry(event_index);
		cut_hist->Fill(val_ADC[channel]);
		if(val_ADC[channel] > lower && val_ADC[channel] < upper) {
			filtered_events.push_back(event_index);
		}
	}

	//Set canvas to display cut.
	TCanvas *cut_can = new TCanvas(((branch + "_cut_can").data()), ((branch + "_cut_can").data()), config::canvas_x, config::canvas_y);

	TLine *cut_line_lower = new TLine(lower, 0, lower, cut_hist->GetMaximum());
	cut_line_lower->SetLineColor(kRed);
	TLine *cut_line_upper = new TLine(upper, 0, upper, cut_hist->GetMaximum());
	cut_line_upper->SetLineColor(kRed);
	cut_hist->Draw();
	cut_line_lower->Draw();
	cut_line_upper->Draw();
	cut_can->Write();

	//Delete heap objects.
	delete cut_can;
	delete cut_line_lower;
	delete cut_line_upper;
	delete cut_hist;

	return(filtered_events);
}
