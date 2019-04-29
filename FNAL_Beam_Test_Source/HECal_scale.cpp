/*
 * EHCal_scale.cpp
 *
 *  Created on: Apr 16, 2019
 *      Author: dylan
 */


#include <vector>

#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TLine.h>
#include <TSpectrum.h>

#include "config.h"

using namespace std;


TF1 *get_HECal_scale(TTree *tree, vector<int> events, vector<vector<float>> ped_fits, double eCal_scale_factor, double low_rms, double high_rms) {
	//Set branches from tree.
	float hCal[config::hCal_channels];
	float eCal[config::eCal_channels];
	TBranch *b_hCal = tree->GetBranch("HCalRawADC");
	TBranch *b_eCal = tree->GetBranch("ECalRawADC");
	b_hCal->SetAddress(hCal);
	b_eCal->SetAddress(eCal);

	//Set Histograms
	TH1D *hadron_sum = new TH1D(("Hadron_Sum"+to_string(eCal_scale_factor)).data(), ("HCal + k*ECal Sum Energy"+to_string(eCal_scale_factor)).data(), config::bins, config::low_bin, config::high_bin*1.5);
	TH2D *ECal_vs_HCal = new TH2D(("ECal_vs_HCal"+to_string(eCal_scale_factor)).data(), ("ECal_vs_HCal"+to_string(eCal_scale_factor)).data(), 100, config::low_bin, config::high_bin, 100, config::low_bin, config::high_bin);

	//Fill channel histograms.
	double hCal_sum_hadron, eCal_sum_hadron, sum_hadron;
	for(auto event_index:events) {
		b_hCal->GetEntry(event_index);
		b_eCal->GetEntry(event_index);

		hCal_sum_hadron = 0;
		for(int i=0; i<config::hCal_channels; i++) {
			if(hCal[i] > ped_fits[0][i+config::hCal_start_channel] + ped_fits[1][i+config::hCal_start_channel] * config::eh_scale_hCal_signal_sigma) {
				hCal_sum_hadron += (hCal[i] - ped_fits[0][i+config::hCal_start_channel]) / config::hCal_weights[i];
			}
		}
		eCal_sum_hadron = 0;
		for(int i=0; i<config::eCal_channels; i++) {
			if(eCal[i] > ped_fits[0][i+config::eCal_start_channel] + ped_fits[1][i+config::eCal_start_channel] * config::eh_scale_eCal_signal_sigma) {
				eCal_sum_hadron += (eCal[i] - ped_fits[0][i+config::eCal_start_channel]) / config::eCal_weights[i];
			}
		}
		ECal_vs_HCal->Fill(eCal_sum_hadron, hCal_sum_hadron);

		sum_hadron = hCal_sum_hadron + eCal_scale_factor * eCal_sum_hadron;
		hadron_sum->Fill(sum_hadron);
	}

	TCanvas *ECal_vs_HCal_can = new TCanvas(("ECal_vs_HCal_can"+to_string(eCal_scale_factor)).data(), ("ECal vs HCal Sum of Hadron Energy"+to_string(eCal_scale_factor)).data(), config::canvas_x, config::canvas_y);
	gPad->SetLogz();
	ECal_vs_HCal->Draw("COLZ");
	ECal_vs_HCal_can->Write();
	ECal_vs_HCal->Write();

	TCanvas *hadron_sum_can = new TCanvas(("Hadron_Sum_can"+to_string(eCal_scale_factor)).data(), ("ECal HCal Hadron Sum Energy"+to_string(eCal_scale_factor)).data(), config::canvas_x, config::canvas_y);

	TF1 *gauss = new TF1(("hadron_sum_gauss "+to_string(eCal_scale_factor)).data(), "gaus"); //Hardcode fix
	TF1 *gauss_bkg = new TF1(("hadron_sum_gauss2 "+to_string(eCal_scale_factor)).data(), "gaus"); //Hardcode fix
	hadron_sum->Fit(gauss, "Q");
	gauss_bkg->SetRange(gauss->GetParameter(1)-low_rms*gauss->GetParameter(2), gauss->GetParameter(1)+high_rms*gauss->GetParameter(2));
	hadron_sum->Fit(gauss_bkg, "QR");
	hadron_sum->Draw();
	hadron_sum_can->Write();
	hadron_sum->Write();

	delete ECal_vs_HCal_can;
	delete ECal_vs_HCal;
	delete hadron_sum_can;
	delete hadron_sum;
	delete gauss;

	return(gauss_bkg);
}



