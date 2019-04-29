/*
 * tree_pass_test.cpp
 *
 *  Created on: Apr 6, 2019
 *      Author: dylan
 */


#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"

TH1D* pass_tree(TTree *tree) {
	float Sc1ADC;
	TBranch *event_num = tree->GetBranch("eventno");
	TBranch *sc1 = tree->GetBranch("Sc1RawADC");
	sc1->SetAddress(&Sc1ADC);

	TH1D *sc1_hist = new TH1D("sc1_hist", "sc1 Hist", 4096, 0.5, 4096.5);

	int eventindex = 0;
	while(event_num->GetEntry(eventindex)) {
		sc1->GetEntry(eventindex);
		sc1_hist->Fill(Sc1ADC);
		eventindex++;
	}

	return(sc1_hist);
}
