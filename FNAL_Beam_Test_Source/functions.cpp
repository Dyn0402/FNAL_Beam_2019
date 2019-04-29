/*
 * functions.cpp
 *
 *  Created on: Apr 12, 2019
 *      Author: dylan
 */


#include <vector>

#include <TTree.h>

using namespace std;

vector<int> get_events_complement(vector<int> events, int num_events) {
	vector<int> complement;
	int event_index = 0;
	int next_event = -1;
	if(events.size() > 0) {next_event = events[event_index];}
	for(int i=0; i<num_events; i++) {
		if(i==next_event) {
			next_event = events[++event_index];
		} else {
			complement.push_back(i);
		}
	}

	return(complement);
}


vector<int> get_overlap_events(vector<int> events1, vector<int> events2, int num_events) {
	vector<int> overlap;
	int event1_index = 0;
	int next_event1 = -1;
	int event2_index = 0;
	int next_event2 = -1;
	if(events1.size() > 0) {next_event1 = events1[event1_index];}
	if(events2.size() > 0) {next_event2 = events2[event2_index];}
	for(int i = 0; i<num_events+1; i++) {
		if(i==next_event1 && i==next_event2) {
			overlap.push_back(i);
		}
		if(i==next_event1) {next_event1 = events1[++event1_index];}
		if(i==next_event2) {next_event2 = events2[++event2_index];}
	}

	return(overlap);
}


int get_num_events(TTree *tree) {
	TBranch *b_event_num = tree->GetBranch("eventno");
	return(b_event_num->GetEntries());
}


vector<double> get_moving_average(vector<double> list, unsigned int n) {
	vector<double> moving_average;

	if(list.size() >= n && n>0) {//Check to make sure there are enough elements (n) in list to get at least one average point.
		//Initialize moving average.
		double sum = 0;
		for(unsigned int i=0; i<n; i++) {
			sum+=list[i];
		}
		moving_average.push_back(sum/n);

		//Fill rest of moving average points by adding new value and subtracting oldest value from sum.
		for(unsigned int i=n; i<list.size(); i++) {
			sum += list[i] - list[i-n];
			moving_average.push_back(sum/n);
		}
	}

	return(moving_average);
}


double get_mean(vector<double> list) {
	double sum = 0;
	for(auto element:list) {
		sum += element;
	}
	return(sum/(int)list.size());
}

