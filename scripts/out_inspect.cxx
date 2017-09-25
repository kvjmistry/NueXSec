#include "../xsecAna/TpcObjectContainer.h"
#include "../xsecAna/ParticleContainer.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TInterpreter.h"
#include "TROOT.h"

#include <iostream>
#include <vector>

#include "../xsecAna/LinkDef.h"

int out_inspect()
{

	const char * _file1 = "../nue_xsec_extraction.root";
	std::cout << "File Path: " << _file1 << std::endl;
	//first we need to open the root file
	TFile * f = new TFile(_file1);
	if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; return 1; }
	TTree * mytree = (TTree*)f->Get("AnalyzeTPCO/tree");

	std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v = nullptr;
	mytree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

	const int total_entries = mytree->GetEntries();

	for(int i = 0; i < total_entries; i++)
	{
		mytree->GetEntry(i);
		const int n_tpc_obj = tpc_object_container_v->size();
		std::cout << "Number of TPC Objects: " << n_tpc_obj << std::endl;
		for(auto const tpc_obj : *tpc_object_container_v)
		{
			const int n_pfp = tpc_obj.NumPFParticles();
			std::cout << " \t Number of PFParticles: " << n_pfp << std::endl;
		}
		//std::cout << tpc_object_container_v->at(0).RunNumber() << std::endl;
	}//end looping events

	return 0;
}

#ifndef __ROOTCLING__


int main(){

	return out_inspect();
}

#endif
