
// V4: using variable binning
// V5: muon jet filtering, using HydJet_x.root files
// V6: using new weigthing method, pthat.  No longer defining binning in this program, need to move to plotting macros.


#include "myProcesses/hiforest/plugin/eventMap_hiForest.h"
#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TLatex.h>
#include <TCut.h>
#include "TDatime.h"
#include <vector>
#include "TCanvas.h"
 #include <dirent.h>  
 #include <stdio.h> 
 #include <string.h> 
 #include <stdlib.h>



void etaPtHistoV6(){

	vector <float> *pf_jtpt=0, *pf_jteta=0, *pf_jtphi=0, *pf_partonFlavor=0; 
	vector <float> *muPt=0, *muEta=0, *muPhi=0, *muChi2NDF=0, *muInnerD0=0, *muInnerDz=0;
	int nMu; //,muIsGlobal,muIsTracker, muMuonHits, muStations, muTrkLayers, muPixelHits;
	Float_t weight;
	vector <int> *muIsGlobal=0, *muIsTracker=0, *muMuonHits=0, *muStations=0, *muTrkLayers=0, *muPixelHits=0;

	/// DEFINE BINNING ///


	// create eta edges array
	const Int_t Nbins = 1000;
	Float_t etaEdges[Nbins+1];
	Float_t deltaEta = 0.2; // width of bin
	Float_t etamin=-1.5;
	Float_t etamax=1.5;
	Float_t jteta_check=etamin;
	etaEdges[0]=etamin;
	int numUsedEtaBins = 0;
	for(int zz=1; jteta_check<etamax; zz++){
		
		jteta_check=jteta_check+deltaEta;
		etaEdges[zz]=etaEdges[zz-1]+deltaEta;
		numUsedEtaBins=zz;
		
		
	}
	

	Float_t usedEtaEdges[numUsedEtaBins+1];
	for(int jj=0; jj<numUsedEtaBins+1; jj++){
		usedEtaEdges[jj]=etaEdges[jj];
	}
	

	// create pT edges array
	
	Float_t ptEdges[Nbins+1];
	Float_t deltaPt = 10.0; // initial width of bin (pT) // 10 GeV up to 100, 20 GeV up to 250, 250 out = 50 GeV, 400-500 = one bin
	Float_t ptmin=50;
	Float_t ptmax=600;
	Float_t jtpt_check=ptmin;
	ptEdges[0]=ptmin;
	int numUsedPtBins = 0;
	for(int yy=1; jtpt_check<ptmax; yy++){

		if(jtpt_check<100){
			deltaPt=10.;	
			ptEdges[yy]=ptEdges[yy-1]+deltaPt;
			numUsedPtBins=numUsedPtBins+1;
		}
		else if (jtpt_check<250){
			deltaPt=20.;
			ptEdges[yy]=ptEdges[yy-1]+deltaPt;
			numUsedPtBins=numUsedPtBins+1;
		}
		else if (jtpt_check<400){
			deltaPt=50;
			ptEdges[yy]=ptEdges[yy-1]+deltaPt;
			numUsedPtBins=numUsedPtBins+1;
			if(ptEdges[yy]>=400){
				ptEdges[yy]=400;
				ptEdges[yy+1]=ptmax;
				numUsedPtBins=numUsedPtBins+1;
				break;
			}
		}
		
		
		
		jtpt_check=jtpt_check+deltaPt;
		
	}
	

	Float_t usedPtEdges[numUsedPtBins+1];
	for(int ii=0; ii<numUsedPtBins+1; ii++){
		usedPtEdges[ii]=ptEdges[ii];
		cout << usedPtEdges[ii] << endl;
	}


	// DEFINE 2D HISTOGRAMS 
		TH1D *h_jetpt = new TH1D("h_jetpt","All jets",numUsedPtBins,usedPtEdges);
		TH1D *h_jeteta = new TH1D("h_jeteta","All jets",numUsedEtaBins,usedEtaEdges);

		TH2D *h2 = new TH2D("h2","All jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		// parse jet by species of particles
		TH2D *h2_d = new TH2D("h2_d","d jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_u = new TH2D("h2_u","u jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_s = new TH2D("h2_s","s jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_c = new TH2D("h2_c","c jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_b = new TH2D("h2_b","b jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_dbar = new TH2D("h2_dbar","#bar{d} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ubar = new TH2D("h2_ubar","#bar{u} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_sbar = new TH2D("h2_sbar","#bar{s} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_cbar = new TH2D("h2_cbar","#bar{c} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_bbar = new TH2D("h2_bbar","#bar{b} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_g = new TH2D("h2_g","gluon jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_hq = new TH2D("h2_hq","heavy quark jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_lq = new TH2D("h2_lq","light quark jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_qall = new TH2D("h2_qall","All quark jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_dall = new TH2D("h2_dall","d & #bar{d} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_uall = new TH2D("h2_uall","u & #bar{u} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_sall = new TH2D("h2_sall","s & #bar{s} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_call = new TH2D("h2_call","c & #bar{c} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ball = new TH2D("h2_ball","b & #bar{b} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ud = new TH2D("h2_ud","u & d jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ubardbar = new TH2D("h2_ubardbar","#bar{u} & #bar{d} jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_other = new TH2D("h2_other","ghost jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);


		// muon-tagged jets
		TH2D *h2_MJ = new TH2D("h2_MJ","All muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		// parse jet by species of particles
		TH2D *h2_d_MJ = new TH2D("h2_d_MJ","d muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_u_MJ = new TH2D("h2_u_MJ","u muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_s_MJ = new TH2D("h2_s_MJ","s muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_c_MJ = new TH2D("h2_c_MJ","c muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_b_MJ = new TH2D("h2_b_MJ","b muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_dbar_MJ = new TH2D("h2_dbar_MJ","#bar{d} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ubar_MJ = new TH2D("h2_ubar_MJ","#bar{u} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_sbar_MJ = new TH2D("h2_sbar_MJ","#bar{s} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_cbar_MJ = new TH2D("h2_cbar_MJ","#bar{c} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_bbar_MJ = new TH2D("h2_bbar_MJ","#bar{b} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_g_MJ = new TH2D("h2_g_MJ","gluon muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_hq_MJ = new TH2D("h2_hq_MJ","heavy quark muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_lq_MJ = new TH2D("h2_lq_MJ","light quark muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_qall_MJ = new TH2D("h2_qall_MJ","All quark muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_dall_MJ = new TH2D("h2_dall_MJ","d & #bar{d} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_uall_MJ = new TH2D("h2_uall_MJ","u & #bar{u} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_sall_MJ = new TH2D("h2_sall_MJ","s & #bar{s} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_call_MJ = new TH2D("h2_call_MJ","c & #bar{c} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ball_MJ = new TH2D("h2_ball_MJ","b & #bar{b} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ud_MJ = new TH2D("h2_ud_MJ","u & d muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_ubardbar_MJ = new TH2D("h2_ubardbar_MJ","#bar{u} & #bar{d} muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);
		TH2D *h2_other_MJ = new TH2D("h2_other_MJ","ghost muon-tagged-jets",numUsedEtaBins,usedEtaEdges,numUsedPtBins,usedPtEdges);

		TH1D *deltaR = new TH1D("deltaR","#Delta r",10,0,1);

h2->Sumw2();
h2_d->Sumw2();
h2_u->Sumw2();
h2_s->Sumw2();
h2_c->Sumw2();
h2_b->Sumw2();
h2_dbar->Sumw2();
h2_ubar->Sumw2();
h2_sbar->Sumw2();
h2_cbar->Sumw2();
h2_bbar->Sumw2();
h2_hq->Sumw2();
h2_lq->Sumw2();
h2_qall->Sumw2();
h2_dall->Sumw2();
h2_uall->Sumw2();
h2_sall->Sumw2();
h2_call->Sumw2();
h2_ball->Sumw2();
h2_g->Sumw2();
h2_ud->Sumw2();
h2_ubardbar->Sumw2();
h2_other -> Sumw2();
h2_MJ->Sumw2();
h2_d_MJ->Sumw2();
h2_u_MJ->Sumw2();
h2_s_MJ->Sumw2();
h2_c_MJ->Sumw2();
h2_b_MJ->Sumw2();
h2_dbar_MJ->Sumw2();
h2_ubar_MJ->Sumw2();
h2_sbar_MJ->Sumw2();
h2_cbar_MJ->Sumw2();
h2_bbar_MJ->Sumw2();
h2_hq_MJ->Sumw2();
h2_lq_MJ->Sumw2();
h2_qall_MJ->Sumw2();
h2_dall_MJ->Sumw2();
h2_uall_MJ->Sumw2();
h2_sall_MJ->Sumw2();
h2_call_MJ->Sumw2();
h2_ball_MJ->Sumw2();
h2_g_MJ->Sumw2();
h2_ud_MJ->Sumw2();
h2_ubardbar_MJ->Sumw2();
h2_other_MJ -> Sumw2();
deltaR -> Sumw2();
	
	


	//// CUTS ///////
	const double etamaxcut = 1.5; // default = 1.5
	const double etamincut = 0.;
	const double pTmincut = 50.; // default = 120.
	const double pTmaxcut = 600.;
	


	/// LOAD DATA ///
	int NFiles = 0;
	
	struct dirent *de;  // Pointer for directory entry 
	DIR *dr = opendir("/home/clayton/Analysis/skims4/0000"); 

   	if (dr == NULL)  // opendir returns NULL if couldn't open directory 
    	{ 
        	printf("Could not open current directory" ); 
        	return 0; 
    	}
		
    	if(dr){ 
		 
		while((de=readdir(dr)) != NULL) {
		if(strcmp(de->d_name,".")!=0 && strcmp(de->d_name,"..")!=0){
		  	NFiles++;
		}
    	    	}
 	closedir(dr);
	}

	cout <<"number of files ="<< NFiles << endl;




for(int file = 1; file < NFiles+1; file++){
	if(file==126){continue;} // missing this root file
	cout << "Processing file " << file << "/"<<NFiles<< endl;
	auto f = TFile::Open(Form("/home/clayton/Analysis/skims4/0000/HydJet_%d.root",file));

    TTree *inp_tree = (TTree*)f->Get("jet_tree;1");
	Long64_t n_evts = inp_tree->GetEntriesFast();

	
	inp_tree->SetBranchAddress("flowpf_jtpt",&pf_jtpt);
	inp_tree->SetBranchAddress("flowpf_jteta",&pf_jteta);
	inp_tree->SetBranchAddress("flowpf_jtphi",&pf_jtphi);
	inp_tree->SetBranchAddress("flowpf_partonFlavor",&pf_partonFlavor);
	inp_tree->SetBranchAddress("muPt",&muPt);
	inp_tree->SetBranchAddress("muEta",&muEta);
	inp_tree->SetBranchAddress("muPhi",&muPhi);
	inp_tree->SetBranchAddress("muChi2NDF",&muChi2NDF);
	inp_tree->SetBranchAddress("muInnerD0",&muInnerD0);
	inp_tree->SetBranchAddress("muInnerDz",&muInnerDz);
	inp_tree->SetBranchAddress("nMu",&nMu);
	inp_tree->SetBranchAddress("muIsGlobal",&muIsGlobal);
	inp_tree->SetBranchAddress("muIsTracker",&muIsTracker);
	inp_tree->SetBranchAddress("muMuonHits",&muMuonHits);
	inp_tree->SetBranchAddress("muStations",&muStations);
	inp_tree->SetBranchAddress("muTrkLayers",&muTrkLayers);
	inp_tree->SetBranchAddress("muPixelHits",&muPixelHits);
	inp_tree->SetBranchAddress("weight",&weight);

	
    
	
	


	// loop variables //	
	int jeti_frac=0;
	int evi = 0;
	int evi_frac = 0;
	double w=0;
	
	int nMuCount=0;
	// Event loop
	for (evi = 0; evi < n_evts; evi++){

		
        	inp_tree->GetEntry(evi);
               if((100*evi/n_evts)%5==0 && 100*evi/n_evts > evi_frac){
               	//cout<<"evt frac: "<<evi_frac<<"%"<<endl;
                }
        evi_frac = 100*evi/n_evts;
        
        w=weight; 
		// Jet loop
		Double_t muPtCut = 10;
		for(int jetj=0; jetj < (int)pf_jteta->size(); jetj++){
			
			Double_t jetPtj = pf_jtpt->at(jetj);
			Double_t jetEtaj = pf_jteta->at(jetj);
			Double_t jetPhij = pf_jtphi->at(jetj);
			
			double x = pf_jteta->at(jetj);
			double u = pf_jtpt->at(jetj);
			
			if(fabs(x)>etamaxcut || u < pTmincut || u>pTmaxcut){continue;}
			int y = pf_partonFlavor->at(jetj);
			/*
			double w=1; // weight.  will equal inverse area
			double de = 0.2;
			double dp = 10.0;
			if(u<100){dp=10.;}
			else if(u<260){dp=20.;}
			else if(u<360){dp=50.;}
			else if(u<400){dp=40.;}
			else{dp=200;}
			double a = dp*de;
			w=1/a;
			*/
			//Float_t w=0;
			


			if(y==0){h2_other -> Fill(x,u,w);}
			if(y==1){h2_d->Fill(x,u,w); h2_lq->Fill(x,u,w); h2_dall->Fill(x,u,w); h2_ud->Fill(x,u,w);}
			if(y==2){h2_u->Fill(x,u,w); h2_lq->Fill(x,u,w); h2_uall->Fill(x,u,w); h2_ud->Fill(x,u,w);}
			if(y==3){h2_s->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_sall->Fill(x,u,w);}
			if(y==4){h2_c->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_call->Fill(x,u,w);}
			if(y==5){h2_b->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_ball->Fill(x,u,w);}
			if(y==-1){h2_dbar->Fill(x,u,w); h2_lq->Fill(x,u,w); h2_dall->Fill(x,u,w); h2_ubardbar->Fill(x,u,w);}
			if(y==-2){h2_ubar->Fill(x,u,w); h2_lq->Fill(x,u,w); h2_uall->Fill(x,u,w); h2_ubardbar->Fill(x,u,w);}
			if(y==-3){h2_sbar->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_sall->Fill(x,u,w);}
			if(y==-4){h2_cbar->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_call->Fill(x,u,w);}
			if(y==-5){h2_bbar->Fill(x,u,w); h2_hq->Fill(x,u,w); h2_ball->Fill(x,u,w);}
			if(y==21){h2_g->Fill(x,u,w);}


			h2->Fill(x,u,w);
			

			if(nMu==0){continue;}
			Double_t deltaRmin=100;
			for(int mui=0; mui<(int) muPt->size();mui++){
				if(muIsTracker->at(mui)==0 || TMath::Abs(muEta->at(mui))>2.4 || muChi2NDF->at(mui)==-99 || muChi2NDF->at(mui)>10
					||TMath::Abs(muInnerD0->at(mui))>0.2 || TMath::Abs(muInnerDz->at(mui))>0.5 || muMuonHits->at(mui)<= 0
					|| muStations->at(mui)<= 1 || muTrkLayers->at(mui)<=5 || muPixelHits->at(mui)<=0 || muPt->at(mui) < muPtCut){continue;}
			
				Double_t muPti = muPt->at(mui);
				Double_t muEtai = muEta->at(mui);
				Double_t muPhii = muPhi->at(mui);
							
				Double_t deltaEtaij = muEtai-jetEtaj;
				Double_t deltaPhiij = acos(cos(muPhii-jetPhij));
				Double_t deltaRij = sqrt(pow(deltaEtaij,2)+pow(deltaPhiij,2));
				
				if(deltaRij<deltaRmin){
					deltaRmin=deltaRij;
				}

			}  // end muon loop
	
			deltaR->Fill(deltaRmin);
			if(deltaRmin<0.4){
				if(y==0){h2_other_MJ -> Fill(x,u,w);}
				if(y==1){h2_d_MJ->Fill(x,u,w); h2_lq_MJ->Fill(x,u,w); h2_dall_MJ->Fill(x,u,w); h2_ud_MJ->Fill(x,u,w);}
				if(y==2){h2_u_MJ->Fill(x,u,w); h2_lq_MJ->Fill(x,u,w); h2_uall_MJ->Fill(x,u,w); h2_ud_MJ->Fill(x,u,w);}
				if(y==3){h2_s_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_sall_MJ->Fill(x,u,w);}
				if(y==4){h2_c_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_call_MJ->Fill(x,u,w);}
				if(y==5){h2_b_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_ball_MJ->Fill(x,u,w);}
				if(y==-1){h2_dbar_MJ->Fill(x,u,w); h2_lq_MJ->Fill(x,u,w); h2_dall_MJ->Fill(x,u,w); h2_ubardbar_MJ->Fill(x,u,w);}
				if(y==-2){h2_ubar_MJ->Fill(x,u,w); h2_lq_MJ->Fill(x,u,w); h2_uall_MJ->Fill(x,u,w); h2_ubardbar_MJ->Fill(x,u,w);}
				if(y==-3){h2_sbar_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_sall_MJ->Fill(x,u,w);}
				if(y==-4){h2_cbar_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_call_MJ->Fill(x,u,w);}
				if(y==-5){h2_bbar_MJ->Fill(x,u,w); h2_hq_MJ->Fill(x,u,w); h2_ball_MJ->Fill(x,u,w);}
				if(y==21){h2_g_MJ->Fill(x,u,w);}

				h2_MJ->Fill(x,u,w);

			}
		} // end jet loop

	} // end event loop


} // end file loop





//auto wf = TFile::Open("etaPtHistoV6.root", "recreate");
//auto wf = TFile::Open("etaPtHistoV6_mutptcut_5.root","recreate");
auto wf = TFile::Open("etaPtHistoV6_mutptcut_10.root","recreate");
//auto wf = TFile::Open("etaPtHistoV6_mutptcut_15.root","recreate");
h2->Write();
h2_d->Write();
h2_u->Write();
h2_s->Write();
h2_c->Write();
h2_b->Write();
h2_dbar->Write();
h2_ubar->Write();
h2_sbar->Write();
h2_cbar->Write();
h2_bbar->Write();
h2_hq->Write();
h2_lq->Write();
h2_qall->Write();
h2_dall->Write();
h2_uall->Write();
h2_sall->Write();
h2_call->Write();
h2_ball->Write();
h2_g->Write();
h2_ud->Write();
h2_ubardbar->Write();
h2_other->Write();
h2_MJ->Write();
h2_d_MJ->Write();
h2_u_MJ->Write();
h2_s_MJ->Write();
h2_c_MJ->Write();
h2_b_MJ->Write();
h2_dbar_MJ->Write();
h2_ubar_MJ->Write();
h2_sbar_MJ->Write();
h2_cbar_MJ->Write();
h2_bbar_MJ->Write();
h2_hq_MJ->Write();
h2_lq_MJ->Write();
h2_qall_MJ->Write();
h2_dall_MJ->Write();
h2_uall_MJ->Write();
h2_sall_MJ->Write();
h2_call_MJ->Write();
h2_ball_MJ->Write();
h2_g_MJ->Write();
h2_ud_MJ->Write();
h2_ubardbar_MJ->Write();
h2_other_MJ->Write();
deltaR->Write();
wf->Close();



}


 // end program
