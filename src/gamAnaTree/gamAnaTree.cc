////////////////////////////////////////////////
// Run run_gamAnaTree.C with condor
// Use int process argument to loop through files of inlist
//////////////////////////////////////////////////

#include "gamAnaTree.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include <ctime>

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>

#include <TSystem.h>
#include <TApplication.h>

#include "../CentralTrackData.h"
#include "ChargedTrack.h"
#include "CentralArmTrack.h"
#include "TrackCollectionBase.h"
#include "TrackCollection.h"
#include "TriggerParticle.h"
#include "../gamConst.h"

#ifdef __MAKECINT__
#pragma link C++ class CentralTrackData+;
#pragma link C++ class ChargedTrack+;
#pragma link C++ class CentralArmTrack+;
#pragma link C++ class TriggerParticle+;
#pragma link C++ class std::vector<CentralArmTrack>+;
#pragma link C++ class TrackCollectionBase+;
#pragma link C++ class TrackCollection<CentralArmTrack>+;
#pragma link C++ class std::map<std::string,cut>+;
#pragma link C++ namespace GA+;
#endif

using namespace std;

const float PI = (float) TMath::Pi();

const float dphibw = (maxdphi - mindphi)/ndphibin;

void gamAnaTree(
    const int iprocess = -1,
    const bool dochargeveto = false,
    const bool usepc3map = true,
    const bool use3x3pc3map = true,
    const bool dofidcut = true,
    const bool dopsimix = true,
    const bool dozvmix = true,
    const int onlyarm = -1,
    const char* inlistname = "inlist_CNT.txt",
    const char* outdir = "wrk/"
    )
{
    time_t now = time(NULL);
    cout << asctime(localtime(&now)) << endl;

    cutmap["chi2ndf"] = CHI2NDF;
    cutmap["dcaz"] = DCAZ;
    cutmap["dca2d"] = DCA2D;
    cutmap["pc3dphi"] = PC3DPHI;
    cutmap["pc3dz"] = PC3DZ;
    cutmap["emcdphi"] = EMCDPHI;
    cutmap["emcdz"] = EMCDZ;

    memset(emcwu,0,sizeof(emcwu));
    memset(pc3bad,0,sizeof(pc3bad));
    memset(pc3bad3x3,0,sizeof(pc3bad3x3));
    memset(pc3dphispar,0,sizeof(pc3dphispar));
    memset(pc3dzpar,0,sizeof(pc3dzpar));
    memset(h2commap,0,sizeof(h2commap));
    memset(h2runmap,0,sizeof(h2runmap));


    verbosity = 1;//default 1
    bool dofulltrackana = false;
    bool dotrackeventana = false;
    bool skippi0 = false;

    cout<<"Applying charged particle veto? "<<dochargeveto<<endl;
    cout<<"Use pc3 map? "<<usepc3map<<endl;
    cout<<"Do pc3 3x3 map? "<<use3x3pc3map<<endl;
    cout<<"Applying fiducial cut? "<<dofidcut<<endl;
    cout<<"Mixing in psi? "<<dopsimix<<endl;
    cout<<"Mixing in zv? "<<dozvmix<<endl;
    cout<<"Which arm? "<<onlyarm<<endl;

    cout<<"Do full track analysis? "<<dofulltrackana<<endl;
    cout<<"Do track event analysis? "<<dotrackeventana<<endl;
    cout<<"Skipping pi0s? "<<skippi0<<endl;

    //==================================================
    // Getting input file name
    //==================================================
    string inname;
    if (iprocess >= 0)
    {
	int filectr = 0;
	ifstream inlist(inlistname,ifstream::in);
	if (inlist.is_open())
	{
	    while (getline(inlist,inname))
	    {
		if (filectr==iprocess) break;
		filectr++;
	    }
	    inlist.close();
	}
	else
	{
	    cerr<<"Cannot open "<<inlistname<<endl;
	    return;
	}
    }
    else
    {

	inname = "test_input.root";
    }
    //==================================================
    // Input file name has been set
    //==================================================


    //==================================================
    // Set up output file name
    //==================================================
    ostringstream outoss;
    if (iprocess >= 0)
    {
	int dotpos = inname.find(".");
	string str1 = inname.substr(dotpos-6,6);
	outoss<<outdir<<str1;
	if (onlyarm >= 0)
	    outoss<<"_arm"<<onlyarm;
	outoss<<"_analyzed.root";
    }
    else
    {
	outoss<<outdir<<"test_analyzed.root";
    }
    //==================================================
    // Output filename has been set
    //==================================================


    //========== cout some file information ==========
    if (verbosity > 0)
    {
	cout<<"Running gamAnaTree"<<endl;
	cout<<"Input = "<<inname<<endl;
	cout<<"Output = "<<outoss.str()<<endl;
    }
    //================================================

    //==================================================
    // Open input file
    //==================================================
    infile = TFile::Open(inname.c_str());
    if (!infile)
    {
	cerr<<"Line "<<__LINE__<<": Unable to open "<<inname <<endl;
	cerr<<"Exiting."<<endl;

	return;
    }
    if (infile->IsZombie())
    {
	cerr<<"Line "<<__LINE__<<": "<<inname<<" is zombie!"<<endl;
	cerr<<"Exiting."<<endl;

	return;
    }
    else if (verbosity > 0)
    {
	cout<<"Input opened successfully"<<endl;
    }
    //==================================================
    // Input file has been opened
    //==================================================

    //==================================================
    // Create output file and initialize histograms
    //==================================================
    outfile = new TFile(outoss.str().c_str(),"RECREATE");
    if (!outfile)
    {
	cerr<<"Line "<<__LINE__<<": Unable to open "<<outoss.str() <<endl;
	cerr<<"Exiting."<<endl;

	return;
    }
    else if (verbosity > 0)
    {
	cout<<"Output opened successfully"<<endl;
    }

    InitHistos();

    if (GetEmcMaps() < 0)
    {
	cerr<<"Could not retrieve emc maps. Exiting."<<endl;
	return;
    }
    if (GetPc3Maps() < 0)
    {
	cerr<<"Could not retrieve pc3 maps. Exiting."<<endl;
	return;
    }
    if (GetPc3Pars() < 0)
    {
	cerr<<"Could not retrieve pc3 matching parameters. Exiting."<<endl;
	return;
    }

    const char* loc_mis = "mistag.root";
    const char* loc_tag = "tageff.root";

    TF1* ftageff[ncent_fun];
    TF1* fdeceff[ncent_fun];
    TF1* fmistag[ncent_fun];

    if (GetFunctions(loc_mis,"fmistag_",fmistag) < 0)
    {
	cerr<<"Could not retrieve mistag functions. Exiting."<<endl;
	return;
    }

    if (GetFunctions(loc_tag,"ftageff_",ftageff) < 0)
    {
	cerr<<"Could not retrieve tageff functions. Exiting."<<endl;
	return;
    }

    if (GetFunctions(loc_tag,"fdeceff_",fdeceff) < 0)
    {
	cerr<<"Could not retrieve deceff functions. Exiting."<<endl;
	return;
    }

    //============= Output ready =======================

    //introducing a method to study all saved tracks
    if (dofulltrackana)
    {
	if (verbosity > 0) cout<<"Conducting a full analysis on tracks."<<endl;
	TrackAnalysis();
	return End();
    }

    //introducing method to study track events
    if (dotrackeventana)
    {
	if (verbosity > 0) cout<<"Conducting an analysis on the track event trees."<<endl;
	TrackEventAnalysis();

	return End();

    }

    //==================================================
    // Getting trigger tree and setting all branches
    //==================================================
    TTree* trigtree = (TTree*) infile->Get("TrigTree");
    if (!trigtree)
    {
	cerr<<"Line "<<__LINE__<<": Cannot find trigtree in file "<< inname <<endl;
	cerr<<"Exiting."<<endl;

	return;
    }
    else if (verbosity > 0)
    {
	cout<<"Found tree "<<trigtree->GetName()<<endl;
    }

    //set trigtree branches
    TriggerParticle* trg = NULL;
    TBranch* trigBranch = trigtree->GetBranch("TrigPart");
    trigBranch->SetAddress(&trg);

    int runno;
    TBranch* runBranch = trigtree->GetBranch("RunNumber");
    runBranch->SetAddress(&runno);

    int evtno;
    TBranch* evtBranch = trigtree->GetBranch("EventNumber");
    evtBranch->SetAddress(&evtno);

    float psi;
    TBranch* psiBranch = trigtree->GetBranch("Psi");
    psiBranch->SetAddress(&psi);

    float zv;
    TBranch* z0Branch = trigtree->GetBranch("Zvertex");
    z0Branch->SetAddress(&zv);

    float cent;
    TBranch* centBranch = trigtree->GetBranch("Centrality");
    centBranch->SetAddress(&cent);

    if (verbosity > 0)
	cout<<"Set branches for "<<trigtree->GetName()<<endl;

    //==================================================
    // Got trigger tree and set all branches
    //==================================================


    //==================================================
    //loop through triggers in trigger tree
    //==================================================
    unsigned int ntrig = trigtree->GetEntries();
    if (verbosity > 0)
	cout<<"Begin looping through "<<ntrig<<" triggers"<<endl;

    int tmprunno = -999;

    for (unsigned int itrig = 0; itrig < ntrig; itrig++)
    {

	trigtree->GetEntry(itrig);

	if (verbosity > 1)
	{
	    cout<<"\n=================================================="<<endl;
	    cout<<"Beginning analysis on trigger "<<itrig<<"/"<<ntrig-1<<endl;
	    cout<<"Event Info: Run #"<<runno<<"; Event #"<<evtno<<endl;
	    trg->Print();
	    now = time(NULL);
	    cout << asctime(localtime(&now));
	    cout<<"==================================================\n"<<endl;
	}
	else if (verbosity > 0 && itrig%50==0)
	{
	    cout<<"Analyzing triggered event "<<itrig<<"/"<<ntrig-1<<endl;
	    now = time(NULL);
	    cout << asctime(localtime(&now)) << endl;
	}

	trigevtno.push_back(evtno);

	if (maxcent != limcent[ncent])
	{
	    cerr<<"maxcent "<<maxcent<<" != limcent[ncent] "<<limcent[ncent]<<endl;
	    return;
	}
	
	if (cent > maxcent)
	{
	    if (verbosity > 1)
		cout<<"Centrality is "<<cent<<" which is greater than "<<maxcent<<". Skipping."<<endl;
	    continue;
	}

	if (onlyarm >= 0)
	{
	    float tmpphi = trg->GetPhi();
	    while (tmpphi < 0)
		tmpphi += 2.*PI;
	    while (tmpphi >= 2.*PI)
		tmpphi -= 2.*PI;

	    if (tmpphi < PI/2. || tmpphi > 3.*PI/2.)
	    {
		if (onlyarm != 0) continue;
	    }
	    else if (onlyarm != 1) continue;
	}

	float ptt = trg->GetPt();

	//=============== Skip if pt out of range ===============

	if (ptt < limptt[0])
	{
	    if (verbosity > 1)
		cout<<"Skipping because trigger pt lower than minimum."<<endl;
	    continue;
	}
	if (ptt > limptt[nptt])
	{
	    if (verbosity > 1)
		cout<<"Skipping because trigger pt higher than maximum."<<endl;
	    continue;
	}

	//get run specific emc hot map
	if (tmprunno != runno)
	{
	    if (GetEmcRunMaps(runno) < 0)
	    {
	    	cerr<<"Could not retrieve emc run specific maps. Exiting."<<endl;
	    	return;
	    }
	    else if (verbosity > 0)
	    {
		cout<<"Emc Run Specific Hot Maps retrieved successfully."<<endl;
	    }
	    tmprunno = runno;
	}


	if (trg->type == PI0)
	{
	    if (skippi0)
	    {
		if (verbosity > 1)
		    cout<<"Skipping because trigger is pi0."<<endl;
		continue;
	    }
	}
	else
	{
	    int sec = -999;
	    int ypos = -999;
	    int zpos = -999;
	    GA::TowerToSYZ(trg->towid,sec,ypos,zpos);

	    if (emcwu[sec][ypos][zpos])
	    {
		if (verbosity > 1)
		    cout<<"Trigger failed emcwu. Skipping."<<endl;

		continue;
	    }

	    if (h2commap[sec]->GetBinContent(zpos+1,ypos+1) > GA::maxhotrun)
	    {
		if (verbosity > 1)
		    cout<<"Trigger failed emc common hot map. Skipping."<<endl;

		continue;
	    }

	    if (h2runmap[sec]->GetBinContent(zpos+1,ypos+1) > 0)
	    {
		if (verbosity > 1)
		    cout<<"Trigger failed emc run specific hot map. Skipping."<<endl;
		continue;
	    }


	    if (dofidcut)
	    {
		bool fidok = GA::PassFidCutTowID(trg->GetMom(),trg->towid);
		if (!fidok)
		{
		    if (verbosity > 1)
			cout<<"Skipping because trigger failed fiducial cut."<<endl;
		    continue;
		}
	    }

	    if (dochargeveto)
	    {
		if (usepc3map)
		{
		    if (pc3bad[sec][ypos][zpos] || (use3x3pc3map && pc3bad3x3[sec][ypos][zpos]))
		    {
			if (verbosity > 1)
			    cout<<"Skipping because trigger failed pc3 dead map."<<endl;
			continue;
		    }
		}

		float pc3dphis = trg->pc3dphis;
		float pc3dz = trg->pc3dz;

		int pc3zbin = GA::GetTowZBin(sec,zpos);
		int pc3ybin = GA::GetTowYBin(sec,ypos);

		float pc3dphis_mean = pc3dphispar[sec][pc3ybin][GA::IMEAN_PC3];
		float pc3dphis_sigma = pc3dphispar[sec][pc3ybin][GA::ISIGMA_PC3];
		float pc3sdphis = (pc3dphis - pc3dphis_mean)/pc3dphis_sigma;
		h1pc3sdphis_trg->Fill(pc3sdphis);

		float pc3dz_mean = pc3dzpar[sec][pc3zbin][GA::IMEAN_PC3];
		float pc3dz_sigma = pc3dzpar[sec][pc3zbin][GA::ISIGMA_PC3];
		float pc3sdz = (pc3dz - pc3dz_mean)/pc3dz_sigma;
		h1pc3sdz_trg->Fill(pc3sdz);

		float pc3sdr = sqrt(pc3sdz*pc3sdz + pc3sdphis*pc3sdphis);
		h1pc3sdr_trg->Fill(pc3sdr);

		if (pc3sdr < GA::NSIGMA_PC3DR)
		{
		    if (verbosity > 1)
			cout<<"Skipping because trigger failed charge veto."<<endl;
		    continue;
		}
	    }

	    h2zypos[sec]->Fill(zpos,ypos);
	}

	int izv_trg = GA::GetZvBin(zv);
	if (verbosity > 2)
	{
	    cout<<"izv_trg = "<<izv_trg<<endl;
	}
	if (izv_trg < 0) return;

	vector<int> izv_tree;
	if (dozvmix)
	{
	    if (verbosity > 2)
	    {
		cout<<"Doing zv mixing. izv_tree.push_back("<<izv_trg<<");"<<endl;
	    }
	    izv_tree.push_back(izv_trg);
	}
	else
	{
	    if (verbosity > 2)
	    {
		cout<<"Not doing zv mixing.  Adding all zvbin to vector."<<endl;
	    }
	    for (int izv = 0; izv < GA::nzv_tree; izv++)
	    {
		if (verbosity > 2)
		{
		    cout<<"izv_tree.push_back("<<izv_trg<<");"<<endl;
		}
		izv_tree.push_back(izv);
	    }
	}

	//======================================================================

	icent_ana = GetCentBin(cent);
	if (verbosity > 2)
	{
	    cout<<"icent_ana = GetCentBin("<<cent<<") = "<<icent_ana<<endl;
	}


	vector<int> icent_tree;

	for (int icent = 0; icent < GA::ncent_tree; icent++)
	{
	    for (float cent = GA::bwcent_tree*icent; cent < GA::bwcent_tree*(icent+1); cent++)
	    {
		if (GetCentBin(cent) == icent_ana)
		{
		    icent_tree.push_back(icent);
		    if (verbosity > 2)
			cout<<"cent tree bin "<<icent<<" pushed to back of vector"<<endl;
		    break;
		}
	    }
	}

	
	int ipsi_trg = GA::GetPsiBin(psi);
	int ipsinear_trg = GA::GetPsiBinNear(psi);
	if (verbosity > 2)
	{
	    cout<<"ipsi_trg = "<<ipsi_trg<<" ipsinear_trg = "<<ipsinear_trg<<endl;
	}
	if (ipsi_trg < 0 || ipsinear_trg < 0) return;

	vector<int> ipsi_tree;
	if (dopsimix)
	{
	    if (verbosity > 2)
	    {
		cout<<"Doing psi mixing."<<endl;
		cout<<"ipsi_tree.push_back("<<ipsi_trg<<");"<<endl;
		cout<<"ipsi_tree.push_back("<<ipsinear_trg<<");"<<endl;
	    }

	    ipsi_tree.push_back(ipsi_trg);
	    ipsi_tree.push_back(ipsinear_trg);
	}
	else
	{
	    if (verbosity > 2)
	    {
		cout<<"Not doing psi mixing.  Adding all psibin to vector."<<endl;
	    }
	    for (int ipsi = 0; ipsi < GA::npsi_tree; ipsi++)
	    {
		if (verbosity > 2)
		{
		    cout<<"ipsi_tree.push_back("<<ipsi_trg<<");"<<endl;
		}
		ipsi_tree.push_back(ipsi);
	    }
	}
    
	//check if roughly enough events for mixing
	{
	    TTree* trktree = NULL;

	    trktree = (TTree*) infile->Get(Form("trktree_C%i_P%i_Z%i",icent_tree.front(),ipsi_trg,izv_trg));
	    if (!trktree)
	    {
		cerr<<"Line "<<__LINE__<<": Unable to find trktree ";
		fprintf(stderr,"trktree_C%i_P%i_Z%i",icent_tree.front(),ipsi_trg,izv_trg);
		cerr<<" in "<<inname <<endl;
		cerr<<"Exiting."<<endl;
		return;
	    }

	    unsigned int nevt = trktree->GetEntries();
	    h1nmixevt_raw->Fill(nevt);

	    h1nmixavail->Fill(nevt/(2*npsimixbin/GA::npsi_tree));
	    unsigned int nmixmin_tmp = nmixmin*2*npsimixbin/GA::npsi_tree;

	    if (nevt-1 < nmixmin_tmp)
	    {
		if (verbosity > 1)
		{
		    cout<<"There are "<<nevt<<" entries in track tree which does not meet the required minimum of "<<nmixmin_tmp+1<<endl;
		}
		continue;
	    }
	    else if (verbosity > 2)
	    {
		cout<<"There are "<<nevt<<" entries in track tree which meets the required minimum of "<<nmixmin_tmp+1<<endl;
	    }
	}


	iptt_ana = GetPtTrigBin(ptt);
	if (verbosity > 2)
	{
	    cout<<"iptt_ana = GetPtTrigBin("<<ptt<<") = "<<iptt_ana<<endl;
	}
	if (iptt_ana < 0) continue;

	idpsi_ana = GetDPsiBin(psi,trg->GetPhi());
	if (verbosity > 2)
	{
	    cout<<"idpsi_ana = GetDPsiBin("<<psi<<","<<trg->GetPhi()<<") = "<<idpsi_ana<<endl;
	}

	//--------------------------------------------------
	//passed any trigger cuts. lets fill histos
	//--------------------------------------------------
	h1cent_trg[trg->type]->Fill(cent);
	h1psi_trg[trg->type]->Fill(psi);
	h1zv_trg[trg->type]->Fill(zv);

	FillTrigHistos(trg);

	if (verbosity > 2)
	{
	    cout<<"h2cent_ptt["<<trg->type<<"]["<<idpsi_ana<<"]->Fill("<<cent<<","<<ptt<<")"<<endl;
	    cout<<"Filling "<<h2cent_ptt[trg->type][idpsi_ana]->GetName()<<endl;
	}
	h2cent_ptt[trg->type][idpsi_ana]->Fill(cent,ptt);
	if (verbosity > 2)
	{
	    cout<<"Filled "<<h2cent_ptt[trg->type][idpsi_ana]->GetName()<<endl;
	}

	// Get dir and dec weights
	double dirwt = -999.;
	double decwt = -999.;
	double dirwtY = -999.;
	double decwtY = -999.;

	if (trg->type != PI0)
	{
	    if (verbosity > 2)
	    {
		cout<<"Getting weights..."<<endl;
	    }
	    if(GetWeights(dirwt,decwt,dirwtY,decwtY,ftageff,fmistag,trg->type,cent,ptt) < 0)
	    {
		cerr<<"Error getting weights. Exiting."<<endl;
		return;
	    }
	    if (verbosity > 2)
	    {
		cout<<"Got weights dirwt = "<<dirwt<<"  decwt = "<<decwt<<"  dirwtY = "<<dirwtY<<"  decwtY = "<<decwtY<<endl;
	    }

	    h2cent_ptt_adj[DIR][idpsi_ana]->Fill(cent,ptt,dirwt);
	    h2cent_ptt_adj[DEC][idpsi_ana]->Fill(cent,ptt,decwt);

	    if (verbosity > 2)
	    {
		cout<<"Evaluating decay efficiency..."<<endl;
	    }
	    double deceff = -999;
	    deceff = fdeceff[ncent_fun-1]->Eval(ptt);
	    if (deceff < 0. || deceff > 2.)
	    {
		cerr<<"Decay efficiency evaluated to "<<deceff<<" is out of bounds."<<endl;
		cerr<<fdeceff[ncent_fun-1]->GetName()<<endl;
		cerr<<"Exiting."<<endl;
		return;
	    }
	    if (verbosity > 2)
	    {
		cout<<"Decay efficiency, decwt, evaluated to "<<decwt<<"."<<endl;
	    }

	    double decwt_eff = decwt/deceff;
	    h2cent_ptt_deceff[idpsi_ana]->Fill(cent,ptt,decwt_eff);

	}
	
	TH1F* h1dphi_tmp[2][npta];
	TH1F* h1dphi_xi_tmp[2][nxi];
	TH1F* h1dphiacc_xi_tmp[nxi];

	if (verbosity > 2)
	{
	    cout<<"Creating temporary dphi histograms..."<<endl;
	}

	for (int imix = 0; imix < 2; imix++)
	{
	    for (int ipta = 0; ipta < npta; ipta++)
	    {
		h1dphi_tmp[imix][ipta] = new TH1F(
		    Form("dphi_tmp_%i_%i_%i",itrig,imix,ipta),
		    Form("dphi_tmp_%i_%i_%i",itrig,imix,ipta),
		    ndphibin,mindphi,maxdphi);

		h1dphi_tmp[imix][ipta]->Sumw2();
	    }

	    for (int ixi = 0; ixi < nxi; ixi++)
	    {
		h1dphi_xi_tmp[imix][ixi] = new TH1F(
		    Form("dphi_xi_tmp_%i_%i_%i",itrig,imix,ixi),
		    Form("dphi_xi_tmp_%i_%i_%i",itrig,imix,ixi),
		    ndphibin,mindphi,maxdphi);

		h1dphi_xi_tmp[imix][ixi]->Sumw2();
	    }
	}

	for (int ixi = 0; ixi < nxi; ixi++)
	{
	    h1dphiacc_xi_tmp[ixi] = new TH1F(
		Form("dphiacc_xi_tmp_%i_%i",itrig,ixi),
		Form("dphiacc_xi_tmp_%i_%i",itrig,ixi),
		ndphibin,mindphi,maxdphi);

	    h1dphiacc_xi_tmp[ixi]->Sumw2();
	}


	if (verbosity > 2)
	{
	    cout<<"Created temporary dphi histograms!"<<endl;
	}

	int runno_tmp = -999;

	float ntrkfound = 0.0;
	float ngoodtrkfound = 0.0;

	unsigned int nmixevt = 0;

	// ========================================
	// Start getting track events
	// ========================================
	for (vector<int>::iterator it_zv = izv_tree.begin(); it_zv != izv_tree.end(); it_zv++)
	{
	    for (vector<int>::iterator it_psi = ipsi_tree.begin(); it_psi != ipsi_tree.end(); it_psi++)
	    {
		for (vector<int>::iterator it_cent = icent_tree.begin(); it_cent != icent_tree.end(); it_cent++)
		{
		    
		    if (verbosity > 1)
			cout<<"Tree indices: *it_zv = "<<*it_zv<<"; icent_tree = "<<*it_cent<<"; *it_psi = "<<*it_psi<<endl;

		    TTree* trktree = NULL;

		    trktree = (TTree*) infile->Get(Form("trktree_C%i_P%i_Z%i",*it_cent,*it_psi,*it_zv));
		    if (!trktree)
		    {
			cerr<<"Line "<<__LINE__<<": Unable to find trktree ";
			fprintf(stderr,"trktree_C%i_P%i_Z%i",*it_cent,*it_psi,*it_zv);
			cerr<<" in "<<inname <<endl;
			cerr<<"Exiting."<<endl;
			return;
		    }
		    else if (verbosity > 1)
			cout<<"Found trktree "<<trktree->GetName()<<endl;


		    TrackCollectionBase* trkhold = NULL;
		    TBranch* trkBranch = trktree->GetBranch("TrackBranch");
		    trkBranch->SetAddress(&trkhold);

		    int evtno_trk;
		    TBranch* evtBranch_trk = trktree->GetBranch("EventNumber");
		    evtBranch_trk->SetAddress(&evtno_trk);

		    int runno_trk;
		    TBranch* runBranch_trk = trktree->GetBranch("RunNumber");
		    runBranch_trk->SetAddress(&runno_trk);

		    float psi_trk;
		    TBranch* psiBranch_trk = trktree->GetBranch("Psi");
		    psiBranch_trk->SetAddress(&psi_trk);

		    float zv_trk;
		    TBranch* z0Branch_trk = trktree->GetBranch("Zvertex");
		    z0Branch_trk->SetAddress(&zv_trk);

		    float cent_trk;
		    TBranch* centBranch_trk = trktree->GetBranch("Centrality");
		    centBranch_trk->SetAddress(&cent_trk);

		    unsigned int nevt = trktree->GetEntries();
		
		    if (verbosity > 1)
			cout<<"Total events in trktree = "<<nevt<<endl;

		    unsigned int nmixevt_curtree = 0;

		    for (unsigned int ievt = 0; ievt < nevt; ievt++)
		    {
			trktree->GetEntry(ievt);

			if (icent_ana != GetCentBin(cent_trk))
			{
			    if (verbosity > 3)
			    {
				cout<<"MISMATCH: Centrality of potential mixing event, "<<cent_trk<<", is not in same bin as the triggered event centrality, "<<cent<<".  Continue."<<endl;
			    }
			    delete trkhold; trkhold = NULL;
			    continue;
			}
			if (verbosity > 2)
			{
			    cout<<"MATCH: Centrality of potential mixing event, "<<cent_trk<<", in same bin as triggered event centrality, "<<cent<<"!"<<endl;
			}

			//check if within psi range
			float psidiff = psi-psi_trk;

			while (psidiff < -PI/2.)
			    psidiff+=PI;

			while (psidiff >= PI/2.)
			    psidiff-=PI;

			if (dopsimix)
			{

			    if (fabs(psidiff) > (psimixbw*PI)/2.)
			    {
				if (verbosity > 2)
				{
				    cout<<"MISMATCH: psi track event ("<<psi_trk<<") - psi trigger event ("<<psi<<") is greater than threshold ("<<psimixbw/2.<<"). Continue."<<endl;
				}
				delete trkhold; trkhold = NULL;
				continue;
			    }

			    if (verbosity > 2)
			    {
				cout<<"MATCH: psi track event ("<<psi_trk<<") - psi trigger event ("<<psi<<") is within threshold ("<<psimixbw/2.<<")!"<<endl;
			    }
			}

			if (verbosity > 1)
			{
			    if (runno_trk != runno_tmp)
			    {
				cout<<"Looping through track events from Run #"<<runno_trk<<"; Events remaining: "<<nevt-ievt<<endl;
				runno_tmp = runno_trk;
			    }
			}

			bool same_evt = runno == runno_trk && evtno == evtno_trk;

			if (!same_evt)
			{
			    if (!trkhold->IsMixable())
			    {
				nomixevtno.push_back(evtno_trk);

				if (verbosity > 1)
				    cout<<"Event: "<<runno_trk<<"."<<evtno_trk<<" was saved for another trigger.  Skip."<<endl;
				delete trkhold; trkhold = NULL;
				continue;
			    }
			    nmixevt++;
			    nmixevt_curtree++;
			}
			else if (verbosity > 1)
			{
			    cout<<"Same event found! Trig Run #"<<runno<<" & Event #"<<evtno<<" = Track Run #"<<runno_trk<<" & Event #"<<evtno_trk<<endl;
			}

			//loop through tracks in TrackCollectionBase trkhold
			unsigned int ntrk = trkhold->GetEntries();
			if (verbosity > 2)
			    cout<<"Looping through "<<ntrk<<" tracks from counted event "<<ievt<<" -> Event #"<<evtno_trk<<endl;
			h2cent_ntrk_trkevt[0]->Fill(cent_trk,ntrk);
			ntrkfound+=ntrk;

			unsigned int ngoodtrk = 0;

			ChargedTrack* trk = NULL;
			for (unsigned int itrk = 0; itrk < ntrk; itrk++)
			{
			    if (itrk % 100 == 0 && verbosity > 2)
			    {
				cout<<"Analyzing track "<<itrk<<"/"<<ntrk-1<<endl;
			    }
			
			    bool trkfound = trkhold->GetEntry(itrk,&trk);
			    if (!trkfound)
			    {
				cerr<<"Could not find track in trkhold!  Exiting."<<endl;
				return;
			    }


			    FillTrkHistos(trk,0,*it_cent);

			    if (!IsTrackGood(trk)) {
				trk = NULL;
				continue;
			    }

			    FillTrkHistos(trk,1,*it_cent);

			    ngoodtrk++;

			    FillCorrels(trg,trk,same_evt,h1dphi_tmp,h1dphi_xi_tmp,h1dphiacc_xi_tmp);
			    trk = NULL;
			
			}//closes track loop

			if (verbosity > 2)
			    cout<<"Found "<<ngoodtrk<<" good tracks out of the possible "<<ntrk<<"."<<endl;

			h2cent_ntrk_trkevt[1]->Fill(cent_trk,ngoodtrk);
			ngoodtrkfound+=ngoodtrk;

			delete trkhold; trkhold = NULL;
		    
		    }//closes event loop (or track tree loop)

		    if (verbosity > 1)
			cout<<"Mixed trigger with "<<nmixevt_curtree<<" events from this track tree."<<endl;
		
		}//end cent loop
	    }//end psi bin loop
	}//end zv bin loop

	if (verbosity > 1)
	    cout<<"Done analyzing events for this trigger. Analyzed "<<nmixevt<<" mixed events.  "<<ntrkfound<<" tracks and "<<ngoodtrkfound<<" good tracks found."<<endl;

	h1nmixevt->Fill(nmixevt);
	

	for (int ipta = 0; ipta < npta; ipta++)
	{
	    double tmpscl = 1./((double)nmixevt);
	    h1dphi_tmp[1][ipta]->Scale(tmpscl);

	    for (int imix = 0; imix < 2; imix++)
	    {

		h1dphi[trg->type][imix][iptt_ana][icent_ana][idpsi_ana][ipta]->Add(h1dphi_tmp[imix][ipta]);

		if (trg->type != PI0)
		{
		    h1dphi_adj[DIR][imix][iptt_ana][icent_ana][idpsi_ana][ipta]->Add(h1dphi_tmp[imix][ipta],dirwtY);
		    h1dphi_adj[DEC][imix][iptt_ana][icent_ana][idpsi_ana][ipta]->Add(h1dphi_tmp[imix][ipta],decwtY);
		}

		delete h1dphi_tmp[imix][ipta];
	    }
	}

	for (int ixi = 0; ixi < nxi; ixi++)
	{
	    double tmpscl = 1./((double)nmixevt);
	    h1dphi_xi_tmp[1][ixi]->Scale(tmpscl);
	    h1dphiacc_xi_tmp[ixi]->Scale(tmpscl);

	    for (int imix = 0; imix < 2; imix++)
	    {
		h1dphi_xi[trg->type][imix][iptt_ana][icent_ana][idpsi_ana][ixi]->Add(h1dphi_xi_tmp[imix][ixi]);

		if (trg->type != PI0)
		{
		    h1dphi_xi_adj[DIR][imix][iptt_ana][icent_ana][idpsi_ana][ixi]->Add(h1dphi_xi_tmp[imix][ixi],dirwtY);
		    h1dphi_xi_adj[DEC][imix][iptt_ana][icent_ana][idpsi_ana][ixi]->Add(h1dphi_xi_tmp[imix][ixi],decwtY);
			
		}

		delete h1dphi_xi_tmp[imix][ixi];
	    }


	    h1dphiacc_xi[trg->type][iptt_ana][icent_ana][idpsi_ana][ixi]->Add(h1dphiacc_xi_tmp[ixi]);

	    if (trg->type != PI0)
	    {
		h1dphiacc_xi_adj[DIR][iptt_ana][icent_ana][idpsi_ana][ixi]->Add(h1dphiacc_xi_tmp[ixi],dirwtY);
		h1dphiacc_xi_adj[DEC][iptt_ana][icent_ana][idpsi_ana][ixi]->Add(h1dphiacc_xi_tmp[ixi],decwtY);
			
	    }

	    delete h1dphiacc_xi_tmp[ixi];
	}

	delete trg; trg = NULL;
    }//closes trigger loop

    return End();
}

void End()
{
    for (vector<int>::iterator it_nomix = nomixevtno.begin(); it_nomix != nomixevtno.end(); it_nomix++)
    {
	bool match = false;
	for (vector<int>::iterator it_trig = trigevtno.begin(); it_trig != trigevtno.end(); it_trig++)
	{
	    if (*it_nomix == *it_trig) match = true;
	}
	if (!match) cerr<<"Triggered event "<<*it_nomix<<" was in the mix pool."<<endl;
	else if (verbosity > 2)
	    cout<<"The TrackCollectionBase for event "<<*it_nomix<<" was correctly marked as nonmixable."<<endl;
    }

    if (verbosity > 0)
	cout<<"Closing input file..."<<endl;
    infile->Close();

    if (verbosity > 0)
	cout<<"Writing output file "<<outfile->GetName()<<endl;
    outfile->Write();

    if (verbosity > 0)
	cout<<"Closing output file..."<<endl;
    outfile->Close();

    if (verbosity > 0)
	cout<<"Successful completion!"<<endl;

    return;
}

int ReadCutVars()
{
    ifstream parfile("fitpar.txt", ifstream::in);
    if (!parfile)
    {
	cerr<<"Could not open fitpar.txt.  Exiting"<<endl;
	return -1;
    }

    string line;
    while (getline(parfile,line))
    {
	istringstream iss_line(line.data());
	string scut;
	int icent, ipt;
	float mean_tmp, sigma_tmp;
	iss_line>>scut>>icent>>ipt>>mean_tmp>>sigma_tmp;
	int icut = cutmap.find(scut)->second;

	sigma_trkcut[icut][icent][ipt] = sigma_tmp;
	mean_trkcut[icut][icent][ipt] = mean_tmp;

	nsig_trkcut[icut][icent][ipt] = (icut == DCAZ || icut == DCA2D) ? 6. : 3.;
    }

    return 1;

}

int GetCentBin_TrkCut(int centrality)
{
    for (int icent = 0; icent < ncent_trkcut; icent++)
    {
	if (centrality >= limcent_trkcut[icent] && centrality < limcent_trkcut[icent+1])
	    return icent;
    }
    return -999;
}

int GetPtBin_TrkCut(float pt)
{
    for (int ipt = 0; ipt < npt_trkcut; ipt++)
    {
	if (pt >= limpt_trkcut[ipt] && pt < limpt_trkcut[ipt+1])
	    return ipt;
    }
    return -999;
}

void TrackEventAnalysis()
{
    if (verbosity > 0)
	cout<<"Beginning Track Event Analysis."<<endl;

    for (int icent_tree = 0; icent_tree < GA::ncent_tree; icent_tree++)
    {
	for (int ipsi_tree = 0; ipsi_tree < GA::npsi_tree; ipsi_tree++)
	{
	    for (int izv_tree = 0; izv_tree < GA::nzv_tree; izv_tree++)
	    {
		TTree* trktree = NULL;
		trktree = (TTree*) infile->Get(Form("trktree_C%i_P%i_Z%i",icent_tree,ipsi_tree,izv_tree));

		if (!trktree)
		{
		    cerr<<"Line "<<__LINE__<<": Unable to find 'trktree"<<endl;
		    cerr<<"Exiting."<<endl;
		    return;
		}

		else if (verbosity > 1)
		    cout<<"Found trktree "<<trktree->GetName()<<endl;

		TrackCollectionBase* trkhold = NULL;
		TBranch* trkBranch = trktree->GetBranch("TrackBranch");
		trkBranch->SetAddress(&trkhold);

		float zv_trk;
		TBranch* z0Branch_trk = trktree->GetBranch("Zvertex");
		z0Branch_trk->SetAddress(&zv_trk);

		float cent_trk;
		TBranch* centBranch_trk = trktree->GetBranch("Centrality");
		centBranch_trk->SetAddress(&cent_trk);

		unsigned int nevt = trktree->GetEntries();

		if (verbosity > 1)
		    cout<<"Total events in trktree = "<<nevt<<endl;

		h1binsize_trkevt->Fill(nevt);

		for (unsigned int ievt = 0; ievt < nevt; ievt++)
		{
		    trktree->GetEntry(ievt);

		    h2cent_zv_trkevt->Fill(cent_trk,zv_trk);

		    unsigned int ntrk = trkhold->GetEntries();
		    h2cent_ntrk_trkevt[0]->Fill(cent_trk,ntrk);

		    delete trkhold; trkhold = NULL;
		}
	    }
	}
    }
}

void TrackAnalysis()
{
    if (verbosity > 0)
	cout<<"Beginning Track Analysis."<<endl;

    for (int icent_tree = 0; icent_tree < GA::ncent_tree; icent_tree++)
    {
	for (int ipsi_tree = 0; ipsi_tree < GA::npsi_tree; ipsi_tree++)
	{
	    for (int izv_tree = 0; izv_tree < GA::nzv_tree; izv_tree++)
	    {
		TTree* trktree = NULL;
		trktree = (TTree*) infile->Get(Form("trktree_C%i_P%i_Z%i",icent_tree,ipsi_tree,izv_tree));

		if (!trktree)
		{
		    cerr<<"Line "<<__LINE__<<": Unable to find 'trktree"<<endl;
		    cerr<<"Exiting."<<endl;
		    return;
		}

		else if (verbosity > 1)
		    cout<<"Found trktree "<<trktree->GetName()<<endl;

		TrackCollectionBase* trkhold = NULL;
		TBranch* trkBranch = trktree->GetBranch("TrackBranch");
		trkBranch->SetAddress(&trkhold);

		int runno_trk;
		TBranch* runBranch_trk = trktree->GetBranch("RunNumber");
		runBranch_trk->SetAddress(&runno_trk);

		int evtno_trk;
		TBranch* evtBranch_trk = trktree->GetBranch("EventNumber");
		evtBranch_trk->SetAddress(&evtno_trk);

		float psi_trk;
		TBranch* psiBranch_trk = trktree->GetBranch("Psi");
		psiBranch_trk->SetAddress(&psi_trk);

		float zv_trk;
		TBranch* z0Branch_trk = trktree->GetBranch("Zvertex");
		z0Branch_trk->SetAddress(&zv_trk);

		float cent_trk;
		TBranch* centBranch_trk = trktree->GetBranch("Centrality");
		centBranch_trk->SetAddress(&cent_trk);

		unsigned int nevt = trktree->GetEntries();

		if (verbosity > 1)
		    cout<<"Total events in trktree = "<<nevt<<endl;

		h1binsize_trkevt->Fill(nevt);

		for (unsigned int ievt = 0; ievt < nevt; ievt++)
		{
		    trktree->GetEntry(ievt);

		    h2cent_zv_trkevt->Fill(cent_trk,zv_trk);

		    unsigned int ntrk = trkhold->GetEntries();
		    h2cent_ntrk_trkevt[0]->Fill(cent_trk,ntrk);

		    unsigned int ngoodtrk = 0;

		    ChargedTrack* trk = NULL;
		    for (unsigned int itrk = 0; itrk < ntrk; itrk++)
		    {
			bool trkfound = trkhold->GetEntry(itrk,&trk);
			if (!trkfound)
			{
			    cerr<<"Could not find trk in trkhold!  Exiting."<<endl;
			    return;
			}

			FillTrkHistos(trk,0,icent_tree);

			if (!IsTrackGood(trk)) {
			    trk = NULL;
			    continue;
			}
			ngoodtrk++;
			
			FillTrkHistos(trk,1,icent_tree);

			h2cent_pta->Fill(cent_trk,trk->GetPt());

			trk = NULL;
		    }

		    h2cent_ntrk_trkevt[1]->Fill(cent_trk,ngoodtrk);

		    delete trkhold; trkhold = NULL;
		}
	    }

	}
    }

    if (verbosity > 0)
	cout<<"Track Analysis Complete!"<<endl;
}


bool FillCorrels(TriggerParticle* trg, ChargedTrack* trk, bool same_evt, TH1F* h1dphi_ptr[][npta], TH1F* h1dphi_xi_ptr[][nxi], TH1F* h1dphiacc_xi_ptr[nxi])
{
    if (verbosity > 4)
    {
	cout<<"Filling correlation functions"<<endl;
    }
    bool ismix = !same_evt;

    float dphi = trk->GetPhi() - trg->GetPhi();
    h1dphi_test_raw->Fill(dphi/PI);

    while (dphi < 0.)
    	dphi+=2.*PI;

    while (dphi >= 2.*PI)
    	dphi-=2.*PI;

    if (dphi > PI)
    	dphi = 2.*PI - dphi;


    h1dphi_test->Fill(dphi/PI);

    float pt_trg = trg->GetPt();
    float pt_trk = trk->GetPt();

    int ipta = GetPtAssocBin(pt_trk);
    if (ipta >= 0)
	h1dphi_ptr[ismix][ipta]->Fill(dphi/PI);

    int ixi = GetXiBin(pt_trg,pt_trk);
    if (ixi >= 0)
    {
    	h1dphi_xi_ptr[ismix][ixi]->Fill(dphi/PI);
    }

    //fill acc histos if pta <= 3
    if (ismix && pt_trk <= ptamixmax1)
    {
	int ixi = GetXiBin(pt_trg,pt_trk);
	if (ixi >= 0)
	{
	    h1dphiacc_xi_ptr[ixi]->Fill(dphi/PI);
	}
	if (pt_trk > ptamixmax0)
	{
	    while (--ixi >= 0)
	    {
		h1dphiacc_xi_ptr[ixi]->Fill(dphi/PI);
	    }
	}
    }

    return true;
}

bool IsTrackGood(ChargedTrack* trk)
{
    float trkpt = trk->GetPt();

    if (trkpt < GA::minpta) return false;
    if (trkpt > GA::maxpta) return false;    

    GA::tracktype_t track_t = trk->GetTrackType();

    CentralTrackData* mycnt = NULL;

    if (track_t == GA::CNT)
    {
	mycnt = trk->GetCentralTrackData();
    }
    else if (track_t == GA::SVXCNT)
    {
	cerr<<"Svx Tracks no longer supported!"<<endl;
	return false;
    }
    else
    {
    	cerr<<"Track type is unknown!!!"<<endl;
    	return false;
    }

    if (!mycnt->EmcMatch(3.0))
    {
	return false;
    }

    if (!mycnt->Pc3Match(3.0))
    {
	return false;
    }

    if (trkpt < 5. && !mycnt->isn0neg)
    {
	return false;
    }

    if (!mycnt->CheckDEP(-2.0))
    {
	return false;
    }
    
    return true;

}

void InitHistos()
{
    if (verbosity > 0)
	cout<<"Initializing histograms"<<endl;
    
    h1nmixavail = new TH1F("nmixavail","Mixed Events Available for Trigger (Estimated);Mixed Events;Triggers",1000,0,10000.);

    h1nmixevt_raw = new TH1F("nmixevt_raw","Mixed Events (In Same Bin);Mixed Events;N Triggers",1000,0,1000000.);

    h1nmixevt = new TH1F("nmixevt","Mixed Events;Mixed Events;N Triggers",10000,0,10000.);


    h2cent_zv_trkevt = new TH2F("cent_zv_trkevt","#splitline{Saved Track Events}{cent_zv_trkevt};centrality;z vertex",101,-0.5,100.5,2*GA::maxzv,-GA::maxzv,GA::maxzv);


    h1binsize_trkevt = new TH1F("binsize_trkevt","Saved Track Events in Tree Bins;N Bin;N Track Events",1000,0,1000);

    h2cent_pta = new TH2F("cent_pta","Good Tracks;centrality;p_{T}",
			  GA::maxcent,0.5,GA::maxcent+0.5,
			  (GA::maxpta-GA::minpta)*10.,GA::minpta,GA::maxpta);

    for (int itype = 0; itype < ntype_trk; itype++)
    {
	h2cent_ntrk_trkevt[itype] = new TH2F(Form("cent_ntrk_type%i_trkevt",itype),
					     Form("#splitline{Saved Track Events}{cent_ntrk_type%i_trkevt};centrality;tracks",itype),
					     101,-0.5,100.5,400,0,800);

	h1phi0[itype] = new TH1F(Form("phi0_type%i",itype),
				 Form("phi0_type%i;phi0",itype),
				 nphibin,minphi,maxphi);

	h1eta[itype] = new TH1F(Form("eta_type%i",itype),
				Form("eta_type%i;eta",itype),
				netabin,mineta,maxeta);

	h1px[itype] = new TH1F(Form("px_type%i",itype),
			       Form("px_type%i",itype),
			       2*npbin,-maxp,maxp);

	h1py[itype] = new TH1F(Form("py_type%i",itype),
			       Form("py_type%i",itype),
			       2*npbin,-maxp,maxp);

	h1pz[itype] = new TH1F(Form("pz_type%i",itype),
			       Form("pz_type%i",itype),
			       npzbin,-maxpz,maxpz);

	h1pt[itype] = new TH1F(Form("pt_type%i",itype),
			       Form("pt_type%i;pt",itype),
			       npbin,0.,maxp);

	h14bitpat[itype] = new TH1F(Form("4bitpat_type%i",itype),
				    Form("4bitpat_type%i",itype),
				    20,0.,20.);

	h1zed[itype] = new TH1F(Form("zed_type%i",itype),
				Form("zed_type%i",itype),
				201,-100.5,100.5);

	h1n0[itype] = new TH1F(Form("n0_type%i",itype),
			       Form("track n0_type%i;n0",itype),
			       100,-50.,50);

	h1ecore[itype] = new TH1F(Form("ecore_type%i",itype),
				  Form("track ecore_type%i;ecore",itype),
				  npbin,-maxp,maxp);

	h1chi2ndf[itype] = new TH1F(Form("chi2ndf_type%i",itype),
				    Form("chi2/ndf_type%i",itype),
				    25,0.,10.);

	h1pc3dphi[itype] = new TH1F(Form("pc3dphi_type%i",itype),
				    Form("pc3dphi_type%i",itype),
				    nphimatchbin,minphimatchbin,maxphimatchbin);
	h1pc3dz[itype] = new TH1F(Form("pc3dz_type%i",itype),
				  Form("pc3dz_type%i",itype),
				  nzmatchbin,minzmatchbin,maxzmatchbin);
	h1pc3sdphi[itype] = new TH1F(Form("pc3sdphi_type%i",itype),
				     Form("pc3sdphi_type%i",itype),
				     nsmatchbin,minsmatchbin,maxsmatchbin);
	h1pc3sdz[itype] = new TH1F(Form("pc3sdz_type%i",itype),
				   Form("pc3sdz_type%i",itype),
				   nsmatchbin,minsmatchbin,maxsmatchbin);

	h1pc3sig[itype] = new TH1F(Form("pc3sig_type%i",itype),
				   Form("pc3sig_type%i;pc3dr*%.0f",itype,GA::CFF),
				   (int)GA::MAXUC + 1,0,(int)GA::MAXUC + 1);

	h1emcdphi[itype] = new TH1F(Form("emcdphi_type%i",itype),
				    Form("emcdphi_type%i",itype),
				    nphimatchbin,minphimatchbin,maxphimatchbin);
	h1emcdz[itype] = new TH1F(Form("emcdz_type%i",itype),
				  Form("emcdz_type%i",itype),
				  nzmatchbin,minzmatchbin,maxzmatchbin);
	h1emcsdphi[itype] = new TH1F(Form("emcsdphi_type%i",itype),
				     Form("emcsdphi_type%i",itype),
				     nsmatchbin,minsmatchbin,maxsmatchbin);
	h1emcsdz[itype] = new TH1F(Form("emcsdz_type%i",itype),
				   Form("emcsdz_type%i",itype),
				   nsmatchbin,minsmatchbin,maxsmatchbin);

	h1emcsig[itype] = new TH1F(Form("emcsig_type%i",itype),
				   Form("emcsig_type%i;emcdr*%.0f",itype,GA::CFF),
				   (int)GA::MAXUC + 1,0,(int)GA::MAXUC + 1);

	h1dep[itype] = new TH1F(Form("dep_type%i",itype),
				Form("dep_type%i;dep*%.0f",itype,GA::CFF),
				(int)GA::MAXC - (int)GA::MINC + 1,(int)GA::MINC,(int)GA::MAXC + 1);

	h1dca2d[itype] = new TH1F(Form("dca2d_type%i",itype),
				  Form("dca2d_type%i",itype),
				  ndcabin,mindcabin,maxdcabin);

	h1dcaz[itype] = new TH1F(Form("dcaz_type%i",itype),
				 Form("dcaz_type%i",itype),
				 ndcabin,mindcabin,maxdcabin);

    }


    for (int isec = 0; isec < GA::nsec; isec++)
    {
	int ny = -999;
	int nz = -999;
	GA::GetNzNy(nz,ny,isec);

	h2zypos[isec] =
	    new TH2F(Form("zypos_sec%i",isec),
		     Form("EMC isec = %i;izpos;iypos",isec),
		     nz,0.,(float)nz,ny,0.,(float)ny);
    }

    for (int itype = 0; itype < ntype; ++itype)
    {
	

	h1cent_trg[itype]  = new TH1F(Form("cent_trg_type%i",itype),
				      Form("cent_trg_type%i",itype),
				      GA::maxcent,0.5,GA::maxcent+0.5);

	h1zv_trg[itype]  = new TH1F(Form("zv_trg_type%i",itype),
				    Form("zv_trg_type%i",itype),
				    20*GA::maxzv,-GA::maxzv,GA::maxzv);

	h1psi_trg[itype]  = new TH1F(Form("psi_trg_type%i",itype),
				     Form("psi_trg_type%i",itype),
				     60,0.,PI);
	
	h2etaphi_trg[itype] = new TH2F(Form("etaphi_trg_type%i",itype),
				       Form("etaphi_trg_type%i",itype),
				       netabin,mineta,maxeta,
				       nphibin,minphi,maxphi);

	h1px_trg[itype] = new TH1F(Form("px_trg_type%i",itype),
				   Form("px_trg_type%i",itype),
				   2*npbin,-maxp,maxp);

	h1py_trg[itype] = new TH1F(Form("py_trg_type%i",itype),
				   Form("py_trg_type%i",itype),
				   2*npbin,-maxp,maxp);

	h1pz_trg[itype] = new TH1F(Form("pz_trg_type%i",itype),
				   Form("pz_trg_type%i",itype),
				   npzbin,-maxpz,maxpz);


	h1pc3dphis_trg[itype] = new TH1F(Form("pc3dphis_trg_type%i",itype),
					 Form("pc3dphis_trg_type%i",itype),
					 200,-100.,100.);
	h1pc3dz_trg[itype] = new TH1F(Form("pc3dz_trg_type%i",itype),
				      Form("pc3dz_trg_type%i",itype),
				      200,-100.,100.);
	h1pc3dr_trg[itype] = new TH1F(Form("pc3dr_trg_type%i",itype),
				      Form("pc3dr_trg_type%i",itype),
				      200,0.,100.);
	
	h1trkdphis_trg[itype] = new TH1F(Form("trkdphis_trg_type%i",itype),
					 Form("trkdphis_trg_type%i",itype),
					 200,-100.,100.);
	h1trkdz_trg[itype] = new TH1F(Form("trkdz_trg_type%i",itype),
				      Form("trkdz_trg_type%i",itype),
				      200,-100.,100.);
	h1trkdr_trg[itype] = new TH1F(Form("trkdr_trg_type%i",itype),
				      Form("trkdr_trg_type%i",itype),
				      200,0.,100.);

	if (verbosity > 1)
	{
	    cout<<"Creating h2cent_ptt histograms..."<<endl;
	}
	for (int idpsi = 0; idpsi < ndpsi; idpsi++)
	{
	    if (verbosity > 2)
	    {
		cout<<"Creating h2cent_ptt["<<itype<<"]["<<idpsi<<"]"<<endl;
	    }
	    h2cent_ptt[itype][idpsi]  = new TH2F(Form("cent_ptt_%s_dpsi%i",spec[itype],idpsi),
						     Form("cent_ptt_%s_dpsi%i;centrality;pt",spec[itype],idpsi),
						     GA::maxcent,0.5,GA::maxcent+0.5,
						     (limptt[nptt]-limptt[0])*10.,limptt[0],limptt[nptt]);
	    h2cent_ptt[itype][idpsi]->Sumw2();
	    if (verbosity > 2)
	    {
		cout<<"Created "<<h2cent_ptt[itype][idpsi]->GetName()<<endl;
	    }

	    if (itype == DIR || itype == DEC)
	    {
		if (verbosity > 2)
		{
		    cout<<"Creating h2cent_ptt_adj["<<itype<<"]["<<idpsi<<"]"<<endl;
		}
		h2cent_ptt_adj[itype][idpsi]  = new TH2F(Form("cent_ptt_adj_%s_dpsi%i",spec[itype],idpsi),
							     Form("cent_ptt_adj_%s_dpsi%i;centrality;pt",spec[itype],idpsi),
							     GA::maxcent,0.5,GA::maxcent+0.5,
							     (limptt[nptt]-limptt[0])*10.,limptt[0],limptt[nptt]);
		h2cent_ptt_adj[itype][idpsi]->Sumw2();
		if (verbosity > 2)
		{
		    cout<<"Created "<<h2cent_ptt_adj[itype][idpsi]->GetName()<<endl;
		}

		if (itype == DEC)
		{
		    if (verbosity > 2)
		    {
			cout<<"Creating h2cent_ptt_deceff["<<idpsi<<"]"<<endl;
		    }
		    h2cent_ptt_deceff[idpsi]  = new TH2F(Form("cent_ptt_deceff_dpsi%i",idpsi),
							     Form("cent_ptt_deceff_dpsi%i;centrality;pt",idpsi),
							     GA::maxcent,0.5,GA::maxcent+0.5,
							     (limptt[nptt]-limptt[0])*10.,limptt[0],limptt[nptt]);
		    h2cent_ptt_deceff[idpsi]->Sumw2();
		    if (verbosity > 2)
		    {
			cout<<"Created "<<h2cent_ptt_deceff[idpsi]->GetName()<<endl;
		    }
		}

	    }
	


	    for (int icent = 0; icent < ncent; icent++)
	    {

		for (int iptt = 0; iptt < nptt; iptt++)
		{

		    for (int imix = 0; imix < 2; imix++)
		    {

			for (int ixi = 0; ixi < nxi; ixi++)
			{

			    h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi] = new TH1F (
				Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],imix,iptt,icent,idpsi,ixi),
				Form("#splitline{dphi_%s_mix%i}{Ptt: [%.2f,%.2f)  Centrality: %i-%i%%  d#Psi: %s  Xi: [%.2f,%.2f)};#Delta#phi;((1/N_{trig}) (dN_{ch}/d#Delta#phi)",spec[itype],imix,limptt[iptt],limptt[iptt+1],limcent[icent],limcent[icent+1],(dpsistr[idpsi]),limxi[ixi],limxi[ixi+1]),
				ndphibin,mindphi,maxdphi);
			    h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi]->Sumw2();			    
			    
			    if (itype == DIR || itype == DEC)
			    {

			    	h1dphi_xi_adj[itype][imix][iptt][icent][idpsi][ixi] = new TH1F (
			    	    Form("dphi_adj_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],imix,iptt,icent,idpsi,ixi),
			    	    Form("#splitline{dphi_adj_%s_mix%i}{Ptt: [%.2f,%.2f)  Centrality: %i-%i%%  d#Psi: %s  Xi: [%.2f,%.2f)};#Delta#phi;((1/N_{trig}) (dN_{ch}/d#Delta#phi)",spec[itype],imix,limptt[iptt],limptt[iptt+1],limcent[icent],limcent[icent+1],(dpsistr[idpsi]),limxi[ixi],limxi[ixi+1]),
			    	    ndphibin,mindphi,maxdphi);

			    	h1dphi_xi_adj[itype][imix][iptt][icent][idpsi][ixi]->Sumw2();

			    }
			}

			for (int ipta = 0; ipta < npta; ipta++)
			{
			    
			    h1dphi[itype][imix][iptt][icent][idpsi][ipta] = new TH1F (
				Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",spec[itype],imix,iptt,icent,idpsi,ipta),
				Form("#splitline{dphi_%s_mix%i}{Ptt: [%.2f,%.2f)  Centrality: %i-%i%%  d#Psi: %s  Pta: [%.2f,%.2f)};#Delta#phi;((1/N_{trig}) (dN_{ch}/d#Delta#phi)",spec[itype],imix,limptt[iptt],limptt[iptt+1],limcent[icent],limcent[icent+1],(dpsistr[idpsi]),limpta[ipta],limpta[ipta+1]),
				ndphibin,mindphi,maxdphi);
			    h1dphi[itype][imix][iptt][icent][idpsi][ipta]->Sumw2();

			    if (itype == DIR || itype == DEC)
			    {
				h1dphi_adj[itype][imix][iptt][icent][idpsi][ipta] = new TH1F (
				    Form("dphi_adj_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",spec[itype],imix,iptt,icent,idpsi,ipta),
				    Form("#splitline{dphi_adj_%s_mix%i}{Ptt: [%.2f,%.2f)  Centrality: %i-%i%%  d#Psi: %s  Pta: [%.2f,%.2f)};#Delta#phi;((1/N_{trig}) (dN_{ch}/d#Delta#phi)",spec[itype],imix,limptt[iptt],limptt[iptt+1],limcent[icent],limcent[icent+1],(dpsistr[idpsi]),limpta[ipta],limpta[ipta+1]),
				    ndphibin,mindphi,maxdphi);
				h1dphi_adj[itype][imix][iptt][icent][idpsi][ipta]->Sumw2();
			    }

			}
			
		    }

		    for (int ixi = 0; ixi < nxi; ixi++)
		    {

			h1dphiacc_xi[itype][iptt][icent][idpsi][ixi] = new TH1F (
			    Form("dphiacc_%s_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],iptt,icent,idpsi,ixi),
			    Form("#splitline{dphiacc_%s}{Ptt: [%.2f,%.2f)  Centrality: %i-%i%%  d#Psi: %s  Xi: [%.2f,%.2f)};#Delta#phi;((1/N_{trig}) (dN_{ch}/d#Delta#phi)",spec[itype],limptt[iptt],limptt[iptt+1],limcent[icent],limcent[icent+1],(dpsistr[idpsi]),limxi[ixi],limxi[ixi+1]),
			    ndphibin,mindphi,maxdphi);
			h1dphiacc_xi[itype][iptt][icent][idpsi][ixi]->Sumw2();			    
			    
			if (itype == DIR || itype == DEC)
			{
				
			    h1dphiacc_xi_adj[itype][iptt][icent][idpsi][ixi] = new TH1F (
				Form("dphiacc_adj_%s_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],iptt,icent,idpsi,ixi),
				Form("#splitline{dphiacc_adj_%s}{Ptt: [%.2f,%.2f)  Centrality: %i-%i%%  d#Psi: %s  Xi: [%.2f,%.2f)};#Delta#phi;((1/N_{trig}) (dN_{ch}/d#Delta#phi)",spec[itype],limptt[iptt],limptt[iptt+1],limcent[icent],limcent[icent+1],(dpsistr[idpsi]),limxi[ixi],limxi[ixi+1]),
				ndphibin,mindphi,maxdphi);

			    h1dphiacc_xi_adj[itype][iptt][icent][idpsi][ixi]->Sumw2();

			}
		    }
		    
		}
		
	    }
	    
	}
	
    }

    h1pc3sdz_trg = new TH1F("pc3sdz_trg","Photon Sigmalized pc3dz;pc3sdz",200,-10.,10);
    h1pc3sdphis_trg = new TH1F("pc3sdphis_trg","Photon Sigmalized pc3dphis;pc3sdphis",200,-10.,10);
    h1pc3sdr_trg = new TH1F("pc3sdr_trg","Photon Sigmalized pc3dr;pc3sdr",200,-10.,10);

    h1dphi_test = new TH1F("dphi_all","dphi_all",8*ndphibin,-4.,4.);

    h1dphi_test_raw = new TH1F("dphi_all_raw","dphi_all_raw",8*ndphibin,-4.,4.);

    if (verbosity > 0)
	cout<<"Histograms initialized successfully"<<endl;

    return;
}


void FillTrkHistos(ChargedTrack* trk, int type, int icent)
{
    if (type != 0 && type != 1)
    {
	cerr<<"array index type is out of range: "<<type<<endl;
	return;
    }

    float phi = trk->GetPhi()/PI;
    while (phi < minphi)
	phi+=2.;
    while (phi >= maxphi)
	phi-=2.;

    h1phi0[type]->Fill(phi);
    h1eta[type]->Fill(trk->GetEta());
    h1pt[type]->Fill(trk->GetPt());
    h1px[type]->Fill(trk->px);
    h1py[type]->Fill(trk->py);
    h1pz[type]->Fill(trk->pz);

    CentralTrackData* mycnt = trk->GetCentralTrackData();
    if (mycnt)
    {
	h1pc3sig[type]->Fill(mycnt->GetPc3SigChar());
	h1emcsig[type]->Fill(mycnt->GetEmcSigChar());

	h1dep[type]->Fill(mycnt->GetDEPChar());
    }
    else if (verbosity > 1)
    {
	cerr<<__LINE__<<": could not find CentralTrackData"<<endl;
    }

    return;
}

void FillTrigHistos(TriggerParticle* trg)
{
    if (verbosity > 2)
	cout<<"Filling trigger histograms"<<endl;
    
    float phi = trg->GetPhi()/PI;
    while (phi < minphi)
	phi+=2.;

    while (phi >= maxphi)
	phi-=2.;

    h2etaphi_trg[trg->type]->Fill(trg->GetEta(),phi);

    h1px_trg[trg->type]->Fill(trg->px);
    h1py_trg[trg->type]->Fill(trg->py);
    h1pz_trg[trg->type]->Fill(trg->pz);

    h1pc3dphis_trg[trg->type]->Fill(trg->pc3dphis);
    h1pc3dz_trg[trg->type]->Fill(trg->pc3dz);
    h1pc3dr_trg[trg->type]->Fill(trg->GetPc3dr());

    h1trkdphis_trg[trg->type]->Fill(trg->trkdphis);
    h1trkdz_trg[trg->type]->Fill(trg->trkdz);
    h1trkdr_trg[trg->type]->Fill(trg->GetTrkdr());

    if (verbosity > 2)
	cout<<"Filled trigger histograms"<<endl;

    return;
}

int GetCentBin(float cent)
{
    for (int i = 0; i < ncent; i++)
    {
	if (cent > limcent[i] && cent <= limcent[i+1])
    	{
    	    return i;
    	}
    }

    return -1;

}

int GetDPsiBin(float psi,float phi)
{
    float dpsi = (psi - phi)/PI;

    while (dpsi < 0.) dpsi+=1.;
    while (dpsi > 1.) dpsi-=1.;

    if (dpsi > .75 || dpsi <= .25) return 0;

    if (dpsi > .25 && dpsi <= .75) return 1;

    cerr<<"Should not be at this line.  Something wrong with dpsi binning scheme."<<endl;
    return -1;
}

int GetPtAssocBin(float pt)
{
    for (int i = 0; i < npta; i++)
    {
	if (limpta[i] <= pt && pt < limpta[i+1])
	{
	    return i;
	}
    }

    if (verbosity > 1)
	cout<<"WARNING!!  Could not get pta bin! pt = "<<pt<<endl;
    return -1;
}

int GetXiBin(float pt_trg, float pt_trk)
{
    float zt = pt_trk/pt_trg;
    float xi = log(1/zt);

    for (int i = 0; i < nxi; i++)
    {
	if (limxi[i] <= xi && xi < limxi[i+1])
	{
	    return i;
	}
    }

    if (verbosity > 1)
	cout<<"WARNING: Could not get xi bin. pt_trg = "<<pt_trg<<" pt_trk = "<<pt_trk<<" zt = "<<zt<<" xi = "<<xi<<endl;
    return -1;
}

int GetPtTrigBin(float pt)
{
    for (int i = 0; i < nptt; i++)
    {
	if (limptt[i] <= pt && pt < limptt[i+1])
	{
	    return i;
	}
    }

    cerr<<"WARNING!!  Could not get ptt bin!"<<endl;
    return -1;
}

int GetEmcRunMaps(int runno)
{
    if (verbosity > 0) cout<<"Reading EMC Run Specific Hot Maps"<<endl;

    for (int isec = 0; isec < GA::nsec; isec++)
    {
    	if (h2runmap[isec] != NULL)
    	{
    	    delete h2runmap[isec];
    	    h2runmap[isec] = NULL;
    	}
    }

    const char* loc_run = "run14_emchotmaps.root";

    TFile* file_run = new TFile(loc_run,"READ");

    if (!file_run)
    {
	cerr<<"Could not open input "<<loc_run<<endl;
	return -1;
    }
    if (file_run->IsZombie())
    {
	cerr<<"Input "<<loc_run<<" is zombie."<<endl;
	return -1;
    }
    if (verbosity > 0)
    {
	cout<<"Input "<<loc_run<<" opened successfully."<<endl;
    }

    for (int isec = 0; isec < GA::nsec; isec++)
    {
	h2runmap[isec] = (TH2I*)file_run->Get(Form("hotmap_run%d_sec%d",runno,isec));
	if (h2runmap[isec] == NULL)
	{
	    printf("Could not find hotmap_run%d_sec%d\n",runno,isec);
	    cerr<<"ABORT!!"<<endl;

	    file_run->Close();
	    delete file_run;

	    return -1;
	}

	h2runmap[isec]->SetDirectory(0);

	if (verbosity > 1)
	    printf("Retrieved hotmap_run%d_sec%d\n",runno,isec);
    }

    file_run->Close();
    delete file_run;

    return 0;
}

int GetPc3Pars()
{
    //==================================================
    // Pc3 Vars
    //==================================================

    if (verbosity > 0) cout<<"Reading PC3 Matching Parameters"<<endl;

    string line;

    const char* dphis_loc = "pc3dphis.txt";
    ifstream ifs_dphis(dphis_loc, ifstream::in);

    if (!ifs_dphis.is_open())
    {
	cerr<<"Error opening pc3dphis file "<<dphis_loc<<endl;
	return -1;
    }
    else
    {
	if (verbosity > 0) cout<<dphis_loc<<" opened successfully"<<endl;
    }
    while ( getline(ifs_dphis, line) )
    {
	if (!isdigit(line.at(0))) continue;
	istringstream iss_line(line.data());
	int isec,iybin;
	double mean,sigma;
	iss_line >> isec >> iybin >> mean >> sigma;
	pc3dphispar[isec][iybin][0] = mean;
	pc3dphispar[isec][iybin][1] = sigma;
    }

    const char* dz_loc = "pc3dz.txt";
    ifstream ifs_dz(dz_loc, ifstream::in);

    if (!ifs_dz.is_open())
    {
	cerr<<"Error opening pc3dz file "<<dz_loc<<endl;
	return -1;
    }
    else
    {
	if (verbosity > 0) cout<<dz_loc<<" opened successfully"<<endl;
    }
    while ( getline(ifs_dz, line) )
    {
	if (!isdigit(line.at(0))) continue;
	istringstream iss_line(line.data());
	int isec,izbin;
	double mean,sigma;
	iss_line >> isec >> izbin >> mean >> sigma;
	pc3dzpar[isec][izbin][0] = mean;
	pc3dzpar[isec][izbin][1] = sigma;
    }
    return 0;
}



int GetPc3Maps()
{
    if (verbosity > 0) cout<<"Reading PC3 Maps"<<endl;

    const char* loc = "pc3bad.txt";

    ifstream ifs(loc, ifstream::in);
    if (!ifs.is_open())
    {
	cerr<<"Error opening file "<<loc<<endl;
	return -1;
    }
    else
    {
	if (verbosity > 0) cout<<loc<<" opened successfully"<<endl;
    }
    string line;
    int sector, iy, iz, flag;
    while ( getline(ifs, line) )
    {
	istringstream iss_line(line.data());
	iss_line >> sector >> iy >> iz >> flag;
	if (verbosity > 4)
	{
	    cout << sector <<" "<< iy <<" "<< iz <<" "<< flag << endl;
	}
	if (flag != 0)
	{
	    pc3bad[sector][iy][iz] = true;
	}
    }


    //3x3 pc3 dead cut
    for (int isec = 0; isec < GA::nsec; isec++)
    {
	int nz = -999;
	int ny = -999;
	GA::GetNzNy(nz,ny,isec);

	for (int iy = 0; iy < ny; iy++)
	{
	    for (int iz = 0; iz < nz; iz++)
	    {

		for (int jy = iy-1; jy <= iy+1; jy++)
		{
		    if (jy < 0 || jy >= ny) continue;

		    for (int jz = iz-1; jz <= iz+1; jz++)
		    {
			if (jz < 0 || jz >= nz) continue;

			if (pc3bad[isec][iy][iz])
			    pc3bad3x3[isec][jy][jz]=true;
		    }
		}
	    }
	}
    }

    return 0;
}

int GetEmcMaps()
{
    if (verbosity > 0) cout<<"Reading EMC Maps"<<endl;

    const char* loc_wu = "emcwarnuncal.txt";

    ifstream ifs_wu(loc_wu, ifstream::in);
    if (!ifs_wu.is_open())
    {
	cerr<<"Error opening file "<<loc_wu<<endl;
	return -1;
    }
    else
    {
	if (verbosity > 0) cout<<loc_wu<<" opened successfully"<<endl;
    }
    string line;
    int sector, iy, iz, flag;
    while ( getline(ifs_wu, line) )
    {
	istringstream iss_line(line.data());
	iss_line >> sector >> iy >> iz >> flag;
	if (flag != 0)
	{
	    emcwu[sector][iy][iz] = true;
	}
    }

    //-------------------- EMC Hot Maps --------------------

    if (verbosity > 0) cout<<"Reading EMC Common Hot Maps"<<endl;

    const char* loc_com = "run14_emchotmaps_common.root";

    TFile* file_com = new TFile(loc_com,"READ");

    if (!file_com)
    {
	cerr<<"Could not open input "<<loc_com<<endl;
	return -1;
    }
    if (file_com->IsZombie())
    {
	cerr<<"Input "<<loc_com<<" is zombie."<<endl;
	return -1;
    }
    if (verbosity > 0)
    {
	cout<<"Input "<<loc_com<<" opened successfully."<<endl;
    }

    for (int isec = 0; isec < GA::nsec; isec++)
    {
	h2commap[isec] = (TH2I*)file_com->Get(Form("hotmap_sec%d",isec));
	if (h2commap[isec] == NULL)
	{
	    cerr<<"ABORT!!"<<endl;
	    return -1;
	}
	if (verbosity > 1)
	    printf("Retrieved hotmap_sec%d\n",isec);

	h2commap[isec]->SetDirectory(0);
    }

    file_com->Close();
    delete file_com;

    return 0;
}

int GetFunctions(const char* loc, const char* name, TF1* f1[ncent_fun])
{
    if (verbosity > 0) cout<<"Getting functions "<<name<<"? from "<<loc<<endl;

    memset(f1,0,sizeof(TF1*)*ncent_fun);

    TFile* file = new TFile(loc,"READ");

    if (!file)
    {
	cerr<<"Could not open input "<<loc<<endl;
	return -1;
    }
    if (file->IsZombie())
    {
	cerr<<"Input "<<loc<<" is zombie."<<endl;
	return -1;
    }
    if (verbosity > 0)
    {
	cout<<"Input "<<loc<<" opened successfully."<<endl;
    }

    for (int icent = 0; icent < ncent_fun; icent++)
    {
	ostringstream oss;
	oss<<name<<icent;

	TF1* f1tmp = (TF1*)file->Get(oss.str().c_str());
	if (f1tmp == NULL)
	{
	    cerr<<"ABORT!!"<<endl;
	    return -1;
	}
	if (verbosity > 1)
	    cout<<"Retrieved "<<f1tmp->GetName()<<endl;

	f1[icent] = (TF1*) f1tmp->Clone();
    }

    file->Close();
    delete file;

    return 0;
}

int GetWeights(double& dirwt, double& decwt, double& dirwtY, double& decwtY, TF1* ftageff[ncent_fun], TF1* fmistag[ncent_fun], int type, float cent, float pt)
{
    int icent = cent - 0.5;

    double tageff = -999;

    tageff = ftageff[icent]->Eval(pt);

    if (tageff < 0. || tageff > 1.)
    {
	cerr<<"Tagging efficiency evaluated to "<<tageff<<" is out of bounds."<<endl;
	cerr<<ftageff[icent]->GetName()<<endl;
	return -1;
    }

    double mistag = -999;

    mistag = fmistag[icent]->Eval(pt);

    if (mistag > 1.)
    {
	cerr<<"Mistag probability out of bounds."<<endl;
	cerr<<fmistag[icent]->GetName()<<"->Eval("<<pt<<") = "<<mistag<<endl;	    
	return -1;
    }

    if (mistag < 0.)
    {
	if (mistag < -.01)
	{
	    cerr<<"Mistag probability out of bounds."<<endl;
	    cerr<<fmistag[icent]->GetName()<<"->Eval("<<pt<<") = "<<mistag<<endl;	    
	    return -1;
	}
	cout<<"Mistag probability less than 0. but greater than -0.01"<<endl;
	cout<<fmistag[icent]->GetName()<<"->Eval("<<pt<<") = "<<mistag<<endl;
	cout<<"Setting mistag prob to 0."<<endl;
	mistag = 0.;
    }

    double r = asr, t = tageff, m = mistag;
    double d = (-r*t*m + r*t + t*m - m);
	
    if (type == DIR)
    {
	dirwt = t/(t-m);
	decwt = m/(m-t);

	dirwtY = r*t/d;
	decwtY = m*(1-t+t*r)/-d;
    }
    else if (type == DEC)
    {
	dirwt = (1-t)/(m-t);
	decwt = (1-m)/(t-m);

	dirwtY = (1-t)/-d;
	decwtY = (1-m)*(1-t*t*r)/d;
    }

    return 0;
}


# ifndef __CINT__
int main()
{
    gamAnaTree();
    return 0;
}
# endif
