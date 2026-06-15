#include "gamAna.h"

#include "ChargedTrack.h"
#include "CentralTrackData.h"
#include "CentralArmTrack.h"
#include "TrackCollectionBase.h"
#include "TrackCollection.h"
#include "TriggerParticle.h"
#include "EmcCluster.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <gsl/gsl_rng.h>

#include <PHCompositeNode.h>
#include <Fun4AllReturnCodes.h>
#include <recoConsts.h>
#include <getClass.h>
#include <TriggerHelper.h>
#include <ErtUtils.h>
#include <PHPoint.h>
#include <PHAngle.h>

#include <SvxClusterList.h>
#include <SvxCluster.h>
#include <SvxClusterInfo.h>
#include <SvxCentralTrackList.h>
#include <SvxCentralTrack.h>
#include <SvxSegmentList.h>
#include <SvxCentralTrackRecalList.h>
#include <SvxCentralTrackRecal.h>
#include <RunHeader.h>
#include <EventHeader.h>
#include <PHGlobal.h>
#include <PHCentralTrack.h>
#include <PHSnglCentralTrack.h>
#include <VtxOut.h>
#include <ErtOut.h>

#include <svxDetectorGeo.hh>
#include <SvxSensor.h>

#include <emcClusterContainer.h>
#include <emcClusterContent.h>
#include <compactCNT/PadHitMap.h>

#include <RpConst.h>
#include <PreviousEvent.h>
#include <RpSumXYObject.h>
#include <ReactionPlaneCalibv1.h>
#include <ReactionPlaneObject.h>
#include <ReactionPlaneSngl.h>

#include <McEvalSingleList.h>


#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TRandom.h>

#include <TSQLServer.h>
#include <TSQLStatement.h>

#include <TTimeStamp.h>

#include "TOAD.h"

using namespace std;
using namespace findNode;

const float PI = (float) TMath::Pi();

const float gamAna::pi0pars[GA::nsec][2] =
{
    {0.138,0.01},
    {0.138,0.01},
    {0.138,0.01},
    {0.138,0.01},
    {0.138,0.01},
    {0.138,0.01},
    {0.138,0.01},
    {0.138,0.01}
};

const char* gamAna::spec[nspec] = {"dir","dec","pi0"};


const std::string gamAna::trkcut_name[ntrkcut] = {"RAW","ZED","QUAL","PC3","EMC","N0","DEP"};

//==============================================================
gamAna::gamAna(std::string &filename) :
    SubsysReco("GAMANA"),
    outfile(NULL),
    outfilename(filename),
    d_runhead(NULL),
    d_eventhead(NULL),
    d_global(NULL),
    d_trx(NULL),
    d_vtxout(NULL),
    d_peve(NULL),
    d_svxtrx(NULL),
    d_rp(NULL),
    d_emccontain(NULL),
    dofullana(true),
    dotrackana(false),
    domistagana(false),
    readeventcount(true),
    useemc3x3(false),
    usemyemc(true),
    saveallevt(false),
    tracktype(GA::CNT),
    rpdet(GA::FVTX)
{

    if (verbosity > 0)
    {
	cout<<"This is a gamAna"<<endl;
	if (!dofullana)
	    cout<<"Doing quick analysis on EMC clusters."<<endl;
	TTimeStamp* tts = new TTimeStamp();
	tts->Print();
	delete tts;
    }

    init();

}
//==============================================================

int gamAna::Init(PHCompositeNode *topNode) {

    evtno=0;
    goodevtno=0;

    if (verbosity > 0)
    {
	cout<<"gamana::begin Init()"<<endl;
	TTimeStamp* tts = new TTimeStamp();
	tts->Print();
	delete tts;
	if (!dofullana)
	    cout<<"Doing quick analysis on EMC clusters."<<endl;
	cout<<"domistagana = "<<domistagana<<"; dotrackana = "<<dotrackana<<endl;
	cout<<"Track type = "<<tracktype<<endl;
    }

    if (tracktype == GA::SVX)
    {
	cout<<"WARNING! tracktype SVX not yet supported. Return ABORTRUN"<<endl;
	return ABORTRUN;
    }

    if (ptmin_pi0trig > ptmin_gamtrig)
    {
	cout<<PHWHERE<<"pi0 trigger minpt ("<<ptmin_pi0trig<<") is greater than gamma trigger minpt ("<<ptmin_gamtrig<<").  Certain decisions are based on the assumption that this is not true.  Returning ABORTRUN."<<endl;
	return ABORTRUN;
    }

    if (dofullana || domistagana)
    {
	if (verbosity > 0) cout<<PHWHERE<<"Reading Default EMC Maps with TOAD"<<endl;

	TOAD toad_loader("zrowan");

	string emcwarnuncal_location = toad_loader.location("emcwarnuncal.txt");
	ifstream ifs_warnuncal(emcwarnuncal_location.c_str(), ifstream::in);
	if (!ifs_warnuncal.is_open())
	{
	    cerr<<"Error opening file "<<emcwarnuncal_location<<endl;
	    return ABORTRUN;
	}
	else
	{
	    if (verbosity > 0) cout<<emcwarnuncal_location<<" opened successfully"<<endl;
	}
	string line;
	int sector, iy, iz, flag;
	while ( getline(ifs_warnuncal, line) )
	{
	    istringstream iss_line(line.data());
	    iss_line >> sector >> iy >> iz >> flag;
	    if (flag != 0)
	    {
		if (((flag&64)==64)||((flag&16)==16)||((flag&8)==8))
		    emcbad[sector][iy][iz] = true;

		else if ((flag&32)==32 || (flag&4)==4 || (flag&2)==2)
		    emcnearbad[sector][iy][iz] = true;
	    }
	}

	//-------------------- EMC Hot Maps (All Runs) --------------------

	if (verbosity > 0) cout<<PHWHERE<<"Reading EMC Common Hot Maps with TOAD"<<endl;
	string commonmap_location = toad_loader.location("run14_emchotmaps_common.root");
	TFile* commonmap_file = new TFile(commonmap_location.c_str(),"READ");

	if (!commonmap_file)
	{
	    cerr<<"Could not open input "<<commonmap_location<<endl;
	    return ABORTRUN;
	}
	if (commonmap_file->IsZombie())
	{
	    cerr<<"Input "<<commonmap_location<<" is zombie."<<endl;
	    return ABORTRUN;
	}
	if (verbosity > 0)
	{
	    cout<<"Input "<<commonmap_location<<" opened successfully."<<endl;
	}

	for (int isec = 0; isec < GA::nsec; isec++)
	{
	    h2commonmap[isec] = (TH2I*)commonmap_file->Get(Form("hotmap_sec%d",isec));
	    if (h2commonmap[isec] == NULL)
	    {
		cout<<"ABORT!!"<<endl;
		return ABORTRUN;
	    }
	    if (verbosity > 0)
		printf("Retrieved hotmap_sec%d\n",isec);

	    h2commonmap[isec]->SetDirectory(0);
	}

	commonmap_file->Close();
	delete commonmap_file;

	GetPc3Pars(toad_loader);
	GetPc3Maps(toad_loader);
    }

    outfile = new TFile(outfilename.c_str(),"RECREATE");

    initHistos();

    if (dofullana)
    {

	trig = new TriggerParticle();

	trigtree = new TTree("TrigTree","TriggerParticleTree");
	trigtree->Branch("TrigPart","TriggerParticle",&trig);
	trigtree->Branch("RunNumber",&runno,"runno/I");
	trigtree->Branch("EventNumber",&evtno,"evtno/I");
	trigtree->Branch("Psi",&psi,"psi/F");
	trigtree->Branch("Zvertex",&z0,"z0/F");
	trigtree->Branch("Centrality",&centrality,"centrality/F");

	delete trig;
    }

    if (dofullana || dotrackana)
    {
	for (int icent = 0; icent < GA::ncent_tree; icent++)
	{
	    for (int ipsi = 0; ipsi < GA::npsi_tree; ipsi++)
	    {
		for (int iz = 0; iz < GA::nzv_tree; iz++)
		{
		    trktree[icent][ipsi][iz] =
			new TTree(Form("trktree_C%i_P%i_Z%i",icent,ipsi,iz),
				  Form("trktree_C%i_P%i_Z%i",icent,ipsi,iz));


		    if (tracktype == GA::SVXCNT)
		    {
			cout<<"Svx Cnt Tracks no longer supported! Aborting Run!"<<endl;
			return ABORTRUN;
		    }
		    else if (tracktype == GA::CNT)
		    {

		    	trkhold = new TrackCollection<CentralArmTrack>();
		    	trktree[icent][ipsi][iz]->Branch("TrackBranch","TrackCollection<CentralArmTrack>",&trkhold);
		    }


		    trktree[icent][ipsi][iz]->Branch("RunNumber",&runno,"runno/I");
		    trktree[icent][ipsi][iz]->Branch("EventNumber",&evtno,"evtno/I");
		    trktree[icent][ipsi][iz]->Branch("Psi",&psi,"psi/F");
		    trktree[icent][ipsi][iz]->Branch("Zvertex",&z0,"z0/F");
		    trktree[icent][ipsi][iz]->Branch("Centrality",&centrality,"centrality/F");

		    delete trkhold;

		}
	    }
	}
    }


    //-------------------------------------------------------------------
    if (verbosity > 0)
    {
	cout<<"ptmin_pi0trig = "<<ptmin_pi0trig<<"; ptmin_gamtrig = "<<ptmin_gamtrig<<"; asymmax = "<<asymmax<<endl;
	cout<<"gamana::end Init()"<<endl;
	TTimeStamp* tts = new TTimeStamp();
	tts->Print();
	delete tts;
    }
    return 0;
}

//==============================================================

void gamAna::initHistos(TFile* myemchdmap)
{

    //------------------- Global Histograms -----------------------

    h1tracktype = new TH1F("tracktype","tracktype",3,-0.5,2.5);

    h1evtraw = new TH1F("evtraw","evtraw",1,-0.5,0.5);
    h1evtmod = new TH1F("evtmod","evtmod",100,0.,100.);
    h1skipevtmod = new TH1F("skipevtmod","skipevtmod",100,0.,100.);
    h1evtinrun = new TH1F("evtinrun","evtinrun",1,-0.5,0.5);

    h1tickcut = new TH1F("tickcut","tickcut",2,-0.5,1.5);

    h1z0 = new TH1F("z0","d_vtxout->get_ZVertex();z0;Events",600,-30.,30.);
    h1cent_raw = new TH1F("cent_raw","before cutsd_global->getCentrality();centrality;Events",101,-0.5,100.5);

    h1bcent = new TH1F("bcent","centbin",GA::ncent_tree,0,GA::ncent_tree);

    h2bz0cent = new TH2F("bz0cent","bz0cent;z vertex bin;centrality",GA::nzv_tree,0,GA::nzv_tree,(int)GA::maxcent,0.5,GA::maxcent+0.5);

    h2z0zbbc = new TH2F("z0zbbc","z0zbbc",40,-20,20,40,-20,20);
    h2z0zbbc->SetXTitle("default z vertex [cm]");
    h2z0zbbc->SetYTitle("bbc z vertex [cm]");

    h1bpsi = new TH1F("bpsi","psibin",GA::npsi_tree,0,GA::npsi_tree);
    h1psi = new TH1F("psi","psi",npsibin,psimin_h,psimax_h);

    h1psi_bbc_raw = new TH1F("psi_bbc_raw","psi bbc raw",npsibin,psimin_h-0.5,psimax_h+0.5);
    h1psi_bbcs = new TH1F("psi_bbcs","psi_bbcs",npsibin,psimin_h,psimax_h);
    h1psi_bbcn = new TH1F("psi_bbcn","psi_bbcn",npsibin,psimin_h,psimax_h);

    h2centpsidiff_bbcns = new TH2F("centpsidiff_bbcns","centpsidiff_bbcns",(int)GA::maxcent,0.5,GA::maxcent+0.5,npsibin,-psimax_h,psimax_h);
    h2centpsidiff_bbcns->SetXTitle("centrality");
    h2centpsidiff_bbcns->SetYTitle("(psi_bbcn-psi_bbcs)/PI");

    h1psi_fvtx_raw = new TH1F("psi_fvtx_raw","psi fvtx raw",npsibin,psimin_h-0.5,psimax_h+0.5);
    h1psi_fvtxs = new TH1F("psi_fvtxs","psi_fvtxs",npsibin,psimin_h,psimax_h);
    h1psi_fvtxn = new TH1F("psi_fvtxn","psi_fvtxn",npsibin,psimin_h,psimax_h);

    h2centpsidiff_fvtxns = new TH2F("centpsidiff_fvtxns","centpsidiff_fvtxns",(int)GA::maxcent,0.5,GA::maxcent+0.5,npsibin,-psimax_h,psimax_h);
    h2centpsidiff_fvtxns->SetXTitle("centrality");
    h2centpsidiff_fvtxns->SetYTitle("(psi_fvtxn-psi_fvtxs)/PI");

    h2psi_bbc_bbcs = new TH2F("psi_bbc_bbcs","psi_bbc_bbcs;psi bbc;psi bbcs",npsibin_2d,psimin_h_2d,psimax_h_2d,npsibin_2d,psimin_h_2d,psimax_h_2d);
    h2psi_bbc_bbcn = new TH2F("psi_bbc_bbcn","psi_bbc_bbcn;psi bbc;psi bbcn",npsibin_2d,psimin_h_2d,psimax_h_2d,npsibin_2d,psimin_h_2d,psimax_h_2d);
    h2psi_bbcn_bbcs = new TH2F("psi_bbcn_bbcs","psi_bbcn_bbcs;psi bbcn;psi bbcs",npsibin_2d,psimin_h_2d,psimax_h_2d,npsibin_2d,psimin_h_2d,psimax_h_2d);

    h2psi_fvtx_fvtxs = new TH2F("psi_fvtx_fvtxs","psi_fvtx_fvtxs;psi fvtx;psi fvtxs",npsibin_2d,psimin_h_2d,psimax_h_2d,npsibin_2d,psimin_h_2d,psimax_h_2d);
    h2psi_fvtx_fvtxn = new TH2F("psi_fvtx_fvtxn","psi_fvtx_fvtxn;psi fvtx;psi fvtxn",npsibin_2d,psimin_h_2d,psimax_h_2d,npsibin_2d,psimin_h_2d,psimax_h_2d);
    h2psi_fvtxn_fvtxs = new TH2F("psi_fvtxn_fvtxs","psi_fvtxn_fvtxs;psi fvtxn;psi fvtxs",npsibin_2d,psimin_h_2d,psimax_h_2d,npsibin_2d,psimin_h_2d,psimax_h_2d);

    if (dofullana || domistagana)
    {
	h2centnclus = new TH2F("centnclus","centnclus",20,0,100,120,0,1200);
	h2centnclus->SetXTitle("centrality");
	h2centnclus->SetYTitle("N EMC Clusters");
    }

    if (dofullana || dotrackana)
    {
	h2centnpart = new TH2F("centnpart","centnpart",20,0,100,100,0,1000);
	h2centnpart->SetXTitle("centrality");
	h2centnpart->SetYTitle("N Central Tracks");
    }

    h1BBC = new TH1F("BBC","BBC",2,-0.5,1.5);
    h1BBCall = new TH1F("BBCall","BBCall",2,-0.5,1.5);
    h1BBCnar = new TH1F("BBCnar","BBCnar",2,-0.5,1.5);
    h1BBCnov = new TH1F("BBCnov","BBCnov",2,-0.5,1.5);

    h1ERT = new TH1F("ERT","ERT",2,-0.5,1.5);
    h1ERTall = new TH1F("ERTall","ERTall",2,-0.5,1.5);
    h1ERTBBC = new TH1F("ERTBBC","ERTBBC",2,-0.5,1.5);

    //--------------- EMCal Sector Dependent Histograms ---------------

    if (dofullana || domistagana)
    {

	for (int isec = 0; isec < GA::nsec + 6; isec++)
	{
	    if (isec >= GA::nsec)
	    {
		int tmpsec = isec - GA::nsec;
		int arm = tmpsec/3;
		tmpsec -= arm*3;
		int sector = tmpsec + arm*4;
		h2mpt[isec] = new TH2F(Form("mpt_%i",isec),
				       Form("mpt_%i sector %i and sector %i",isec,sector,sector+1),
				       500,0.,1.,120,0.,60.);
	    }

	    else
		h2mpt[isec] = new TH2F(Form("mpt_%i",isec),
				       Form("mpt_%i",isec),
				       500,0.,1.,120,0.,60.);

	    h2mpt[isec]->SetYTitle("p_{T} [GeV/c]");
	    h2mpt[isec]->SetXTitle("m [GeV/c^{2}]");
	}

	for (int isec = 0; isec < GA::nsec; isec++)
	{
	    h1r[isec] = new TH1F(Form("r_sec%i",isec),
				 Form("r_sec%i",isec),
				 300,450.,600.);

	    h1s[isec] = new TH1F(Form("s_sec%i",isec),
				 Form("s_sec%i",isec),
				 300,450.,600.);


	    if (myemchdmap!=NULL)
	    {
		h2emchdmap[isec] = new TH2I();

		if (verbosity > 0)
		{
		    cout<<"WARNING!!!!! emchdmap.root has opposite sector number convention."
			<<" This might no longer be the case.  Be careful!"<<endl;
		}
		*h2emchdmap[isec] = *(TH2I*)myemchdmap->Get(Form("emchdmap_%i",(int)fabs(isec-7)));
	    }

	    double why = GA::nySc;
	    double zee = GA::nzSc;
	    if (isec < 2) {
		why = GA::nyGl;
		zee = GA::nzGl;
	    }

	    h2iziy_emcbad[isec] =
		new TH2I(Form("iziy_emcbad_sec%i",isec),
			 Form("EMC Warn Uncal Map.  Sector %i",isec), zee,0,zee,why,0,why);
	    h2iziy_emcbad[isec]->SetYTitle("iypos");
	    h2iziy_emcbad[isec]->SetXTitle("izpos");

	    h2iziy_nomapcut[isec] =
		new TH2I(Form("iziy_nomapcut%i",isec),
			 Form("EMC Cluster Position.  No Map Cuts.  Sector %i",isec), zee,0,zee,why,0,why);
	    h2iziy_nomapcut[isec]->SetYTitle("iypos");
	    h2iziy_nomapcut[isec]->SetXTitle("izpos");

	    h2iziy_bordwu[isec] =
		new TH2I(Form("iziy_bordwu%i",isec),
			 Form("EMC Cluster Position.  One Tower Border & EMC Warn Uncal Map Cuts.  Sector %i",isec), zee,0,zee,why,0,why);
	    h2iziy_bordwu[isec]->SetYTitle("iypos");
	    h2iziy_bordwu[isec]->SetXTitle("izpos");

	    if (domistagana)
	    {
		h2iziy_fid[isec] =
		    new TH2I(Form("iziy_fid%i",isec),
			     Form("EMC Clusters Outside Fiducial Region.  Sector %i",isec), zee,0,zee,why,0,why);
		h2iziy_fid[isec]->SetYTitle("iypos");
		h2iziy_fid[isec]->SetXTitle("izpos");
	    }

	}

	h2zphi_arm[1] = new TH2F("zphi_arm1","zphi_arm1",200,-200,200,200,-.3,.4);
	h2zphi_arm[1]->SetXTitle("clus->z()");
	h2zphi_arm[1]->SetYTitle("clus->phi()");

	h2zphi_arm[0] = new TH2F("zphi_arm0","zphi_arm0",200,-200,200,200,.6,1.3);
	h2zphi_arm[0]->SetXTitle("clus->z()");
	h2zphi_arm[0]->SetYTitle("clus->phi()");
    }

    //--------------- Histograms for pi0s ------------------

    if (domistagana)
    {
	string stradd[ntrigcat] = {"","_emc3x3","_myemc","_cv","_pc3","3x3"};
	string strsum = "";
	for (int icat = 0; icat < ntrigcat; icat++)
	{
	    strsum += stradd[icat];

	    h2ptcent_trigcat[icat] = new TH2F(Form("ptcent%s",strsum.c_str()),Form("ptcent%s;pt;centrality",strsum.c_str()),nebin,0.,emax_h,100,0.,100.);
	    h2ecent_trigcat[icat] = new TH2F(Form("ecent%s",strsum.c_str()),Form("ecent%s;ecore;centrality",strsum.c_str()),nebin,0.,emax_h,100,0.,100.);

	    for (int icent = 0; icent < ninvmcentbin; icent++)
	    {
		for (int ipt = 0; ipt < ninvmptbin; ipt++)
		{
		    h1invm_trigcat[icat][icent][ipt] = new TH1F(
			Form("invm%s_%i_%i",strsum.c_str(),icent,ipt),
			Form("invm%s_%i_%i  integer pt = %i",strsum.c_str(),icent,ipt,ipt + ((int)ptmin_gamtrig)),
			125,.05,.3);
		}
	    }
	}
    }

    if (dofullana || domistagana)
    {
	h1asym_pi0 = new TH1F("asym_pi0","asym_pi0",100,0.,1.);
    }

    //--------------- Histograms for all clusters ------------------

    if (dofullana || domistagana)
    {

	h1histosec = new TH1F("histosec","histosec",17,-2.5,14.5);

	h2ptcent = new TH2F("ptcent_notrigcat","ptcent OK Clusters;pt;centrality",nebin,0.,emax_h,100,0.,100.);
	h2ecent = new TH2F("ecent_notrigcat","ecent OK Clusters;ecore;centrality",nebin,0.,emax_h,100,0.,100.);

	h1e_nocut = new TH1F("e_nocut","e_nocut",nebin,0.,emax_h);
	h1ex = new TH1F("ex","ex",nebin,-emax_h,emax_h);
	h1ey = new TH1F("ey","ey",nebin,-emax_h,emax_h);
	h1ez = new TH1F("ez","ez",nezbin,-ezmax_h,ezmax_h);

	h1prob_nocut = new TH1F("prob_nocut","prob_nocut",100,0.,1.);

	h2armphi_emclus = new TH2F("armphi_emclus","EMCal Cluster; Arm; phi",2,-0.5,1.5,200,-2.,2.);
    }

    //--------------- Histograms for candidates only -----------------

    if (dofullana)
    {

	h1lowasym = new TH1F("lowasym","lowest asym that tags photon as decay;lowest asym;tagged decays",1000,0.,1.);
	h1ncand = new TH1F("ncand","ncand",1,0,1);

	//--------------- Histograms for all frag func anas -----------------


	for (int ispec = 0; ispec < nspec; ispec++)
	{
	    h1e_spec[ispec] = new TH1F(Form("e_%s",spec[ispec]),
				       Form("e_%s",spec[ispec]),
				       nebin,0.,emax_h);

	    h1ntrigperevt[ispec] = new TH1F(Form("ntrigperevt_%s",spec[ispec]),
					    Form("ntrigperevt_%s",spec[ispec]),
					    100,0,100);
	}
    }

    //--------------- Histograms for central tracks -----------------

    if (dofullana || dotrackana)
    {

	string tmpstr = "";
	for (unsigned int icut = 0; icut < ntrkcut; icut++)
	{

	    tmpstr+=trkcut_name[icut];

	    h2ptcent_trk[icut] = new TH2F(Form("ptcent_trk_%i",icut),
					  Form("Central Track: %s;p_{T};centrality",tmpstr.c_str()),
					  (GA::maxpta-GA::minpta)*4,GA::minpta,GA::maxpta,
					  (int)GA::maxcent,0.5,GA::maxcent+.5);

	    if (icut == 0)
		tmpstr = "";
	    else
		tmpstr+=" ";
	}

	h1zed_trk_raw =  new TH1F("zed_trk_raw","zed_trk_raw",201,-100.5,100.5);

	h1dcarm = new TH1F("dcarm","dcarm",2,-.5,1.5);

	h1px = new TH1F("px","track px;px",nebin,-emax_h,emax_h);
	h1py = new TH1F("py","track py;py",nebin,-emax_h,emax_h);
	h1pz = new TH1F("pz","track pz;pz",nezbin,-ezmax_h,ezmax_h);

	h1n0 = new TH1F("n0","track n0;n0",100,-50.,50);
	h1ecore = new TH1F("ecore","track ecore;ecore",nebin,-emax_h,emax_h);
	h1eop = new TH1F("eop","eop",100,0.,2.);
	h1dep = new TH1F("dep","dep",200,-40.,10.);

	h1svxnclus = new TH1F("svxnclus","svxnclus",20,0.,20.);

	h1chi2ndf = new TH1F("chi2ndf","chi2/ndf",25,0.,10.);

	h1quality = new TH1F("quality","quality",102,-1.5,100.5);
	h1charge = new TH1F("charge","charge",5,-2.5,2.5);
	h1pattern = new TH1F("pattern","pattern",18,-1.5,16.5);

	h1pc3dphi = new TH1F("pc3dphi","pc3dphi",nphimatchbin,phimatchmin_h,phimatchmax_h);
	h1pc3dz = new TH1F("pc3dz","pc3dz",nzmatchbin,zmatchmin_h,zmatchmax_h);

	h1pc3sdphi = new TH1F("pc3sdphi","pc3sdphi",nmatchbin,matchmin_h,matchmax_h);
	h1pc3sdz = new TH1F("pc3sdz","pc3sdz",nmatchbin,matchmin_h,matchmax_h);
	h1pc3sig = new TH1F("pc3sig","pc3sig",nmatchbin,0.,matchmax_h);

	h1emcdphi = new TH1F("emcdphi","emcdphi",nphimatchbin,phimatchmin_h,phimatchmax_h);
	h1emcdz = new TH1F("emcdz","emcdz",nzmatchbin,zmatchmin_h,zmatchmax_h);

	h1emcsdphi = new TH1F("emcsdphi","emcsdphi",nmatchbin,matchmin_h,matchmax_h);
	h1emcsdz = new TH1F("emcsdz","emcsdz",nmatchbin,matchmin_h,matchmax_h);
	h1emcsig = new TH1F("emcsig","emcsig",nmatchbin,0.,matchmax_h);

	h1dca2d = new TH1F("dca2d","dca2d",ndcabin,dcamin_h,dcamax_h);
	h1dcaz = new TH1F("dcaz","dcaz",ndcabin,dcamin_h,dcamax_h);
    }
}
//End initHistos()

//==============================================================

int gamAna::GetNodes(PHCompositeNode *topNode)
{

    //-------------------EventHeader----------------------------------
    d_eventhead = findNode::getClass<EventHeader>(topNode, "EventHeader");
    if(d_eventhead == NULL) { cerr << PHWHERE << " EventHeader node not found." << endl; return DISCARDEVENT; }
    //------------------------------------------------------------------

    //-------------------PreviousEvent----------------------------------
    d_peve = NULL;
    d_peve    = getClass<PreviousEvent>(topNode, "PreviousEvent"); // Tick
    if(!d_peve){ cerr << PHWHERE << " No PreviousEvent object !" << endl;   return DISCARDEVENT; }
    //-------------------------------------------------------

    //-------------------PHGlobal----------------------------------
    d_global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
    if (!d_global) { cerr << PHWHERE << " PHGlobal node not found." << endl; return DISCARDEVENT; }
    //-------------------------------------------------------

    if (dofullana || dotrackana)
    {
	//-------------------Central arm tracks-----------------
	d_trx = getClass<PHCentralTrack>(topNode,"PHCentralTrack");
	if(d_trx == NULL) {cerr <<"No CNT"<<endl;}
	//------------------------------------------------------
    }

    //------------------Vertex------------------------------
    d_vtxout = findNode::getClass<VtxOut>(topNode,"VtxOut");
    if(d_vtxout == NULL) { cerr << "VtxOut node not found." << endl; return DISCARDEVENT;}
    //-----------------------------------------------------

    if ((dofullana || dotrackana) && tracktype != GA::CNT)
    {
	d_svxtrx = findNode::getClass<SvxCentralTrackList>(topNode, "SvxCentralTrackList");
	if(d_svxtrx == NULL) { cerr << PHWHERE << " SvxCentralTrackList node not found." << endl; return DISCARDEVENT; }
    }

    d_rp = getClass<ReactionPlaneObject>(topNode, "ReactionPlaneObject");
    if(d_rp == NULL) { cerr << PHWHERE << " ReactionPlaneObject node not found." << endl; return DISCARDEVENT; }

    if (dofullana || domistagana)
    {
	d_emccontain  = findNode::getClass<emcClusterContainer>(topNode, "emcClusterContainer");
	if(d_emccontain==NULL) { cerr<<PHWHERE<<" emcClusterContainer node not found"<<endl; return DISCARDEVENT; }
    }

    return 0;
}

//==============================================================

int gamAna::InitRun(PHCompositeNode *topNode) {

    if (verbosity > 0)
    {
	cout<<PHWHERE<<" "<<endl;
    }

    //-------------------RunHeader----------------------------------
    d_runhead = findNode::getClass<RunHeader>(topNode, "RunHeader");
    if(d_runhead == NULL) { cerr << PHWHERE << " RunHeader node not found." << endl; return DISCARDEVENT; }
    //-------------------------------------------------------

    h1tracktype->Fill((int)tracktype);

    runno=d_runhead->get_RunNumber();

    char cmd[200];
    sprintf(cmd, "select eventsinrun from run where runnumber=%d", runno);
    TSQLServer *db = TSQLServer::Connect("odbcn://daq", "phnxrc", "");
    if (verbosity > 0)
        cout << PHWHERE <<" DB: " << db->ServerInfo() << endl;
    TSQLStatement *stmt = db->Statement(cmd, 1);
    unsigned int events_in_run = 0;
    if (stmt->Process())
    {
	stmt->StoreResult();
	while (stmt->NextResultRow())
	{
	    events_in_run = stmt->GetInt(0);
	}
    }
    delete stmt;
    delete db;

    h1evtinrun->Fill(0.,(double)events_in_run);

    TOAD toad_loader("zrowan");

    if (readeventcount)
    {

	if (verbosity > 0) cout<<PHWHERE<<"Reading Events per Bin with TOAD"<<endl;
	string evtct_loc = toad_loader.location("EventCount.root");
	TFile* evtct_file = new TFile(evtct_loc.c_str(),"READ");

	if (!evtct_file)
	{
	    cerr<<"Could not open input "<<evtct_loc<<endl;
	    return ABORTRUN;
	}
	if (evtct_file->IsZombie())
	{
	    cerr<<"Input "<<evtct_loc<<" is zombie."<<endl;
	    return ABORTRUN;
	}
	if (verbosity > 0)
	{
	    cout<<"Input "<<evtct_loc<<" opened successfully."<<endl;
	}


	h2bz0cent_evtmod = (TH2F*)evtct_file->Get(Form("bz0cent_%d",runno));
	if (h2bz0cent_evtmod == NULL)
	{
	    printf("Could not find bz0cent_%d\n",runno);
	    cout<<"ABORTRUN!"<<endl;
	    return ABORTRUN;
	}
	if (verbosity > 0)
	    printf("bz0cent_%d",runno);

	h2bz0cent_evtmod->SetDirectory(0);

	evtct_file->Close();
	delete evtct_file;
    }

    //-------------------- EMC Hot Maps (Run Specific) --------------------

    if (verbosity > 0) cout<<PHWHERE<<"Reading EMC Run Hot Maps with TOAD"<<endl;
    string runmap_location = toad_loader.location("run14_emchotmaps.root");
    TFile* runmap_file = new TFile(runmap_location.c_str(),"READ");

    if (!runmap_file)
    {
	cerr<<"Could not open input "<<runmap_location<<endl;
	return ABORTRUN;
    }
    if (runmap_file->IsZombie())
    {
	cerr<<"Input "<<runmap_location<<" is zombie."<<endl;
	return ABORTRUN;
    }

    if (verbosity > 0)
    {
	cout<<"Input "<<runmap_location<<" opened successfully."<<endl;
    }

    for (int isec = 0; isec < GA::nsec; isec++)
    {
	h2runmap[isec] = (TH2I*)runmap_file->Get(Form("hotmap_run%d_sec%d",runno,isec));
	if (h2runmap[isec] == NULL)
	{
	    cerr<<"Could not find runmap. Will be returning ABORTEVENT in process_event"<<endl;
	    continue;
	}
	if (verbosity > 0)
	    printf("Retrieved hotmap_run%d_sec%d\n",runno,isec);

	h2runmap[isec]->SetDirectory(0);
    }

    runmap_file->Close();
    delete runmap_file;

    return 0;
}
//==============================================================

int gamAna::process_event(PHCompositeNode *topNode) {

    if((evtno%1000)==0 && verbosity > 0)
    {
	cout<<"event processed# "<< evtno<<"  "<<"run ="<<" "<<runno<<endl;
	TTimeStamp* tts = new TTimeStamp();
	tts->Print();
	delete tts;
    }

    evtno++;
    h1evtraw->Fill(0);

    for (int isec = 0; isec < GA::nsec; isec++)
    {
	if (h2runmap[isec] == NULL)
	{
	    return ABORTEVENT;
	}
    }

    if( GetNodes(topNode) == DISCARDEVENT) return DISCARDEVENT;

    for( int i = 0; i < 3; i++ ) pticks[i] = (d_peve!=NULL) ? d_peve->get_clockticks(i) : -999;
    if( TickCut() != 0 )
    {
	h1tickcut->Fill(1);
	return EVENT_OK;
    }
    else
	h1tickcut->Fill(0);

    TriggerHelper trighelp(topNode);

    isBBC(&trighelp);
    isERT(&trighelp);

    // ====================================================================
    // --------------- Determine event type for tree buffer ---------------
    // ====================================================================

    // ------------------------------
    // Z Vertex
    // ------------------------------
    float zbbc = d_vtxout->get_ZVertex("BBC");
    z0 = d_vtxout->get_ZVertex();
    h2z0zbbc->Fill(z0,zbbc);

    if (fabs(z0) > GA::maxzv) return 0;

    if (tracktype == GA::SVX || tracktype == GA::SVXCNT)
    {
	cout<<"Svx Cnt Tracks no longer supported! Aborting Run!"<<endl;
	return ABORTRUN;
    }


    bz0 = GA::GetZvBin(z0);
    if (bz0 < 0) return 0;

    // ------------------------------
    // Centrality
    // ------------------------------
    centrality=d_global->getCentrality();
    h1cent_raw->Fill(centrality);

    if (centrality < 1. || centrality > GA::maxcent)
    {
	float bbcq = d_global->getBbcChargeN()+d_global->getBbcChargeS();
	if (verbosity > 1)
	    cout<<PHWHERE<<"Unknown centrality: "<<centrality<<" bbcq: " << bbcq << endl;
	return 0;
    }

    bcent = GA::GetCentBin(centrality);
    if(bcent < 0) return 0;

    // ------------------------------
    // Reaction Plane
    // ------------------------------

    // if BBC
    if (rpdet == GA::BBC)
    {
	ReactionPlaneSngl *rpsngl = NULL;

	// BBC South
	rpsngl = d_rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 0, 1));
	if (!rpsngl)
	{
	    if (verbosity>0)
	    {
		cerr<<"Could not find bbc south rp."<<endl;
	    }
	    return DISCARDEVENT;
	}
	float psi_bbcs = rpsngl->GetPsi();
	if (psi_bbcs < 0.) psi_bbcs+=PI;
	if (psi_bbcs >= PI) psi_bbcs-=PI;
	h1psi_bbcs->Fill(psi_bbcs/PI);

	// BBC North
	rpsngl = d_rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 1, 1));
	if (!rpsngl)
	{
	    if (verbosity>0)
	    {
		cerr<<"Could not find bbc north rp."<<endl;
	    }
	    return DISCARDEVENT;
	}
	float psi_bbcn = rpsngl->GetPsi();
	if (psi_bbcn < 0.) psi_bbcn+=PI;
	if (psi_bbcn >= PI) psi_bbcn-=PI;
	h1psi_bbcn->Fill(psi_bbcn/PI);
	
	// BBC Sum
	rpsngl = d_rp->getReactionPlane(RP::calcIdCode(RP::ID_BBC, 2, 1));
	if (!rpsngl)
	{
	    if (verbosity>0)
	    {
		cerr<<"Could not find bbc sum rp."<<endl;
	    }
	    return DISCARDEVENT;
	}
	float psi_bbc = rpsngl->GetPsi();
	h1psi_bbc_raw->Fill(psi_bbc/PI);
	if (psi_bbc < 0.) psi_bbc+=PI;
	if (psi_bbc >= PI) psi_bbc-=PI;
	
	if (verbosity > 2)
	    rpsngl->Print();

	h2centpsidiff_bbcns->Fill(centrality,(psi_bbcn-psi_bbcs)/PI);
	h2psi_bbc_bbcn->Fill(psi_bbc/PI,psi_bbcn/PI);
	h2psi_bbc_bbcs->Fill(psi_bbc/PI,psi_bbcs/PI);
	h2psi_bbcn_bbcs->Fill(psi_bbcn/PI,psi_bbcs/PI);

	psi = psi_bbc;
    }

    // if FVTX
    else if (rpdet == GA::FVTX)
    {
	ReactionPlaneSngl *rpsngl = NULL;

	// FVTX South
	rpsngl = d_rp->getReactionPlane(RP::calcIdCode(RP::ID_FVT, 40, 1));
	if (!rpsngl)
	{
	    if (verbosity>0)
	    {
		cerr<<"Could not find fvtx south rp."<<endl;
	    }
	    return DISCARDEVENT;
	}
	float psi_fvtxs = rpsngl->GetPsi();
	if (psi_fvtxs < 0.) psi_fvtxs+=PI;
	if (psi_fvtxs >= PI) psi_fvtxs-=PI;
	h1psi_fvtxs->Fill(psi_fvtxs/PI);

	// FVTX North
	rpsngl = d_rp->getReactionPlane(RP::calcIdCode(RP::ID_FVT, 41, 1));
	if (!rpsngl)
	{
	    if (verbosity>0)
	    {
		cerr<<"Could not find fvtx north rp."<<endl;
	    }
	    return DISCARDEVENT;
	}
	float psi_fvtxn = rpsngl->GetPsi();
	if (psi_fvtxn < 0.) psi_fvtxn+=PI;
	if (psi_fvtxn >= PI) psi_fvtxn-=PI;
	h1psi_fvtxn->Fill(psi_fvtxn/PI);
	
	// FVTX Sum
	rpsngl = d_rp->getReactionPlane(RP::calcIdCode(RP::ID_FVT, 42, 1));
	if (!rpsngl)
	{
	    if (verbosity>0)
	    {
		cerr<<"Could not find fvtx sum rp."<<endl;
	    }
	    return DISCARDEVENT;
	}
	float psi_fvtx = rpsngl->GetPsi();
	h1psi_fvtx_raw->Fill(psi_fvtx/PI);
	if (psi_fvtx < 0.) psi_fvtx+=PI;
	if (psi_fvtx >= PI) psi_fvtx-=PI;
	
	if (verbosity > 2)
	    rpsngl->Print();

	h2centpsidiff_fvtxns->Fill(centrality,(psi_fvtxn-psi_fvtxs)/PI);
	h2psi_fvtx_fvtxn->Fill(psi_fvtx/PI,psi_fvtxn/PI);
	h2psi_fvtx_fvtxs->Fill(psi_fvtx/PI,psi_fvtxs/PI);
	h2psi_fvtxn_fvtxs->Fill(psi_fvtxn/PI,psi_fvtxs/PI);

	psi = psi_fvtx;
    }

    else
    {
	cerr<<"Invalid rp detector id"<<endl;
	return ABORTRUN;
    }

    h1psi->Fill(psi/PI);
	
    if (psi < 0. || psi >= PI)
    {
    	if (verbosity > 0)
    	    cout<<PHWHERE<<"Psi OOR: "<<psi<<endl;
    	return 0;
    }

    bpsi = GA::GetPsiBin(psi);
    h1bpsi->Fill(bpsi);
    if (bpsi < 0) return 0;

    //-------------------------------------------------------------
    //  Done Determining Centrality and RP Classes
    //-------------------------------------------------------------

    goodevtno++;

    h1z0->Fill(z0);
    h2bz0cent->Fill(bz0,centrality);
    h1bcent->Fill(bcent);

    if (dofullana || domistagana)
    {
	nclus = d_emccontain->size();
	h2centnclus->Fill(centrality,nclus);
    }

    if (domistagana)
    {
	return MistagAna(true);
    }

    if (dofullana || dotrackana)
    {
	float ncnt = d_trx->get_npart();
	h2centnpart->Fill(centrality,ncnt);

	if (tracktype == GA::CNT){
	    ntrk = ncnt;
	}
	if (tracktype == GA::SVXCNT){
	    cout<<"Svx Cnt Tracks no longer supported! Aborting Run!"<<endl;
	    return ABORTRUN;
	}
    }

    bool saveformix = false;
    bool saveevt = false;
    if (dofullana)
    {
	vector<int> savedIDs;

	float emin_pi0gam = (ptmin_pi0trig/2.) * (1. - asymmax);
	float emin_gampal = ptmin_gamtrig * (1. - asymmax)/(1. + asymmax);
	float emin = (emin_pi0gam < emin_gampal) ? emin_pi0gam : emin_gampal;
	if (emin < emin_clus) emin = emin_clus;

	int ntrigperevt[3] = {0,0,0};

	std::vector <EmcCluster*> myclusters;

	for (int iclus = 0; iclus < nclus; iclus++)
	{
	    EmcCluster* myclus_i = new EmcCluster();
	    myclusters.push_back(myclus_i);
	    emcClusterContent* clus_i = d_emccontain->getCluster(iclus);

	    myclus_i->isok = isClusterOK(clus_i,myclus_i,emin,true);
	    if (!myclus_i->isok) continue;

	    float e_i = clus_i->ecore();
	    myclus_i->CalcVars(e_i,clus_i->x(),clus_i->y(),clus_i->z(),z0);
	    float ex_i = myclus_i->ex;
	    float ey_i = myclus_i->ey;
	    float ez_i = myclus_i->ez;
	    float pt_i = myclus_i->pt;

	    float tmpPhi = clus_i->phi()/PI;
	    h2armphi_emclus->Fill(myclus_i->arm,tmpPhi);

	    if (tmpPhi < -.5) tmpPhi+=2.;
	    h2zphi_arm[myclus_i->arm]->Fill(clus_i->z(),tmpPhi);

	    h1r[myclus_i->sector]->Fill(myclus_i->r);
	    h1s[myclus_i->sector]->Fill(myclus_i->s);
	    h1ex->Fill(ex_i);
	    h1ey->Fill(ey_i);
	    h1ez->Fill(ez_i);
	    h2ptcent->Fill(pt_i,centrality);
	    h2ecent->Fill(e_i,centrality);

	    if (!dofullana) continue;

	    if (pt_i >= ptmin_gamtrig)
	    {
		if (verbosity > 1)
		{
		    cout<<"Saving photon "<<iclus<<" with pt = "<<pt_i<<" from Event #"<<evtno<<endl;
		}
		savedIDs.push_back(iclus);
	    }

	    for (int jclus = 0; jclus < iclus; jclus++)
	    {
		EmcCluster* myclus_j = myclusters[jclus];
		if(!myclus_j->isok) continue;

		emcClusterContent* clus_j = d_emccontain->getCluster(jclus);

		if (myclus_i->arm != myclus_j->arm) continue;

		float e_j = clus_j->ecore();

		if ((e_i + e_j) < ptmin_pi0trig) continue;

		float asym = fabs(e_i - e_j) / (e_i + e_j);
		if (asym > asymmax) continue;

		float ex_j = myclus_j->ex;
		float ey_j = myclus_j->ey;
		float ez_j = myclus_j->ez;

		float px_pair = ex_i+ex_j;
		float py_pair = ey_i+ey_j;
		float pz_pair = ez_i+ez_j;
		float pt_pair = sqrt(px_pair*px_pair + py_pair*py_pair);

		float invm = sqrt(2.0*e_i*e_j - 2.0*(ex_i*ex_j+ey_i*ey_j+ez_i*ez_j));

		int histosec = GetHistoSector(myclus_i->sector,myclus_j->sector);
		h1histosec->Fill(histosec);
		if (histosec < 0 || histosec >= GA::nsec + 6)
		{
		    if (verbosity > 0)
		    {
			cout<<PHWHERE<<"WARNING!! GetHistoSector returned invalid sector number = "<<histosec<<"; This may be due to low pair energy. e_i = "<<e_i<<". e_j = "<<+e_j<<". asym = "<<asym<<endl;
		    }

		}
		else
		    h2mpt[histosec]->Fill(invm,pt_pair);

		float mu = (pi0pars[myclus_i->sector][0]+pi0pars[myclus_j->sector][0])/2.;
		float sigma = (pi0pars[myclus_i->sector][1]+pi0pars[myclus_j->sector][1])/2.;

		if (fabs(invm - mu) < npi0sigma*sigma)
		{
		    if (asym < myclus_i->lowasym) myclus_i->lowasym = asym;
		    if (asym < myclus_j->lowasym) myclus_j->lowasym = asym;

		    myclus_i->isdecay = true;
		    myclus_j->isdecay = true;

		    if (verbosity > 1)
		    {
			cout<<"Photon "<<iclus<<" with E = "<<e_i
			    <<" and photon "<<jclus<<" with E = "<<e_j
			    <<" from Event #"<<evtno<<" are decays."<<endl;
		    }

		    if (pt_pair < ptmin_pi0trig) continue;

		    ntrigperevt[2]++;
		    h1asym_pi0->Fill(asym);
		    h1e_spec[2]->Fill(e_j+e_i);

		    trig = new TriggerParticle(2,px_pair,py_pair,pz_pair);
		    trigtree->Fill();
		    delete trig;

		}
	    }
	}

	if (!dofullana) return EVENT_OK;

	h1ncand->Fill(0.,savedIDs.size());

	for (vector<int>::iterator id = savedIDs.begin(); id != savedIDs.end(); ++id)
	{
	    emcClusterContent* clus = d_emccontain->getCluster(*id);
	    EmcCluster* myclus = myclusters[*id];

	    if (myclus->isdecay)
	    {
		h1lowasym->Fill(myclus->lowasym);
	    }

	    float e = clus->ecore();
	    float ex = myclus->ex;
	    float ey = myclus->ey;
	    float ez = myclus->ez;
	    float s = myclus->s;

	    ntrigperevt[myclus->isdecay]++;

	    h1e_spec[myclus->isdecay]->Fill(e);

	    trig = new TriggerParticle(myclus->isdecay,ex,ey,ez);

	    float emcpc3dphis = clus->emcpc3dphi() * s;

	    float emctrkdphis = clus->emctrkdphi() * s;

	    float emcpc3dz = clus->emcpc3dz();

	    float emctrkdz = clus->emctrkdz();

	    trig->SetPc3Veto(emcpc3dphis,emcpc3dz);
	    trig->SetTrkVeto(emctrkdphis,emctrkdz);

	    int towid = -9999;
	    GA::TowerFromSYZ(towid,myclus->sector,clus->iypos(),clus->izpos());
	    trig->towid = towid;

	    trigtree->Fill();
	    delete trig;

	}

	for (std::vector<EmcCluster*>::iterator it = myclusters.begin(); it != myclusters.end(); ++it)
	    delete *it;

	bool trigfound = false;
	for (int ispec = 0; ispec < nspec; ispec++)
	{
	    if (ntrigperevt[ispec] > 0) trigfound = true;
	    h1ntrigperevt[ispec]->Fill(ntrigperevt[ispec]);
	}

	int nevtest = h2bz0cent_evtmod->GetBinContent(bz0+1,(int)centrality+1)/((float)npsimixbin) + 0.5;
	int evtmod = nevtest/nevtdesire;
	if (evtmod<1) evtmod=1;
	h1evtmod->Fill(evtmod);
	saveformix = goodevtno%evtmod==0;
	if (saveformix)
	{
	    int nsaveest = nevtest/evtmod;
	    if (nsaveest > nevtdesire)
	    {
		int skipevtmod = nsaveest/(nsaveest-nevtdesire) + 1;
		h1skipevtmod->Fill(skipevtmod);
		saveformix = goodevtno%(evtmod*skipevtmod)!=0;
	    }
	}
	saveevt = trigfound || saveformix;
    }

    if (saveallevt) saveevt = true;

    if (saveevt)
    {
	if (tracktype == GA::SVXCNT) {
	    cout<<"Svx Cnt Tracks no longer supported! Aborting Run!"<<endl;
	    return ABORTRUN;
	}
	if (tracktype == GA::CNT) trkhold = new TrackCollection<CentralArmTrack>(saveformix);
    }

    if (dotrackana || saveevt)
    {
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
	    if (verbosity > 2)
	    {
		cout<<"Making my track "<<itrk<<" of "<<ntrk-1<<endl;
	    }
	    bool filltrkhists = true;
	    ChargedTrack* mytrk = NULL;
	    bool trkok = makeChargedTrack(&mytrk,itrk,filltrkhists);

	    if (!trkok)
	    {
		if (verbosity > 2)
		{
		    cout<<"Track is not OK. Deleting."<<endl;
		}
		delete mytrk;
		continue;
	    }

	    if (saveevt)
	    {
		if (verbosity > 2)
		{
		    cout<<"Adding track to track holder."<<endl;
		}
		trkhold->AddEntry(mytrk);
	    }

	    delete mytrk;

	}

	if (saveevt)
	{
	    trktree[bcent][bpsi][bz0]->Fill();
	    delete trkhold;
	}
    }

    return EVENT_OK;
}

int gamAna::End(PHCompositeNode *topNode) {
    if (verbosity > 0) cout << "This is the end..." << endl;

    for (int isec = 0; isec < GA::nsec; isec++)
    {
	if (h2emchdmap[isec] != NULL)
	    delete h2emchdmap[isec];

	if (h2runmap[isec] != NULL)
	    delete h2runmap[isec];

	if (h2commonmap[isec] != NULL)
	    delete h2commonmap[isec];

    }

    delete h2bz0cent_evtmod;

    if (verbosity > 0) cout << "Writing out..." << endl;
    outfile->Write();
    if (verbosity > 0) cout << "Closing output file..." << endl;
    outfile->Close();
    delete outfile;

    return 0;
}

//==============================================================

int gamAna::GetHistoSector(int sector_i, int sector_j)
{
    if (sector_i == sector_j) return sector_i;

    int arm_i = sector_i/4;
    int sec_i = sector_i - arm_i*4;

    int arm_j = sector_j/4;
    int sec_j = sector_j - arm_j*4;

    if (arm_i == arm_j)
    {
	if (sec_i == sec_j - 1)
	    return sec_i + GA::nsec;
	else if (sec_j == sec_i - 1)
	    return sec_j + GA::nsec;
    }

    if (verbosity > 0)
	cout<<PHWHERE<<"GetHistoSector: returning -1; sector_i = "<<sector_i<<" arm_i = "<<arm_i<<" sec_i = "<<sec_i<<";  sector_j = "<<sector_j<<" arm_j = "<<arm_j<<" sec_j = "<<sec_j<<endl;

    return -1;
}

//==============================================================
// Consolidating Track Cuts
//==============================================================
bool gamAna::makeChargedTrack(ChargedTrack** mytrk, int itrk, bool fillhists)
{
    if (tracktype == GA::CNT)
    {
	*mytrk = new CentralArmTrack();
	CentralArmTrack* mycnttrk = dynamic_cast<CentralArmTrack*>(*mytrk);

	float mom = d_trx->get_mom(itrk);
	float the0 = d_trx->get_the0(itrk);
	float phi0 = d_trx->get_phi0(itrk);
	(*mytrk)->Set3MomFromMomThePhi(mom,the0,phi0);

	if (!(*mytrk)->MomOk())
	{
	    if (verbosity > 2)
	    {
		cout<<"Track momentum out of range"<<endl;
		(*mytrk)->Print();
	    }
	    return false;
	}

	float pt = (*mytrk)->GetPt();

	CentralTrackData* mycnt = mycnttrk->GetCentralTrackData();
	bool is_svxcnt = false;
	bool is_mycnt_ok = makeCentralTrackData(mycnt,itrk,is_svxcnt,fillhists,pt);
	if (!is_mycnt_ok) return false;

	if (verbosity > 2)
	{
	    cout<<"Cnt Track has passed all online cuts"<<endl;
	    (*mytrk)->Print();
	}
    }
    else if (tracktype == GA::SVXCNT)
    {
	cerr<<"Svx Cnt Tracks no longer supported! Aborting Run!"<<endl;
	return false;
    }
    else
    {
	cerr<<PHWHERE<<"Invalid track type. Returning false. (SVX not yet supported.)"<<endl;
	return false;
    }

    if (fillhists)
    {
	if (verbosity > 2)
	{
	    cout<<"Filling track momentum histos."<<endl;
	}

	h1px->Fill((*mytrk)->px);
	h1py->Fill((*mytrk)->py);
	h1pz->Fill((*mytrk)->pz);
    }

    return true;
}

bool gamAna::makeCentralTrackData(CentralTrackData* mycnt, int icnt, bool is_svxcnt, bool fillhists, float pt)
{
    if (verbosity > 2)
    {
	cout<<"Making My Cnt."<<endl;
    }

    short dcarm = d_trx->get_dcarm(icnt);
    if (dcarm != 0 && dcarm != 1)
    {
	if (verbosity > 2)
	{
	    cout<<"dcarm ="<<dcarm<<" not ok. Will not save track."<<endl;
	}
	return false;
    }
    mycnt->arm = (bool)dcarm;


    if (fillhists)
	h2ptcent_trk[trkcut_map["RAW"]]->Fill(pt,centrality);


    float zed = d_trx->get_zed(icnt);
    if (fillhists)
    {
	h1zed_trk_raw->Fill(zed);
    }
    if (!mycnt->ZedOk(zed))
    {
	if (verbosity > 2)
	{
	    cout<<"zed = "<<zed<<" not ok. Will not save track."<<endl;
	}
	return false;
    }

    if (fillhists)
	h2ptcent_trk[trkcut_map["ZED"]]->Fill(pt,centrality);


    int quality = d_trx->get_quality(icnt);
    if (!mycnt->QualOk(quality))
    {
	if (verbosity > 2)
	{
	    cout<<"quality ="<<quality<<" not ok. Will not save track."<<endl;
	}
	return false;
    }

    if (fillhists)
	h2ptcent_trk[trkcut_map["QUAL"]]->Fill(pt,centrality);


    //PC match
    float pc3dphi = d_trx->get_pc3dphi(icnt);
    float pc3dz = d_trx->get_pc3dz(icnt);

    float pc3sdphi = d_trx->get_pc3sdphi(icnt);
    float pc3sdz = d_trx->get_pc3sdz(icnt);
    float pc3sig_f = sqrt(pow(pc3sdphi,2) + pow(pc3sdz,2));
    mycnt->SetPc3Sig(pc3sig_f);

    if (verbosity > 2)
    {
	printf("Pc3 Info: pc3dphi %.3f; pc3dz %.3f; pc3sdphi %.3f; pc3sdz %.3f; pc3sig_f %.3f; pc3sigchar %d;\n",pc3dphi,pc3dz,pc3sdphi,pc3sdz,pc3sig_f,mycnt->GetPc3SigChar());
    }

    if (!is_svxcnt && !mycnt->Pc3Ok())
    {
	if (verbosity > 2)
	{
	    cout<<"PC3 not ok. Will not save track."<<endl;
	}
	return false;
    }

    if (fillhists)
	h2ptcent_trk[trkcut_map["PC3"]]->Fill(pt,centrality);

    //EMC match
    float emcdphi = d_trx->get_emcdphi(icnt);
    float emcdz = d_trx->get_emcdz(icnt);

    float emcsdphi = d_trx->get_emcsdphi(icnt);
    float emcsdz = d_trx->get_emcsdz(icnt);
    float emcsig_f = sqrt(pow(emcsdphi,2) + pow(emcsdz,2));
    mycnt->SetEmcSig(emcsig_f);

    if (verbosity > 2)
    {
	printf("Emc Info: emcdphi %.3f; emcdz %.3f; emcsdphi %.3f; emcsdz %.3f; emcsig_f %.3f; emcsigchar %d;\n",emcdphi,emcdz,emcsdphi,emcsdz,emcsig_f,mycnt->GetEmcSigChar());
    }

    if (!is_svxcnt && !mycnt->EmcOk())
    {
	if (verbosity > 2)
	{
	    cout<<"EMC not ok. Will not save track."<<endl;
	}
	return false;
    }

    if (fillhists)
	h2ptcent_trk[trkcut_map["EMC"]]->Fill(pt,centrality);

    float n0 = d_trx->get_n0(icnt);
    mycnt->isn0neg = (n0 < 0.);

    if (verbosity > 2)
    {
	cout<<"n0 = "<<n0<<"; isn0neg = "<<mycnt->isn0neg<<endl;
    }

    float dep = d_trx->get_dep(icnt);
    mycnt->SetDEP(dep);

    if (verbosity > 2)
    {
	cout<<"dep = "<<dep<<"; depchar = "<<(int)mycnt->GetDEPChar()<<endl;
    }

    if (fillhists)
    {
	if (mycnt->isn0neg && pt<5.)
	{
	    h2ptcent_trk[trkcut_map["N0"]]->Fill(pt,centrality);

	    if (mycnt->CheckDEP(-3.))
	    {
		h2ptcent_trk[trkcut_map["DEP"]]->Fill(pt,centrality);

		if (verbosity > 2)
		{
		    cout<<"mycnt->CheckDEP(-3.) = true"<<endl;
		}
	    }
	    else if (verbosity > 2)
	    {
		cout<<"mycnt->CheckDEP(-3.) = false"<<endl;
	    }

	    if (verbosity > 2)
	    {
		cout<<"dep = "<<dep<<"; depchar = "<<(int)mycnt->GetDEPChar()<<endl;
	    }
	}
    }

    float ecore = d_trx->get_ecore(icnt);
    float eop = ecore/d_trx->get_mom(icnt);

    if (fillhists)
    {
	if (verbosity > 2)
	{
	    cout<<"Filling track histos."<<endl;
	}

	h1dcarm->Fill(dcarm);
	h1quality->Fill(quality);
	h1charge->Fill(d_trx->get_charge(icnt));
	h1pc3dphi->Fill(pc3dphi);
	h1pc3dz->Fill(pc3dz);
	h1pc3sdphi->Fill(pc3sdphi);
	h1pc3sdz->Fill(pc3sdz);
	h1pc3sig->Fill(pc3sig_f);
	h1emcdphi->Fill(emcdphi);
	h1emcdz->Fill(emcdz);
	h1emcsdphi->Fill(emcsdphi);
	h1emcsdz->Fill(emcsdz);
	h1emcsig->Fill(emcsig_f);
	h1n0->Fill(n0);
	h1ecore->Fill(ecore);
	h1eop->Fill(eop);
	h1dep->Fill(dep);

	if (verbosity > 2)
	{
	    cout<<"Track histos filled."<<endl;
	}
    }

    if (verbosity > 2)
    {
	mycnt->Print();
    }

    return true;
}

//==============================================================
// Consolidating EMC Cluster Cuts
//==============================================================
bool gamAna::isClusterOK(emcClusterContent* clus, EmcCluster* myclus, float emin, bool fillhists)
{

    float e = clus->ecore();
    float prob = clus->prob_photon();

    if (fillhists)
    {
	h1e_nocut->Fill(e);
	h1prob_nocut->Fill(prob);
    }

    if (e < emin) return false;
    if (prob < .02) return false;

    int sector = clus->sector();
    if (sector < 0 || sector > 3)
    {
	cerr<<PHWHERE<<"Event " << evtno << " WARNING! emc sector is out of range -> "<<sector<<endl;
	return false;
    }

    int arm = clus->arm();
    if (arm < 0 || arm > 1)
    {
	cerr<<PHWHERE<<"Event " << evtno << " WARNING! emc arm is out of range -> "<<arm<<endl;
	return false;
    }

    myclus->SetSector(arm,sector);
    int sec = myclus->sector;
    int iypos = clus->iypos();
    int izpos = clus->izpos();

    if (fillhists)
    {
	h2iziy_nomapcut[sec]->SetBinContent(izpos+1,iypos+1,1);
    }

    if (sec > 1)
    {
	if (iypos >= GA::nySc || iypos < 0)
	{
	    cerr<<PHWHERE<<"Event " << evtno << " WARNING! emc Sc iypos is out of range -> "<<iypos<<endl;
	    return false;
	}

	if (iypos == GA::nySc - 1 || iypos == 0) return false;

	if (izpos >= GA::nzSc || izpos < 0)
	{
	    cerr<<PHWHERE<<"Event " << evtno << " WARNING! emc Sc izpos is out of range -> "<<izpos<<endl;
	    return false;
	}

	if (izpos == GA::nzSc - 1 || izpos == 0) return false;
    }
    else
    {
	if (iypos >= GA::nyGl || iypos < 0)
	{
	    cerr<<PHWHERE<<"Event " << evtno << " WARNING! emc Gl iypos is out of range -> "<<iypos<<endl;
	    return false;
	}

	if (iypos == GA::nyGl - 1 || iypos == 0) return false;

	if (izpos >= GA::nzGl || izpos < 0)
	{
	    cerr<<PHWHERE<<"Event " << evtno << " WARNING! emc Gl izpos is out of range -> "<<izpos<<endl;
	    return false;
	}

	if (izpos == GA::nzGl - 1 || izpos == 0) return false;
    }
    if (fillhists)
    {
	h2iziy_bordwu[sec]->SetBinContent(izpos+1,iypos+1,1);
    }

    // --------------- warn and uncal maps ---------------
    if (emcbad[sec][iypos][izpos])
    {
	if (verbosity > 2)
	{
	    cout<<"EMCal Cluster from sec = "<<sec<<"; iy = "<<iypos<<"; iz = "<<izpos
		<<"; did not pass emc warn uncal"<<endl;
	}
	h2iziy_emcbad[sec]->SetBinContent(izpos+1,iypos+1,1);
	return false;
    }
    if (useemc3x3 && emcnearbad[sec][iypos][izpos])
    {
	if (verbosity > 2)
	{
	    cout<<"EMCal Cluster from sec = "<<sec<<"; iy = "<<iypos<<"; iz = "<<izpos
		<<"; did not pass emc warn uncal 3x3"<<endl;
	}
	return false;
    }
    if (usemyemc)
    {
	if (h2commonmap[sec]->GetBinContent(izpos+1,iypos+1) > GA::maxhotrun)
	{
	    if (verbosity > 2)
	    {
		cout<<"EMCal Cluster from sec = "<<sec<<"; iy = "<<iypos<<"; iz = "<<izpos
		    <<"; did not pass emc common hot map"<<endl;
	    }
	    return false;
	}
	if (h2runmap[sec]->GetBinContent(izpos+1,iypos+1) > 0)
	{
	    if (verbosity > 2)
	    {
		cout<<"EMCal Cluster from sec = "<<sec<<"; iy = "<<iypos<<"; iz = "<<izpos
		    <<"; did not pass emc run specific hot map"<<endl;
	    }
	    return false;
	}
    }
    if (fillhists)
    {
	h2iziy_bordwu[sec]->SetBinContent(izpos+1,iypos+1,1);
    }

    return true;
}

bool gamAna::TickCut()
{
    return ( ( 50<pticks[0]&&pticks[0]<120)||
    	     (700<pticks[1]&&pticks[1]<780) );
}

bool gamAna::isBBC(TriggerHelper* trighelp)
{
    bool bbc = trighelp->trigScaled( "BBCLL1(>1 tubes)" );

    bool bbcnar = (
	trighelp->trigScaled( "BBCLL1(>1 tubes) narrowvtx" ) ||
	trighelp->trigScaled( "BBCLL1(>1 tubes) narrowvtx CopyA" ) ||
	trighelp->trigScaled( "BBCLL1(>1 tubes) narrowvtx CopyB" )
	);

    bool bbcnov = trighelp->trigScaled( "BBCLL1(>1 tubes) novertex" );

    bool bbcall = bbc || bbcnov || bbcnar;

    h1BBCall->Fill(bbcall);
    h1BBC->Fill(bbc);
    h1BBCnar->Fill(bbcnar);
    h1BBCnov->Fill(bbcnov);

    return bbcall;
}

bool gamAna::isERT(TriggerHelper* trighelp)
{
    bool ertbbc = (
	trighelp->trigScaled( "ERT_4x4a&BBCLL1" ) ||
	trighelp->trigScaled( "ERT_4x4b&BBCLL1" ) ||
	trighelp->trigScaled( "ERT_4x4c&BBCLL1" )
	);

    bool ert = (
	trighelp->trigScaled( "ERT_2x2" ) ||
	trighelp->trigScaled( "ERT_4x4b" )
	);

    bool ertall = ertbbc || ert;

    h1ERTall->Fill(ertall);
    h1ERT->Fill(ert);
    h1ERTBBC->Fill(ertbbc);

    return ertall;
}

int gamAna::MistagAna(bool plotlowephotons)
{
    float emin_i = plotlowephotons ? emin_clus : ptmin_gamtrig;

    for (int iclus = 0; iclus < nclus; iclus++)
    {
	EmcCluster* myclus_i = new EmcCluster();
	emcClusterContent* clus_i = d_emccontain->getCluster(iclus);

	myclus_i->isok = isClusterOK(clus_i,myclus_i,emin_i,true);
	if (!myclus_i->isok)
	{
	    delete myclus_i;
	    continue;
	}

	float e_i = clus_i->ecore();
	myclus_i->CalcVars(e_i,clus_i->x(),clus_i->y(),clus_i->z(),z0);
	float ex_i = myclus_i->ex;
	float ey_i = myclus_i->ey;
	float ez_i = myclus_i->ez;
	float pt_i = myclus_i->pt;

	int sec_i = myclus_i->sector;
	int iy_i = clus_i->iypos();
	int iz_i = clus_i->izpos();

	h2ptcent->Fill(pt_i,centrality);
	h2ecent->Fill(e_i,centrality);
	h1ex->Fill(ex_i);
	h1ey->Fill(ey_i);
	h1ez->Fill(ez_i);
	h1r[sec_i]->Fill(myclus_i->r);
	h1s[sec_i]->Fill(myclus_i->s);

	if (pt_i < ptmin_gamtrig)
	{
	    delete myclus_i;
	    continue;
	}
	if (!GA::PassFidCut(e_i,sec_i,iy_i,iz_i))
	{
	    delete myclus_i;
	    continue;
	}

	int trigcat = 0;
	if (!emcnearbad[sec_i][iy_i][iz_i])
	{
	    trigcat++;

	    if (h2runmap[sec_i]->GetBinContent(iz_i+1,iy_i+1) == 0 && h2commonmap[sec_i]->GetBinContent(iz_i+1,iy_i+1) <= GA::maxhotrun)
	    {
		trigcat++;

		//==============================
		//        PC3 MATCHING
		//==============================
		float pc3dphis = clus_i->emcpc3dphi() * myclus_i->s;
		float pc3dz = clus_i->emcpc3dz();

		int pc3zbin = GA::GetTowZBin(sec_i,iz_i);
		int pc3ybin = GA::GetTowYBin(sec_i,iy_i);

		float pc3dphis_mean = pc3dphispar[sec_i][pc3ybin][GA::IMEAN_PC3];
		float pc3dphis_sigma = pc3dphispar[sec_i][pc3ybin][GA::ISIGMA_PC3];
		float pc3dphis_s = (pc3dphis - pc3dphis_mean)/pc3dphis_sigma;

		float pc3dz_mean = pc3dzpar[sec_i][pc3zbin][GA::IMEAN_PC3];
		float pc3dz_sigma = pc3dzpar[sec_i][pc3zbin][GA::ISIGMA_PC3];
		float pc3dz_s = (pc3dz - pc3dz_mean)/pc3dz_sigma;

		float pc3dr_s = sqrt(pc3dz_s*pc3dz_s + pc3dphis_s*pc3dphis_s);

		if (pc3dr_s >= GA::NSIGMA_PC3DR)
		{
		    trigcat++;

		    if (!pc3bad[sec_i][iy_i][iz_i])
		    {
			trigcat++;

			if (!pc3nearbad[sec_i][iy_i][iz_i])
			{
			    trigcat++;
			}
		    }
		}
	    }
	}


	for (int icat = 0; icat <= trigcat; icat++)
	{
	    h2ptcent_trigcat[icat]->Fill(pt_i,centrality);
	    h2ecent_trigcat[icat]->Fill(e_i,centrality);
	}

	float emin_j = e_i * (1. - asymmax)/(1. + asymmax);
	if (emin_j < emin_clus) emin_j = emin_clus;

	for (int jclus = 0; jclus < nclus; jclus++)
	{
	    if (iclus == jclus) continue;

	    EmcCluster* myclus_j = new EmcCluster();
	    emcClusterContent* clus_j = d_emccontain->getCluster(jclus);

	    myclus_j->isok = isClusterOK(clus_j,myclus_j,emin_j,false);
	    if (!myclus_j->isok)
	    {
		delete myclus_j;
		continue;
	    }
	    if (myclus_i->arm != myclus_j->arm)
	    {
		delete myclus_j;
		continue;
	    }

	    float e_j = clus_j->ecore();
	    myclus_j->CalcVars(e_j,clus_j->x(),clus_j->y(),clus_j->z(),z0);
	    float ex_j = myclus_j->ex;
	    float ey_j = myclus_j->ey;
	    float ez_j = myclus_j->ez;

	    float invm = sqrt(2.0*e_i*e_j - 2.0*(ex_i*ex_j+ey_i*ey_j+ez_i*ez_j));

	    int centbin = centrality;
	    int ptbin = pt_i - ptmin_gamtrig;
	    if (ptbin >= 0 && ptbin < ninvmptbin && centbin >= 0 && centbin < ninvmcentbin)
	    {
		for (int icat = 0; icat <= trigcat; icat++)
		{
		    h1invm_trigcat[icat][centbin][ptbin]->Fill(invm);
		}
	    }

	    float px_pair = ex_i+ex_j;
	    float py_pair = ey_i+ey_j;
	    float pt_pair = sqrt(px_pair*px_pair + py_pair*py_pair);

	    int histosec = GetHistoSector(sec_i,myclus_j->sector);
	    h1histosec->Fill(histosec);
	    if (histosec < 0 || histosec >= GA::nsec + 6)
	    {
		if (verbosity > 0)
		    cout<<PHWHERE<<"WARNING!! GetHistoSector returned invalid sector number = "<<histosec<<"; This may be due to low pi0 energy = "<<e_i+e_j<<endl;
	    }
	    else
		h2mpt[histosec]->Fill(invm,pt_pair);

	    delete myclus_j;
	}
	delete myclus_i;
    }
    return EVENT_OK;
}

void gamAna::init() {

    ntrk = -999;

    h2bz0cent_evtmod = NULL;
    runno   = -999;

    goodevtno = 0;
    evtno = 0;

    pticks[0]=pticks[1]=pticks[2]=0;

    memset(emcbad,0,sizeof(emcbad));
    memset(emcnearbad,0,sizeof(emcnearbad));

    memset(h2commonmap,0,sizeof(h2commonmap));
    memset(h2runmap,0,sizeof(h2runmap));

    memset(pc3bad,0,sizeof(pc3bad));
    memset(pc3nearbad,0,sizeof(pc3nearbad));

    memset(pc3dzpar,0,sizeof(pc3dzpar));
    memset(pc3dphispar,0,sizeof(pc3dphispar));

    z0 = -9999.;
    centrality = -9999.;
    psi = -9999.;

    bz0 = -999;
    bcent = -999;
    bpsi = -999;

    nclus = -999;

    for (unsigned int i = 0; i < ntrkcut; i++)
    {
    	trkcut_map[trkcut_name[i]] = i;
    }

    trig = NULL;
    trigtree = NULL;

    trkhold = NULL;
    memset(trktree,0,sizeof(trktree));

    h1tracktype = NULL;
    h1evtraw = NULL;
    h1evtmod = NULL;
    h1skipevtmod = NULL;
    h1evtinrun = NULL;
    h1tickcut = NULL;
    h1z0 = NULL;
    h1cent_raw = NULL;
    h1bcent = NULL;
    h2bz0cent = NULL;

    h2z0zbbc = NULL;

    h1bpsi = NULL;
    h1psi = NULL;

    h1psi_bbc_raw = NULL;
    h1psi_bbcs = NULL;
    h1psi_bbcn = NULL;

    h2centpsidiff_bbcns = NULL;
    h2psi_bbc_bbcn = NULL;
    h2psi_bbc_bbcs = NULL;
    h2psi_bbcn_bbcs = NULL;

    h1psi_fvtx_raw = NULL;
    h1psi_fvtxs = NULL;
    h1psi_fvtxn = NULL;

    h2centpsidiff_fvtxns = NULL;
    h2psi_fvtx_fvtxn = NULL;
    h2psi_fvtx_fvtxs = NULL;
    h2psi_fvtxn_fvtxs = NULL;

    memset(h2iziy_emcbad,0,sizeof(h2iziy_emcbad));
    memset(h2iziy_bordwu,0,sizeof(h2iziy_bordwu));
    memset(h2iziy_nomapcut,0,sizeof(h2iziy_nomapcut));
    memset(h2iziy_fid,0,sizeof(h2iziy_fid));
    memset(h2emchdmap,0,sizeof(h2emchdmap));
    memset(h2zphi_arm,0,sizeof(h2zphi_arm));

    h2centnclus = NULL;
    h2centnpart = NULL;

    h1BBC = NULL;
    h1BBCall = NULL;
    h1BBCnar = NULL;
    h1BBCnov = NULL;

    h1ERT = NULL;
    h1ERTall = NULL;
    h1ERTBBC = NULL;

    h1asym_pi0 = NULL;

    memset(h2ecent_trigcat,0,sizeof(h2ecent_trigcat));
    memset(h2ptcent_trigcat,0,sizeof(h2ptcent_trigcat));

    memset(h1invm_trigcat,0,sizeof(h1invm_trigcat));

    memset(h2mpt,0,sizeof(h2mpt));

    h1histosec = NULL;

    h2ecent = NULL;
    h2ptcent = NULL;

    h1e_nocut = NULL;
    h1ex = NULL;
    h1ey = NULL;
    h1ez = NULL;

    memset(h1r,0,sizeof(h1r));
    memset(h1s,0,sizeof(h1s));

    h1prob_nocut = NULL;

    h2armphi_emclus = NULL;

    h1ncand = NULL;
    h1lowasym = NULL;

    memset(h2ptcent_trk,0,sizeof(h2ptcent_trk));

    h1zed_trk_raw = NULL;

    h1dcarm = NULL;

    h1quality = NULL;
    h1charge = NULL;

    h1pc3dphi = NULL;
    h1pc3dz = NULL;
    h1pc3sdphi = NULL;
    h1pc3sdz = NULL;
    h1pc3sig = NULL;

    h1emcdphi = NULL;
    h1emcdz = NULL;
    h1emcsdphi = NULL;
    h1emcsdz = NULL;
    h1emcsig = NULL;

    h1px = NULL;
    h1py = NULL;
    h1pz = NULL;

    h1n0 = NULL;
    h1ecore = NULL;
    h1eop = NULL;
    h1dep = NULL;

    h1svxnclus = NULL;

    h1chi2ndf = NULL;

    h1pattern = NULL;

    h1dca2d = NULL;
    h1dcaz = NULL;

    memset(h1ntrigperevt,0,sizeof(h1ntrigperevt));
    memset(h1e_spec,0,sizeof(h1e_spec));

    return;
}

void gamAna::setTrackType(GA::tracktype_t track_t)
{
    try
    {
	if ( !(track_t == GA::SVX || track_t == GA::CNT || track_t == GA::SVXCNT))
	    throw invalid_argument("integer does not correspond to valid tracktype");
	tracktype = track_t;
    }
    catch (const invalid_argument& e)
    {
	cerr << e.what() << endl;
	cerr << "Default track type of " << tracktype << " will be used." << endl;
    }
}

void gamAna::setRpDet(GA::detid_t detid)
{
    try
    {
	if ( !(detid == GA::FVTX || detid == GA::BBC))
	    throw invalid_argument("integer does not correspond to valid detid");
	rpdet = detid;
    }
    catch (const invalid_argument& e)
    {
	cerr << e.what() << endl;
	cerr << "Default detid of " << rpdet << " will be used." << endl;
    }
}

int gamAna::GetPc3Pars(TOAD toad_loader)
{

    if (verbosity > 0) cout<<"Reading PC3 Matching Parameters"<<endl;

    string dphis_loc = toad_loader.location("pc3dphis.txt");
    ifstream ifs_dphis(dphis_loc.c_str(), ifstream::in);

    if (!ifs_dphis.is_open())
    {
	cerr<<"Error opening pc3dphis file "<<dphis_loc<<endl;
	return -1;
    }

    if (verbosity > 0)
	cout<<dphis_loc<<" opened successfully"<<endl;

    string line;
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

    string dz_loc = toad_loader.location("pc3dz.txt");
    ifstream ifs_dz(dz_loc.c_str(), ifstream::in);

    if (!ifs_dz.is_open())
    {
	cerr<<"Error opening pc3dz file "<<dz_loc<<endl;
	return -1;
    }

    if (verbosity > 0)
	cout<<dz_loc<<" opened successfully"<<endl;

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


int gamAna::GetPc3Maps(TOAD toad_loader)
{
    if (verbosity > 0) cout<<"Reading PC3 Maps"<<endl;

    string loc = toad_loader.location("pc3bad.txt");
    ifstream ifs(loc.c_str(), ifstream::in);

    if (!ifs.is_open())
    {
	cerr<<"Error opening file "<<loc<<endl;
	return -1;
    }

    if (verbosity > 0)
	cout<<loc<<" opened successfully"<<endl;

    string line;
    while ( getline(ifs, line) )
    {
	if (!isdigit(line.at(0))) continue;
	istringstream iss_line(line.data());
	int isec, iy, iz, flag;
	iss_line >> isec >> iy >> iz >> flag;
	if (verbosity > 2)
	{
	    cout << isec <<" "<< iy <<" "<< iz <<" "<< flag << endl;
	}
	if (flag != 0)
	{
	    pc3bad[isec][iy][iz] = true;
	}
    }

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
			    pc3nearbad[isec][jy][jz]=true;
		    }
		}
	    }
	}
    }

    return 0;
}
