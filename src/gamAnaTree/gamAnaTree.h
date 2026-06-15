#ifndef __GAMANATREE_H__
#define __GAMANATREE_H__

// NOTE: nptt is changed to easily allow artist and hadronaccept code to deal with ptt rebinned gamAnaTree output.  When running gamAnaTree the finer bins may be used (nptt = 4)
// keyword: REBIN
// be careful. nxi is transitioning to behave like nptt. no longer rebin in artist
// but again finer bins use on train!

#include "../gamConst.h"

#include <map>
#include <vector>

#define ncent_fun 94

//---------- classes needed ----------

//histos
class TH1F;
class TH2F;
class TH2I;

class TF1;

class TFile;

//my classes
class ChargedTrack;
class TriggerParticle;

bool emcwu[GA::nsec][GA::nyGl][GA::nzGl];
bool pc3bad[GA::nsec][GA::nyGl][GA::nzGl];
bool pc3bad3x3[GA::nsec][GA::nyGl][GA::nzGl];

double pc3dzpar[GA::nsec][GA::nztowbin][2];
double pc3dphispar[GA::nsec][GA::nytowbin][2];

const unsigned int nmixmin = 200;

const int npt_trkcut = 9;
const float limpt_trkcut[npt_trkcut+1] = {.75,1.,1.25,1.75,2.25,3.,4.,5.,7.,10.};

const int ncent_trkcut = 5;
const int limcent_trkcut[ncent_trkcut+1] = {0,10,20,40,60,100};

const static int ntrkcut = 7;

float nsig_trkcut[ntrkcut][ncent_trkcut][npt_trkcut];
float sigma_trkcut[ntrkcut][ncent_trkcut][npt_trkcut];
float mean_trkcut[ntrkcut][ncent_trkcut][npt_trkcut];

enum cut
{
    CHI2NDF,
    DCAZ,
    DCA2D,
    PC3DPHI,
    PC3DZ,
    EMCDPHI,
    EMCDZ
};

std::map<std::string,cut> cutmap;

int GetCentBin_TrkCut(int centrality);
int GetPtBin_TrkCut(float pt);

const int npta = 5;
const float limpta[npta+1] = {0.5,1.,2.,3.,4.,5.};

const float ptamixmax0 = 2.0;
const float ptamixmax1 = 3.0;

//keyword: REBIN
const int nxi = 3;
const float limxi[nxi+1] = {0.0,0.8,1.6,2.4};

const int maxcent = 60;

//centrality binning
const int ncent = 8;
const int limcent[ncent+1] = {0,5,10,15,20,30,40,50,60};

//pt binning
// keyword: REBIN
const int nptt = 4;
const float limptt[nptt+1] = {5.,7.,9.,12.,15.};

//psi binning
const int ndpsi = 2;
const char* dpsistr[ndpsi+1] =
{
    "[0,1/4], (3/4,0]",
    "(1/4,3/4]",
    "[0,1]"
};

const float psimixbw = 0.025;//limit for psi mixing
const int npsimixbin = 40;// = 1/psimixbw

//triggers
const int ntype = 3;
const char* spec[ntype+1] = {"dir","dec","pi0","inc"};
enum TRIG
{
    DIR=0,
    DEC=1,
    PI0=2,
    INC=3
};
const int ntype_adj = 2;

//away side ratio
//integrated away side yield ratio for tagged and missed decay photons
const double asr = 1.00;

//track types: before and after cuts
const int ntype_trk = 2;

bool dofidcut;

bool dochargeveto;

//histo constants
const float maxp = 50.;
const int npbin = 1000;

const float maxpz = 10.;
const int npzbin = 1000;

const int nphibin = 200;
const float minphi = -.5;
const float maxphi = 1.5;

const int ndphibin = 12;
const float mindphi = 0.;
const float maxdphi = 1.;

const int netabin = 400;
const float mineta = -.5;
const float maxeta = .5;

const float nphimatchbin = 300;
const float minphimatchbin = -.3;
const float maxphimatchbin = .3;

const float nzmatchbin = 100;
const float minzmatchbin = -10.;
const float maxzmatchbin = 10.;

const float nsmatchbin = 80;
const float minsmatchbin = -4.;
const float maxsmatchbin = 4.;

const int ndcabin = 2000;
const float mindcabin = -1.;
const float maxdcabin = 1.;

const int nfragbin = 100;
const float minfrag = 0.;
const float maxfrag = 2.;

//cut consts
const float dphicut = 0.2;

//---------- variables ----------

int verbosity;

int iptt_ana;
int icent_ana;
int idpsi_ana;

//test to ensure ismixable tag set properly
std::vector<int> trigevtno;
std::vector<int> nomixevtno;

TFile* infile;
TFile* outfile;

//---------- some methods ----------

int GetEmcMaps();
int GetPc3Maps();
int GetPc3Pars();
int GetEmcRunMaps(int runno);
int GetFunctions(const char* loc, const char* name, TF1* f1[ncent_fun]);
int GetWeights(double& dirwt, double& decwt, double& dirwtY, double& decwtY, TF1* ftageff[ncent_fun], TF1* fmistag[ncent_fun], int type, float cent, float pt);

//binning methods. these indices will be used for made plots
int GetCentBin(float cent);
int GetDPsiBin(float psi,float phi);
int GetPtTrigBin(float pt);
int GetPtAssocBin(float pt);
int GetXiBin(float pt_trg, float pt_trk);

//trigger methods
void FillTrigHistos(TriggerParticle* trggr);

//track methods
void FillTrkHistos(ChargedTrack* trk, int type, int icent);
void TrackAnalysis();
void TrackEventAnalysis();
bool IsTrackGood(ChargedTrack* trk);

void End();
void InitHistos();
int ReadCutVars();

bool FillCorrels(TriggerParticle* trg, ChargedTrack* trk, bool same_evt, TH1F* h1dphi_ptr[][npta], TH1F* h1dphi_xi_ptr[][nxi], TH1F* h1dphiacc_xi_ptr[nxi]);

//---------- histograms ----------

TH2I* h2commap[GA::nsec];
TH2I* h2runmap[GA::nsec];

TH2F* h2zypos[GA::nsec];

TH1F* h1nmixavail;

TH1F* h1nmixevt_raw;
TH1F* h1nmixevt;

TH1F* h1binsize_trkevt;

//===== track event histos =====

TH2F* h2cent_zv_trkevt;
TH2F* h2cent_ntrk_trkevt[ntype_trk];

//===== track histos =====

TH2F* h2cent_pta;

TH1F* h1quality[ntype_trk];
TH1F* h1charge[ntype_trk];

//pc3 matching
TH1F* h1pc3dphi[ntype_trk];
TH1F* h1pc3dz[ntype_trk];
TH1F* h1pc3sdphi[ntype_trk];
TH1F* h1pc3sdz[ntype_trk];

TH1F* h1pc3sig[ntype_trk];

//emc matching
TH1F* h1emcdphi[ntype_trk];
TH1F* h1emcdz[ntype_trk];
TH1F* h1emcsdphi[ntype_trk];
TH1F* h1emcsdz[ntype_trk];

TH1F* h1emcsig[ntype_trk];

TH1F* h1phi0[ntype_trk];
TH1F* h1eta[ntype_trk];
TH1F* h1pt[ntype_trk];
TH1F* h1px[ntype_trk];
TH1F* h1py[ntype_trk];
TH1F* h1pz[ntype_trk];

//4bit hit pattern
TH1F* h14bitpat[ntype_trk];
TH1F* h1zed[ntype_trk];
TH1F* h1n0[ntype_trk];
TH1F* h1ecore[ntype_trk];

TH1F* h1dep[ntype_trk];

//chi^2 / ndf
TH1F* h1chi2ndf[ntype_trk];

TH1F* h1pattern[ntype_trk];

//DCA
TH1F* h1dca2d[ntype_trk];
TH1F* h1dcaz[ntype_trk];

//===== trigger histos =====
TH1F* h1cent_trg[ntype];
TH1F* h1psi_trg[ntype];
TH1F* h1zv_trg[ntype];
TH1F* h1phi_trg[ntype];
TH2F* h2etaphi_trg[ntype];
TH1F* h1py_trg[ntype];
TH1F* h1px_trg[ntype];
TH1F* h1pz_trg[ntype];

TH1F* h1pc3dphis_trg[ntype];
TH1F* h1pc3dz_trg[ntype];
TH1F* h1pc3dr_trg[ntype];

TH1F* h1trkdphis_trg[ntype];
TH1F* h1trkdz_trg[ntype];
TH1F* h1trkdr_trg[ntype];

TH1F* h1pc3sdz_trg;
TH1F* h1pc3sdphis_trg;
TH1F* h1pc3sdr_trg;

TH2F* h2cent_ptt[ntype][ndpsi];
TH2F* h2cent_ptt_adj[2][ndpsi];
TH2F* h2cent_ptt_deceff[ndpsi];

//track+trigger histos
vector<vector<vector<vector<vector<vector<TH1F*> > > > > > h1dphi(ntype+1, vector<vector<vector<vector<vector<TH1F*> > > > >(3, vector<vector<vector<vector<TH1F*> > > >(nptt, vector<vector<vector<TH1F*> > >(ncent+1, vector<vector<TH1F*> >(ndpsi+1, vector<TH1F*>(npta))))));

vector<vector<vector<vector<vector<vector<TH1F*> > > > > > h1dphi_adj(ntype+1, vector<vector<vector<vector<vector<TH1F*> > > > >(3, vector<vector<vector<vector<TH1F*> > > >(nptt, vector<vector<vector<TH1F*> > >(ncent+1, vector<vector<TH1F*> >(ndpsi+1, vector<TH1F*>(npta))))));

vector<vector<vector<vector<vector<vector<TH1F*> > > > > > h1dphi_xi(ntype+1, vector<vector<vector<vector<vector<TH1F*> > > > >(3, vector<vector<vector<vector<TH1F*> > > >(nptt, vector<vector<vector<TH1F*> > >(ncent+1, vector<vector<TH1F*> >(ndpsi+1, vector<TH1F*>(nxi))))));

vector<vector<vector<vector<vector<vector<TH1F*> > > > > > h1dphi_xi_adj(ntype+1, vector<vector<vector<vector<vector<TH1F*> > > > >(3, vector<vector<vector<vector<TH1F*> > > >(nptt, vector<vector<vector<TH1F*> > >(ncent+1, vector<vector<TH1F*> >(ndpsi+1, vector<TH1F*>(nxi))))));

vector<vector<vector<vector<vector<TH1F*> > > > > h1dphiacc_xi(ntype+1, vector<vector<vector<vector<TH1F*> > > >(nptt, vector<vector<vector<TH1F*> > >(ncent+1, vector<vector<TH1F*> >(ndpsi+1, vector<TH1F*>(nxi)))));

vector<vector<vector<vector<vector<TH1F*> > > > > h1dphiacc_xi_adj(ntype+1, vector<vector<vector<vector<TH1F*> > > >(nptt, vector<vector<vector<TH1F*> > >(ncent+1, vector<vector<TH1F*> >(ndpsi+1, vector<TH1F*>(nxi)))));

TH1F* h1xi[ntype][nptt][ncent+1][ndpsi+1];
TH1F* h1xi_adj[ntype][nptt][ncent+1][ndpsi+1];

TH1F* h1dphi_test;
TH1F* h1dphi_test_raw;

#endif
