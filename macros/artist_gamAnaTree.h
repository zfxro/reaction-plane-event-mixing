#include "../src/gamAnaTree/gamAnaTree.h"
#include "../src/gamConst.h"
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TFile.h>
#include <TROOT.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <fstream>
#include <sstream>


void getCounters(bool isadj);
void MixSubNorm();
void GetRemBg();
void SubRemBg();
void Decompose(const char* filename_tageff, const char* filename_mistag);
void GetNhExp();
void GetNhExpXi(int);
void dPhiHistos(TFile* input, TFile* input_nopsimix, bool isadj);
void dPhiHistosXi(TFile* input, TFile* input_nopsimix, bool isadj);


TCanvas* cansolo;

const static int rebinint = 2;
const bool useksic = false;

bool skipsave = false;
bool skippi0 = false;
bool savesolo = false;//takes while so keeping false until needed

float ntrg_pttcentdpsi[ntype][nptt+1][ncent+1][ndpsi+1];//+1 for all.
float nmix_pttcentdpsi[ntype][nptt+1][ncent+1][ndpsi+1];//+1 for all.

const static int NCAT = 4;
enum category
{
    CENT=0,
    DPSI=1,
    PTT=2,
    PTA=3
};
const char* catstr[NCAT] =
{
    "cent",
    "dpsi",
    "ptt",
    "pta"
};
const int NDRAW[NCAT] =
{
    ncent+1,
    ndpsi+1,
    nptt,
    npta
};

const static int NMIXTYPE = 3;
enum mixtype
{
    RAW=0,
    MIX=1,
    MSN=2
};
const char* mixstr[NMIXTYPE] =
{
    "raw",
    "mix",
    "msn"
};

const static int nxi_new = 2;

//FVTX resolution
static const float reso[ncent+1] = {
    0.674302,
    0.745883,
    0.279465,
    0.68742
};

static const float PI = TMath::Pi();
const static int nYpt = 18;
const static float minYpt = 5.25;
const static float maxYpt = 14.25;
const static float Ybw = .5;

const static float detah = .7;

vector<vector<float>> expnh;

vector<vector<vector<vector<vector<float>>>>> expnhxi;


static const float pt_expect[2][nYpt][2] = {
    {{5.5, 0.0000459108}, {6., 0.0000284513}, {6.5, 0.0000185074}, {7., 0.0000122213}, {7.5, 8.26532e-6}, {8., 5.75881e-6}, {8.5, 4.04008e-6}, {9., 2.98215e-6}, {9.5, 2.25102e-6}, {10., 1.66023e-6}, {10.5, 1.21401e-6}, {11., 8.80615e-7}, {11.5, 6.68967e-7}, {12., 5.19039e-7}, {12.5, 4.13354e-7}, {13., 3.38868e-7}, {13.5, 2.69615e-7}, {14., 2.14348e-7}},
    {{5.5, 0.0000223201}, {6., 0.0000120278}, {6.5, 6.77668e-6}, {7., 4.00986e-6}, {7.5, 2.4668e-6}, {8., 1.56493e-6}, {8.5, 1.02312e-6}, {9., 6.80489e-7}, {9.5, 4.58995e-7}, {10., 3.16338e-7}, {10.5, 2.23017e-7}, {11., 1.60029e-7}, {11.5, 1.17542e-7}, {12., 8.87529e-8}, {12.5, 6.75602e-8}, {13., 5.12616e-8}, {13.5, 3.8708e-8}, {14., 2.90775e-8}}
};

TH1F* h1dphi_nopsimix[ntype+1][nptt][ncent+1][ndpsi+1][npta];
TH1F* h1dphi_bg[ntype+1][nptt][ncent+1][ndpsi+1][npta];

TH1F* h1dphi_xi_nopsimix[ntype+1][nptt][ncent+1][ndpsi+1][nxi];
TH1F* h1dphi_xi_nopsimix_ksic[ntype+1][nptt][ncent+1][ndpsi+1][nxi];
TH1F* h1dphi_xi_bg[ntype+1][nptt][ncent+1][ndpsi+1][nxi];
