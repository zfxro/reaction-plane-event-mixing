#ifndef __GAMANA_H__
#define __GAMANA_H__

#include <SubsysReco.h>
#include <vector>
#include <map>

#include "gamConst.h"

class TOAD;

class PHCompositeNode;

class RunHeader;
class EventHeader;
class PreviousEvent;
class VtxOut;
class PHGlobal;
class PHCentralTrack;

class TriggerHelper;

class SvxCentralTrackList;

class ReactionPlaneObject;

class emcClusterContainer;
class emcClusterContent;

// ROOT
class TFile;
class TNtuple;
class TTree;
class TH1I;
class TH2I;
class TH1F;
class TH1D;
class TH2F;
class TH2D;
class TH3F;
class TH3D;
class TProfile;

class CentralTrackData;
class ChargedTrack;
class CentralArmTrack;
class TrackCollectionBase;
class TriggerParticle;
class EmcCluster;

class gamAna: public SubsysReco {
public:

    gamAna(std::string &filename);
    virtual ~gamAna() {}

    virtual int  Init(PHCompositeNode *topNode);
    virtual int  process_event(PHCompositeNode *topNode);
    virtual int  End(PHCompositeNode *topNode);

    virtual int  InitRun(PHCompositeNode *topNode);

    void doFullAna(bool boolean) {dofullana = boolean;}
    void doTrackAna(bool boolean) {dotrackana = boolean;}
    void doMistagAna(bool boolean) {domistagana = boolean;}
    void readEventCount(bool boolean) {readeventcount = boolean;}
    void useEmc3x3(bool boolean) {useemc3x3 = boolean;}
    void useMyEmc(bool boolean) {usemyemc = boolean;}
    void saveAllEvents(bool boolean) {saveallevt = boolean;}

    void setTrackType(GA::tracktype_t track_t);
    void setRpDet(GA::detid_t detid);

private:
    void init();

    int MistagAna(bool plotlowephotons = false);

    bool makeChargedTrack(ChargedTrack** tmptrk, int itrk, bool fillhists = false);
    bool makeCentralTrackData(CentralTrackData* mycnt, int icnt, bool is_svxcnt, bool fillhists, float pt);
    bool isClusterOK(emcClusterContent* clus, EmcCluster* myclus, float emin, bool fillhists);

    int  GetNodes(PHCompositeNode *topNode);
    void initHistos(TFile* map = NULL);

    bool TickCut();

    bool isBBC(TriggerHelper* trig);
    bool isERT(TriggerHelper* trig);

    void remapCntSvx();

    int GetHistoSector(int sector_i, int sector_j);

    int GetPc3Pars(TOAD toad_loader);
    int GetPc3Maps(TOAD toad_loader);

private:

    // Output histogram file
    TFile*           outfile;
    std::string      outfilename;

    // Data nodes
    RunHeader*                d_runhead;
    EventHeader*              d_eventhead;
    PHGlobal*                 d_global;
    PHCentralTrack*           d_trx;
    VtxOut*                   d_vtxout;
    PreviousEvent*            d_peve;

    SvxCentralTrackList*      d_svxtrx;

    ReactionPlaneObject*      d_rp;

    emcClusterContainer*      d_emccontain;

    bool dofullana;
    bool dotrackana;
    bool domistagana;
    bool readeventcount;
    bool useemc3x3;
    bool usemyemc;
    bool saveallevt;

    GA::tracktype_t tracktype;
    GA::detid_t rpdet;

    const static int ncut = 25;

    // Histogram settings
    const static int nebin = 160;
    const static float emax_h = 40.;

    const static int nezbin = 100;
    const static float ezmax_h = 5.;

    const static int nphimatchbin = 300;
    const static float phimatchmin_h = -.3;
    const static float phimatchmax_h = .3;

    const static int nzmatchbin = 100;
    const static float zmatchmin_h = -10.;
    const static float zmatchmax_h = 10.;

    const static int nmatchbin = 80;
    const static float matchmin_h = -10.;
    const static float matchmax_h = 10.;

    const static int ndcabin = 200;
    const static float dcamin_h = -1.;
    const static float dcamax_h = 1.;

    const static int nsvxphibin = 200;
    const static float svxphimin_h = -.5;
    const static float svxphimax_h = 1.5;

    const static int netabin = 400;
    const static float etamin_h = -.5;
    const static float etamax_h = .5;

    const static int ninvmptbin = 10.;
    const static int ninvmcentbin = 95.;

    const static int ntrigcat = 6;

    const static int npsibin = 200.;
    const static int psimin_h = 0.0;
    const static int psimax_h = 1.0;

    const static int npsibin_2d = 50.;
    const static int psimin_h_2d = 0.0;
    const static int psimax_h_2d = 1.0;

    const static int nspec = 3;
    const static char* spec[nspec];

    const static float pi0pars[GA::nsec][2];

    const static float emin_clus = 0.5;
    const static float ptmin_pi0trig = 5.0;
    const static float ptmin_gamtrig = 5.0;
    const static float asymmax = .8;
    const static float npi0sigma = 3.;

    // EMC hot/dead tower maps
    bool emcbad[GA::nsec][GA::nyGl][GA::nzGl];
    bool emcnearbad[GA::nsec][GA::nyGl][GA::nzGl];

    bool pc3bad[GA::nsec][GA::nyGl][GA::nzGl];
    bool pc3nearbad[GA::nsec][GA::nyGl][GA::nzGl];

    double pc3dzpar[GA::nsec][GA::nztowbin][2];
    double pc3dphispar[GA::nsec][GA::nytowbin][2];

    TH2I* h2commonmap[GA::nsec];
    TH2I* h2runmap[GA::nsec];

    TH2F* h2bz0cent_evtmod;
    const static int nevtdesire = 500;
    const static int npsimixbin = 40;

    // Event variables
    float z0;
    float centrality;
    float psi;

    int bz0;
    int bcent;
    int bpsi;

    int runno;
    int evtno;
    int goodevtno;
    int pticks[3];

    int ntrk;
    int nclus;

    const static unsigned int ntrkcut = 7;
    const static std::string trkcut_name[ntrkcut];
    std::map<std::string,unsigned int> trkcut_map;

    TrackCollectionBase* trkhold;
    TriggerParticle* trig;

    // Trees: [centrality][reaction plane][z-vertex]
    TTree* trktree[GA::ncent_tree][GA::npsi_tree][GA::nzv_tree];
    TTree* trigtree;

    //------------------- Global Histograms -----------------------

    TH1F* h1tracktype;

    TH1F* h1evtraw;
    TH1F* h1evtmod;
    TH1F* h1skipevtmod;
    TH1F* h1evtinrun;

    TH1F* h1tickcut;

    TH1F* h1z0;
    TH1F* h1cent_raw;
    TH1F* h1bcent;
    TH2F* h2bz0cent;

    TH2F* h2z0zbbc;

    TH1F* h1bpsi;
    TH1F* h1psi;

    TH1F* h1psi_bbc_raw;
    TH1F* h1psi_bbcs;
    TH1F* h1psi_bbcn;

    TH2F* h2centpsidiff_bbcns;
    TH2F* h2psi_bbc_bbcn;
    TH2F* h2psi_bbc_bbcs;
    TH2F* h2psi_bbcn_bbcs;

    TH1F* h1psi_fvtx_raw;
    TH1F* h1psi_fvtxs;
    TH1F* h1psi_fvtxn;

    TH2F* h2centpsidiff_fvtxns;
    TH2F* h2psi_fvtx_fvtxn;
    TH2F* h2psi_fvtx_fvtxs;
    TH2F* h2psi_fvtxn_fvtxs;

    // EMC map histograms
    TH2I* h2iziy_emcbad[GA::nsec];
    TH2I* h2iziy_bordwu[GA::nsec];
    TH2I* h2iziy_nomapcut[GA::nsec];
    TH2I* h2emchdmap[GA::nsec];
    TH2I* h2iziy_fid[GA::nsec];
    TH2F* h2zphi_arm[GA::narm];

    // Centrality
    TH2F* h2centnclus;
    TH2F* h2centnpart;

    // Trigger histograms
    TH1F* h1BBC;
    TH1F* h1BBCall;
    TH1F* h1BBCnar;
    TH1F* h1BBCnov;

    TH1F* h1ERT;
    TH1F* h1ERTall;
    TH1F* h1ERTBBC;

    //--------------- Pi0 histograms ------------------

    TH1F* h1asym_pi0;

    TH1F* h1invm_trigcat[ntrigcat][ninvmcentbin][ninvmptbin];

    TH2F* h2ecent_trigcat[ntrigcat];
    TH2F* h2ptcent_trigcat[ntrigcat];

    TH2F* h2mpt[GA::nsec+6];

    //--------------- Cluster histograms ------------------

    TH1F* h1histosec;

    TH2F* h2ecent;
    TH2F* h2ptcent;

    TH1F* h1e_nocut;
    TH1F* h1ex;
    TH1F* h1ey;
    TH1F* h1ez;

    TH1F* h1r[GA::nsec];
    TH1F* h1s[GA::nsec];

    TH1F* h1prob_nocut;

    TH2F* h2armphi_emclus;

    //--------------- Candidate histograms -----------------

    TH1F* h1ncand;
    TH1F* h1lowasym;

    //--------------- Central track histograms -----------------

    TH2F* h2ptcent_trk[ntrkcut];

    TH1F* h1zed_trk_raw;
    TH1F* h1dcarm;

    TH1F* h1quality;
    TH1F* h1charge;

    // PC3 matching
    TH1F* h1pc3dphi;
    TH1F* h1pc3dz;
    TH1F* h1pc3sdphi;
    TH1F* h1pc3sdz;
    TH1F* h1pc3sig;

    // EMC matching
    TH1F* h1emcdphi;
    TH1F* h1emcdz;
    TH1F* h1emcsdphi;
    TH1F* h1emcsdz;
    TH1F* h1emcsig;

    TH1F* h1px;
    TH1F* h1py;
    TH1F* h1pz;

    TH1F* h1ecore;
    TH1F* h1n0;
    TH1F* h1eop;
    TH1F* h1dep;

    TH1F* h1svxnclus;
    TH1F* h1chi2ndf;
    TH1F* h1pattern;

    // DCA
    TH1F* h1dca2d;
    TH1F* h1dcaz;

    //--------------- Fragmentation function histograms -----------------

    TH1F* h1ntrigperevt[nspec];
    TH1F* h1e_spec[nspec];

};

#endif
