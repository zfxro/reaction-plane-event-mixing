#ifndef __GAMCONST_H__
#define __GAMCONST_H__

#include <string>

namespace GA
{
    int GetCentBin(float centrality);
    float GetPsiBinWidth();
    int GetPsiBin(float psi);
    int GetPsiBinNear(float psi);
    int GetZvBin(float zv);
    int GetNFidTow(float ecore, std::string type);
    bool PassFidCutTowID(float energy, int towid);
    bool PassFidCut(float energy, int sec, int ypos, int zpos);
    bool TowerFromSYZ(int& towid, int sec, int ypos, int zpos);
    bool TowerToSYZ(int towid, int& sec, int& ypos, int& zpos);
    int GetTowZBin(int isec, int iz);
    int GetTowYBin(int isec, int iy);
    int GetTowZBinTowID(int towid);
    int GetTowYBinTowID(int towid);
    int GetTowBin(int itow, int ntow, int nbin);
    int GetEcoreBin(float ecore);
    bool IsSectorGl(int isec);
    bool IsSectorSc(int isec);
    int GetNzNy(int& nz, int& ny, int isec);

    enum tracktype_t
    {
        SVXCNT = 0,
        CNT = 1,
        SVX = 2
    };

    enum detid_t
    {
        BBC = 0,
        FVTX = 1
    };

    const static int narm = 2;
    const static int npsi_tree = 10;

    const static int ncent_tree = 9;
    const static float bwcent_tree = 10.;
    const static float maxcent = 93.;

    const static int nzv_tree = 6;
    const static float maxzv = 15.;

    const float minpta = 0.5;
    const float maxpta = 10.;

    const int nsec = 8;
    const int nyGl = 48;
    const int nzGl = 96;
    const int nySc = 36;
    const int nzSc = 72;
    const int ntower = 24768;

    const int nztowbin = 12;
    const int nytowbin = 6;

    const int neGl = 13;
    const float nfidtoweGl[neGl][2] = {{2, 61.2354}, {3, 30.6183}, {4, 20.4129}, {5, 15.3104}, {6, 12.249}, {7, 10.2082}, {8, 8.75067}, {9, 7.65758}, {10, 6.8075}, {11, 6.12751}, {12, 5.57123}, {13, 5.10772}, {14, 4.71559}};

    const int neSc = 9;
    const float nfidtoweSc[neSc][2] = {{2, 41.8528}, {3, 20.9273}, {4, 13.9525}, {5, 10.4654}, {6, 8.37338}, {7, 6.97889}, {8, 5.98299}, {9, 5.23622}, {10, 4.65552}};

    const int necore = 21;
    const float limecore[necore] = {0.5,0.75,1.,1.25,1.5,2.,2.5,3.,3.5,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.};

    const float emcsig_absmax = 3.;
    const float pc3sig_absmax = 3.;
    const float CFF = 10.;
    const float zed_absmax = 75.0;
    const unsigned char MAXUC = 255;
    const char MAXC = 127;
    const char MINC = -127;

    const int maxhotrun = 10;

    const int NSIGMA_PC3DR = 3;
    const int IMEAN_PC3 = 0;
    const int ISIGMA_PC3 = 1;

};

#endif
