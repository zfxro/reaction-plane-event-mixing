#ifndef __EMCCLUSTER_H_
#define __EMCCLUSTER_H_

#include <Rtypes.h>

class EmcCluster
{
public:
    float lowasym;
    float ex;
    float ey;
    float ez;
    float pt;
    float r;
    float s;
    int sector;
    int arm;
    bool isok;
    bool isdecay;

    EmcCluster();
    virtual ~EmcCluster();

    void CalcVars(float e, float x, float y, float z, float z0);
    void SetSector(int emcarm, int emcsec);
    void Print();

};

#endif
