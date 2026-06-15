#ifndef __CHARGEDTRACK_H_
#define __CHARGEDTRACK_H_

#include "gamConst.h"
#include <Rtypes.h>

class CentralTrackData;

class ChargedTrack
{
public:

    virtual GA::tracktype_t GetTrackType() = 0;

    void Set3MomFromMomThePhi(float, float, float);

    bool MomOk();

    float GetPt();
    float GetMom();
    float GetPhi();
    float GetEta();
    float GetThe();

    virtual CentralTrackData* GetCentralTrackData() {return NULL;};

    virtual void Print();

    ChargedTrack();
    virtual ~ChargedTrack();

public:

    float px;
    float py;
    float pz;

public:

    ClassDef(ChargedTrack,1)
};

#endif
