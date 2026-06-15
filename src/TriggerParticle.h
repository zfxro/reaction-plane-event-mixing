#ifndef __TRIGGERPARTICLE_H_
#define __TRIGGERPARTICLE_H_

#include <Rtypes.h>

class TriggerParticle
{
public:
    int type;
    float px;
    float py;
    float pz;

    float pc3dphis;
    float pc3dz;
    float trkdphis;
    float trkdz;

    int towid;

    TriggerParticle();
    virtual ~TriggerParticle();

    TriggerParticle(int, float, float, float);

    void SetTrkVeto(float, float);
    void SetPc3Veto(float, float);

    float GetPt();
    float GetMom();
    float GetPhi();
    float GetEta();

    float GetPc3dr();
    float GetTrkdr();

    void Print();

    ClassDef(TriggerParticle,1)
};

#endif
