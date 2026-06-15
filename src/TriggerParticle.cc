#include "TriggerParticle.h"
#include "gamConst.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <TMath.h>

ClassImp(TriggerParticle);

TriggerParticle::TriggerParticle()
{
    type = towid = -9999;
    px = py = pz = -9999.;
    pc3dphis = pc3dz = trkdphis = trkdz = -9999.;
}

TriggerParticle::TriggerParticle(int t, float momx, float momy, float momz)
    : type(t), px(momx), py(momy), pz(momz)
{
    towid = -9999;
    pc3dphis = pc3dz = trkdphis = trkdz = -9999.;
}

TriggerParticle::~TriggerParticle()
{}

void TriggerParticle::SetTrkVeto(float dcu, float dcz)
{
    trkdphis = dcu;
    trkdz = dcz;
}

void TriggerParticle::SetPc3Veto(float padu, float padz)
{
    pc3dphis = padu;
    pc3dz = padz;
}

void TriggerParticle::Print()
{
    std::cout<<"Trigger: type = "<<type<<"  px = "<<px<<"  py = "<<py
             <<"  pt = "<<GetPt()<<"  pz = "<<pz<<std::endl;
}

float TriggerParticle::GetPt()
{
    return sqrt(px*px + py*py);
}

float TriggerParticle::GetMom()
{
    return sqrt(px*px + py*py + pz*pz);
}

float TriggerParticle::GetPhi()
{
    return atan2(py,px);
}

float TriggerParticle::GetEta()
{
    float mom = GetMom();
    assert(mom > 0. && "TriggerParticle::GetEta: momentum is zero");
    return atanhf(pz/mom);
}

float TriggerParticle::GetPc3dr()
{
    if (type == 2) return -9999.;
    if (pc3dphis == -9999. || pc3dz == -9999.) return -9999.;
    return sqrt(pc3dphis*pc3dphis + pc3dz*pc3dz);
}

float TriggerParticle::GetTrkdr()
{
    if (type == 2) return -9999.;
    if (trkdphis == -9999. || trkdz == -9999.) return -9999.;
    return sqrt(trkdphis*trkdphis + trkdz*trkdz);
}
