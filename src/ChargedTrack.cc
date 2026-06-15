#include "ChargedTrack.h"
#include <iostream>
#include <cmath>
#include <cassert>

ClassImp(ChargedTrack);

ChargedTrack::ChargedTrack()
{
    px = py = pz = -9999.;
}

ChargedTrack::~ChargedTrack()
{}

void ChargedTrack::Set3MomFromMomThePhi(float mom, float the, float phi)
{
    float sphi = sin(phi);
    float cphi = cos(phi);
    float sthe = sin(the);
    float cthe = cos(the);

    pz = mom * cthe;
    float pt = mom * sthe;
    px = pt * cphi;
    py = pt * sphi;
}

bool ChargedTrack::MomOk()
{
    if (fabs(px) > 900.) return false;
    if (fabs(py) > 900.) return false;
    if (fabs(pz) > 900.) return false;
    if (GetPt() < GA::minpta) return false;
    if (GetPt() > GA::maxpta) return false;
    return true;
}

float ChargedTrack::GetPt()
{
    return sqrt(px*px + py*py);
}

float ChargedTrack::GetMom()
{
    return sqrt(px*px + py*py + pz*pz);
}

float ChargedTrack::GetPhi()
{
    return atan2(py,px);
}

float ChargedTrack::GetEta()
{
    float mom = GetMom();
    assert(mom > 0. && "ChargedTrack::GetEta: momentum is zero");
    return atanhf(pz/mom);
}

float ChargedTrack::GetThe()
{
    float mom = GetMom();
    assert(mom > 0. && "ChargedTrack::GetThe: momentum is zero");
    return acos(pz/mom);
}

void ChargedTrack::Print()
{
    std::cout<<"ChargedTrack: px = "<<px<<", py = "<<py<<", pz = "<<pz
             <<", pt = "<<GetPt()<<", p = "<<GetMom()<<std::endl;
}
