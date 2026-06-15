#include "CentralTrackData.h"
#include "gamConst.h"
#include <iostream>
#include <math.h>

ClassImp(CentralTrackData);

CentralTrackData::CentralTrackData()
{
    arm = isn0neg = 0;
    pc3sig = emcsig = GA::MAXUC;
    dep = GA::MAXC;
}

CentralTrackData::~CentralTrackData()
{}

bool CentralTrackData::QualOk(int quality)
{
    if (!(quality == 31 || quality == 63)) return false;
    return true;
}

bool CentralTrackData::ZedOk(float zed)
{
    if (fabs(zed) > GA::zed_absmax) return false;
    return true;
}

bool CentralTrackData::EmcOk()
{
    unsigned char max = GA::emcsig_absmax*GA::CFF;
    return (emcsig < max);
}

bool CentralTrackData::EmcMatch(float nsigma)
{
    float max_f = roundf(nsigma * GA::CFF);
    if (max_f <= 0. || max_f > (float)GA::MAXC)
    {
        std::cerr<<"CentralTrackData::EmcMatch: nsigma "<<nsigma<<" is out of range!"<<std::endl;
        return false;
    }
    unsigned char max = max_f;
    return (emcsig < max);
}

void CentralTrackData::SetEmcSig(float emcsig_f)
{
    emcsig_f = emcsig_f*GA::CFF;
    if (emcsig_f > (float)GA::MAXUC)
    {
        emcsig = GA::MAXUC;
        return;
    }
    else if (emcsig_f < 0.)
    {
        std::cerr<<"CentralTrackData::SetEmcSig: negative value, setting to MAXUC"<<std::endl;
        emcsig = GA::MAXUC;
        return;
    }
    else emcsig = emcsig_f;
}

bool CentralTrackData::CheckDEP(float sigmax)
{
    float max_f = roundf(sigmax * GA::CFF);
    if (max_f < (float)GA::MINC || max_f > (float)GA::MAXC)
    {
        std::cerr<<"CentralTrackData::CheckDEP: sigmax "<<sigmax<<" is out of range!"<<std::endl;
        return false;
    }
    char max = max_f;
    if (max < 0)
        return (dep <= max);

    return (dep < max);
}

int CentralTrackData::SetDEP(float dep_f)
{
    if (dep_f < -900.)
    {
        dep = GA::MAXC;
        return -1;
    }
    dep_f = dep_f*GA::CFF;
    if (dep_f > (float)GA::MAXC)
    {
        dep = GA::MAXC;
        return 2;
    }
    else if (dep_f < (float)GA::MINC)
    {
        dep = GA::MINC;
        return 1;
    }
    else dep = dep_f;
    return 0;
}

bool CentralTrackData::Pc3Ok()
{
    unsigned char max = GA::pc3sig_absmax*GA::CFF;
    return (pc3sig < max);
}

bool CentralTrackData::Pc3Match(float nsigma)
{
    float max_f = roundf(nsigma * GA::CFF);
    if (max_f <= 0 || max_f > (float)GA::MAXUC)
    {
        std::cerr<<"CentralTrackData::Pc3Match: nsigma "<<nsigma<<" is out of range!"<<std::endl;
        return false;
    }
    unsigned char max = max_f;
    return (pc3sig < max);
}

void CentralTrackData::SetPc3Sig(float pc3sig_f)
{
    pc3sig_f = pc3sig_f*GA::CFF;
    if (pc3sig_f > (float)GA::MAXUC)
    {
        pc3sig = GA::MAXUC;
        return;
    }
    else if (pc3sig_f < 0.)
    {
        std::cerr<<"CentralTrackData::SetPc3Sig: negative value, setting to MAXUC"<<std::endl;
        pc3sig = GA::MAXUC;
        return;
    }
    else pc3sig = pc3sig_f;
}

void CentralTrackData::Print()
{
    std::cout<<"CentralTrackData: pc3sig = "<<(int)pc3sig<<", emcsig = "<<(int)emcsig
             <<", isn0neg = "<<isn0neg<<", arm = "<<arm<<", dep = "<<(int)dep<<std::endl;
}
