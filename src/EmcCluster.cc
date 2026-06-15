#include "EmcCluster.h"
#include <iostream>
#include <math.h>
#include <TMath.h>

EmcCluster::EmcCluster()
{
    lowasym = 999.;
    r = s = -999.;
    ex = ey = ez = pt = -999.;
    arm = sector = -1;
    isok = isdecay = false;
}

EmcCluster::~EmcCluster()
{}

void EmcCluster::CalcVars(float e, float x, float y, float z, float z0)
{
    float zd = z - z0;
    s = sqrt(x*x + y*y);
    r = sqrt(x*x + y*y + zd*zd);
    ex = e*(x/r);
    ey = e*(y/r);
    ez = e*(zd/r);
    pt = sqrt(ex*ex + ey*ey);
}

void EmcCluster::SetSector(int emcarm, int emcsec)
{
    arm = fabs(emcarm - 1);
    sector = arm*4 + emcsec;
}

void EmcCluster::Print()
{
    std::cout<<"Cluster: ex = "<<ex<<"  ey = "<<ey<<"  ez = "<<ez<<std::endl;
}
