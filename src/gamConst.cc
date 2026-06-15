#include "gamConst.h"

#include <TMath.h>

#include <iostream>

namespace GA
{

    int GetZvBin(float zv)
    {
        float binwidth = (2*maxzv)/((float)nzv_tree);
        int bin = (zv+maxzv)/binwidth;

        if (bin < 0 || bin >= nzv_tree)
        {
            if (fabs(zv - maxzv) < .001*maxzv)
            {
                bin = nzv_tree - 1;
                return bin;
            }
            else
            {
                std::cerr<<"gamConst.cc : zv = "<<zv<<" and bin = "<<bin<<" which is < 0 || >= "<<nzv_tree<<std::endl;
                return -999;
            }
        }
        else return bin;
    }

    float GetPsiBinWidth()
    {
        return ((float) TMath::Pi())/((float)npsi_tree);
    }

    int GetPsiBin(float psi)
    {
        float binwidth = GetPsiBinWidth();
        int bin = psi/binwidth;

        if (bin < 0 || bin >= npsi_tree)
        {
            if (fabs(psi - TMath::Pi()) < .001*TMath::Pi())
            {
                bin = npsi_tree - 1;
                return bin;
            }
            else
            {
                std::cerr<<"gamConst.cc : psi bin out of bounds: "<<bin<<std::endl;
                return -999;
            }
        }
        else return bin;
    }

    int GetPsiBinNear(float psi)
    {
        float binwidth = GetPsiBinWidth();
        int bin = psi/binwidth;

        if (bin < 0 || bin >= npsi_tree)
        {
            if (fabs(psi - TMath::Pi()) < .001*TMath::Pi())
            {
                bin = npsi_tree - 1;
            }
            else
            {
                std::cerr<<"gamConst.cc : psi bin out of bounds: "<<bin<<std::endl;
                return -999;
            }
        }

        float psibinmid = (bin + 0.5) * binwidth;
        int binnear = -999;
        if (psi >= psibinmid)
        {
            binnear = bin + 1;
            if (binnear == npsi_tree)
                binnear = 0;
        }
        else
        {
            binnear = bin - 1;
            if (binnear == -1)
                binnear = npsi_tree - 1;
        }
        return binnear;
    }


    int GetCentBin(float centrality)
    {
        int bin = centrality/bwcent_tree;

        if (bin >= ncent_tree)
            bin = ncent_tree - 1;

        if (bin < 0)
        {
            std::cerr<<"gamConst.cc : cent bin out of bounds: "<<bin<<std::endl;
            return -999;
        }
        else return bin;
    }

    int GetNFidTow(float ecore, std::string type)
    {
        const float (*nfidtowe)[2];
        int ne;
        if (type == "Sc")
        {
            nfidtowe = &nfidtoweSc[0];
            ne = neSc;
        }
        else if (type == "Gl")
        {
            nfidtowe = &nfidtoweGl[0];
            ne = neGl;
        }
        else
        {
            std::cerr<<"Invalid type "<<type<<std::endl;
            return -1;
        }

        for (int ie = 0; ie < ne; ie++)
        {
            if (ecore > nfidtowe[ie][1])
            {
                return nfidtowe[ie][0];
            }
        }

        return -1;
    }

    bool PassFidCut(float energy, int sec, int ypos, int zpos)
    {
        int nz;
        int ny;
        int nfidtow;
        if (sec > 1)
        {
            nfidtow = GetNFidTow(energy,"Sc");
            nz = nzSc;
            ny = nySc;
        }
        else
        {
            nfidtow = GetNFidTow(energy,"Gl");
            nz = nzGl;
            ny = nyGl;
        }

        if (nfidtow < 0)
        {
            std::cerr<<"GA::PassFidCut() returned "<<nfidtow<<" which is out of bounds."<<std::endl;
            return false;
        }

        if (zpos < nfidtow || zpos >= (nz - nfidtow)) return false;

        if (sec == 3 || sec == 7)
        {
            if (ypos >= (ny - nfidtow)) return false;
        }

        if (sec == 0 || sec == 4)
        {
            if (ypos < nfidtow) return false;
        }

        return true;
    }

    bool PassFidCutTowID(float energy, int towid)
    {
        int sec = -9999;
        int ypos = -9999;
        int zpos = -9999;

        bool success = TowerToSYZ(towid,sec,ypos,zpos);

        if (!success) return false;

        return PassFidCut(energy,sec,ypos,zpos);
    }

    bool TowerFromSYZ(int& towid, int sec, int ypos, int zpos)
    {
        if (sec >= 0 && sec < 2)
        {
            if (ypos < 0 || ypos >= GA::nyGl)
            {
                std::cerr<<"gamConst::TowerFromSYZ: out of range ypos = "<<ypos<<std::endl;
                return false;
            }
            if (zpos < 0 || zpos >= GA::nzGl)
            {
                std::cerr<<"gamConst::TowerFromSYZ: out of range zpos = "<<zpos<<std::endl;
                return false;
            }
            towid = (sec * GA::nzGl * GA::nyGl) + (ypos * GA::nzGl) + zpos;
            return true;
        }
        else if (sec >=2 && sec < GA::nsec)
        {
            if (ypos < 0 || ypos >= GA::nySc)
            {
                std::cerr<<"gamConst::TowerFromSYZ: out of range ypos = "<<ypos<<std::endl;
                return false;
            }
            if (zpos < 0 || zpos >= GA::nzSc)
            {
                std::cerr<<"gamConst::TowerFromSYZ: out of range zpos = "<<zpos<<std::endl;
                return false;
            }
            towid = (2 * GA::nzGl * GA::nyGl) + ((sec - 2) * GA::nzSc * GA::nySc) + (ypos * GA::nzSc) + zpos;
            return true;
        }
        else
        {
            std::cerr<<"gamConst::TowerFromSYZ: out of range sector = "<<sec<<std::endl;
            return false;
        }
    }

    bool TowerToSYZ(int towid, int& sec, int& ypos, int& zpos)
    {
        if (towid < 0 || towid >= 2*GA::nzGl*GA::nyGl + 6*GA::nzSc*GA::nySc)
        {
            std::cerr<<"gamConst::TowerToSYZ: out of range tower id = "<<towid<<std::endl;
            sec = -9999;
            ypos = -9999;
            zpos = -9999;
            return false;
        }

        int id = towid;

        if (id >= 2*GA::nzGl*GA::nyGl)
        {
            id -= 2 * GA::nzGl * GA::nyGl;
            sec = id/(GA::nzSc * GA::nySc) + 2;
            ypos = (id - (sec - 2) * GA::nzSc * GA::nySc)/GA::nzSc;
            zpos = (id - (sec - 2) * GA::nzSc * GA::nySc) - ypos * GA::nzSc;
        }
        else
        {
            sec = id/(GA::nzGl * GA::nyGl);
            ypos = (id - sec * GA::nzGl * GA::nyGl)/GA::nzGl;
            zpos = (id - sec * GA::nzGl * GA::nyGl) - ypos * GA::nzGl;
        }

        return true;
    }

    bool IsSectorGl(int isec)
    {
        if (isec < 0 || isec >= nsec)
        {
            std::cerr<<"gamConst.cc : isec out of bounds: "<<isec<<std::endl;
            return false;
        }
        if (isec < 2) return true;
        else return false;
    }

    bool IsSectorSc(int isec)
    {
        if (isec < 0 || isec >= nsec)
        {
            std::cerr<<"gamConst.cc : isec out of bounds: "<<isec<<std::endl;
            return false;
        }
        if (isec < 2) return false;
        else return true;
    }

    int GetNzNy(int& nz, int& ny, int isec)
    {
        if (IsSectorGl(isec))
        {
            nz = GA::nzGl;
            ny = GA::nyGl;
            return 0;
        }
        else if (IsSectorSc(isec))
        {
            nz = GA::nzSc;
            ny = GA::nySc;
            return 1;
        }
        return -1;
    }

    int GetTowZBin(int isec, int iz)
    {
        if (IsSectorSc(isec)) return GetTowBin(iz,nzSc,nztowbin);
        if (IsSectorGl(isec)) return GetTowBin(iz,nzGl,nztowbin);
        return -999;
    }

    int GetTowYBin(int isec, int iy)
    {
        if (IsSectorSc(isec)) return GetTowBin(iy,nySc,nytowbin);
        if (IsSectorGl(isec)) return GetTowBin(iy,nyGl,nytowbin);
        return -999;
    }

    int GetTowZBinTowID(int towid)
    {
        int isec = -999;
        int iy = -999;
        int iz = -999;
        bool success = TowerToSYZ(towid,isec,iy,iz);
        if (!success) return -999;
        return GetTowZBin(isec,iz);
    }

    int GetTowYBinTowID(int towid)
    {
        int isec = -999;
        int iy = -999;
        int iz = -999;
        bool success = TowerToSYZ(towid,isec,iy,iz);
        if (!success) return -999;
        return GetTowYBin(isec,iy);
    }

    int GetTowBin(int itow, int ntow, int nbin)
    {
        int binsize = ntow/nbin;
        int ibin = itow/binsize;
        if (ibin < 0 || ibin >= nbin)
        {
            std::cerr<<"gamConst.cc : ibin out of bounds: "<<ibin<<std::endl;
            return -999;
        }
        return ibin;
    }

    int GetEcoreBin(float ecore)
    {
        for (int ibin = 0; ibin < necore; ibin++)
        {
            if (ecore < limecore[ibin])
            {
                return ibin-1;
            }
        }

        return (necore-1);
    }
};
