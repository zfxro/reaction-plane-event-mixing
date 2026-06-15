// flags:
// -2 :: canvas for each
// -1 :: draw same
// >= 0 :: use flag as index (no looping)

// REBIN: change nptt if using rebinned ptt

#define ntype 3
// keyword: REBIN
#define nptt 2
#define ncent 4
#define ndpsi 2
#define nxi 3

#define SAME -1
#define CAN -2

//warning. pi0 at end. unlike inc at end elsewhere
string typestr[ntype+1] = {"dir","dec","inc","pi0"};

// keyword: REBIN
const float limptt[nptt+1] = {5.,9.,15.};
const int limcent[ncent] = {0,20,60,100};
const float limxi[nxi+1] = {0.0,0.8,1.6,2.4};

int ican;

int getindex(int flag, int ican, int isam);
void makename(ostringstream& oss, int flag, string str);

void aspiringArtist(int ftype = 0, int fptt = -1, int fcent = -2, int fdpsi = 2, bool isadj = true, const char* inname = "artist_out.root")
{
    gStyle->SetOptStat(0);

    TFile* input = new TFile(inname,"READ");

    vector<int> flags;
    flags.push_back(ftype);
    flags.push_back(fptt);
    flags.push_back(fcent);
    flags.push_back(fdpsi);

    vector<int> nvec;
    nvec.push_back(ntype);
    nvec.push_back(nptt);
    nvec.push_back(ncent);
    nvec.push_back(ndpsi);

    int ncan = -999;
    int nsam = -999;

    for (vector<int>::iterator it = flags.begin(); it < flags.end(); it++)
    {
        if (*it == CAN)
            ncan = nvec[it - flags.begin()];
        else if (*it == SAME)
            nsam = nvec[it - flags.begin()];
    }

    TLine* tl = new TLine(0.0,0.0,1.0,0.0);
    tl->SetLineStyle(9);
    tl->SetLineWidth(2);

    for (ican = 0; ican < ncan; ican++)
    {
        ostringstream osscan;
        makename(osscan,ftype,"type");
        makename(osscan,fptt,"ptt");
        makename(osscan,fcent,"cent");
        makename(osscan,fdpsi,"dpsi");

        cout<<osscan.str()<<endl;

        TCanvas* tcan = new TCanvas(osscan.str().c_str(),"",1300,1500);
        tcan->Divide(1,3,0.001,0.001);

        for (int ixi = 0; ixi < nxi; ixi++)
        {
            THStack* stk = new THStack(Form("stk%i%i",ican,ixi),Form("%.1f < #xi < %.1f",limxi[ixi],limxi[ixi+1]));

            TPad* pad = (TPad*) tcan->cd(fabs(ixi-nxi));

            for (int isam = 0; isam < nsam; isam++)
            {
                int itype = getindex(ftype,ican,isam);
                int iptt = getindex(fptt,ican,isam);
                int icent = getindex(fcent,ican,isam);
                int idpsi = getindex(fdpsi,ican,isam);

                ostringstream ossh;
                ossh<<"dphi_";
                if (isadj)
                    ossh<<"adj_";
                ossh<<typestr[itype]<<"_mix2_ptt"<<iptt<<"_cent"<<icent<<"_dpsi"<<idpsi<<"_xi"<<ixi;

                TH1F* h1 = (TH1F*) gROOT->FindObject(ossh.str().c_str());
                if (!h1)
                {
                    cerr<<"Could not find "<<ossh.str()<<endl;
                    return;
                }

                cout<<"Found "<<ossh.str()<<endl;
                h1->Print();

                h1->SetMarkerStyle(24+isam);
                h1->SetMarkerColor(2+isam);
                h1->SetLineColor(2+isam);
                h1->SetMarkerSize(1.5);

                ostringstream osshn;
                if (ftype == SAME)
                {
                    osshn<<typestr[isam];
                }
                else if (fptt == SAME)
                {
                    osshn<<limptt[isam]<<" < p^{#gamma}_{T} < "<<limptt[isam+1];
                }
                else if (fcent == SAME)
                {
                    if (isam == ncent-1)
                        osshn<<"0 - 100%";
                    else
                        osshn<<limcent[isam]<<" - "<<limcent[isam+1]<<"%";
                }
                else if (fdpsi == SAME)
                {
                    if (isam == 0)
                        osshn<<"#frac{-#pi}{4}<#Delta#Psi<#frac{-#pi}{4}";
                    if (isam == 1)
                        osshn<<"#frac{#pi}{2}<#Delta#Psi<#frac{3#pi}{4}";
                    if (isam == 2)
                        osshn<<"0<#Delta#Psi<#pi";
                }

                h1->SetTitle(osshn.str().c_str());

                stk->Add(h1);
            }

            float x0 = 0.4;

            stk->Draw("nostack");
            gPad->BuildLegend(x0,0.7,x0+0.2,0.9,"");
            tl->Draw();
            pad->Update();
        }
        ostringstream osspic;
        osspic<<osscan.str()<<".png";
        tcan->SaveAs(osspic.str().c_str());
    }
}

int getindex(int flag, int ican, int isam)
{
    if (flag == -2) return ican;
    if (flag == -1) return isam;
    return flag;
}

void makename(ostringstream& oss, int flag, string str)
{
    if (flag == -1)
    {
        std::transform(str.begin(), str.end(),str.begin(), ::toupper);
        oss<<str;
        return;
    }

    else if (flag == -2)
    {
        oss<<str<<ican;
    }

    else
    {
        oss<<str<<flag;
    }
}
