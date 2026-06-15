// flags:
// -2 :: canvas for each
// -1 :: draw same
// >= 0 :: use flag as index (no looping)

// REBIN: change nptt if using rebinned ptt

#define ntype 4
#define ndpsi 3
// keyword: REBIN
#define nptt 2
#define ncent 4 //includes MB
#define nxi 2

#define SAME -1
#define CAN -2

string typestr[ntype] = {"dir","dec","pi0","inc"};

const float limptt[nptt+1] = {5.,9.,15.};
const int limcent[ncent] = {0,20,60,100};
const float limxi[nxi+1] = {0.0,1.2,2.4};

int ican;

void makename(ostringstream& oss, int flag, string str);

void aspiringArtist_NoLoop(int itype = 0, int iptt = 0, int icent = 0, int idpsi = 2, bool isadj = true, const char* inname = "artist_out.root")
{
    gStyle->SetOptStat(0);

    TFile* input = new TFile(inname,"READ");

    TLine* tl = new TLine(0.0,0.0,1.0,0.0);
    tl->SetLineStyle(9);
    tl->SetLineWidth(2);

    ostringstream osscan;

    makename(osscan,itype,"type");
    makename(osscan,iptt,"ptt");
    makename(osscan,icent,"cent");
    makename(osscan,idpsi,"dpsi");

    cout<<osscan.str()<<endl;

    TCanvas* tcan = new TCanvas(osscan.str().c_str(),"",1300,1500);
    tcan->Divide(2,3,0.001,0.001);

    for (int ixi = 0; ixi < nxi; ixi++)
    {
        TPad* pad = (TPad*) tcan->cd(fabs(ixi-nxi));

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

        h1->SetTitle(Form("%.1f < #xi < %.1f",limxi[ixi],limxi[ixi+1]));

        h1->SetMarkerStyle(24);
        h1->SetMarkerColor(2);
        h1->SetLineColor(2);
        h1->SetMarkerSize(1.5);

        h1->Draw();
        tl->Draw();

        pad->Update();

    }

        ostringstream osspic;
        osspic<<osscan.str()<<".png";
        tcan->SaveAs(osspic.str().c_str());
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
