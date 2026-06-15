// flags:
// -2 :: canvas for each
// -1 :: draw same
// >= 0 :: use flag as index (no looping)

// REBIN: change nptt if using rebinned ptt

#define ntype 3
// keyword: REBIN
#define nptt 4
#define ncent 3
#define ndpsi 3
#define nxi 2

#define SAME -1
#define CAN -2

string typestr[ntype] = {"dir","dec","inc"};

const float limptt[nptt+1] = {5.,7.,9.,12.,15.};
const int limcent[ncent] = {0,40,100};
const float limxi[nxi+1] = {0.0,0.8,1.6,2.4};

int ican;

int getindex(int flag, int ican, int isam);
void makename(ostringstream& oss, int flag, string str);

void aspiringArtist_RemBG(bool isadj = true)
{
    gStyle->SetOptStat(0);

    TFile* input = new TFile("artist_out.root","READ");

    TLine* tl = new TLine(0.0,0.0,1.0,0.0);
    tl->SetLineStyle(9);
    tl->SetLineWidth(2);

    for (int itype = 0; itype < ntype; itype++)
    {
        for (int iptt = 0; iptt < nptt; iptt++)
        {
            for (int icent = 0; icent < ncent; icent++)
            {
                for (int idpsi = 0; idpsi < ndpsi; idpsi++)
                {

                    ostringstream osscan;

                    makename(osscan,itype,"type");
                    makename(osscan,iptt,"ptt");
                    makename(osscan,icent,"cent");
                    makename(osscan,idpsi,"dpsi");

                    cout<<osscan.str()<<endl;

                    TCanvas* tcan = new TCanvas(osscan.str().c_str(),"",1600,1000);
                    tcan->Divide(1,2,0.001,0.001);

                    for (int ixi = 0; ixi < nxi; ixi++)
                    {
                        int ipad = fabs(ixi-nxi);
                        TPad* pad = (TPad*) tcan->cd(ipad);

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
                        h1->SetMarkerColor(4);
                        h1->SetLineColor(4);
                        h1->SetMarkerSize(1.5);

                        h1->Draw();
                        tl->Draw();

                        ossh.str("");
                        ossh<<"dphi_xi_bg_"<<typestr[itype]<<"_ptt"<<iptt<<"_cent"<<icent<<"_dpsi"<<idpsi<<"_xi"<<ixi;

                        h1 = NULL;
                        h1 = (TH1F*) gROOT->FindObject(ossh.str().c_str());
                        if (!h1)
                        {
                            cerr<<"Could not find "<<ossh.str()<<endl;
                            return;
                        }

                        cout<<"Found "<<ossh.str()<<endl;
                        h1->Print();

                        h1->SetTitle("Subtracted Remaining Background");

                        h1->SetMarkerStyle(24);
                        h1->SetMarkerColor(3);
                        h1->SetLineColor(3);
                        h1->SetMarkerSize(1.5);

                        h1->Draw("same");

                    }

                    ostringstream osspic;
                    osspic<<"pics/"<<osscan.str()<<".png";
                    tcan->SaveAs(osspic.str().c_str());
                }
            }
        }
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
