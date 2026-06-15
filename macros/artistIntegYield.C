// Input: list of root files from artist_gamAnaTree.C
//        these files differ in their dphi binning
//
// dPhi Binning: dphi histos are binned such that the last bin is the only one of interest
//               integrated yields of different dphi regions will be plotted together
//
// Comparisons: a) Peripheral to Central
//                      - all-plane histos used here
//              b) In vs Out of Plane
//                      - 20-60% centrality used here (where we expect largest effect)
//
// Graphs:
//          X-axis = Xi (0-1.2,1.2,2.4)
//          Y-axis = above ratios
//          New canvas for each trigger pt and ratio (a or b) configuration
//          dPhi integration regions plotted together (THStack)

namespace GAT {
    #include "gamAnaTree.h"
}

#define nfile 4
#define nptt 2
#define ntype 4
#define nxi 2
#define ncent 3

const int midcent = 1;

const float ximin = 0.;
const float ximax = 2.4;

const float pttlim[nptt+1] = {5.0,9.0,15.0};

const char* file[nfile] = {
    "artist_out_ndphi6.root",
    "artist_out_ndphi4.root",
    "artist_out_ndphi3.root",
    "artist_out_ndphi2.root"
};

const char* name[nfile] = {
    "|#delta#Phi - #pi| < #pi/6",
    "|#delta#Phi - #pi| < #pi/4",
    "|#delta#Phi - #pi| < #pi/3",
    "|#delta#Phi - #pi| < #pi/2",
};


void artistIntegYield(const char* outname = "integYield.root")
{
    TFile* output = new TFile(outname,"RECREATE");

    TH1F* hratcent[ntype][nptt][nfile];
    TH1F* hratdpsi[ntype][nptt][nfile];

    for (int ifile = 0; ifile < nfile; ifile++)
    {
        cout<<"Reading file "<<file[ifile]<<endl;
        TFile* input = new TFile(file[ifile]);
        assert(input);

        for (int itype = 0; itype < ntype; itype++)
        {
            for (int iptt = 0; iptt < nptt; iptt++)
            {
                hratcent[itype][iptt][ifile] = new TH1F(Form("cent_%s_ptt%i_file%i",GAT::spec[itype],iptt,ifile),name[ifile],nxi,ximin,ximax);
                hratcent[itype][iptt][ifile]->SetDirectory(0);

                hratdpsi[itype][iptt][ifile] = new TH1F(Form("dpsi_%s_ptt%i_file%i",GAT::spec[itype],iptt,ifile),name[ifile],nxi,ximin,ximax);
                hratdpsi[itype][iptt][ifile]->SetDirectory(0);

                for (int ixi = 0; ixi < nxi; ixi++)
                {

                    ostringstream pref;
                    pref<<"dphi_";
                    if (itype != 2)
                        pref<<"adj_";
                    pref<<GAT::spec[itype]<<"_mix2_ptt"<<iptt<<"_cent";

                    ostringstream suf;
                    suf<<"_dpsi2_xi"<<ixi;

                    TH1F* hcent = (TH1F*) gROOT->FindObject(Form("%s0%s",pref.str().c_str(),suf.str().c_str()));
                    assert(hcent);

                    TH1F* hperi = (TH1F*) gROOT->FindObject(Form("%s%d%s",pref.str().c_str(),ncent-1,suf.str().c_str()));
                    assert(hperi);

                    hcent->Divide(hperi);

                    int nbin = hcent->GetNbinsX();

                    hratcent[itype][iptt][ifile]->SetBinContent(ixi+1,hcent->GetBinContent(nbin));
                    hratcent[itype][iptt][ifile]->SetBinError(ixi+1,hcent->GetBinError(nbin));

                    pref.str("");
                    pref<<"dphi_";
                    if (itype != 2)
                        pref<<"adj_";
                    pref<<GAT::spec[itype]<<"_mix2_ptt"<<iptt<<"_cent"<<midcent<<"_dpsi";

                    suf.str("");
                    suf<<"_xi"<<ixi;

                    TH1F* hin = (TH1F*) gROOT->FindObject(Form("%s0%s",pref.str().c_str(),suf.str().c_str()));
                    TH1F* hout = (TH1F*) gROOT->FindObject(Form("%s1%s",pref.str().c_str(),suf.str().c_str()));

                    hin->Divide(hout);

                    nbin = hin->GetNbinsX();

                    hratdpsi[itype][iptt][ifile]->SetBinContent(ixi+1,hin->GetBinContent(nbin));
                    hratdpsi[itype][iptt][ifile]->SetBinError(ixi+1,hin->GetBinError(nbin));
                }
            }
        }

        input->Close();
        delete input;
    }

    output->cd();

    THStack* stk_cent[ntype][nptt];
    THStack* stk_dpsi[ntype][nptt];

    for (int itype = 0; itype < ntype; itype++)
    {
        string typestr = "";
        if (itype == 0) typestr = "Direct Photons";
        if (itype == 1) typestr = "Decay Photons";
        if (itype == 2) typestr = "Neutral Pions";
        if (itype == 3) typestr = "Inclusive Photons";

        for (int iptt = 0; iptt < nptt; iptt++)
        {
            stk_cent[itype][iptt] = new THStack(Form("cent_%s_ptt%i",GAT::spec[itype],iptt),Form("Central over Peripheral  %s  %.1f < ptt < %.1f  (All-Plane)",typestr.c_str(),pttlim[iptt],pttlim[iptt+1]));

            stk_dpsi[itype][iptt] = new THStack(Form("dpsi_%s_ptt%i",GAT::spec[itype],iptt),Form("In over Out of Plane  %s  %.1f < ptt < %.1f  (%i-%i%%)",typestr.c_str(),pttlim[iptt],pttlim[iptt+1],GAT::limcent[midcent],GAT::limcent[midcent+1]));

            for (int ifile = 0; ifile < nfile; ifile++)
            {
                hratcent[itype][iptt][ifile]->SetMarkerStyle(24+ifile);
                hratcent[itype][iptt][ifile]->SetMarkerColor(2+ifile);
                hratcent[itype][iptt][ifile]->SetLineColor(2+ifile);
                hratcent[itype][iptt][ifile]->SetMarkerSize(1.5);
                stk_cent[itype][iptt]->Add(hratcent[itype][iptt][ifile]);

                hratdpsi[itype][iptt][ifile]->SetMarkerStyle(24+ifile);
                hratdpsi[itype][iptt][ifile]->SetMarkerColor(2+ifile);
                hratdpsi[itype][iptt][ifile]->SetLineColor(2+ifile);
                hratdpsi[itype][iptt][ifile]->SetMarkerSize(1.5);
                stk_dpsi[itype][iptt]->Add(hratdpsi[itype][iptt][ifile]);

            }

            stk_cent[itype][iptt]->Write();
            stk_dpsi[itype][iptt]->Write();
        }
    }

    cout<<"Writing to "<<output->GetName()<<endl;
    output->Write();
    output->Close();

}
