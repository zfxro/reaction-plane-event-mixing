// reads two output of artist_gamanatree.C
//     1) artist input = merged gamanatree no psimix
//     2) artist input = merged gamanatree yes psimix

// produces histograms that will be used (in conjunction with rp reso) to subtract remaining flow background

#include "gamAnaTree.h"
#include "artist_gamAnaTree.h"

void getFlow(
    const char* inname_nopsimix = "",
    const char* inname_psimix = "",
    const char* outname = "getFlow.root"
    )
{

    const int mixtype = MIX;

    TFile* output = new TFile(outname,"RECREATE");

    TH1F* h1flow[ntype+1][nptt][ncent+1][npta];

    //only use msn and dphi: [0,1] histos
    for (int itype = 0; itype < ntype + 1; itype ++)
    {
        if (itype==PI0) continue;

        for (int iptt = 0; iptt < nptt; iptt++)
        {
            for (int ipta = 0; ipta < npta; ipta++)
            {
                for (int icent = 0; icent < ncent+1; icent++)
                {
                    h1flow[itype][iptt][icent][ipta] = new TH1F(
                        Form("flow_%s_ptt%i_cent%i_pta%i",spec[itype],iptt,icent,ipta),
                        Form("flow_%s_ptt%i_cent%i_pta%i",spec[itype],iptt,icent,ipta),
                        ndphibin,mindphi,maxdphi);
                }
            }
        }
    }

    TFile* input_nopsimix = new TFile(inname_nopsimix,"READ");
    if (!input_nopsimix)
    {
        cout<<"Could not open file!"<<endl;
        return;
    }

    TH1F* h1dphi_nopsimix[ntype+1][nptt][ncent+1][npta];

    //only use msn and dphi: [0,1] histos
    for (int itype = 0; itype < ntype + 1; itype ++)
    {
        if (itype==PI0) continue;

        for (int iptt = 0; iptt < nptt; iptt++)
        {
            for (int ipta = 0; ipta < npta; ipta++)
            {
                for (int icent = 0; icent < ncent+1; icent++)
                {
                    TH1F* h1tmp = (TH1F*) gROOT->FindObject(Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",spec[itype],mixtype,iptt,icent,ndpsi,ipta));

                    if (!h1tmp)
                    {
                        printf("Could not find dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i\n",spec[itype],mixtype,iptt,icent,ndpsi,ipta);
                        return;
                    }

                    h1dphi_nopsimix[itype][iptt][icent][ipta] = (TH1F*) h1tmp->Clone(Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i_nopsimix",spec[itype],mixtype,iptt,icent,ndpsi,ipta));
                }
            }
        }
    }

    TFile* input_psimix = new TFile(inname_psimix,"READ");
    if (!input_psimix)
    {
        cout<<"Could not open file!"<<endl;
        return;
    }

    TH1F* h1dphi_psimix[ntype+1][nptt][ncent+1][npta];

    //only use msn and dphi: [0,1] histos
    for (int itype = 0; itype < ntype + 1; itype ++)
    {
        if (itype==PI0) continue;

        for (int iptt = 0; iptt < nptt; iptt++)
        {
            for (int ipta = 0; ipta < npta; ipta++)
            {
                for (int icent = 0; icent < ncent+1; icent++)
                {
                    TH1F* h1tmp = (TH1F*) gROOT->FindObject(Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",spec[itype],mixtype,iptt,icent,ndpsi,ipta));

                    if (!h1tmp)
                    {
                        printf("Could not find dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i\n",spec[itype],mixtype,iptt,icent,ndpsi,ipta);
                        return;
                    }

                    h1dphi_psimix[itype][iptt][icent][ipta] = (TH1F*) h1tmp->Clone(Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i_psimix",spec[itype],mixtype,iptt,icent,ndpsi,ipta));
                }
            }
        }
    }

    TCanvas* canvas[ntype+1][nptt][ncent+1][npta];
    gStyle->SetOptStat(0);

    //do the subtraction
    for (int itype = 0; itype < ntype + 1; itype ++)
    {
        if (itype==PI0) continue;

        for (int iptt = 0; iptt < nptt; iptt++)
        {
            for (int ipta = 0; ipta < npta; ipta++)
            {
                for (int icent = 0; icent < ncent+1; icent++)
                {
                    TH1F* h1tmp = h1flow[itype][iptt][icent][ipta];
                    TH1F* h1no = h1dphi_nopsimix[itype][iptt][icent][ipta];
                    TH1F* h1yes = h1dphi_psimix[itype][iptt][icent][ipta];

                    if (mixtype==MSN)
                    {
                        h1tmp->Add(h1no);
                        h1tmp->Add(h1yes,-1.);
                    }
                    else if (mixtype==MIX)
                    {
                        h1tmp->Add(h1yes);
                        h1tmp->Add(h1no,-1.);
                    }

                    TCanvas* can = canvas[itype][iptt][icent][ipta];
                    can = new TCanvas(Form("can_%s_ptt%i_cent%i_pta%i",spec[itype],iptt,icent,ipta),"",1500,800);
                    can->Divide(3,1);

                    can->cd(1);
                    h1no->Draw();

                    can->cd(2);
                    h1yes->Draw();

                    can->cd(3);
                    h1tmp->Draw();

                    can->SaveAs(Form("getFlow_%s_ptt%i_cent%i_pta%i.png",spec[itype],iptt,icent,ipta));

                }
            }
        }
    }

    output->Write();
    output->Close();

    return;
}
