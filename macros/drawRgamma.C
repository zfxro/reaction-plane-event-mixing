void drawRgamma(const char* inname = "artist_out.root")
{
    gStyle->SetTitleFontSize(0.08);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetPadBottomMargin(0.2);
    gStyle->SetPadTopMargin(0.1);
        
    TCanvas* can = new TCanvas("can","",1300,600);
    can->Divide(2,2,0.001,0.001);

    TFile* input = new TFile(inname,"READ");

    const int ncent = 3;

    int limcent[4] = {0,20,60,92};
        
    for (int icent = 0; icent <= ncent; icent++)
    {
        TH1F* hinc = (TH1F*) gROOT->FindObject(Form("pt_adj_inc_cent%i_dpsi2",icent));
        TH1F* hdec = (TH1F*) gROOT->FindObject(Form("pt_adj_dec_cent%i_dpsi2",icent));

        hinc->Divide(hdec);

        int locent;
        int hicent;
        if (icent == ncent)
        {
            locent = limcent[0];
            hicent = limcent[ncent];
        }
        else
        {
            locent = limcent[icent];
            hicent = limcent[icent+1];
        }

        hinc->SetTitle(Form("%i - %i %%;p_{T};R_{#gamma}",locent,hicent));
        hinc->GetYaxis()->SetTitleOffset(0.5);
        hinc->GetYaxis()->SetRangeUser(0.,5.);
        hinc->SetTitleSize(0.08,"X");
        hinc->SetTitleSize(0.08,"Y");
        hinc->SetLabelSize(0.06,"X");
        hinc->SetLabelSize(0.06,"Y");

        can->cd(icent+1);

        hinc->DrawCopy();

    }

    can->SaveAs("Rgamma.eps");
}
