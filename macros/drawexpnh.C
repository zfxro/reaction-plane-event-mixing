void drawexpnh()
{
    int ncent = 4;

    TCanvas* can = new TCanvas("can","can",1600,900);
    can->Divide(3,2);

    for (int icent = 0; icent <= ncent; icent++)
    {
        can->cd(icent+1)->SetLogy();
        TH1F* h1 = (TH1F*) gROOT->FindObject(Form("expnh%i",icent));
        h1->Draw();
    }
}
