void draw_chargedvetoplots()
{
    TCanvas* c1 = new TCanvas("c1","c1",1200,675);
    c1->Divide(3,2);

    gStyle->SetOptStat(1111111);

    int i = 0;

    i++;
    c1->cd(i);
    gPad->SetLogy();
    pc3dphis_trg_type0->Rebin(10);
    pc3dphis_trg_type0->Draw();

    i++;
    c1->cd(i);
    gPad->SetLogy();
    pc3dz_trg_type0->Rebin(10);
    pc3dz_trg_type0->Draw();

    i++;
    c1->cd(i);
    gPad->SetLogy();
    pc3dr_trg_type0->Rebin(10);
    pc3dr_trg_type0->Draw();

    i++;
    c1->cd(i);
    gPad->SetLogy();
    trkdphis_trg_type0->Rebin(10);
    trkdphis_trg_type0->Draw();

    i++;
    c1->cd(i);
    gPad->SetLogy();
    trkdz_trg_type0->Rebin(10);
    trkdz_trg_type0->Draw();

    i++;
    c1->cd(i);
    gPad->SetLogy();
    trkdr_trg_type0->Rebin(10);
    trkdr_trg_type0->Draw();
}
