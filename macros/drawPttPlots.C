//open gamAnaTree merged file and run

void drawPttPlots()
{

    int ndpsi = 3;
    int ncent = 3;

    THStack* st_dir = new THStack("st_dir","Tagged Direct");
    THStack* st_dec = new THStack("st_dec","Tagged Decay");
    THStack* st_adj_dir  = new THStack("st_adj_dir","Corrected Direct");
    THStack* st_adj_dec  = new THStack("st_adj_dec","Corrected Decay");

    for (int idpsi = ndpsi-1; idpsi >= 0; idpsi--)
    {
        for (int icent = ncent-1; icent >= 0; icent--)
        {
            TH1F* pt_dir = (TH1F*) gROOT->FindObject(Form("pt_trg_dir_centb%i_dpsib%i",icent,idpsi));
            float scl_dir = pt_dir->GetBinContent(1);
            pt_dir->Scale(1./scl_dir);
            st_dir->Add(pt_dir);

            TH1F* pt_dec = (TH1F*) gROOT->FindObject(Form("pt_trg_dec_centb%i_dpsib%i",icent,idpsi));
            float scl_dec = pt_dec->GetBinContent(1);
            pt_dec->Scale(1./scl_dec);
            st_dec->Add(pt_dec);

            TH1F* pt_adj_dir = (TH1F*) gROOT->FindObject(Form("pt_trg_adj_dir_centb%i_dpsib%i",icent,idpsi));
            float scl_adj_dir = pt_adj_dir->GetBinContent(1);
            pt_adj_dir->Scale(1./scl_adj_dir);
            st_adj_dir->Add(pt_adj_dir);

            TH1F* pt_adj_dec = (TH1F*) gROOT->FindObject(Form("pt_trg_adj_dec_centb%i_dpsib%i",icent,idpsi));
            float scl_adj_dec = pt_adj_dec->GetBinContent(1);
            pt_adj_dec->Scale(1./scl_adj_dec);
            st_adj_dec->Add(pt_adj_dec);
        }
    }

    TCanvas* can = new TCanvas("can","",1300,600);
    can->Divide(2,2);

    TPad* pad = NULL;

    pad = (TPad*) can->cd(1);
    pad->SetLogy();
    st_dir->Draw("plc nostack");
    gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

    pad = (TPad*) can->cd(2);
    pad->SetLogy();
    st_dec->Draw("plc nostack");
    gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

    pad = (TPad*) can->cd(3);
    pad->SetLogy();
    st_adj_dir->Draw("plc nostack");
    gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

    pad = (TPad*) can->cd(4);
    pad->SetLogy();
    st_adj_dec->Draw("plc nostack");
    gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

}
