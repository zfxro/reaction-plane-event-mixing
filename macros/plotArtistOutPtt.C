void plotArtistOutPtt()
{
    int ncent = 3;
    for (int icent = 0; icent < ncent+1; icent++)
    {
        TCanvas* c1 = new TCanvas(Form("c%i",icent),"",1500,875);
        c1->Divide(3,2);

        c1->cd(1);
        TH1F* pt_tag_dir = (TH1F*) gROOT->FindObject(Form("pt_tag_dir_cent%i_dpsi2",icent));
        pt_tag_dir->DrawCopy();

        c1->cd(2);
        TH1F* pt_tag_dec = (TH1F*) gROOT->FindObject(Form("pt_tag_dec_cent%i_dpsi2",icent));
        pt_tag_dec->DrawCopy();

        c1->cd(3);
        TH1F* pt_tag_inc = (TH1F*) gROOT->FindObject(Form("pt_tag_inc_cent%i_dpsi2",icent));
        pt_tag_inc->DrawCopy();

        //==================================================
        c1->cd(4);
        TH1F* pt_adj_dir = (TH1F*) gROOT->FindObject(Form("pt_adj_dir_cent%i_dpsi2",icent));
        pt_adj_dir->DrawCopy();

        c1->cd(5);
        TH1F* pt_adj_dec = (TH1F*) gROOT->FindObject(Form("pt_adj_dec_cent%i_dpsi2",icent));
        pt_adj_dec->DrawCopy();

        //==================================================
        c1->cd(6);
        pt_tag_inc->Divide(pt_adj_dec);
        pt_tag_inc->DrawCopy();
        //==================================================
    }
}
