void projectandnormalize(const char* hname = "centNpart")
{
    TH2F* h2 = (TH2F*) gROOT->FindObject(hname);
    
    int nx = h2->GetNbinsX();
    int ny = h2->GetNbinsY();

    TH1F* h1 = new TH1F("h1","h1",nx,0,100);

    cout<<"Nbins: x = "<<nx<<"; y = "<<ny<<endl;

    for (int ix = 1; ix <= nx; ix++)
    {
        float nentry = 0.;
        float bc = 0.;
        for (int iy = 1; iy <= ny; iy++)
        {       
            float tmpentry = h2->GetBinContent(ix,iy);
            nentry += tmpentry;
            float ycenter = h2->GetYaxis()->GetBinCenter(iy);
            bc += tmpentry*ycenter;
        }
        if (nentry != 0) bc = bc/nentry;
        h1->SetBinContent(ix,bc);
    }
    h1->Draw();
}
