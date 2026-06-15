void CalcIncGamPtSpectra(
    const char* filename =
    "EMCalClusAna.root",
    const char* outfilename =
    "IncGamSpectra.root"
    )
{

    //check bin widths, limits etc match with what is done below.

    const int ncent = 4;
    const int limcent[ncent+1] = {0,20,40,60,100};

    TFile* file = new TFile(filename,"READ");
    TH1F* h1cent = (TH1F*) gROOT->FindObject("cent");
    TH2F* h2ptcent = (TH2F*) gROOT->FindObject("ptcent");
    TH2F* h2ecent = (TH2F*) gROOT->FindObject("ecent");

    h1cent->SetDirectory(0);
    h2ptcent->SetDirectory(0);
    h2ecent->SetDirectory(0);

    file->Close();

    TFile* outfile = new TFile(outfilename,"RECREATE");
    TH1F* h1pt[ncent];
    TH1F* h1e[ncent];

    h1cent->SetDirectory(outfile);
    h2ptcent->SetDirectory(outfile);
    h2ecent->SetDirectory(outfile);

    int ncentbin = h2ptcent->GetNbinsY();
    int centbw = 100/ncentbin;

    cout<<"ncentbin = "<<ncentbin<<"; centbw = "<<centbw<<endl;

    for (int icent = 0; icent < ncent; icent++)
    {
        int firstbin = limcent[icent]/centbw + 1;
        int lastbin = limcent[icent+1]/centbw;
        printf("icent = %i; centrality: {%i,%i}; firstbin = %i; lastbin = %i;\n",icent,limcent[icent],limcent[icent+1],firstbin,lastbin);

        h1pt[icent] = (TH1F*) h2ptcent->ProjectionX(Form("pt%i",icent),firstbin,lastbin,"e");
        h1e[icent] = (TH1F*) h2ecent->ProjectionX(Form("e%i",icent),firstbin,lastbin,"e");

        //bin width for h1cent is 1; range is -0.5, 100.5
        double nevt = h1cent->Integral(limcent[icent]+1,limcent[icent+1]);
        printf("nevt = %f\n",nevt);

        h1pt[icent]->Scale(1./nevt);
        h1e[icent]->Scale(1./nevt);
    }
    outfile->Write();
    outfile->Close();
}
