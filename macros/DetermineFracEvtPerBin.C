void DetermineFracEvtPerBin(
    const char* inname = "gamAna_histos.root",
    const char* outname = "fracevtperbin_requireprecise.txt"
    )
{
    TFile* infile = new TFile(inname);

    //it seems evt raw is less than the number of events recorded in database
    TH1F* h1evtraw = (TH1F*) gROOT->FindObject("evtraw");
    float nevtraw = h1evtraw->GetEntries();

    TH1F* h1evtinrun = (TH1F*) gROOT->FindObject("evtinrun");
    float nevtinrun = h1evtinrun->GetBinContent(1);
    nevtinrun = nevtinrun/(h1evtinrun->GetEntries());//divide by number of segments
    nevtinrun = nevtinrun * 331.;//multiply by number of runs

    float nevt = nevtinrun;

    //should we use as is or take projection for increased stats?
    TH2F* h2bz0cent = (TH2F*) gROOT->FindObject("bz0cent");

    TH1F* h1bz0 = (TH1F*) h2bz0cent->ProjectionX("bz0");
    TH1F* h1cent = (TH1F*) h2bz0cent->ProjectionY("cent");

    int nz = h1bz0->GetNbinsX();
    int ncent = h1cent->GetNbinsX();

    float nent_z = h1bz0->GetEntries();

    ofstream* ofs_pro = new ofstream("fracperevt_pro.txt");
    for (int iz = 0; iz < nz; iz++)
    {
        float f_z = h1bz0->GetBinContent(iz+1);
        f_z = f_z/nent_z;
        //"good" = passed event cuts (require precise vertex +- 10)
        for (int icent = 0; icent < ncent; icent++)
        {
            float f_c = h1cent->GetBinContent(icent+1);
            f_c = f_c/nevt;

            float f_evt = f_z * f_c;

            *ofs_pro<<iz<<" "<<icent<<" "<<f_evt<<endl;
        }
    }

    ofstream* ofs = new ofstream("fracperevt.txt");
    for (int iz = 0; iz < nz; iz++)
    {
        for (int icent = 0; icent < ncent; icent++)
        {
            int binno = h2bz0cent->GetBin(iz+1,icent+1);
            float f_evt = h2bz0cent->GetBinContent(binno);
            f_evt = f_evt/nevt;

            *ofs<<iz<<" "<<icent<<" "<<f_evt<<endl;
        }
    }
}
