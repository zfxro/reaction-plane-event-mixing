#include "gamAnaTree.cc"

void mergeinvm(
    const char* inname,
    const char* outname,
    const char* tag
    )
{
    cout<<"WARNING!! make sure all cent and pt bin limits are integer!"<<endl;

    TFile* outfile = new TFile(outname,"RECREATE");
    TH2F* h2ptcent_rebin = new TH2F("ptcent_rebin","ptcent_rebin",nptt,0,nptt,ncent+1,0,ncent+1);
    TH1F* h1invm_rebin[ncent+1][nptt];
    for (int icent = 0; icent < ncent+1; icent++)
    {
        for (int ipt = 0; ipt < nptt; ipt++)
        {
            h1invm_rebin[icent][ipt] = new TH1F(Form("invm_rebin_cent%i_pt%i",icent,ipt),Form("invm_rebin_cent%i_pt%i",icent,ipt),125,.05,.3);
        }
    }

    TFile* infile = new TFile(inname,"READ");
    if (!infile || infile->IsZombie())
    {
        cerr<<"Difficulty reading input."<<endl;
        return;
    }


    for (int icent = 0; icent < 95; icent++)
    {
        for (int ipt = 0; ipt < 10; ipt++)
        {
            TH1F* h1invm = (TH1F*) gROOT->FindObject(Form("invm%s_%i_%i",tag,icent,ipt));

            float cent = icent;
            float pt = 5. + ipt;

            int icentrebin = GetCentBin(cent);
            int iptrebin = GetPtTrigBin(pt);

            printf("icent %i   cent %.2f   icentrebin %i   ipt %i   pt %.2f   iptrebin %i\n",icent,cent,icentrebin,ipt,pt,iptrebin);

            h1invm_rebin[icentrebin][iptrebin]->Add(h1invm);
            h1invm_rebin[ncent][iptrebin]->Add(h1invm);

            delete h1invm;
        }
    }

    TH2F* h2ptcent = (TH2F*) gROOT->FindObject(Form("ptcent%s",tag));
    int ncentbin = h2ptcent->GetNbinsY();
    int nptbin = h2ptcent->GetNbinsX();
    for (int icent = 0; icent < ncentbin; icent++)
    {
        for (int ipt = 0; ipt < nptbin; ipt++)
        {
            float pt = h2ptcent->GetXaxis()->GetBinCenter(ipt+1);
            if (pt < limptt[0] || pt >= limptt[nptt]) continue;
            float cent = h2ptcent->GetYaxis()->GetBinCenter(icent+1);
            float binvalue = h2ptcent->GetBinContent(ipt+1,icent+1);

            int icentrebin = GetCentBin(cent);
            int iptrebin = GetPtTrigBin(pt);

            printf("icent %i   cent %.2f   icentrebin %i   ipt %i   pt %.2f   iptrebin %i\n",icent,cent,icentrebin,ipt,pt,iptrebin);

            h2ptcent_rebin->Fill(iptrebin,icentrebin,binvalue);
            h2ptcent_rebin->Fill(iptrebin,ncent,binvalue);
        }
    }


    cout<<"Closing input file..."<<endl;
    infile->Close();


    cout<<"Writing output file "<<outfile->GetName()<<endl;
    outfile->Write();


    cout<<"Closing output file..."<<endl;
    outfile->Close();


    cout<<"Successful completion!"<<endl;
}
