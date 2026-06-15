// tag efficiency (pt) for each integer centrality
// saving fits to file
// also print fits

// icent = 94 is all cent
#define ncent 4
#define ntype 2
#define ndpsi 2
#define nfile_ppg 10

const char* name_ppg[nfile_ppg] = {
    "Rg_0_5",
    "Rg_5_10",
    "Rg_10_15",
    "Rg_15_20",
    "Rg_20_30",
    "Rg_30_40",
    "Rg_40_50",
    "Rg_50_60",
    "Rg_60_93",
    "Rg_0_93"
};

void GetHistos(const char* name, TH1F* h1[ncent]);
void PlotHistos(const char* name, TH1F* h1dir[ncent], TH1F* h1dec[ncent], TH1F* h1rat[ncent], TH1F* h1rg[ncent]);
void GetRgPPG139(const char* indir, TGraphErrors* tgr[nfile_ppg]);
void PlotRg(const char* name, TH1F* h1rg[ncent], TGraphErrors* tgr[nfile_ppg]);

int centlim[ncent][2] = {
    {0,20},
    {20,60},
    {60,93},
    {0,93}
};


void plot_anaTrigAna(
    const char* inname = "out_anaTrigAna.root"
    )
{

    TFile* input = new TFile(inname,"READ");
    assert(input!=NULL);
    input->Print();

    TH1F* h1pt_dir[ncent];
    TH1F* h1pt_dir_adj[ncent];

    TH1F* h1pt_dec[ncent];
    TH1F* h1pt_dec_adj[ncent];
    TH1F* h1pt_dec_eff[ncent];

    GetHistos("pt_trg_dir_",h1pt_dir);
    GetHistos("pt_trg_dir_adj_",h1pt_dir_adj);

    GetHistos("pt_trg_dec_",h1pt_dec);
    GetHistos("pt_trg_adj_",h1pt_dec_adj);
    GetHistos("ptt_deceff_",h1pt_dec_eff);

    //==================================================

    TH1F* h1rat[ncent];
    TH1F* h1rat_adj[ncent];
    TH1F* h1rat_eff[ncent];

    TH1F* h1rg[ncent];
    TH1F* h1rg_adj[ncent];
    TH1F* h1rg_eff[ncent];

    PlotHistos("Raw",h1pt_dir,h1pt_dec,h1rat,h1rg);
    PlotHistos("Adj",h1pt_dir_adj,h1pt_dec_adj,h1rat_adj,h1rg_adj);
    PlotHistos("Eff",h1pt_dir_adj,h1pt_dec_eff,h1rat_eff,h1rg_eff);

    TGraphErrors* rg_ppg[nfile_ppg];
    GetRgPPG139("ppg139_data/PbSc/",rg_ppg);
    for (int ifile = 0; ifile < nfile_ppg; ifile++)
    {
        cout<<"Printing ifile "<<ifile<<endl;
        rg_ppg[ifile]->Print();
    }

    PlotRg("Raw",h1rg,rg_ppg);
    PlotRg("Adj",h1rg_adj,rg_ppg);
    PlotRg("Eff",h1rg_eff,rg_ppg);

    cout<<"Complete!"<<endl;

}

void PlotHistos(const char* name, TH1F* h1dir[ncent], TH1F* h1dec[ncent], TH1F* h1rat[ncent], TH1F* h1rg[ncent])
{
    TCanvas* can = new TCanvas("can","",1400,800);
    can->Divide(1,2);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);

    float ymax = 10.;
    if (strcmp(name,"Raw") == 0)
        ymax = 30.;

    for (int icent = 0; icent < ncent; icent++)
    {
        TPad* pad = (TPad*) can->cd(1);
        pad->SetLogy();
        THStack* hs = new THStack("hs",Form("%s Photons (%i-%i%%)",name,centlim[icent][0],centlim[icent][1]));

        h1dir[icent]->SetMarkerColor(kBlue);
        h1dir[icent]->SetMarkerStyle(kCircle);
        h1dir[icent]->SetLineColor(kBlue);

        TF1* fdir = new TF1(Form("fdir_%i",icent),("[0]*x^[1]"),5.,15.);
        fdir->SetParameter(0,1.);
        fdir->SetParameter(1,-7.);
        fdir->SetLineColor(kBlue);
        fdir->SetLineWidth(1);
        h1dir[icent]->Fit(fdir);


        hs->Add(h1dir[icent]);

        h1dec[icent]->SetMarkerColor(kRed);
        h1dec[icent]->SetMarkerStyle(kCircle);
        h1dec[icent]->SetLineColor(kRed);

        TF1* fdec = new TF1(Form("fdec_%i",icent),("[0]*x^[1]"),5.,15.);
        fdec->SetParameter(0,1.);
        fdec->SetParameter(1,-7.);
        fdec->SetLineColor(kRed);
        fdec->SetLineWidth(1);
        h1dec[icent]->Fit(fdec);

        hs->Add(h1dec[icent]);

        hs->Draw("nostack");
        hs->GetXaxis()->SetRangeUser(5.,15.);
        pad->Update();

        TLegend* tleg = new TLegend(0.7,0.65,0.9,0.9);
        tleg->AddEntry(h1dir[icent],Form("direct %.1e p_{T}^{%.2f #pm %.2f} ",fdir->GetParameter(0),fdir->GetParameter(1),fdir->GetParError(1)));
        tleg->AddEntry(h1dec[icent],Form("decay %.1e p_{T}^{%.2f #pm %.2f} ",fdec->GetParameter(0),fdec->GetParameter(1),fdec->GetParError(1)));
        tleg->Draw();

        pad = (TPad*) can->cd(2);

        h1rat[icent] = (TH1F*) h1dir[icent]->Clone(Form("rat_%s_%i",name,icent));
        h1rat[icent]->Divide(h1dec[icent]);
        h1rat[icent]->SetTitle("Ratio");
        h1rat[icent]->SetMarkerColor(kBlack);
        h1rat[icent]->SetLineColor(kBlack);

        TF1* frat = new TF1(Form("frat_%i",icent),"[0]*x^[1]",5.,15.);
        frat->SetParameter(0,1.);
        frat->SetParameter(1,1.);
        h1rat[icent]->Fit(frat);

        h1rat[icent]->Draw();
        h1rat[icent]->GetYaxis()->SetRangeUser(0.,ymax);

        TLine* tl = new TLine(5.,1.,15.,1.);
        tl->SetLineStyle(7);
        tl->SetLineColor(kBlack);
        tl->Draw();

        pad->Update();
        can->SaveAs(Form("%s_%i.png",name,icent));

        h1rg[icent] = (TH1F*) h1dir[icent]->Clone(Form("rg_%s_%i",name,icent));
        h1rg[icent]->SetTitle("R_{#gamma}");
        h1rg[icent]->Add(h1dec[icent]);
        h1rg[icent]->Divide(h1dec[icent]);

        TF1* frg = new TF1(Form("frg_%i",icent),"1 + [0]*x^[1]",5.,15.);
        frg->SetParameter(0,1.);
        frg->SetParameter(1,1.);
        h1rg[icent]->Fit(frg);

        h1rg[icent]->Draw();
        h1rg[icent]->GetYaxis()->SetRangeUser(0.,ymax);

        pad->Update();
        can->SaveAs(Form("Rg_%s_%i.png",name,icent));

        //avoid drawing later
        h1rg[icent]->GetFunction(Form("frg_%i",icent))->SetBit(TF1::kNotDraw);

        delete hs;
    }

    delete can;
}

void GetHistos(const char* name, TH1F* h1[ncent])
{

    for (int icent = 0; icent < ncent-1; icent++)
    {
        ostringstream oss;
        oss << name << "centb"<<icent<<"_dpsib2"<<endl;;
        cout<<oss.str()<<endl;

        cout<<"(TH1F*) gROOT->FindObject("<<oss.str().c_str()<<");"<<endl;
        cout<<oss.str()<<endl;
        h1[icent] = (TH1F*) gROOT->FindObject(oss.str().c_str());
        assert(h1[icent] != NULL);
        h1[icent]->SetDirectory(0);
        h1[icent]->Rebin(2);
        h1[icent]->GetXaxis()->SetRangeUser(5.,15.);

        if (icent == 0)
        {
            oss.str("");
            oss<<name<<ncent-1;
            h1[ncent-1] = (TH1F*) h1[icent]->Clone(oss.str().c_str());
        }
        else
        {
            h1[ncent-1]->Add(h1[icent]);
        }
    }
}

void GetRgPPG139(const char* indir, TGraphErrors* tgr[nfile_ppg])
{

    vector<float> vpt[nfile_ppg];
    vector<float> vrg[nfile_ppg];
    vector<float> vex[nfile_ppg];;//error on x needed by tgrapherrors
    vector<float> vstat[nfile_ppg];

    for (int ifile = 0; ifile < nfile_ppg; ifile++)
    {
        ostringstream ossfile;
        ossfile << indir << name_ppg[ifile] << ".txt";

        cout << "Reading file \"" << ossfile.str() << "\"" << endl;

        ifstream ifsfile(ossfile.str().c_str(), ifstream::in);
        if (!ifsfile.is_open())
        {
            cerr << "Error opening file." << endl;
            return;
        }


        string line;

        while (getline(ifsfile,line))
        {
            if (!isdigit(line.at(0))) continue;
            istringstream issline(line.data());

            float pt, rg, stat, sysA, sysB, sysC;
            issline >> pt >> rg >> stat >> sysA >> sysB >> sysC;

            vpt[ifile].push_back(pt);
            vrg[ifile].push_back(rg);
            vstat[ifile].push_back(stat);
        }

        //make some errors for x
        for (vector<float>::iterator it = vpt[ifile].begin(); it != vpt[ifile].end(); ++it)
        {
            float tmpex = -999.;
            if (it == vpt[ifile].begin())
            {
                tmpex = *(it+1) - *it;
                tmpex/=2.;
            }
            else
            {
                tmpex = *it - *(it-1);
                tmpex = *it - *(it-1) - vex[ifile][it - vpt[ifile].begin() - 1];
            }
            vex[ifile].push_back(tmpex);
        }
    }

    for (int ifile = 0; ifile < nfile_ppg; ifile++)
    {
        tgr[ifile] = new TGraphErrors(vpt[ifile].size(),&(vpt[ifile][0]),&(vrg[ifile][0]),&(vex[ifile][0]),&(vstat[ifile][0]));
        tgr[ifile]->SetTitle(name_ppg[ifile]);
    }
}

void PlotRg(const char* name, TH1F* h1rg[ncent], TGraphErrors* tgr[nfile_ppg])
{
    int ppglim[ncent][2] = {
        {0,3},
        {4,7},
        {8,8},
        {9,9}
    };

    TCanvas* can = new TCanvas("can","",1400,800);

    for (int icent = 0; icent < ncent; icent++)
    {
        TLegend* tleg = new TLegend(0.7,0.65,0.9,0.9);

        h1rg[icent]->SetMarkerColor(1);
        h1rg[icent]->SetLineColor(1);
        h1rg[icent]->SetTitle(Form("R_{#gamma} %i-%i%%",centlim[icent][0],centlim[icent][1]));
        h1rg[icent]->Draw();
        can->Update();
        tleg->AddEntry(h1rg[icent]);

        int icolor = 3;
        int imarker = 25;

        for (int ippg = ppglim[icent][0]; ippg <= ppglim[icent][1]; ippg++)
        {
            tgr[ippg]->SetMarkerColor(icolor);
            tgr[ippg]->SetLineColor(icolor);
            tgr[ippg]->SetMarkerStyle(imarker);

            tgr[ippg]->Draw("P");
            can->Update();
            tleg->AddEntry(tgr[ippg]);

            icolor++;
            imarker++;
        }

        tleg->Draw();
        can->Update();
        can->SaveAs(Form("Rg_PPG_%i.png",icent));
    }

}
