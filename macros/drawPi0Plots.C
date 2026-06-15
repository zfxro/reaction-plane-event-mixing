ofstream ofs;
bool watcher= false;

Double_t background(Double_t *x, Double_t *par) {
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t gaussian(Double_t *x,Double_t *par)
{
    Double_t arg = 0;
    if (par[2] != 0) arg = (x[0] - par[1])/par[2];
    Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}

Double_t fitFunction(Double_t *x, Double_t *par) {
    return background(x,par) + gaussian(x,&par[3]);
}

//loc commented out.  no canvas division atm
void Fitting(TH1F *histo, TCanvas* can, int loc = 1)
{  
    TVirtualPad* tvp = can->cd(loc);

    double amp = histo->GetMaximum();
    double mu = histo->GetBinCenter(histo->GetMaximumBin());
    double sig = .01;


    cout<<"amp = "<<amp<<"; mu = "<<mu<<"; sig = "<<sig<<endl;


    int ibin = histo->FindBin(mu-3.*sig);
    int maxbin = histo->FindBin(mu+3.*sig);

    cout<<"ibin = "<<ibin<<"; maxbin = "<<maxbin<<endl;

    TH1F* hclone = (TH1F*)histo->Clone("clone");

    if (watcher)
    {
        hclone->Draw();
        cout<<"Histogram clone has been drawn.  Press return to continue."<<endl;
        tvp->Update();
        gSystem->ProcessEvents();
        cin.ignore();
    }

    //set bins in pi0 peak to 0. fitting just the background now.
    for (;ibin<=maxbin;ibin++)
    {
        hclone->SetBinContent(ibin,0);
        hclone->SetBinError(ibin,0);
    }

    if (watcher)
    {
        hclone->Draw();
        cout<<"pi0 peak bins set to 0.  Press return to continue."<<endl;
        tvp->Update();
        gSystem->ProcessEvents();
        cin.ignore();
    }

    double fitXi = .05;
    double fitXf = .3;
    TF1* iniPoly = new TF1("iniPoly",background,fitXi,fitXf,3);
    iniPoly->SetParameters(amp/100.,10000.,-10000.);
    hclone->Fit("iniPoly","R");

    if (watcher)
    {
        cout<<"Background has been fit with 2nd order polynomial.  Press return to continue."<<endl;
        tvp->Update();
        gSystem->ProcessEvents();
        cin.ignore();
    }

    double polyPars[3];
    iniPoly->GetParameters(polyPars);
    amp -= iniPoly->Eval(mu);
    TF1 *fitFcn = new TF1("fitFcn",fitFunction,fitXi,fitXf,6);
    fitFcn->SetParNames("p0","p1","p2","Constant","Mean","Sigma");
    fitFcn->SetLineColor(kMagenta);
    fitFcn->SetParameters(polyPars[0],polyPars[1],polyPars[2],amp,mu,sig);  


    histo->Fit("fitFcn","R","ep");

    if (watcher)
    {
        cout<<"Background and peak have been fit.  Press return to continue."<<endl;
        tvp->Update();
        gSystem->ProcessEvents();
        cin.ignore();
    }


    // return;
  
    // improve the picture:
    TF1 *backFcn = new TF1("backFcn",background,0,3,3);
    backFcn->SetLineColor(kRed);
  
    TF1 *signalFcn = new TF1("signalFcn",gaussian,0,3,3);
    signalFcn->SetLineColor(kBlue);
    Double_t par[6];
    
    // writes the fit results into the par array
    fitFcn->GetParameters(par);
      
    backFcn->SetParameters(par);
    backFcn->Draw("same");
   
    signalFcn->SetParameters(&par[3]);
    signalFcn->Draw("same"); 
  
    // draw the legend
    TLegend *legend=new TLegend(0.57,0.75,.95,.95);
    legend->SetTextFont(72);
    legend->SetTextSize(0.04);
    legend->SetFillColor(0);
    legend->AddEntry(histo,"Data","lpe");
    legend->AddEntry(backFcn,"Background fit","l");
    legend->AddEntry(signalFcn,"Signal fit","l");
    legend->AddEntry(fitFcn,"Global Fit","l");
    legend->Draw();
  
    double integral = par[3]*sqrt(2*TMath::Pi())*par[5];

    ofs<<par[4]<<endl;
    ofs<<par[5]<<endl;
  
    TPaveStats *ps = new TPaveStats();
    // {
    // }
    // {
    // }

    tvp->Update(); 
    ps = (TPaveStats*)tvp->GetPrimitive("stats");
    ps->SetName("mystats");
  
    TList *list = ps->GetListOfLines();
  
  
    TLatex *hitT = new TLatex(0,0,Form("N#pi^{0} = %.0f",integral));
    hitT->SetTextFont(72);
    hitT->SetTextSize(0.04);
    list->Add(hitT);
    ps->Draw();
  
    histo->SetStats(0);
    tvp->Modified();

    if (watcher)
    {
        cout<<"Fitting method has ended.  Press return to continue."<<endl;
        tvp->Update();
        gSystem->ProcessEvents();
        cin.ignore();
    }  


    // return;

}


void drawPi0Plots(const char* input = "gamAna_13466_histos.root")
{
    gStyle->SetOptStat(0);
    gStyle->SetStatFont(72);
    gStyle->SetStatFontSize(0.02);
    gStyle->SetTitleFontSize(0.04);
    gStyle->SetOptFit(0011);
    gStyle->SetStatW(.25);
    gStyle->SetStatX(.95);
    gStyle->SetStatY(.75);      
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.075);
    gStyle->SetPadTopMargin(0.05);

    TGaxis::SetMaxDigits(3);

    static const int nsec = 8;

    // ofstream ofs;
    ofs.open("pi0fitpar.txt");


    TFile in(input);

    TCanvas *can = new TCanvas("can","",1000,600);
    can->Divide(4,2,0.001,0.001);


    TH1F* h1invm_allsec = NULL;

    for (int isec = 0; isec < nsec+6; isec++)
    {

        printf("\n\nSector %i\n",isec);

        TH2F* h2invm = (TH2F*) gROOT->FindObject(Form("mpt_%i",isec));
        cout<<"WARNING! Reset projection to full pt range if using for mistag analysis!"<<endl;
        int bin0 = h2invm->GetYaxis()->FindBin(5.0001);
        int bin1 = h2invm->GetYaxis()->FindBin(14.9999);
        TH1F* h1invm = (TH1F*) h2invm->ProjectionX(Form("mpt_%i_px",isec),bin0,bin1);
        float binw = h1invm->GetBinWidth(1);
        h1invm->Scale(1./binw);
        if (!h1invm_allsec)
        {
            h1invm_allsec = (TH1F*) h1invm->Clone("invm_allsec");
            h1invm_allsec->SetTitle("");
            h1invm_allsec->GetYaxis()->SetTitle("N_{pair}");
        }
        else
            h1invm_allsec->Add(h1invm);

        //sectors above nsec (photons from neighboring sectors) must be included in allsec
        if (isec >= nsec) continue;
        
        string type = (isec > 1) ? "PbSc" : "PbGl";
        h1invm->SetTitle(Form("Sector %i %s",isec,type.c_str()));
        h1invm->GetYaxis()->SetTitle("N_{pair}");
        h1invm->GetXaxis()->SetRangeUser(.05,.3);
        Fitting(h1invm,can,isec+1);
        // Fitting(h1invm,can);

        //couldnt get pad copying to work. aborting method.

        // can_arm[iarm]->cd(iloc);

    }

    can->SaveAs("pi0fit.pdf");

    // {
    // }

    TCanvas *can_allsec = new TCanvas("can_allsec","",1000,600);
    h1invm_allsec->GetXaxis()->SetRangeUser(.05,.3);
    Fitting(h1invm_allsec,can_allsec);
    can_allsec->SaveAs("pi0fit_allsec.pdf");

    TCanvas *can_arm[2];
    for (int iarm = 0; iarm < 2; iarm++)
    {
        can_arm[iarm] = new TCanvas(Form("can_arm%i",iarm),"",3000,3500);
        can_arm[iarm]->Divide(2,2,0.001,0.001);

        // {
        // }
    };

    //Refitting all because im a bad coder but want a different drawing
    for (int isec = 0; isec < nsec; isec++)
    {
        printf("\n\nSector %i\n",isec);

        TH2F* h2invm = (TH2F*) gROOT->FindObject(Form("mpt_%i",isec));
        cout<<"WARNING! Reset projection to full pt range if using for mistag analysis!"<<endl;
        int bin0 = h2invm->GetYaxis()->FindBin(5.0001);
        int bin1 = h2invm->GetYaxis()->FindBin(14.9999);
        TH1F* h1invm = (TH1F*) h2invm->ProjectionX(Form("mpt_%i_px_extra",isec),bin0,bin1);
        float binw = h1invm->GetBinWidth(1);
        h1invm->Scale(1./binw);
        string type = (isec > 1) ? "PbSc" : "PbGl";
        h1invm->SetTitle(Form("Sector %i %s",isec,type.c_str()));
        h1invm->GetYaxis()->SetTitle("N_{pair}");
        h1invm->GetXaxis()->SetRangeUser(.05,.3);

        int iarm = (isec < 4) ? 0 : 1;
        int iloc = (isec < 4) ? isec+1 : isec-4+1;
        Fitting(h1invm,can_arm[iarm],iloc);
    }

    for (int iarm = 0; iarm < 2; iarm++)
    {
        can_arm[iarm]->SaveAs(Form("pi0fit_arm%i.pdf",iarm));
    }

    cout<<"\nAll Done!\n";
}
