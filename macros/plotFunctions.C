#define NCENT 90

int getFunctions(const char* loc, const char* name, TF1* f1[NCENT]);

void plotFunctions()
{
    TCanvas* can = new TCanvas("can","",1000,1200);

    // get mistag and tageff functions
    TF1* ftageff[NCENT];
    TF1* fmistag[NCENT];
    TF1* ftageff_ksic[NCENT];
    TF1* fmistag_ksic[NCENT];

    const char* loc_mis = "mistag.root";
    const char* loc_tag = "tageff.root";
    const char* loc_mis_ksic = "mistag_ksic.root";
    const char* loc_tag_ksic = "tageff_ksic.root";

    if (getFunctions(loc_mis_ksic,"fmistag_ksic_",fmistag_ksic) < 0)
    {
        cerr<<"Could not retrieve mistag functions from "<<loc_mis_ksic<<endl;
        cerr<<"Exiting..."<<endl;
        return;
    }

    if (getFunctions(loc_tag_ksic,"ftageff_ksic_",ftageff_ksic) < 0)
    {
        cerr<<"Could not retrieve tageff functions from "<<loc_tag_ksic<<endl;
        cerr<<"Exiting..."<<endl;
        return;
    }

    if (getFunctions(loc_mis,"fmistag_",fmistag) < 0)
    {
        cerr<<"Could not retrieve mistag functions from "<<loc_mis<<endl;
        cerr<<"Exiting..."<<endl;
        return;
    }

    if (getFunctions(loc_tag,"ftageff_",ftageff) < 0)
    {
        cerr<<"Could not retrieve tageff functions from "<<loc_tag<<endl;
        cerr<<"Exiting..."<<endl;
        return;
    }

    TH1F* h1mistag = new TH1F("mistag","mistag probability ptt = 6;centrality;mistag",NCENT,0,NCENT);
    TH1F* h1tageff = new TH1F("tageff","tageff probability ptt = 6;centrality;tageff",NCENT,0,NCENT);
    TH1F* h1mistag_ksic = new TH1F("mistag_ksic","mistag probability (ksic applied) ptt = 6;centrality;mistag",NCENT,0,NCENT);
    TH1F* h1tageff_ksic = new TH1F("tageff_ksic","tageff probability (ksic applied) ptt = 6;centrality;tageff",NCENT,0,NCENT);

    for (int ifun = 0; ifun < NCENT; ifun++)
    {
        h1mistag->SetBinContent(ifun+1,fmistag[ifun]->Eval(6.0));

        h1tageff->SetBinContent(ifun+1,ftageff[ifun]->Eval(6.0));

        h1mistag_ksic->SetBinContent(ifun+1,fmistag_ksic[ifun]->Eval(6.0));

        h1tageff_ksic->SetBinContent(ifun+1,ftageff_ksic[ifun]->Eval(6.0));
    }

    h1mistag->Draw();
    can->SaveAs("mistag.png");

    h1tageff->Draw();
    can->SaveAs("tageff.png");

    h1mistag_ksic->Draw();
    can->SaveAs("mistag_ksic.png");

    h1tageff_ksic->Draw();
    can->SaveAs("tageff_ksic.png");

    cout<<"Done!"<<endl;
}

int getFunctions(const char* loc, const char* name, TF1* f1[NCENT])
{
    cout<<"Getting functions "<<name<<"? from "<<loc<<endl;

    memset(f1,0,sizeof(TF1*)*NCENT);

    TFile* file = new TFile(loc,"READ");

    if (!file)
    {
        cerr<<"Could not open input "<<loc<<endl;
        return -1;
    }
    if (file->IsZombie())
    {
        cerr<<"Input "<<loc<<" is zombie."<<endl;
        return -1;
    }

    cout<<"Input "<<loc<<" opened successfully."<<endl;


    for (int icent = 0; icent < NCENT; icent++)
    {
        ostringstream oss;
        oss<<name<<icent;

        TF1* f1tmp = (TF1*)file->Get(oss.str().c_str());
        if (f1tmp == NULL)
        {
            cerr<<"ABORT!!"<<endl;
            return -1;
        }

        oss<<"_clone";
        f1[icent] = (TF1*) f1tmp->Clone(oss.str().c_str());
    }

    file->Close();
    delete file;

    return 0;
}
