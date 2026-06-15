// apply the partner photon ksi correction
// save new mistag and tageff functions

#define ncent 90

int GetFunctions(const char* loc, const char* name, TF1* f1[ncent]);

void ksiCorrectPartner() {

    //------------------------------
    // Get Functions
    //------------------------------

    TF1* ftageff[ncent];
    TF1* fmistag[ncent];

    const char* loc_mis = "mistag.root";
    const char* loc_tag = "tageff.root";

    if (GetFunctions(loc_mis,"fmistag_",fmistag) < 0)
    {
        cerr<<"Could not retrieve mistag functions. Exiting."<<endl;
        return;
    }

    if (GetFunctions(loc_tag,"ftageff_",ftageff) < 0)
    {
        cerr<<"Could not retrieve tageff functions. Exiting."<<endl;
        return;
    }

    //------------------------------
    // Get Ksic for Partner Photons
    //------------------------------

    int ipta_ksic = 0; // use 0.5-5, forget 1.0-5.0
    int iptt_ksic = 2; // trigger 5-15
    double ksic[ncent] = {};

    const char* loc_ksic = "xicorrect_inc.txt";
    ifstream ifs_ksic(loc_ksic,ifstream::in);
    if (!ifs_ksic)
    {
        cerr<<"Could not open "<<loc_ksic<<endl;
        return;
    }
    string line;
    while (getline(ifs_ksic,line))
    {
        if (!isdigit(line.at(0))) continue;

        istringstream iss_line(line.data());
        int iptt, ipta, icent;
        double tmpksic;
        iss_line >> iptt >> ipta >> icent >> tmpksic;

        if (ipta != ipta_ksic) continue;
        if (iptt != iptt_ksic) continue;

        ksic[icent] = tmpksic;
        cout<<icent<<"  "<<tmpksic<<endl;
    }

    //----------------------------------------
    // Scale mistag by ksic
    // Update initial 3 parameters of tageff
    //----------------------------------------

    TCanvas* can = new TCanvas("can","",1600,1000);

    TFile* out_tageff = new TFile("tageff_ksic.root","RECREATE");
    TFile* out_mistag = new TFile("mistag_ksic.root","RECREATE");

    TF1* ftageff_ksic[ncent];
    TF1* fmistag_ksic[ncent];

    // the heavyside business is not supported by root 5?
    // if you wish to use ksic online (in gamanatree) these functions must be remade in similar fashion but creating own function
    // user defined functions (meaning using a method) cannot be written and read

    for (int icent = 0; icent < ncent; icent++)
    {
        out_mistag->cd();
        fmistag_ksic[icent] = new TF1(Form("fmistag_ksic_%i",icent),"[0]*(([1]/x + [2]/x^2 + [3]/x^3) > 0.0 ? ([1]/x + [2]/x^2 + [3]/x^3) : 0.0)",5.,15.);
        fmistag_ksic[icent]->SetParameter(0,ksic[icent]);
        for (int ipar = 0; ipar < 3; ipar++)
        {
            fmistag_ksic[icent]->SetParameter(ipar+1,fmistag[icent]->GetParameter(ipar));
        }
        fmistag_ksic[icent]->Write();

        TH1D* h1 = new TH1D(Form("mistag_ksic_cent%i",icent),Form("mistag_ksic_cent%i",icent),200,5.,15.);
        h1->Add(fmistag_ksic[icent]);
        h1->Sumw2();
        h1->Draw();
        can->SaveAs(Form("mistag_ksic_cent%i.png",icent));

        TH1D* h2 = new TH1D(Form("mistag_cent%i",icent),Form("mistag_cent%i",icent),200,5.,15.);
        h2->Add(fmistag[icent]);
        h2->Sumw2();
        h2->Draw();
        can->SaveAs(Form("mistag_cent%i.png",icent));

        TRatioPlot* rp = new TRatioPlot(h1,h2,"divsym");
        rp->Draw();
        ostringstream namemistag;
        namemistag<<"ratio_mistag_cent"<<icent<<".png";
        can->SaveAs(namemistag.str().c_str());
        delete rp;
        delete h1;
        delete h2;

        out_tageff->cd();
        ftageff_ksic[icent] = new TF1(Form("ftageff_ksic_%i",icent),"[0]*(([1]/x + [2]/x^2 + [3]/x^3) > 0.0 ? ([1]/x + [2]/x^2 + [3]/x^3) : 0.0) + [4] + [5]*x + [6]*x^2",5.,15.);
        ftageff_ksic[icent]->SetParameter(0,ksic[icent]);
        for (int ipar = 0; ipar < 6; ipar++)
        {
            ftageff_ksic[icent]->SetParameter(ipar+1,ftageff[icent]->GetParameter(ipar));
        }
        ftageff_ksic[icent]->Write();

        h1 = new TH1D(Form("tageff_ksic_cent%i",icent),Form("tageff_ksic_cent%i",icent),200,5.,15.);
        h1->Add(ftageff_ksic[icent]);
        h1->Sumw2();
        h1->Draw();
        can->SaveAs(Form("tageff_ksic_cent%i.png",icent));

        h2 = new TH1D(Form("tageff_cent%i",icent),Form("tageff_cent%i",icent),200,5.,15.);
        h2->Add(ftageff[icent]);
        h2->Sumw2();
        h2->Draw();
        can->SaveAs(Form("tageff_cent%i.png",icent));

        rp = new TRatioPlot(h1,h2,"divsym");
        rp->Draw();
        ostringstream nametageff;
        nametageff<<"ratio_tageff_cent"<<icent<<".png";
        can->SaveAs(nametageff.str().c_str());
        delete rp;
        delete h1;
        delete h2;
    }

    out_tageff->Close();
    out_mistag->Close();
}


int GetFunctions(const char* loc, const char* name, TF1* f1[ncent])
{
    cout<<"Getting functions "<<name<<"? from "<<loc<<endl;

    memset(f1,0,sizeof(TF1*)*ncent);

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


    for (int icent = 0; icent < ncent; icent++)
    {
        ostringstream oss;
        oss<<name<<icent;

        TF1* f1tmp = (TF1*)file->Get(oss.str().c_str());
        if (f1tmp == NULL)
        {
            cerr<<"ABORT!!"<<endl;
            return -1;
        }

        f1[icent] = (TF1*) f1tmp->Clone();
    }

    file->Close();
    delete file;

    return 0;
}
