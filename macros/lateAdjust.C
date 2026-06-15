//purpose:
//
//  get weighted averages of mistag and tageff (also allow for away side ratio correction)
//  then apply the adjustment (get adj1_dir stuff) (aka late adjustment as opposed to during gamAnaTree)
//
//input:
//
//  - integawaysideyield.root
//      + iasy (as function of ptt and xi)
//         ipttyield_adj0_[dir,dec]_cent[0-2]_xi[0-2]
//         ixiyield_adj0_[dir,dec]_ptt[0-3]_cent[0-2]
//         ** also the sys errors associated with those... **
//         ipttyield_syserr_adj0_[dir,dec]_cent[0-2]_xi[0-2]
//         ixiyield_syserr_adj0_[dir,dec]_ptt[0-3]_cent[0-2]
//
//    Note: y values are same for syserr plots. only errors differ
//
//  - dphiplotsdrawn.root
//      + v2 (as function of ptt and xi)
//         ipttv2jfromfit_adj0_[dir,dec]_cent[0-2]_xi[0-2]
//         ixiv2jfromfit_adj0_[dir,dec]_ptt[0-3]_cent[0-2]
//         ** also the sys errors associated with those... **
//         ipttv2jfromfit_syserr_adj0_[dir,dec]_cent[0-2]_xi[0-2]
//         ixiv2jfromfit_syserr_adj0_[dir,dec]_ptt[0-3]_cent[0-2]
//
//  - gamtree_psimix.root
//      + centptt histos for adj0 dir and dec (need just inclusive)
//         cent_ptt_[dir,dec]_dpsi[0,1] (will have to add dpsi0 and 1
//         see acceptCorrect for code to retrieve and merge these histos
//
//  - mistag.root
//      + fmistag_[0-59]  Note: cent convention. 0=1%, 59=60%
//
//  - tageff.root
//      + ftageff_[0-59]
//
//    Note: useful to see how functions are grabbed in gamAnaTree.cc
//
//algo:
//
//  ncent = 2; nptt = 4;
//  misavg[ncent][nptt], tagavg[ncent][nptt]
//  loop cent and ptt bins of cent_ptt inc histo (must make inclusive first aka add adj0_dir and adj0_dec)
//  get weighted avg of mis and tag: mis += ngam*fmistag_cent->Eval(ptt); tag += ngam*ftageff_cent->Eval(ptt); totgam += ngam;
//  divide by total photons in that cent/ptt bin: misavg[icent][iptt] = mis/totgam; tagavg[icent][iptt] = tag/totgam;
//
//  then work on making the new graphs using our cool formula
//  allow for asr (or r) factor which will be binned in ptt and xi
//  algo to be continued...
//
//output:
//
//      ipttyield_lateadj_dir_cent[0-2]_xi[0-2]
//      ixiyield_lateadj_dir_ptt[0-3]_cent[0-2]
//      ipttyield_syserr_lateadj_dir_cent[0-2]_xi[0-2]
//      ixiyield_syserr_lateadj_dir_ptt[0-3]_cent[0-2]
//
//      ipttv2jfromfit_lateadj_dir_cent[0-2]_xi[0-2]
//      ixiv2jfromfit_lateadj_dir_ptt[0-3]_cent[0-2]
//      ipttv2jfromfit_syserr_lateadj_dir_cent[0-2]_xi[0-2]
//      ixiv2jfromfit_syserr_lateadj_dir_ptt[0-3]_cent[0-2]

const float PI = TMath::Pi();
const int nspec = 2;
const char* cspec[nspec] = {"dir","dec"};
const int ncent = 3;
const int nxi = 3;
const int nptt = 4;
enum TRIG {DIR=0,DEC=1,INC=2};
enum RPO {IN=0,OUT=1,ALL=2};

//yields
TGraphErrors* ipttyield_adj0[nspec][ncent][nxi];//[dir,dec][cent0-2][xi0-2]
TGraphErrors* ixiyield_adj0[nspec][nptt][ncent];//[dir,dec][ptt0-3][cent0-2]
TGraphErrors* ipttyield_syserr_adj0[nspec][ncent][nxi];
TGraphErrors* ixiyield_syserr_adj0[nspec][nptt][ncent];

//yields output
TGraphErrors* ipttyield_lateadj_dir[ncent][nxi];
TGraphErrors* ixiyield_lateadj_dir[nptt][ncent];
TGraphErrors* ipttyield_syserr_lateadj_dir[ncent][nxi];
TGraphErrors* ixiyield_syserr_lateadj_dir[nptt][ncent];

TGraphErrors* ipttcentrat_lateadj_dir[nxi];
TGraphErrors* ixicentrat_lateadj_dir[nptt];
TGraphErrors* ipttcentrat_syserr_lateadj_dir[nxi];
TGraphErrors* ixicentrat_syserr_lateadj_dir[nptt];

//v2
TGraphErrors* ipttv2jfromfit_adj0[nspec][ncent][nxi];//[dir,dec][cent0-2][xi0-2]
TGraphErrors* ixiv2jfromfit_adj0[nspec][nptt][ncent];//[dir,dec][ptt0-3][cent0-2]
TGraphErrors* ipttv2jfromfit_syserr_adj0[nspec][ncent][nxi];
TGraphErrors* ixiv2jfromfit_syserr_adj0[nspec][nptt][ncent];

//v2j output
TGraphErrors* ipttv2jfromfit_lateadj_dir[ncent][nxi];
TGraphErrors* ixiv2jfromfit_lateadj_dir[nptt][ncent];
TGraphErrors* ipttv2jfromfit_syserr_lateadj_dir[ncent][nxi];
TGraphErrors* ixiv2jfromfit_syserr_lateadj_dir[nptt][ncent];

//inclusive photon cent ptt
TH2F* h2centptt_inc;
TH1F* h1ptt_inc;
//also for unadjusted dir and dec
TH1F* h1ptt_dir[ncent];
TH1F* h1ptt_dec[ncent];

//tageff and mistag functions
const int ncent_fun = 60;//max is 94 but we only need 0-59
TF1* ftageff[ncent_fun];
TF1* fmistag[ncent_fun];
//averages of those in cent and ptt bins
double mis[ncent][nptt] {};
double tag[ncent][nptt] {};
double ninc[ncent][nptt] {};
//also need counters for ndir and ndec unadjusted
double ndir_adj0[ncent][nptt] {};
double ndec_adj0[ncent][nptt] {};

//error widths for sys err bars. copying from integAwaySide.C
double sex_xi = .035;
double xioff[nptt] = {-.12,-.04,.04,.12};
double sex_ptt = .04;
double pttoff[nxi] = {-.09,.0,.09};

TFile* OpenFile(const char* name);
void GetGraphs(TFile* input_yield, TFile* input_v2);
int GetFunctions(const char* loc, const char* name, TF1* f1[ncent_fun]);
void GetWeights(const int icent,const int iptt,double& rrwt,double& rcwt,double& rrwt_asr,double& rcwt_asr,const double r = 1);
int GetPttBin(double ptt);
void GetOldSysErr(TFile* input_yield, TFile* input_v2);

struct errorinfo_t {
    Int_t iadj, itype, iptt, ixi, icent, idpsi;
    Double_t ntrg;
    Double_t flate[2] = {};//flat term error
    Double_t cose[2] = {};//cosine term error
};
errorinfo_t errinfo[2][3][nptt][ncent][3][nxi] = {};
void FillSysError(TFile*);

void MakeRatioGraphs();
    
void lateAdjust(
    const char* inname_yield = "integawaysideyield.root",
    const char* inname_v2 = "dphiplotsdrawn.root",
    const char* inname_ngam = "gamtree_psimix.root",
    const char* inname_mis = "mistag.root",
    const char* inname_tag = "tageff.root",
    const char* inname_bgsub = "weighted_bg_subtracted.root",
    const char* outname = "lateadjust.root"
    )
{
    cout<<"WARNING! Sys errors calculated incorrectly. There are correlations. Cannot just add in quadrature. See integAwaySide.C for proper method including errinfo[]"<<endl;

    TFile* input_bgsub = OpenFile(inname_bgsub);
    FillSysError(input_bgsub);

    TCanvas* can = new TCanvas("can","",4000,3000);

    for (int icent = 0; icent < ncent; icent++)
        for (int iptt = 0; iptt < nptt; iptt++)
            if (mis[icent][iptt] != 0. || tag[icent][iptt] != 0. || ninc[icent][iptt] != 0. || ndir_adj0[icent][iptt] != 0. || ndec_adj0[icent][iptt] != 0.)
                return;

    //---------- Get Away Side Yield and V2 Graphs ----------
    TFile* input_yield = OpenFile(inname_yield);
    TFile* input_v2 = OpenFile(inname_v2);
    GetGraphs(input_yield,input_v2);

    //---------- Get Inclusive Photons Ptt vs Cent 2d Histo ----------
    TFile* input_ngam = OpenFile(inname_ngam);
    {
        TH2F* h2centptt_tmp[nspec][3];
        for (int ispec = 0; ispec < nspec; ispec++)
        {
            for (int idpsi = 0; idpsi < 2; idpsi++)
            {
                input_ngam->GetObject(Form("cent_ptt_%s_dpsi%d",cspec[ispec],idpsi),h2centptt_tmp[ispec][idpsi]);
            }
            h2centptt_tmp[ispec][ALL] = (TH2F*) h2centptt_tmp[ispec][IN]->Clone(Form("cent_ptt_%s_dpsi%d",cspec[ispec],ALL));
            h2centptt_tmp[ispec][ALL]->Add(h2centptt_tmp[ispec][OUT]);
        }

        //----- projections to obtain ptt dists -----
        h1ptt_dir[0] = (TH1F*) h2centptt_tmp[DIR][ALL]->ProjectionY("ptt_dir_cent0",1,20);
        h1ptt_dir[0]->Draw();
        can->SaveAs("ptt_dir_cent0.pdf");
        
        h1ptt_dir[1] = (TH1F*) h2centptt_tmp[DIR][ALL]->ProjectionY("ptt_dir_cent1",21,60);
        h1ptt_dir[1]->Draw();
        can->SaveAs("ptt_dir_cent1.pdf");

        h1ptt_dir[2] = (TH1F*) h2centptt_tmp[DIR][ALL]->ProjectionY("ptt_dir_cent2",1,60);
        h1ptt_dir[2]->Draw();
        can->SaveAs("ptt_dir_cent2.pdf");

        h1ptt_dec[0] = (TH1F*) h2centptt_tmp[DEC][ALL]->ProjectionY("ptt_dec_cent0",1,20);
        h1ptt_dec[0]->Draw();
        can->SaveAs("ptt_dec_cent0.pdf");
        
        h1ptt_dec[1] = (TH1F*) h2centptt_tmp[DEC][ALL]->ProjectionY("ptt_dec_cent1",21,60);
        h1ptt_dec[1]->Draw();
        can->SaveAs("ptt_dec_cent1.pdf");

        h1ptt_dec[2] = (TH1F*) h2centptt_tmp[DEC][ALL]->ProjectionY("ptt_dec_cent2",1,60);
        h1ptt_dec[2]->Draw();
        can->SaveAs("ptt_dec_cent2.pdf");

        //----- fill dir and dec (unadjusted) counters -----
        for (int bptt = h1ptt_dir[0]->GetNbinsX(); bptt >= 1; bptt--)
        {
            int iptt = GetPttBin(h1ptt_dir[0]->GetBinCenter(bptt));
            for (int icent = ncent - 1; icent >= 0; icent--)
            {
                ndir_adj0[icent][iptt] += h1ptt_dir[icent]->GetBinContent(bptt);
                ndec_adj0[icent][iptt] += h1ptt_dec[icent]->GetBinContent(bptt);
            }
        }

        //----- print out counted results -----
        for (int iptt = 0; iptt < nptt; iptt++)
        {
            for (int icent = 0; icent < ncent; icent++)
            {
                cout<<"----- cent"<<icent<<" ptt"<<iptt<<" -----"<<endl;
                cout<<"ndir "<<ndir_adj0[icent][iptt]<<endl;
                cout<<"ndec "<<ndec_adj0[icent][iptt]<<endl;
            }
        }
            

        h2centptt_tmp[DIR][ALL]->Draw("colz");
        can->SaveAs("centptt_dir.pdf");
        
        h2centptt_tmp[DEC][ALL]->Draw("colz");
        can->SaveAs("centptt_dec.pdf");

        h2centptt_inc = (TH2F*) h2centptt_tmp[DIR][ALL]->Clone("centptt_inc");
        h2centptt_inc->Add(h2centptt_tmp[DEC][ALL]);
        h2centptt_inc->Draw("colz");
        can->SaveAs("centptt_inc.pdf");

        h1ptt_inc = (TH1F*) h2centptt_inc->ProjectionY("ptt_inc",1,60);
        h1ptt_inc->Draw();
        can->SaveAs("ptt_inc.pdf");
        
    }

    //---------- Get Mistag Functions ----------
    if (GetFunctions(inname_mis,"fmistag_",fmistag) < 0)
    {
        cerr<<"Could not retrieve mistag functions. Exiting."<<endl;
        return;
    }

    //---------- Get Tageff Functions ----------
    if (GetFunctions(inname_tag,"ftageff_",ftageff) < 0)
    {
        cerr<<"Could not retrieve tageff functions. Exiting."<<endl;
        return;
    }

    //---------- Open output for writing ----------
    TFile* output = new TFile(outname,"RECREATE");

        
    for (int bcent = 1; bcent <= h2centptt_inc->GetNbinsX(); bcent++)
    {
        int ccent = h2centptt_inc->GetXaxis()->GetBinCenter(bcent);
        int icent = -999;
        if (ccent <= 20)
            icent = 0;
        else if (ccent <= 60)
            icent = 1;
        else
            continue;

        
        for (int bptt = 1; bptt <= h2centptt_inc->GetNbinsY(); bptt++)
        {
            double cptt = h2centptt_inc->GetYaxis()->GetBinCenter(bptt);
            int iptt = GetPttBin(cptt);
            if (iptt < 0) continue;
        
            int ngam = h2centptt_inc->GetBinContent(bcent,bptt);

            double m = fmistag[ccent-1]->Eval(cptt);
            if (m < 0)
            {
                m = 0;
            }
            double t = ftageff[ccent-1]->Eval(cptt);
            if (t < 0)
            {
                t = 0;
            }
            
            mis[icent][iptt] += ngam*m;
            tag[icent][iptt] += ngam*t;
            ninc[icent][iptt] += ngam;

            mis[ncent-1][iptt] += ngam*m;
            tag[ncent-1][iptt] += ngam*t;
            ninc[ncent-1][iptt] += ngam;
            
        }
    }

    //normalize
    cout<<"Calculated weighted averages of mistag and tageff."<<endl;
    for (int icent = 0; icent < ncent; icent++)
    {
        for (int iptt = 0; iptt < nptt; iptt++)
        {
            mis[icent][iptt] /= ninc[icent][iptt];
            tag[icent][iptt] /= ninc[icent][iptt];
            printf("icent%d iptt%d averages: mistag = %f, tageff = %f\n",icent,iptt,mis[icent][iptt],tag[icent][iptt]);
        }
    }

    //create adjusted yield and v2
    for (int icent = 0; icent < ncent; icent++)
    {
        for (int iptt = 0; iptt < nptt; iptt++)
        {
            ixiyield_lateadj_dir[iptt][icent] = new TGraphErrors(nxi);
            ixiyield_lateadj_dir[iptt][icent]->SetName(Form("ixiyield_adj1_dir_ptt%d_cent%d_dpsi2",iptt,icent));
            ixiyield_lateadj_dir[iptt][icent]->SetTitle(";#xi;1/N_{#gamma} dN/d#xi");
            
            ixiyield_syserr_lateadj_dir[iptt][icent] = new TGraphErrors(nxi);
            ixiyield_syserr_lateadj_dir[iptt][icent]->SetName(Form("ixiyield_syserr_adj1_dir_ptt%d_cent%d_dpsi2",iptt,icent));
            ixiyield_syserr_lateadj_dir[iptt][icent]->SetTitle(";#xi;1/N_{#gamma} dN/d#xi");
            
            ixiv2jfromfit_lateadj_dir[iptt][icent] = new TGraphErrors(nxi);
            ixiv2jfromfit_lateadj_dir[iptt][icent]->SetName(Form("ixiv2jfromfit_adj1_dir_ptt%d_cent%d",iptt,icent));
            ixiv2jfromfit_lateadj_dir[iptt][icent]->SetTitle(";#xi;1/N_{#gamma} dN/d#xi");
            
            ixiv2jfromfit_syserr_lateadj_dir[iptt][icent] = new TGraphErrors(nxi);
            ixiv2jfromfit_syserr_lateadj_dir[iptt][icent]->SetName(Form("ixiv2jfromfit_syserr_adj1_dir_ptt%d_cent%d",iptt,icent));
            ixiv2jfromfit_syserr_lateadj_dir[iptt][icent]->SetTitle(";#xi;1/N_{#gamma} dN/d#xi");

            double rrwt, rcwt, rrwt_asr, rcwt_asr;
            GetWeights(icent,iptt,rrwt,rcwt,rrwt_asr,rcwt_asr,1.2);

            rrwt *= ndir_adj0[icent][iptt];
            rrwt_asr *= ndir_adj0[icent][iptt];

            rcwt *= ndec_adj0[icent][iptt];
            rcwt_asr *= ndec_adj0[icent][iptt];

            for (int ixi = 0; ixi < nxi; ixi++)
            {
                double xi = ixiyield_adj0[DIR][iptt][icent]->GetX()[ixi];

                //----- yields -----
                double ydir0 = ixiyield_adj0[DIR][iptt][icent]->GetY()[ixi];
                double eydir0 = ixiyield_adj0[DIR][iptt][icent]->GetEY()[ixi];
                double seydir0 = ixiyield_syserr_adj0[DIR][iptt][icent]->GetEY()[ixi];

                double ydec0 = ixiyield_adj0[DEC][iptt][icent]->GetY()[ixi];
                double eydec0 = ixiyield_adj0[DEC][iptt][icent]->GetEY()[ixi];
                double seydec0 = ixiyield_syserr_adj0[DEC][iptt][icent]->GetEY()[ixi];
                
                //----- calc adjusted dir yield -----
                double ydir1 = rrwt_asr*ydir0 + rcwt_asr*ydec0;
                double eydir1 = sqrt(pow(rrwt_asr*eydir0,2.) + pow(rcwt_asr*eydec0,2.));
                double seydir1 = sqrt(pow(rrwt_asr*seydir0,2.) + pow(rcwt_asr*seydec0,2.));

                //----- v2jfromfit -----
                double v2dir0 = ixiv2jfromfit_adj0[DIR][iptt][icent]->GetY()[ixi];
                double ev2dir0 = ixiv2jfromfit_adj0[DIR][iptt][icent]->GetEY()[ixi];
                double sev2dir0 = ixiv2jfromfit_syserr_adj0[DIR][iptt][icent]->GetEY()[ixi];

                double v2dec0 = ixiv2jfromfit_adj0[DEC][iptt][icent]->GetY()[ixi];
                double ev2dec0 = ixiv2jfromfit_adj0[DEC][iptt][icent]->GetEY()[ixi];
                double sev2dec0 = ixiv2jfromfit_syserr_adj0[DEC][iptt][icent]->GetEY()[ixi];
                
                //----- calc adjusted dir v2j -----
                double v2dir1 = rrwt*v2dir0 + rcwt*v2dec0;
                double ev2dir1 = sqrt(pow(rrwt*ev2dir0,2.) + pow(rcwt*ev2dec0,2.));
                double sev2dir1 = sqrt(pow(rrwt*sev2dir0,2.) + pow(rcwt*sev2dec0,2.));

                //----- calc adjusted dir count -----
                double ndir1 = rrwt + rcwt;

                //----- apply normalization -----
                ydir1 /= ndir1;
                eydir1 /= ndir1;
                seydir1 /= ndir1;

                v2dir1 /= ndir1;
                ev2dir1 /= ndir1;
                sev2dir1 /= ndir1;
                
                //----- set tgraphs points -----
                double xibw = .8;
                double x = xi;
                //yields
                ixiyield_lateadj_dir[iptt][icent]->SetPoint(ixi,x,ydir1);
                ixiyield_lateadj_dir[iptt][icent]->SetPointError(ixi,0,eydir1);

                ixiyield_syserr_lateadj_dir[iptt][icent]->SetPoint(ixi,x,ydir1);

                //v2
                ixiv2jfromfit_lateadj_dir[iptt][icent]->SetPoint(ixi,x,v2dir1);
                ixiv2jfromfit_lateadj_dir[iptt][icent]->SetPointError(ixi,0,ev2dir1);

                ixiv2jfromfit_syserr_lateadj_dir[iptt][icent]->SetPoint(ixi,x,v2dir1);
            }

            //----- save our tgraphs -----
            ixiyield_lateadj_dir[iptt][icent]->Write();
            ixiv2jfromfit_lateadj_dir[iptt][icent]->Write();
        }
        for (int ixi = 0; ixi < nxi; ixi++)
        {
            ipttyield_lateadj_dir[icent][ixi] = new TGraphErrors(nptt);
            ipttyield_lateadj_dir[icent][ixi]->SetName(Form("ipttyield_adj1_dir_cent%d_dpsi2_xi%d",icent,ixi));
            ipttyield_lateadj_dir[icent][ixi]->SetTitle(";p_{TT};1/N_{#gamma} dN/dp_{TT}");
                
            ipttyield_syserr_lateadj_dir[icent][ixi] = new TGraphErrors(nptt);
            ipttyield_syserr_lateadj_dir[icent][ixi]->SetName(Form("ipttyield_syserr_adj1_dir_cent%d_dpsi2_xi%d",icent,ixi));
            ipttyield_syserr_lateadj_dir[icent][ixi]->SetTitle(";p_{TT};1/N_{#gamma} dN/dp_{TT}");
            
            ipttv2jfromfit_lateadj_dir[icent][ixi] = new TGraphErrors(nptt);
            ipttv2jfromfit_lateadj_dir[icent][ixi]->SetName(Form("ipttv2jfromfit_adj1_dir_cent%d_xi%d",icent,ixi));
            ipttv2jfromfit_lateadj_dir[icent][ixi]->SetTitle(";p_{TT};1/N_{#gamma} dN/dp_{TT}");
            
            ipttv2jfromfit_syserr_lateadj_dir[icent][ixi] = new TGraphErrors(nptt);
            ipttv2jfromfit_syserr_lateadj_dir[icent][ixi]->SetName(Form("ipttv2jfromfit_syserr_adj1_dir_cent%d_xi%d",icent,ixi));
            ipttv2jfromfit_syserr_lateadj_dir[icent][ixi]->SetTitle(";p_{TT};1/N_{#gamma} dN/dp_{TT}");

            for (int iptt = 0; iptt < nptt; iptt++)
            {
                double rrwt, rcwt, rrwt_asr, rcwt_asr;
                GetWeights(icent,iptt,rrwt,rcwt,rrwt_asr,rcwt_asr,1.2);

                rrwt *= ndir_adj0[icent][iptt];
                rrwt_asr *= ndir_adj0[icent][iptt];

                rcwt *= ndec_adj0[icent][iptt];
                rcwt_asr *= ndec_adj0[icent][iptt];

                double ptt = ipttyield_adj0[DIR][icent][ixi]->GetX()[iptt];

                //----- yields -----
                double ydir0 = ipttyield_adj0[DIR][icent][ixi]->GetY()[iptt];
                double eydir0 = ipttyield_adj0[DIR][icent][ixi]->GetEY()[iptt];
                double seydir0 = ipttyield_syserr_adj0[DIR][icent][ixi]->GetEY()[iptt];

                double ydec0 = ipttyield_adj0[DEC][icent][ixi]->GetY()[iptt];
                double eydec0 = ipttyield_adj0[DEC][icent][ixi]->GetEY()[iptt];
                double seydec0 = ipttyield_syserr_adj0[DEC][icent][ixi]->GetEY()[iptt];
                
                //----- calc adjusted dir yield -----
                double ydir1 = rrwt_asr*ydir0 + rcwt_asr*ydec0;
                double eydir1 = sqrt(pow(rrwt_asr*eydir0,2.) + pow(rcwt_asr*eydec0,2.));
                double seydir1 = sqrt(pow(rrwt_asr*seydir0,2.) + pow(rcwt_asr*seydec0,2.));

                //----- v2jfromfit -----
                double v2dir0 = ipttv2jfromfit_adj0[DIR][icent][ixi]->GetY()[iptt];
                double ev2dir0 = ipttv2jfromfit_adj0[DIR][icent][ixi]->GetEY()[iptt];
                double sev2dir0 = ipttv2jfromfit_syserr_adj0[DIR][icent][ixi]->GetEY()[iptt];

                double v2dec0 = ipttv2jfromfit_adj0[DEC][icent][ixi]->GetY()[iptt];
                double ev2dec0 = ipttv2jfromfit_adj0[DEC][icent][ixi]->GetEY()[iptt];
                double sev2dec0 = ipttv2jfromfit_syserr_adj0[DEC][icent][ixi]->GetEY()[iptt];
                
                //----- calc adjusted dir v2j -----
                double v2dir1 = rrwt*v2dir0 + rcwt*v2dec0;
                double ev2dir1 = sqrt(pow(rrwt*ev2dir0,2.) + pow(rcwt*ev2dec0,2.));
                double sev2dir1 = sqrt(pow(rrwt*sev2dir0,2.) + pow(rcwt*sev2dec0,2.));

                //----- calc adjusted dir count -----
                double ndir1 = rrwt + rcwt;

                //----- apply normalization -----
                ydir1 /= ndir1;
                eydir1 /= ndir1;
                seydir1 /= ndir1;

                v2dir1 /= ndir1;
                ev2dir1 /= ndir1;
                sev2dir1 /= ndir1;
                
                //----- set tgraphs points -----
                double pttbw = 2.;
                if (iptt > 1) pttbw = 3.;                   
                double x = ptt;
                //yields
                ipttyield_lateadj_dir[icent][ixi]->SetPoint(iptt,x,ydir1);
                ipttyield_lateadj_dir[icent][ixi]->SetPointError(iptt,0,eydir1);

                ipttyield_syserr_lateadj_dir[icent][ixi]->SetPoint(iptt,x,ydir1);

                //v2
                ipttv2jfromfit_lateadj_dir[icent][ixi]->SetPoint(iptt,x,v2dir1);
                ipttv2jfromfit_lateadj_dir[icent][ixi]->SetPointError(iptt,0,ev2dir1);

                ipttv2jfromfit_syserr_lateadj_dir[icent][ixi]->SetPoint(iptt,x,v2dir1);
            }

            //----- save our tgraphs -----
            ipttyield_lateadj_dir[icent][ixi]->Write();
            ipttv2jfromfit_lateadj_dir[icent][ixi]->Write();
        }
    }
            
    MakeRatioGraphs();

    //set sys err from oldsyserr then write
    GetOldSysErr(input_yield,input_v2);

    for (int iptt = 0; iptt < nptt; iptt++)
    {
        ixicentrat_syserr_lateadj_dir[iptt]->Write();

        for (int icent = 0; icent < ncent; icent++)
        {
            ixiyield_syserr_lateadj_dir[iptt][icent]->Write();
            ixiv2jfromfit_syserr_lateadj_dir[iptt][icent]->Write();
        }
    }

    for (int ixi = 0; ixi < nxi; ixi++)
    {
        ipttcentrat_syserr_lateadj_dir[ixi]->Write();

        for (int icent = 0; icent < ncent; icent++)
        {
            ipttyield_syserr_lateadj_dir[icent][ixi]->Write();
            ipttv2jfromfit_syserr_lateadj_dir[icent][ixi]->Write();
        }
    }
        
    cout<<"Closing output "<<output->GetName()<<endl;
    output->Write();
    output->Close();
}

TFile* OpenFile(const char* name)
{
    cout<<"Reading file "<<name<<endl;

    TFile* file = new TFile(name,"READ");

    if (!file)
    {
        cerr<<"Could not open "<<name<<endl;
    }
    if (file->IsZombie())
    {
        cerr<<"File "<<name<<" is zombie!"<<endl;
    }
    
    return file;
}

void GetGraphs(TFile* input_yield, TFile* input_v2)
{
    for (int ispec = 0; ispec < nspec; ispec++)
    {
        for (int icent = 0; icent < ncent; icent++)
        {
            for (int iptt = 0; iptt < nptt; iptt++)
            {
                input_yield->GetObject(Form("ixiyield_adj0_%s_ptt%d_cent%d_dpsi2",cspec[ispec],iptt,icent),ixiyield_adj0[ispec][iptt][icent]);
                input_yield->GetObject(Form("ixiyield_syserr_adj0_%s_ptt%d_cent%d_dpsi2",cspec[ispec],iptt,icent),ixiyield_syserr_adj0[ispec][iptt][icent]);

                input_v2->GetObject(Form("ixiv2jfromfit_adj0_%s_ptt%d_cent%d",cspec[ispec],iptt,icent),ixiv2jfromfit_adj0[ispec][iptt][icent]);
                ixiv2jfromfit_adj0[ispec][iptt][icent]->Print();
                input_v2->GetObject(Form("ixiv2jfromfit_syserr_adj0_%s_ptt%d_cent%d",cspec[ispec],iptt,icent),ixiv2jfromfit_syserr_adj0[ispec][iptt][icent]);
            }
            for (int ixi = 0; ixi < nxi; ixi++)
            {
                input_yield->GetObject(Form("ipttyield_adj0_%s_cent%d_dpsi2_xi%d",cspec[ispec],icent,ixi),ipttyield_adj0[ispec][icent][ixi]);
                input_yield->GetObject(Form("ipttyield_syserr_adj0_%s_cent%d_dpsi2_xi%d",cspec[ispec],icent,ixi),ipttyield_syserr_adj0[ispec][icent][ixi]);

                input_v2->GetObject(Form("ipttv2jfromfit_adj0_%s_cent%d_xi%d",cspec[ispec],icent,ixi),ipttv2jfromfit_adj0[ispec][icent][ixi]);
                input_v2->GetObject(Form("ipttv2jfromfit_syserr_adj0_%s_cent%d_xi%d",cspec[ispec],icent,ixi),ipttv2jfromfit_syserr_adj0[ispec][icent][ixi]);
            }
        }
    }
}

int GetFunctions(const char* loc, const char* name, TF1* f1[ncent_fun])
{
    cout<<"Getting functions "<<name<<"[0-"<<ncent_fun-1<<"] from "<<loc<<endl;

    memset(f1,0,sizeof(TF1*)*ncent_fun);

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

    for (int icent = 0; icent < ncent_fun; icent++)
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

void GetWeights(const int icent,const int iptt,double& rrwt,double& rcwt,double& rrwt_asr,double& rcwt_asr,const double r = 1)
{
    double m = mis[icent][iptt];
    double t = tag[icent][iptt];
    double d = (-r*t*m + r*t + t*m - m);
    
    rrwt = t/(t-m);
    rrwt_asr = r*t/d;

    rcwt = (1-t)/(m-t);
    rcwt_asr = (1-t)/-d;

    printf("Weights obtained for icent%d iptt%d with m %f, t %f : rrwt %f, rcwt %f, rrwt_asr %f, rcwt_asr %f\n",icent,iptt,m,t,rrwt,rcwt,rrwt_asr,rcwt_asr);
}

int GetPttBin(double ptt)
{
    if (ptt < 5.)
        return -1.;
    else if (ptt < 7.)
        return 0;
    else if (ptt < 9.)
        return 1;
    else if (ptt < 12.)
        return 2;
    else if (ptt < 15.)
        return 3;
    else
        return -2;
}
    

void FillSysError(TFile* input)
{

    TTree* tree;
    input->GetObject("Terrorinfo",tree);
    errorinfo_t info;
    tree->SetBranchAddress("info",&info);

    int nentries = tree->GetEntries();

    for (int ientry = 0; ientry < nentries; ientry++)
    {
        tree->GetEntry(ientry);
        errinfo[info.iadj][info.itype][info.iptt][info.icent][info.idpsi][info.ixi] = info;
    }                       
}

void MakeRatioGraphs()
{
    for (int iptt = 0; iptt < nptt; iptt++)
    {
        for (int ixi = 0; ixi < nxi; ixi++)
        {
            double centyield = ixiyield_lateadj_dir[0][ixi]->GetY()[iptt];
            double midcentyield = ixiyield_lateadj_dir[1][ixi]->GetY()[iptt];
            double ratio = centyield/midcentyield;

            double e0 = ixiyield_lateadj_dir[0][ixi]->GetEY()[iptt];
            double e1 = ixiyield_lateadj_dir[1][ixi]->GetEY()[iptt];
            e0 /= midcentyield;
            e1 = -e1 * ratio / midcentyield;

            double error = sqrt(e0*e0 + e1*e1);

            double syserror = -9999.;

            if (false)
            {
                // sys error section not yet implemented for adj1 dir
                double a0 = PI/3.;
                double a1 = sqrt(3.)/4.;
                double fe0c = errinfo[0][DIR][iptt][0][ALL][ixi].flate[0];                     
                fe0c *= a0;
                double ce0c = errinfo[0][DIR][iptt][0][ALL][ixi].cose[0];
                ce0c *= a1;

                double ntrgp = errinfo[0][DIR][iptt][1][ALL][ixi].ntrg;
                double fe0p = errinfo[0][DIR][iptt][1][ALL][ixi].flate[0];
                fe0p *= a0;
                double ce0p = errinfo[0][DIR][iptt][1][ALL][ixi].cose[0];
                ce0p *= a1;

                double EC = (fe0c*fe0c + ce0c*ce0c)/(centyield*centyield);
                double EP = (fe0p*fe0p + ce0p*ce0p)/(midcentyield*midcentyield);
                syserror = sqrt(ratio*ratio*(EC+EP));
            }
            else
                syserror = 0;

            //---------- cent ratio vs xi ----------
            if (!ixicentrat_lateadj_dir[iptt])
            {
                ixicentrat_lateadj_dir[iptt] = new TGraphErrors();
                ixicentrat_lateadj_dir[iptt]->SetName(Form("ixicentrat_adj1_dir_ptt%i",iptt));
                ixicentrat_lateadj_dir[iptt]->SetTitle(";#xi;Y_{0-20%}/Y_{20-60%}");

                ixicentrat_syserr_lateadj_dir[iptt] = new TGraphErrors();
                ixicentrat_syserr_lateadj_dir[iptt]->SetName(Form("ixicentrat_syserr_adj1_dir_ptt%i",iptt));
                ixicentrat_syserr_lateadj_dir[iptt]->SetTitle(";#xi;Y_{0-20%}/Y_{20-60%}");
            }
                    
            {
                double x = ixiyield_lateadj_dir[iptt][0]->GetX()[ixi];
                double sex = ixiyield_syserr_lateadj_dir[iptt][0]->GetEX()[ixi];

                ixicentrat_lateadj_dir[iptt]->SetPoint(ixi,x,ratio);
                ixicentrat_lateadj_dir[iptt]->SetPointError(ixi,0,error);

                ixicentrat_syserr_lateadj_dir[iptt]->SetPoint(ixi,x,ratio);
            }

            //---------- cent ratio vs ptt ----------
            if (!ipttcentrat_lateadj_dir[ixi])
            {
                ipttcentrat_lateadj_dir[ixi] = new TGraphErrors();
                ipttcentrat_lateadj_dir[ixi]->SetName(Form("ipttcentrat_adj1_dir_xi%i",ixi));
                ipttcentrat_lateadj_dir[ixi]->SetTitle(";p_{TT};Y_{0-20%}/Y_{20-60%}");

                ipttcentrat_syserr_lateadj_dir[ixi] = new TGraphErrors();
                ipttcentrat_syserr_lateadj_dir[ixi]->SetName(Form("ipttcentrat_syserr_adj1_dir_xi%i",ixi));
                ipttcentrat_syserr_lateadj_dir[ixi]->SetTitle(";p_{TT};Y_{0-20%}/Y_{20-60%}");
            }

            {
                double x = ipttyield_lateadj_dir[0][ixi]->GetX()[iptt];
                double sex = ipttyield_syserr_lateadj_dir[0][ixi]->GetEX()[iptt];

                ipttcentrat_lateadj_dir[ixi]->SetPoint(iptt,x,ratio);
                ipttcentrat_lateadj_dir[ixi]->SetPointError(iptt,0,error);

                ipttcentrat_syserr_lateadj_dir[ixi]->SetPoint(iptt,x,ratio);
            }
        }
    }

    for (int iptt = 0; iptt < nptt; iptt++)
    {
        ixicentrat_lateadj_dir[iptt]->Write();
    }
    for (int ixi = 0; ixi < nxi; ixi++)
    {
        ipttcentrat_lateadj_dir[ixi]->Write();
    }
}

void GetOldSysErr(TFile* input_yield, TFile* input_v2)
{
    for (int iptt = 0; iptt < nptt; iptt++)
    {
        TGraphErrors* ixicentrat;
        input_yield->GetObject(Form("ixicentrat_syserr_adj1_dir_ptt%i",iptt),ixicentrat);

        for (int ixi = 0; ixi < nxi; ixi++)
        {
            ixicentrat_syserr_lateadj_dir[iptt]->SetPointError(ixi,ixicentrat->GetEX()[ixi],ixicentrat->GetEY()[ixi]);
        }

        for (int icent = 0; icent < ncent; icent++)
        {
            TGraphErrors* ixiv2jfromfit;
            input_v2->GetObject(Form("ixiv2jfromfit_syserr_adj1_dir_ptt%i_cent%i",iptt,icent),ixiv2jfromfit);

            TGraphErrors* ixiyield;
            input_yield->GetObject(Form("ixiyield_syserr_adj1_dir_ptt%i_cent%i_dpsi2",iptt,icent),ixiyield);

            for (int ixi = 0; ixi < nxi; ixi++)
            {
                ixiv2jfromfit_syserr_lateadj_dir[iptt][icent]->SetPointError(ixi,ixiv2jfromfit->GetEX()[ixi],ixiv2jfromfit->GetEY()[ixi]);
                ixiyield_syserr_lateadj_dir[iptt][icent]->SetPointError(ixi,ixiyield->GetEX()[ixi],ixiyield->GetEY()[ixi]);
            }
        }
    }

    for (int ixi = 0; ixi < nxi; ixi++)
    {
        TGraphErrors* ipttcentrat;
        input_yield->GetObject(Form("ipttcentrat_syserr_adj1_dir_xi%i",ixi),ipttcentrat);
        for (int iptt = 0; iptt < nptt; iptt++)
        {
            ipttcentrat_syserr_lateadj_dir[ixi]->SetPointError(iptt,ipttcentrat->GetEX()[iptt],ipttcentrat->GetEY()[iptt]);
        }

        for (int icent = 0; icent < ncent; icent++)
        {
            TGraphErrors* ipttv2jfromfit;
            input_v2->GetObject(Form("ipttv2jfromfit_syserr_adj1_dir_cent%i_xi%i",icent,ixi),ipttv2jfromfit);

            TGraphErrors* ipttyield;
            input_yield->GetObject(Form("ipttyield_syserr_adj1_dir_cent%i_dpsi2_xi%i",icent,ixi),ipttyield);
            
            for (int iptt = 0; iptt < nptt; iptt++)
            {
                ipttv2jfromfit_syserr_lateadj_dir[icent][ixi]->SetPointError(iptt,ipttv2jfromfit->GetEX()[iptt],ipttv2jfromfit->GetEY()[iptt]);
                ipttyield_syserr_lateadj_dir[icent][ixi]->SetPointError(iptt,ipttyield->GetEX()[iptt],ipttyield->GetEY()[iptt]);
            }
        }
    }
}
