//recreate adj histos from nonadj. using averages of tageff and mistag

//draw dphi plots resulting from subtractWeightedBG.C

const double asr = 1.21;

const int nadj = 2;
const int ncent = 2;
const int limcent[ncent+1] = {0,20,60};

enum RPO {IN,OUT,ALL};
string strplane[3] = {"in-plane","out-of-plane","all-plane"};

const int markstyle[] = {20,21,33,47,34};

const float PI = TMath::Pi();

const int ncent_fun = 94;

TF1* ftageff[ncent_fun];
TF1* fdeceff[ncent_fun];
TF1* fmistag[ncent_fun];


int GetFunctions(const char* loc, const char* name, TF1* f1[ncent_fun]);

namespace GAT {
    #include "gamAnaTree.h"
}

const double ptavg[GAT::nptt] = {5.67,7.76,10.08,13.15};
const int icentavg[ncent] = {8,34};

int GetWeights(double& dirwt, double& decwt, double& dirwtY, double& decwtY, TF1* ftageff[ncent_fun], TF1* fmistag[ncent_fun], int itype, int icent, int ipt);


//in
// adj0_dir_ptt0 sumcent0 257531
// adj0_dir_ptt0 sumcent1 213041
// adj0_dir_ptt1 sumcent0 45659
// adj0_dir_ptt1 sumcent1 32815
// adj0_dir_ptt2 sumcent0 14776
// adj0_dir_ptt2 sumcent1 9716
// adj0_dir_ptt3 sumcent0 3205
// adj0_dir_ptt3 sumcent1 2012

//out
// adj0_dir_ptt0 sumcent0 242709
// adj0_dir_ptt0 sumcent1 180565
// adj0_dir_ptt1 sumcent0 43664
// adj0_dir_ptt1 sumcent1 29116
// adj0_dir_ptt2 sumcent0 13970
// adj0_dir_ptt2 sumcent1 8729
// adj0_dir_ptt3 sumcent0 3044
// adj0_dir_ptt3 sumcent1 1806

//in
// adj0_dec_ptt0 sumcent0 60786
// adj0_dec_ptt0 sumcent1 51125
// adj0_dec_ptt1 sumcent0 6265
// adj0_dec_ptt1 sumcent1 6044
// adj0_dec_ptt2 sumcent0 1429
// adj0_dec_ptt2 sumcent1 1369
// adj0_dec_ptt3 sumcent0 201
// adj0_dec_ptt3 sumcent1 164


//out
// adj0_dec_ptt0 sumcent0 52877
// adj0_dec_ptt0 sumcent1 40773
// adj0_dec_ptt1 sumcent0 5809
// adj0_dec_ptt1 sumcent1 5010
// adj0_dec_ptt2 sumcent0 1292
// adj0_dec_ptt2 sumcent1 1173
// adj0_dec_ptt3 sumcent0 158
// adj0_dec_ptt3 sumcent1 150


double ndir[ncent][GAT::ndpsi][GAT::nptt] = {
    //sumcent0
    {
        //in
        {257531,45659,14776,3205},
        //out
        {242709,43664,13970,3044}
    },
    //sumcent1
    {
        //in
        {213041,32815,9716,2012},
        //out
        {180565,29116,8729,1806}
    }
};

double ndec[ncent][GAT::ndpsi][GAT::nptt] = {
    //sumcent0
    {
        //in
        {60786,6265,1429,201},
        //out
        {52877,5809,1292,158}
    },
    //sumcent1
    {
        //in
        {51125,6044,1369,164},
        //out
        {40773,5010,1173,150}
    }
};

void drawDphiPlots_LateAdj(
    const char* inname = "weighted_bg_subtracted.root"
    )
{

    // //get Ntrg
    // //-------------------- in-plane case --------------------
    // ntrgin[centbin] += ntrg;

    
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
    
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0.);
    gStyle->SetEndErrorSize(30);
    gStyle->SetMarkerSize(4);
    gStyle->SetHistLineWidth(2);
    gStyle->SetTitleOffset(1.4,"Y");
    gStyle->SetImageScaling(3.);
    gROOT->ForceStyle();

    TFile* input = new TFile(inname,"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Error opening "<<inname<<endl;
        return;
    }

    TH1F* h1dphi[GAT::ntype][ncent][GAT::ndpsi+1][GAT::nptt][GAT::nxi] {};
    TH1F* h1dphi_redux[2][ncent][GAT::ndpsi+1][GAT::nptt][GAT::nxi] {};

    for (int itype = 0; itype < GAT::ntype; itype++)
    {
        if (itype == GAT::PI0) continue;
        for (int icent = 0; icent < ncent; icent++)
        {
            for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        ostringstream oss;
                        oss<<"dphi_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<icent<<"_dpsi"<<idpsi<<"_xi"<<ixi;

                        input->GetObject(oss.str().c_str(),h1dphi[itype][icent][idpsi][iptt][ixi]);
                        if (itype == GAT::DIR)
                        {
                            if (idpsi == GAT::ndpsi)
                            {
                                h1dphi[itype][icent][idpsi][iptt][ixi]->Scale(ndir[icent][0][iptt]+ndir[icent][1][iptt]);
                            }
                            else
                            {
                                h1dphi[itype][icent][idpsi][iptt][ixi]->Scale(ndir[icent][idpsi][iptt]);
                            }
                        }

                        if (itype == GAT::DEC)
                        {
                            if (idpsi == GAT::ndpsi)
                            {
                                h1dphi[itype][icent][idpsi][iptt][ixi]->Scale(ndec[icent][0][iptt]+ndec[icent][1][iptt]);
                            }
                            else
                            {
                                h1dphi[itype][icent][idpsi][iptt][ixi]->Scale(ndec[icent][idpsi][iptt]);
                            }
                        }
                    }
                }
            }
        }
    }

    //calc the redux -- need ntrig
    for (int icent = 0; icent < ncent; icent++)
    {
        for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
        {
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    {
                        ostringstream oss;
                        oss<<"dphi_redux_dir_ptt"<<iptt<<"_sumcent"<<icent<<"_dpsi"<<idpsi<<"_xi"<<ixi;

                        double ntrg = 0.;

                        for (int itype = 0; itype < GAT::ntype; itype++)
                        {
                            if (itype == GAT::PI0) continue;

                            double dirwt,decwt,dirwtY,decwtY;
                            GetWeights(dirwt,decwt,dirwtY,decwtY,ftageff,fmistag,itype,icent,iptt);

                            if (itype == 0)
                            {
                                h1dphi_redux[GAT::DIR][icent][idpsi][iptt][ixi] = (TH1F*) h1dphi[itype][icent][idpsi][iptt][ixi]->Clone(oss.str().c_str());
                                h1dphi_redux[GAT::DIR][icent][idpsi][iptt][ixi]->Scale(dirwtY);
                            }
                            else
                            {
                                h1dphi_redux[GAT::DIR][icent][idpsi][iptt][ixi]->Add(h1dphi[itype][icent][idpsi][iptt][ixi],dirwtY);
                            }

                            if (idpsi == GAT::ndpsi)
                            {
                                if (itype == GAT::DIR)
                                    ntrg += (ndir[icent][0][iptt]+ndir[icent][1][iptt])*dirwt;
                                else
                                    ntrg += (ndec[icent][0][iptt]+ndec[icent][1][iptt])*dirwt;
                            }
                            else
                            {
                                if (itype == GAT::DIR)
                                    ntrg += ndir[icent][idpsi][iptt]*dirwt;
                                else
                                    ntrg += ndec[icent][idpsi][iptt]*dirwt;
                            }
                        }

                        h1dphi_redux[GAT::DIR][icent][idpsi][iptt][ixi]->Scale(1./ntrg);
                    }

                        
                    {
                        ostringstream oss;
                        oss<<"dphi_redux_dec_ptt"<<iptt<<"_sumcent"<<icent<<"_dpsi"<<idpsi<<"_xi"<<ixi;

                        double ntrg = 0.;

                        for (int itype = 0; itype < GAT::ntype; itype++)
                        {
                            if (itype == GAT::PI0) continue;

                            double dirwt,decwt,dirwtY,decwtY;
                            GetWeights(dirwt,decwt,dirwtY,decwtY,ftageff,fmistag,itype,icent,iptt);

                            if (itype == 0)
                            {
                                h1dphi_redux[GAT::DEC][icent][idpsi][iptt][ixi] = (TH1F*) h1dphi[itype][icent][idpsi][iptt][ixi]->Clone(oss.str().c_str());
                                h1dphi_redux[GAT::DEC][icent][idpsi][iptt][ixi]->Scale(decwtY);
                            }
                            else
                            {
                                h1dphi_redux[GAT::DEC][icent][idpsi][iptt][ixi]->Add(h1dphi[itype][icent][idpsi][iptt][ixi],decwtY);
                            }

                            if (idpsi == GAT::ndpsi)
                            {
                                if (itype == GAT::DIR)
                                    ntrg += (ndir[icent][0][iptt]+ndir[icent][1][iptt])*decwt;
                                else
                                    ntrg += (ndec[icent][0][iptt]+ndec[icent][1][iptt])*decwt;
                            }
                            else
                            {
                                if (itype == GAT::DIR)
                                    ntrg += ndir[icent][idpsi][iptt]*decwt;
                                else
                                    ntrg += ndec[icent][idpsi][iptt]*decwt;
                            }
                        }

                        h1dphi_redux[GAT::DEC][icent][idpsi][iptt][ixi]->Scale(1./ntrg);
                    }


                    // {
                    //     h1dphi_redux[itype][icent][idpsi][iptt][ixi]->Scale(1./((ndir[icent][0][iptt]+ndir[icent][1][iptt])*dirwt + (ndec[icent][0][iptt]+ndec[icent][1][iptt])*decwt));
                    // }
                    // {
                    //     h1dphi_redux[itype][icent][idpsi][iptt][ixi]->Scale(1./(ndir[icent][idpsi][iptt]*dirwt + ndec[icent][idpsi][iptt]*decwt));
                    // }

                    // oss<<"dphi_redux_dir_ptt"<<iptt<<"_sumcent"<<icent<<"_dpsi"<<idpsi<<"_xi"<<ixi;

                    // // GetWeights(dirwt,decwt,dirwtY,decwtY,ftageff,fmistag,itype,icent,iptt);

                    // // h1dphi_redux[itype][icent][idpsi][iptt][ixi]->Scale(dirwtY);
                    // // h1dphi_redux[itype][icent][idpsi][iptt][ixi]->Add(h1dphi[GAT::DIR][icent][idpsi][iptt][ixi],decwtY);
                        
                    // // {
                    // //     h1dphi_redux[itype][icent][idpsi][iptt][ixi]->Scale(1./((ndir[icent][0][iptt]+ndir[icent][1][iptt])*dirwt + (ndec[icent][0][iptt]+ndec[icent][1][iptt])*decwt));
                    // // }
                    // // {
                    // //     h1dphi_redux[itype][icent][idpsi][iptt][ixi]->Scale(1./(ndir[icent][idpsi][iptt]*dirwt + ndec[icent][idpsi][iptt]*decwt));
                    // // }

                }
            }
        }
    }


    //in and out same plot
    for (int iredux = 0; iredux < 2; iredux++)
    {
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0) continue;
            for (int icent = 0; icent < ncent; icent++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        TCanvas can ("can","",2500,2000);//2x3
                        TMultiGraph mgr;
                        TLegend leg (0.65,0.75,0.9,0.9);
                        leg.SetTextSize(0.04);
                        
                        for (int idpsi = 0; idpsi < GAT::ndpsi; idpsi++)
                        {
                            TH1F* h1ptr = h1dphi[itype][icent][idpsi][iptt][ixi];
                            if (iredux) h1ptr = h1dphi_redux[itype][icent][idpsi][iptt][ixi];
                        
                            TGraphErrors* tgr = new TGraphErrors(h1ptr);
                            tgr->SetMarkerColor(TColor::GetColorPalette(440/3*idpsi));
                            tgr->SetLineColor(TColor::GetColorPalette(440/3*idpsi));
                            tgr->SetMarkerStyle(markstyle[idpsi]);

                            mgr.Add(tgr,"p");
                            TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                            tgrclone->SetMarkerSize(0);                 
                            mgr.Add(tgrclone,"p");
                            leg.AddEntry(tgr,Form("%s",strplane[idpsi].c_str()),"p");
                        }

                        mgr.SetTitle(";#Delta#phi/#pi;1/N_{#gamma}  dN_{h}/d#Delta#phi");
                        mgr.Draw("A");
                        leg.Draw();

                        TLatex txtmom (0.4,.85,Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                        txtmom.SetTextAlign(22);
                        txtmom.SetNDC();
                        txtmom.SetTextSize(0.04);
                        txtmom.Draw();


                        TLatex txtcent (0.4,.8,Form("Au+Au %i-%i%% 200GeV",limcent[icent],limcent[icent+1]));
                        txtcent.SetTextAlign(22);
                        txtcent.SetNDC();
                        txtcent.SetTextSize(0.04);
                        txtcent.Draw();

                        double ymin = min(0.,mgr.GetHistogram()->GetMinimum());
                        double ymax = max(0.,mgr.GetHistogram()->GetMaximum());

                        ymin = ymin - (ymax-ymin)*.25;
                        ymax = ymax + (ymax-ymin)*.4;
    
                        mgr.GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);
                        mgr.GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                        TLine line(0.,0.,1.,0.);
                        line.SetLineWidth(2);
                        line.SetLineStyle(2);
                        line.Draw();

                        mgr.Draw();

                        can.Update();
                        can.SaveAs(Form("dphiplots_redux/inout_%s_redux%i_cent%i_ptt%i_xi%i.eps",GAT::spec[itype],iredux,icent,iptt,ixi));
                    }
                }
            }
        }
    }    


    //central and mid central same plot

    for (int iredux = 0; iredux < 2; iredux++)
    {
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0) continue;
            
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    TCanvas can ("can","",2500,2000);//2x3
                    TMultiGraph mgr;
                    TLegend leg (0.75,0.75,0.9,0.9);
                    leg.SetTextSize(0.04);
                        
                    for (int icent = 0; icent < ncent; icent++)
                    {
                        TH1F* h1ptr = h1dphi[itype][icent][ALL][iptt][ixi];
                        if (iredux) h1ptr = h1dphi_redux[itype][icent][ALL][iptt][ixi];
                        
                        TGraphErrors* tgr = new TGraphErrors(h1ptr);
                    
                        tgr->SetMarkerColor(TColor::GetColorPalette(440/3*icent));
                        tgr->SetLineColor(TColor::GetColorPalette(440/3*icent));
                        tgr->SetMarkerStyle(markstyle[icent]);

                        mgr.Add(tgr,"p");
                        TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                        tgrclone->SetMarkerSize(0);                     
                        mgr.Add(tgrclone,"p");
                        leg.AddEntry(tgr,Form("%i-%i%%",limcent[icent],limcent[icent+1]),"p");
                    }

                    mgr.SetTitle(";#Delta#phi/#pi;1/N_{#gamma}  dN_{h}/d#Delta#phi");
                    mgr.Draw("A");
                    leg.Draw();

                    TLatex txtmom (0.5,.85,Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                    txtmom.SetTextAlign(22);
                    txtmom.SetNDC();
                    txtmom.SetTextSize(0.04);
                    txtmom.Draw();


                    TLatex txtcent (0.5,.8,"Au+Au 200GeV");
                    txtcent.SetTextAlign(22);
                    txtcent.SetNDC();
                    txtcent.SetTextSize(0.04);
                    txtcent.Draw();

                    double ymin = min(0.,mgr.GetHistogram()->GetMinimum());
                    double ymax = max(0.,mgr.GetHistogram()->GetMaximum());

                    ymin = ymin - (ymax-ymin)*.25;
                    ymax = ymax + (ymax-ymin)*.4;
    
                    mgr.GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);
                    mgr.GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                    TLine line(0.,0.,1.,0.);
                    line.SetLineWidth(2);
                    line.SetLineStyle(2);
                    line.Draw();

                    mgr.Draw();

                    can.Update();
                    can.SaveAs(Form("dphiplots_redux/central_%s_redux%i_ptt%i_xi%i.eps",GAT::spec[itype],iredux,iptt,ixi));
                }
            }
        }
    }
}

int GetFunctions(const char* loc, const char* name, TF1* f1[ncent_fun])
{
    cout<<"Getting functions "<<name<<"? from "<<loc<<endl;

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

    cout<<"Input "<<loc<<" opened successfully."<<endl;

    for (int icent = 0; icent < ncent_fun; icent++)
    {
        ostringstream oss;
        oss<<name<<icent;

        // {
        // }

        // f1[icent]->SetDirectory(0);

        TF1* f1tmp = (TF1*)file->Get(oss.str().c_str());
        if (f1tmp == NULL)
        {
            cerr<<"ABORT!!"<<endl;
            return -1;
        }
        
        cout<<"Retrieved "<<f1tmp->GetName()<<endl;

        f1[icent] = (TF1*) f1tmp->Clone();
    }

    file->Close();
    delete file;

    return 0;
}

// {
//     {
//     }
//     {
//     }
// }

// {
//     {
//     }
// }

int GetWeights(double& dirwt, double& decwt, double& dirwtY, double& decwtY, TF1* ftageff[ncent_fun], TF1* fmistag[ncent_fun], int itype, int icent, int ipt)
{
    double tageff = ftageff[icentavg[icent]]->Eval(ptavg[ipt]);
    double mistag = fmistag[icentavg[icent]]->Eval(ptavg[ipt]);

    double r = asr, t = tageff, m = mistag;
    double d = (-r*t*m + r*t + t*m - m);
        
    if (itype == GAT::DIR)
    {
        dirwt = t/(t-m);
        decwt = m/(m-t);

        dirwtY = r*t/d;
        decwtY = m*(1-t+t*r)/-d;
    }
    else if (itype == GAT::DEC)
    {
        dirwt = (1-t)/(m-t);
        decwt = (1-m)/(t-m);

        dirwtY = (1-t)/-d;
        decwtY = (1-m)*(1-t*t*r)/d;
    }

    return 0;
}
