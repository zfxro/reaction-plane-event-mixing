//same as drawAcceptCorrect but use tgraphs -- this serves as a check

//draw acceptance correction (result of acceptCorrect.C)

const float PI = TMath::Pi();
const int nadj = 2;
const int nmix = 2;
const int ntype = 3;

const int markstyle[] = {20,24,21,25,22,26,23,32,34,47};

namespace GAT {
    #include "../src/gamAnaTree/gamAnaTree.h"
}

int GetTGraphs(const char* pref, const char* inname, TGraphErrors* tgr[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1]);

void drawAcceptCorrect_fromTGraphs(
    const char* inname_raw = "tgraphs_rawnh.root",
    const char* inname_exp = "tgraphs_expnh.root"
    )
{

    TGraphErrors* tgrraw[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1] = {};
    GetTGraphs("tgrhad",inname_raw,tgrraw);

    TGraphErrors* tgrexp[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1] = {};
    GetTGraphs("tgrexpnh",inname_exp,tgrexp);



    TCanvas* can = new TCanvas("can","",1000,1500);
    const int npad = 4;
    TPad* pads[npad] = {};
    
    pads[0] = new TPad("pad0","",0.,.5,.5,1.);
    pads[0]->SetTopMargin(.12);
    pads[0]->SetBottomMargin(0.);
    pads[0]->SetLeftMargin(.18);
    pads[0]->SetRightMargin(0.);
    pads[0]->SetTicks(0,1);
    pads[0]->Draw();

    pads[1] = new TPad("pad1","",.5,.5,1.,1.);
    pads[1]->SetTopMargin(.12);
    pads[1]->SetBottomMargin(0.);
    pads[1]->SetLeftMargin(0.);
    pads[1]->SetRightMargin(.18);
    pads[1]->Draw();
    
    pads[2] = new TPad("pad2","",0.,0.,.5,.5);
    pads[2]->SetTopMargin(0.);
    pads[2]->SetBottomMargin(.12);
    pads[2]->SetLeftMargin(.18);
    pads[2]->SetRightMargin(0.);
    pads[2]->SetTicks(1,1);
    pads[2]->Draw();
    
    pads[3] = new TPad("pad3","",.5,0.,1.,.5);
    pads[3]->SetTopMargin(0.);
    pads[3]->SetBottomMargin(.12);
    pads[3]->SetLeftMargin(0.);
    pads[3]->SetRightMargin(.18);
    pads[3]->SetTicks(1,0);
    pads[3]->Draw();

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;

            for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    pads[iptt]->cd();               

                    TLegend* leg = nullptr;
                    if (iptt == 0)
                    {
                        leg = new TLegend(0.69,0.52,.96,.87);
                    }
                    if (iptt == 1)
                    {
                        leg = new TLegend(0.53,0.64,.80,.87);
                    }

                    TMultiGraph* multi = new TMultiGraph();
                    
                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {
                        vector<double> vacc;
                        vector<double> verr;
                        vector<double> vxi;
                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {                           
                            double raw = (tgrraw[iadj][itype][iptt][icent][idpsi]->GetY())[ixi];
                            double rawerr = (tgrraw[iadj][itype][iptt][icent][idpsi]->GetEY())[ixi];
                            double exp = (tgrexp[iadj][itype][iptt][icent][idpsi]->GetY())[ixi];
                            double experr = (tgrexp[iadj][itype][iptt][icent][idpsi]->GetEY())[ixi];
                            double acc = raw/exp;
                            vacc.push_back(acc);

                            double e0 = rawerr / exp;
                            double e1 = - experr * raw / (exp*exp);
                            double err = sqrt(e0*e0 + e1*e1);
                            verr.push_back(err);

                            vxi.push_back(ixi);
                        }
                            

                        TGraphErrors* tgr = new TGraphErrors(GAT::nxi,&vxi[0],&vacc[0],0,&verr[0]);
                            
                        gStyle->SetPalette(kBird);
                        tgr->SetMarkerColor(TColor::GetColorPalette(254*icent/(GAT::ncent-1)));
                        tgr->SetLineColor(TColor::GetColorPalette(254*icent/(GAT::ncent-1)));

                        tgr->SetMarkerSize(3.0);
                        tgr->SetLineWidth(3);

                        tgr->SetMarkerStyle(markstyle[icent]);

                        multi->Add(tgr,"LP");

                        if ((iptt == 0 && icent < 6) || (iptt == 1 && icent >= 6))
                        {
                            ostringstream oss;
                            oss<<GAT::limcent[icent]<<"-"<<GAT::limcent[icent+1]<<"%"<<endl;
                            leg->AddEntry(tgr,oss.str().c_str(),"lp");
                        }
                    }

                    multi->SetTitle(";#xi;#epsilon");
                    multi->Draw("A");
                    multi->GetHistogram()->GetYaxis()->SetRangeUser(0.065,0.19);
                    multi->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
                    multi->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
                    multi->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
                    multi->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
                    if (leg)
                    {
                        leg->SetMargin(.4);
                        leg->SetTextSize(.05);
                        leg->Draw();
                    }

                    pads[iptt]->Modified();
                    pads[iptt]->Update();

                    TLatex latex;
                    latex.SetTextSize(0.055);
                    latex.DrawLatex(.125,.175,Form("#pi^{0} p_{T} #in [%.1f,%.1f)",GAT::limptt[iptt],GAT::limptt[iptt+1]));


                        
                }

                can->SaveAs(Form("acc_fromtgrs/acc_fromtgrs_adj%i_%s_dpsi%i.png",iadj,GAT::spec[itype],idpsi));
            }
        }
    }
                    

    cout<<"Finished drawing acceptance*efficiency correction."<<endl;
    return;
}

int GetTGraphs(const char* pref, const char* inname, TGraphErrors* tgr[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1])
{
    TFile* input = new TFile(inname,"READ");
    
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {       
                for (int icent = 0; icent < GAT::ncent; icent++)
                {
                    for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                    {
                        TGraphErrors* tgrtmp = (TGraphErrors*) gROOT->FindObject(Form("%s_adj%i_%s_ptt%i_cent%i_dpsi%i",pref,iadj,GAT::spec[itype],iptt,icent,idpsi));
                        tgr[iadj][itype][iptt][icent][idpsi] = (TGraphErrors*) tgrtmp->Clone();
                    }
                }
            }
        }
    }

    input->Close();
    delete input;
    
    return 0;
}
