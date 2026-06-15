//draw acceptance correction (result of acceptCorrect.C)

const float PI = TMath::Pi();
const int nadj = 2;
const int nmix = 2;
const int ntype = 3;

const int markstyle[] = {20,24,21,25,22,26,23,32,34,47};

int GetAccHistos(const char* name);

namespace GAT {
    #include "../src/gamAnaTree/gamAnaTree.h"
}

TH1F* h1acc[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1][GAT::nxi] {};
TGraphErrors* tgracc[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1] {};

void drawAcceptCorrect(
    const char* inname = "gamtree_13466_acceptcorrect.root"
    )
{

    GetAccHistos(inname);

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
                        leg = new TLegend(0.68,0.64,.96,.87);
                    }
                    if (iptt == 1)
                    {
                        leg = new TLegend(0.52,0.64,.80,.87);
                    }
                    
                    TMultiGraph* multi = new TMultiGraph();
                    
                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {
                        TGraphErrors* tgr = tgracc[iadj][itype][iptt][icent][idpsi];

                        gStyle->SetPalette(kBird);
                        tgr->SetMarkerColor(TColor::GetColorPalette(254*icent/(GAT::ncent+0)));
                        tgr->SetLineColor(TColor::GetColorPalette(254*icent/(GAT::ncent+0)));

                        tgr->SetMarkerSize(3.0);
                        tgr->SetLineWidth(1);
                        
                        tgr->SetMarkerStyle(markstyle[icent]);

                        multi->Add(tgr,"LP");

                        if ((iptt == 0 && icent < 4) || (iptt == 1 && icent >= 4))
                        {
                            ostringstream oss;
                            oss<<GAT::limcent[icent]<<"-"<<GAT::limcent[icent+1]<<"%"<<endl;
                            leg->AddEntry(tgr,oss.str().c_str(),"lp");
                        }
                    }

                    multi->SetTitle(";#xi;#epsilon");
                    multi->Draw("A");
                    double ymin = multi->GetHistogram()->GetYaxis()->GetXmin();
                    double ymax = multi->GetHistogram()->GetYaxis()->GetXmax();
                    if (itype == GAT::PI0)
                        multi->GetHistogram()->GetYaxis()->SetRangeUser(0.055,0.195);
                    else
                        multi->GetHistogram()->GetYaxis()->SetRangeUser(0.07,0.25);

                    multi->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
                    multi->GetHistogram()->GetYaxis()->SetLabelSize(0.05);

                    multi->GetHistogram()->GetXaxis()->SetLimits(0.25,2.15);
                    multi->GetHistogram()->GetXaxis()->SetNdivisions(210);
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

                    if (itype == GAT::PI0)
                        latex.DrawLatex(.425,.175,Form("#pi^{0} p_{T} #in [%.1f,%.1f)",GAT::limptt[iptt],GAT::limptt[iptt+1]));
                    else
                        latex.DrawLatex(.425,.23,Form("#gamma_{dir} p_{T} #in [%.1f,%.1f)",GAT::limptt[iptt],GAT::limptt[iptt+1]));

                        
                }

                can->SaveAs(Form("pics/acc_adj%i_%s_dpsi%i.pdf",iadj,GAT::spec[itype],idpsi));
            }
        }
    }
                    

    cout<<"Finished drawing acceptance*efficiency correction."<<endl;
    return;
}

int GetAccHistos(const char* name)
{
    TFile* input = new TFile(name,"READ");
    
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
                        vector<double> vxi;
                        vector<double> vacc;
                        vector<double> verr;
                        
                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {
                            h1acc[iadj][itype][iptt][icent][idpsi][ixi] = (TH1F*) gROOT->FindObject(Form("acc_adj%i_%s_ptt%i_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],iptt,icent,idpsi,ixi));

                            TH1F* h1accptr = h1acc[iadj][itype][iptt][icent][idpsi][ixi];

                            double nbin = h1accptr->GetNbinsX();

                            double error = -9999.;
                            double integral = h1accptr->IntegralAndError(1,nbin,error);

                            double avgacc = integral/nbin;
                            error = error/nbin;

                            vacc.push_back(avgacc);
                            verr.push_back(error);

                            vxi.push_back((GAT::limxi[ixi+1]+GAT::limxi[ixi])/2.);
                        }

                        tgracc[iadj][itype][iptt][icent][idpsi] = new TGraphErrors(GAT::nxi,&vxi[0],&vacc[0],0,&verr[0]);
                        tgracc[iadj][itype][iptt][icent][idpsi]->SetName(Form("tgracc_adj%i_%s_ptt%i_cent%i_dpsi%i",iadj,GAT::spec[itype],iptt,icent,idpsi));
                    }
                }
            }
        }
    }

    input->Close();
    delete input;
    
    return 0;
}
