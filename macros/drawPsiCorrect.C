//draw psi dependent correction to acceptance (from acceptCorrect.C)

const float PI = TMath::Pi();
const int nadj = 2;
const int nmix = 2;
const int ntype = 3;

const int markstyle[] = {20,24,21,25,22,26,23,32,34,47};

int GetFunctions(const char* name);

namespace GAT {
    #include "gamAnaTree.h"
}

TF1* fcor[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi][GAT::nxi] = {};
TGraphErrors* tgrcor[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1] = {};

void drawPsiCorrect(
    const char* inname = "gamtree_13466_acceptcorrect.root"
    )
{

    GetFunctions(inname);

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
                        leg = new TLegend(0.69,0.64,.96,.87);
                    }
                    if (iptt == 1)
                    {
                        leg = new TLegend(0.53,0.64,.80,.87);
                    }

                    TMultiGraph* multi = new TMultiGraph();

                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {
                        TGraphErrors* tgr = tgrcor[iadj][itype][iptt][icent][idpsi];

                        gStyle->SetPalette(kBird);
                        tgr->SetMarkerColor(TColor::GetColorPalette(254*icent/(GAT::ncent-1)));
                        tgr->SetLineColor(TColor::GetColorPalette(254*icent/(GAT::ncent-1)));
                        tgr->SetMarkerSize(3.0);
                        tgr->SetLineWidth(3);
                        tgr->SetMarkerStyle(markstyle[icent]);

                        multi->Add(tgr,"LP");

                        if ((iptt == 0 && icent < 4) || (iptt == 1 && icent >= 4))
                        {
                            ostringstream oss;
                            oss<<GAT::limcent[icent]<<"-"<<GAT::limcent[icent+1]<<"%"<<endl;
                            leg->AddEntry(tgr,oss.str().c_str(),"lp");
                        }
                    }

                    if (idpsi == 0)
                        multi->SetTitle(";#xi;#epsilon^{IN}_{C}(#Delta#phi=#pi)");
                    else if (idpsi == 1)
                        multi->SetTitle(";#xi;#epsilon^{OUT}_{C}(#Delta#phi=#pi)");
                    else
                        multi->SetTitle(";#xi;#epsilon^{OUT}_{C}(#Delta#phi=#pi) - #epsilon^{IN}_{C}(#Delta#phi=#pi)");

                    multi->Draw("A");
                    if (idpsi == GAT::ndpsi)
                        multi->GetHistogram()->GetYaxis()->SetRangeUser(0.0,0.08);
                    else
                        multi->GetHistogram()->GetYaxis()->SetRangeUser(.95,1.05);
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

                can->SaveAs(Form("pics/psicor_adj%i_%s_dpsi%i.png",iadj,GAT::spec[itype],idpsi));
            }
        }
    }


    cout<<"Finished drawing psi correction."<<endl;
    return;
}

int GetFunctions(const char* name)
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
                        vector<double> vcor;
                        vector<double> verr;

                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {
                            double cor = 0;
                            double err = 0;

                            if (idpsi < GAT::ndpsi)
                            {
                                fcor[iadj][itype][iptt][icent][idpsi][ixi] = (TF1*) gROOT->FindObject(Form("fpsicor_adj%i_%s_ptt%i_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],iptt,icent,idpsi,ixi));

                                cor = fcor[iadj][itype][iptt][icent][idpsi][ixi]->Eval(1);
                            }
                            else
                            {
                                cor = fcor[iadj][itype][iptt][icent][1][ixi]->Eval(1);
                                cor -= fcor[iadj][itype][iptt][icent][0][ixi]->Eval(1);
                            }

                            vcor.push_back(cor);
                            verr.push_back(err);

                            vxi.push_back(ixi);
                        }

                        tgrcor[iadj][itype][iptt][icent][idpsi] = new TGraphErrors(GAT::nxi,&vxi[0],&vcor[0],0,&verr[0]);
                        tgrcor[iadj][itype][iptt][icent][idpsi]->SetName(Form("tgrpsicor_adj%i_%s_ptt%i_cent%i_dpsi%i",iadj,GAT::spec[itype],iptt,icent,idpsi));
                    }
                }
            }
        }
    }

    input->Close();
    delete input;

    return 0;
}
