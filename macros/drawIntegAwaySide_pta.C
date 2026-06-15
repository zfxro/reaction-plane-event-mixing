//draw results from integAwaySide.C

enum RPO {IN,OUT,ALL};
string strplane[3] = {"in-plane","out-of-plane","all-plane"};

const float PI = TMath::Pi();
const int nadj = 2;
const int ntype = 3;

const int ncent = 8;
const int nsumcent = 2;

const int limsumcent[nsumcent+1] = {0,20,60};

// // const int markstyle[8] = {20,21,22,23,29,33,34,47};
// const int markstyle[] = {20,24,21,25,22,26,23,32,34,47};
// // const int markstyle[] = {20,24,21,25,22,26,23,32,34,28,47,46};

const int markstyle[] = {20,21,33,47,34};
const int markcol[] = {kBlue,kRed,kGreen,kMagenta,kCyan,kYellow};


namespace GAT {
#include "../src/gamAnaTree/gamAnaTree.h"
}

TGraphErrors* tgriptayield[nadj][ntype][GAT::nptt][nsumcent][GAT::ndpsi+1] = {};
TGraphErrors* tgripttyield[nadj][ntype][GAT::npta][nsumcent][GAT::ndpsi+1] = {};
TGraphErrors* tgriptacentrat[nadj][ntype][GAT::nptt] = {};
//v2jet vs pta
TGraphErrors* tgriptav2j[nadj][ntype][GAT::nptt][nsumcent] = {};
TGraphErrors* tgricentv2j[nadj][ntype][GAT::nptt][GAT::npta] = {};


void drawIntegAwaySide_pta(
    const char* inname = "integawaysideyield_pta.root"
    )
{

    TFile* input = new TFile(inname,"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Error opening "<<inname<<endl;
        return;
    }

    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype != GAT::PI0) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                input->GetObject(Form("iptacentrat_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt),tgriptacentrat[iadj][itype][iptt]);
                
                for (int ipta = 0; ipta < GAT::npta; ipta++)
                {
                    input->GetObject(Form("icentv2j_adj%i_%s_ptt%i_pta%i",iadj,GAT::spec[itype],iptt,ipta),tgricentv2j[iadj][itype][iptt][ipta]);
                }
                for (int icent = 0; icent < nsumcent; icent++)
                {
                    input->GetObject(Form("iptav2j_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent),tgriptav2j[iadj][itype][iptt][icent]);
                    
                    for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                    {
                        input->GetObject(Form("iptayield_adj%i_%s_ptt%i_cent%i_dpsi%i",iadj,GAT::spec[itype],iptt,icent,idpsi),tgriptayield[iadj][itype][iptt][icent][idpsi]);
                        if (iptt==0)
                        {
                            for (int ipta = 0; ipta < GAT::npta; ipta++)
                            {
                                input->GetObject(Form("ipttyield_adj%i_%s_pta%i_cent%i_dpsi%i",iadj,GAT::spec[itype],ipta,icent,idpsi),tgripttyield[iadj][itype][ipta][icent][idpsi]);
                                //normalize by binwidth
                                // tgriptayield[iadj][itype][iptt][icent][idpsi]->GetY()[ipta] /= GAT::limpta[ipta+1]-GAT::limpta[ipta];
                                // tgriptayield[iadj][itype][iptt][icent][idpsi]->SetTitle(";p_{Ta};1/N_{#gamma} dN/dp_{Ta}");
                            
                            }
                        }
                    }
                }
            }
        }
    }

    //set palette
    // gStyle->SetPalette(kBird);
    
    //function of pta
    //same plot for pt
    // TCanvas* cancentrat = new TCanvas("cancentrat","",1000,1000);
    gStyle->SetEndErrorSize(10);
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype != GAT::PI0) continue;

            TCanvas cancentrat ("cancentrat","",1000,1000);         
            TMultiGraph multi;
            // TLegend leg (0.15,0.7,.45,.875);
            TLegend leg (0.3,0.775,0.9,0.9);
            leg.SetTextSize(.035);
            leg.SetTextAlign(22);
            leg.SetNColumns(2);

            for (int iptt = GAT::nptt-1; iptt >= 0; iptt--)
            {
                TGraphErrors* tgr = tgriptacentrat[iadj][itype][iptt];
                
                // tgr->SetMarkerColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));
                // tgr->SetLineColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));

                tgr->SetMarkerColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                tgr->SetLineColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));

                // tgr->SetMarkerColor(markcol[iptt]);
                // tgr->SetLineColor(markcol[iptt]);
                
                tgr->SetMarkerSize(2.);
                tgr->SetLineWidth(3);
                tgr->SetMarkerStyle(markstyle[iptt]);

                multi.Add(tgr,"LP");
                // multi.Add(tgr,"P");

                multi.SetTitle(tgr->GetTitle());
            }

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                TGraphErrors* tgr = tgriptacentrat[iadj][itype][iptt];
                ostringstream oss;
                leg.AddEntry(tgr,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]),"lp");
            }
                
            multi.Draw("A");
            multi.GetHistogram()->GetYaxis()->SetRangeUser(0.25,1.75);
            multi.GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
            multi.GetHistogram()->GetXaxis()->SetLimits(0.,5.0);
            multi.GetHistogram()->GetXaxis()->SetNdivisions(-5);
            leg.Draw();
            double xmin = multi.GetHistogram()->GetXaxis()->GetXmin();
            double xmax = multi.GetHistogram()->GetXaxis()->GetXmax();
            TLine line (xmin,1.0,xmax,1.0);
            line.SetLineStyle(7);
            line.Draw();
            cancentrat.SaveAs(Form("integpics/centrat_adj%i_%s.png",iadj,GAT::spec[itype]));
        }
    }

    //Draw v2 
    //function of pta
    //new plot for each cent
    //same plot for pt
    gStyle->SetEndErrorSize(10);
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype != GAT::PI0) continue;

            TCanvas canptav2 ("canptav2","",1000,1000);

            for (int icent = 0; icent < nsumcent; icent++)
            {

                TMultiGraph multi;
                TLegend leg (0.3,0.775,0.9,0.9);
                leg.SetTextSize(.035);
                leg.SetTextAlign(22);
                leg.SetNColumns(2);
            
                for (int iptt = GAT::nptt-1; iptt >= 0; iptt--)
                {
                    TGraphErrors* tgr = tgriptav2j[iadj][itype][iptt][icent];
                
                    // tgr->SetMarkerColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));
                    // tgr->SetLineColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));

                    tgr->SetMarkerColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                    tgr->SetLineColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));

                    // tgr->SetMarkerColor(markcol[iptt]);
                    // tgr->SetLineColor(markcol[iptt]);
                
                    tgr->SetMarkerSize(2.);
                    tgr->SetLineWidth(3);
                    tgr->SetMarkerStyle(markstyle[iptt]);

                    multi.Add(tgr,"LP");
                    // multi.Add(tgr,"P");

                    multi.SetTitle(tgr->GetTitle());
                }

                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    TGraphErrors* tgr = tgriptav2j[iadj][itype][iptt][icent];
                    ostringstream oss;
                    leg.AddEntry(tgr,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]),"lp");
                }

            
                multi.Draw("A");
                // if (icent == 0)
                //     multi.GetHistogram()->GetYaxis()->SetRangeUser(-1.8,2.9);
                // else if (icent == 1)
                //     multi.GetHistogram()->GetYaxis()->SetRangeUser(-1.3,1.1);
                //for actual v2
                if (icent == 0)
                    multi.GetHistogram()->GetYaxis()->SetRangeUser(-1.8/((4.*sqrt(3))/PI),2.9/((4.*sqrt(3))/PI));
                else if (icent == 1)
                    multi.GetHistogram()->GetYaxis()->SetRangeUser(-1.3/((4.*sqrt(3))/PI),1.1/((4.*sqrt(3))/PI));
                
                multi.GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                multi.GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi.GetHistogram()->GetXaxis()->SetLimits(0.,5.0);
                multi.GetHistogram()->GetXaxis()->SetNdivisions(-5);
                leg.Draw();
                double xmin = multi.GetHistogram()->GetXaxis()->GetXmin();
                double xmax = multi.GetHistogram()->GetXaxis()->GetXmax();
                TLine line (xmin,0.0,xmax,0.0);
                line.SetLineStyle(7);
                line.Draw();
                canptav2.Modified();
                canptav2.SaveAs(Form("integpics/ptav2_adj%i_%s_cent%i.png",iadj,GAT::spec[itype],icent));
            }
        }
    }

    //Draw integrated away side yield
    //function of pta
    //A) all-plane. new plot for each pt. same plot for each cent.
    //B) all-plane. new plot for each cent. same plot for each pt
    //C) in-out-all plane cases same plot. new plot for each pt, cent
    gStyle->SetEndErrorSize(10);
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype != GAT::PI0) continue;
            
            // vector<TMultiGraph> mgra (GAT::nptt,TMultiGraph());
            // vector<TMultiGraph> mgrb (GAT::ncent,TMultiGraph());
            // vector<TMultiGraph> mgrc (GAT::nptt*GAT::ncent,TMultiGraph());

            // vector<TLegend> lega (GAT::nptt,TLegend(0.15,0.7,.45,.875));
            // vector<TLegend> legb (GAT::ncent,TLegend(0.15,0.7,.45,.875));
            // vector<TLegend> legc (GAT::nptt*GAT::ncent,TLegend(0.15,0.7,.45,.875));

            //B) all-plane. new plot for each cent. same plot for each pt
            for (int icent = 0; icent < nsumcent; icent++)
            {
                TCanvas can ("can","",1000,1000);
                can.SetLogy();
                TMultiGraph* mgr = new TMultiGraph();
                // TLegend leg (.15,.7,.45,.875);
                // TLegend leg (.55,.15,.85,.325);
                TLegend leg (0.3,0.775,0.9,0.9);
                // TLegend leg (0.3,0.8,0.9,0.9);
                leg.SetTextSize(.035);
                leg.SetTextAlign(22);
                leg.SetNColumns(2);
                for (int iptt = GAT::nptt-1; iptt >= 0; iptt--)         
                {
                    TGraphErrors* tgr = tgriptayield[iadj][itype][iptt][icent][ALL];

                    tgr->SetMarkerColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                    tgr->SetLineColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                
                    tgr->SetMarkerSize(2.);
                    tgr->SetLineWidth(3);
                    tgr->SetMarkerStyle(markstyle[iptt]);

                    mgr->Add(tgr,"LP");

                    mgr->SetTitle(tgr->GetTitle());
                }

                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    TGraphErrors* tgr = tgriptayield[iadj][itype][iptt][icent][ALL];
                    ostringstream oss;
                    leg.AddEntry(tgr,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]),"lp");
                }
            
                mgr->GetHistogram()->GetYaxis()->SetRangeUser(1.e-3,5.);
                mgr->Draw("A");
                
                // if (icent == 0)
                //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.4,2.2);
                // else if (icent == 1)
                //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.0,1.0);
                mgr->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                mgr->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                mgr->GetHistogram()->GetXaxis()->SetLimits(0.,5.0);
                mgr->GetHistogram()->GetXaxis()->SetNdivisions(-5);
                
                leg.Draw();

                double xmin = mgr->GetHistogram()->GetXaxis()->GetXmin();
                double xmax = mgr->GetHistogram()->GetXaxis()->GetXmax();
                TLine* line = new TLine(xmin,0.0,xmax,0.0);
                line->SetLineStyle(7);
                line->Draw();

                can.Modified();
                can.SaveAs(Form("integpics/y_adj%i_%s_cent%i.png",iadj,GAT::spec[itype],icent));
            }

            //B+) all-plane. new plot for each cent. same plot for each pta
            for (int icent = 0; icent < nsumcent; icent++)
            {
                TCanvas can ("can","",1000,1000);
                can.SetLogy();
                TMultiGraph* mgr = new TMultiGraph();
                // TLegend leg (.15,.7,.45,.875);
                // TLegend leg (.55,.15,.85,.325);
                TLegend leg (0.3,0.775,0.9,0.9);
                // TLegend leg (0.3,0.8,0.9,0.9);
                leg.SetTextSize(.035);
                leg.SetTextAlign(22);
                leg.SetNColumns(2);
                for (int ipta = GAT::npta-1; ipta >= 0; ipta--)         
                {
                    TGraphErrors* tgr = tgripttyield[iadj][itype][ipta][icent][ALL];

                    tgr->SetMarkerColor(TColor::GetColorPalette(220*ipta/(GAT::npta-1)));
                    tgr->SetLineColor(TColor::GetColorPalette(220*ipta/(GAT::npta-1)));
                
                    tgr->SetMarkerSize(2.);
                    tgr->SetLineWidth(3);
                    tgr->SetMarkerStyle(markstyle[ipta]);

                    mgr->Add(tgr,"LP");

                    mgr->SetTitle(tgr->GetTitle());
                }

                for (int ipta = 0; ipta < GAT::npta; ipta++)
                {
                    TGraphErrors* tgr = tgripttyield[iadj][itype][ipta][icent][ALL];
                    ostringstream oss;
                    leg.AddEntry(tgr,Form("%.1f #leq p_{Ta} < %.1f",GAT::limpta[ipta],GAT::limpta[ipta+1]),"lp");
                }
            
                mgr->GetHistogram()->GetYaxis()->SetRangeUser(5.e-4,5.e1);
                mgr->Draw("A");
                
                // if (icent == 0)
                //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.4,2.2);
                // else if (icent == 1)
                //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.0,1.0);
                mgr->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                mgr->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                mgr->GetHistogram()->GetXaxis()->SetLimits(5.,15.);
                mgr->GetHistogram()->GetXaxis()->SetNdivisions(-5);
                
                leg.Draw();

                double xmin = mgr->GetHistogram()->GetXaxis()->GetXmin();
                double xmax = mgr->GetHistogram()->GetXaxis()->GetXmax();
                TLine* line = new TLine(xmin,0.0,xmax,0.0);
                line->SetLineStyle(7);
                line->Draw();

                can.Modified();
                can.SaveAs(Form("integpics/y_adj%i_%s_cent%i_pta.png",iadj,GAT::spec[itype],icent));
            }

            //A) all-plane. new plot for each pt. same plot for each cent.
            for (int iptt = 0; iptt < GAT::nptt; iptt++)                    
            {
                TCanvas can ("can","",1000,1000);
                can.SetLogy();
                TMultiGraph* mgr = new TMultiGraph();
                TLegend leg (.15,.78,.35,.875);

                for (int icent = 0; icent < nsumcent; icent++)
                {
                    TGraphErrors* tgr = tgripttyield[iadj][itype][iptt][icent][ALL];

                    tgr->SetMarkerColor(TColor::GetColorPalette(220*icent/(GAT::nptt-1)));
                    tgr->SetLineColor(TColor::GetColorPalette(220*icent/(GAT::nptt-1)));
                
                    tgr->SetMarkerSize(2.);
                    tgr->SetLineWidth(3);
                    tgr->SetMarkerStyle(markstyle[icent]);

                    mgr->Add(tgr,"LP");

                    ostringstream oss;
                    leg.AddEntry(tgr,Form("%i-%i%%",limsumcent[icent],limsumcent[icent+1]),"lp");

                    mgr->SetTitle(tgr->GetTitle());
                }

            
                mgr->Draw("A");
                // if (icent == 0)
                //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.4,2.2);
                // else if (icent == 1)
                //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.0,1.0);
                mgr->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                mgr->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                mgr->GetHistogram()->GetXaxis()->SetLimits(0.,5.0);
                mgr->GetHistogram()->GetXaxis()->SetNdivisions(-5);
                
                leg.Draw();

                double xmin = mgr->GetHistogram()->GetXaxis()->GetXmin();
                double xmax = mgr->GetHistogram()->GetXaxis()->GetXmax();
                TLine* line = new TLine(xmin,0.0,xmax,0.0);
                line->SetLineStyle(7);
                line->Draw();

                can.Modified();
                can.SaveAs(Form("integpics/y_adj%i_%s_ptt%i.png",iadj,GAT::spec[itype],iptt));
            }

            //C) in-out-all plane cases same plot. new plot for each pt, cent
            for (int iptt = 0; iptt < GAT::nptt; iptt++)                    
            {
                for (int icent = 0; icent < nsumcent; icent++)
                {
                    TCanvas can ("can","",1000,1000);
                    can.SetLogy();
                    TMultiGraph* mgr = new TMultiGraph();
                    TLegend leg (.15,.7875,.40,.875);

                    for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                    {
                        if (idpsi == GAT::ndpsi) continue;//skip all plane case
                        TGraphErrors* tgr = tgriptayield[iadj][itype][iptt][icent][idpsi];
                

                        tgr->SetMarkerColor(TColor::GetColorPalette(220*idpsi/(GAT::nptt-1)));
                        tgr->SetLineColor(TColor::GetColorPalette(220*idpsi/(GAT::nptt-1)));
                
                        tgr->SetMarkerSize(2.);
                        tgr->SetLineWidth(3);
                        tgr->SetMarkerStyle(markstyle[idpsi]);

                        mgr->Add(tgr,"LP");

                        ostringstream oss;
                        leg.AddEntry(tgr,Form("%s",strplane[idpsi].c_str()),"lp");

                        mgr->SetTitle(tgr->GetTitle());
                    }

            
                    mgr->Draw("A");
                    // if (icent == 0)
                    //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.4,2.2);
                    // else if (icent == 1)
                    //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.0,1.0);
                    mgr->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                    mgr->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                    mgr->GetHistogram()->GetXaxis()->SetLimits(0.,5.0);
                    mgr->GetHistogram()->GetXaxis()->SetNdivisions(-5);
                
                    leg.Draw();

                    double xmin = mgr->GetHistogram()->GetXaxis()->GetXmin();
                    double xmax = mgr->GetHistogram()->GetXaxis()->GetXmax();
                    TLine* line = new TLine(xmin,0.0,xmax,0.0);
                    line->SetLineStyle(7);
                    line->Draw();

                    can.Modified();
                    can.SaveAs(Form("integpics/y_adj%i_%s_ptt%i_cent%i.png",iadj,GAT::spec[itype],iptt,icent));

                }
            }

        }
    }

    
}
