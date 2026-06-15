//draw dphi plots resulting from subtractWeightedBG.C

const int nbinppg = 10;
double vy0[nbinppg] = {-0.0128,-0.0363,-0.0437,-0.0479,0.0044,0.0031,0.1130,0.0607,0.0539,0.0148};
double vey0[nbinppg] = {0.0173,0.0197,0.0234,0.0284,0.0413,0.0407,0.0296,0.0222,0.0204,0.0197};
double vy1[nbinppg] = {0.0109,0.0175,0.0120,0.0022,0.0267,0.0258,0.0683,0.0469,0.0519,0.0464};
double vey1[nbinppg] = {0.0079,0.0087,0.0096,0.0137,0.0191,0.0178,0.0131,0.0101,0.0090,0.0085};
double vy2[nbinppg] = {-0.0095,0.0178,0.0170,0.0195,0.0324,0.0263,0.0153,0.0203,0.0246,0.0324};
double vey2[nbinppg] = {0.0051,0.0042,0.0079,0.0088,0.0121,0.0130,0.0084,0.0050,0.0054,0.0038};
double vy3[nbinppg] = {0.00419,0.00529,0.00392,-0.00079,0.00254,0.00473,0.01592,0.00833,0.0169 ,0.0102};
double vey3[nbinppg] = {0.00334,0.00319,0.00349,0.00436,0.00726,0.00727,0.00494,0.00378,0.00319,0.00349};
double vy4[nbinppg] = {0.00335,0.00220,0.00068,-0.00067,0.00162,-0.00025,0.00287,0.00277,0.00278,0.00508};
double vey4[nbinppg] = {0.00077,0.00096,0.00096,0.00211,0.00211,0.00268,0.00191,0.00077,0.00057,0.00153};
double vy5[nbinppg] = {0.001874,0.000088,0.000248,0.000649,0.000235,-0.001043,0.000087,0.000696,0.000868,0.001156};
double vey5[nbinppg] = {0.000614,0.000818,0.000613,0.000409,0.000921,0.000922,0.000819,0.000615,0.000612,0.000616};


const int nadj = 2;
const int ncent = 2;
const int limcent[ncent+1] = {0,20,60};

const double efferr = .076;//from pgg023 (our cutoff of pt < 3? used in efficiency calculation)

enum RPO {IN,OUT,ALL};
string strplane[3] = {"in-plane","out-of-plane","all-plane"};

const int markstyle[] = {20,21,33,47,34};

const float PI = TMath::Pi();

namespace GAT {
    #include "gamAnaTree.h"
}

vector<double> vx;
vector<double> vyxim[GAT::nxi];
vector<double> veyxim[GAT::nxi];

TGraph* gsmearpct[nadj][GAT::ntype+1][GAT::nptt] = {};

struct errorinfo_t {
    Int_t iadj, itype, iptt, ixi, icent, idpsi;
    Double_t ntrg;
    Double_t flate[2] = {};//flat term error
    Double_t cose[2] = {};//cosine term error
};
    
errorinfo_t errinfo[nadj][GAT::ntype+1][GAT::nptt][ncent+1][GAT::ndpsi+1][GAT::nxi] = {};

void FillSysError(TFile* input);

double GetError(TGraphErrors* tgr, TF1* f1, double xmin);
void DrawError(TF1* f1, double err, double xmin, int kColor = kRed, double alpha = .10, int ipar = 0);

void drawDphiPlots(
    const char* inname = "weighted_bg_subtracted.root",
    const char* inname_k = "kComponents/kComponents.root",
    const char* inname_smear = "smeared.root",
    const char* inname_sig = "kComponents/sigmak.root",
    const char* outname = "dphiplotsdrawn.root"    
    )
{

    //Merge ppg210 data
    {
        for (int ibin = 0; ibin < nbinppg; ibin++)
        {
            vyxim[0].push_back((vy4[ibin]+vy5[ibin])/2.);
            veyxim[0].push_back(sqrt(vey4[ibin]*vey4[ibin]+vey5[ibin]*vey5[ibin]));

            vyxim[1].push_back((vy2[ibin]+vy3[ibin])/2.);
            veyxim[1].push_back(sqrt(vey2[ibin]*vey2[ibin]+vey3[ibin]*vey3[ibin]));

            vyxim[2].push_back((vy0[ibin]+vy1[ibin])/2.);
            veyxim[2].push_back(sqrt(vey0[ibin]*vey0[ibin]+vey1[ibin]*vey1[ibin]));

            vx.push_back(.05+.1*ibin);
        }
    }
    
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0.);
    gStyle->SetEndErrorSize(20);
    gStyle->SetMarkerSize(3);
    gStyle->SetHistLineWidth(1);
    gStyle->SetTitleOffset(1.4,"Y");
    gStyle->SetImageScaling(3.);
    gROOT->ForceStyle();

    TFile* input_smear = new TFile(inname_smear,"READ");
    cout<<"Reading file: "<<inname<<endl;
    if (!input_smear || input_smear->IsZombie())
    {
        cerr<<"Error opening "<<inname<<endl;
        return;
    }

    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                input_smear->GetObject(Form("pct_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt),gsmearpct[iadj][itype][iptt]);
            }
        }
    }

    //reaction plane resolution
    TGraphErrors* tgrexpsig = nullptr;
    TGraphErrors* tgrksig[nadj][GAT::ntype+1][GAT::nptt] = {};

    TFile* input_sig = new TFile(inname_sig,"READ");
    cout<<"Reading file: "<<inname<<endl;
    if (!input_sig || input_sig->IsZombie())
    {
        cerr<<"Error opening "<<inname<<endl;
        return;
    }

    input_sig->GetObject("expsig",tgrexpsig);

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                input_sig->GetObject(Form("ksigma_Np_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt),tgrksig[iadj][itype][iptt]);
            }
        }
    }


    TFile* input = new TFile(inname,"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Error opening "<<inname<<endl;
        return;
    }

    FillSysError(input);

    TH1F* h1dphi[nadj][GAT::ntype+1][ncent+1][GAT::ndpsi+1][GAT::nptt][GAT::nxi] {};
    TGraphErrors* tgrdphi[nadj][GAT::ntype+1][ncent+1][GAT::ndpsi+1][GAT::nptt][GAT::nxi] {};
    TGraphErrors* tgrdphi_sys[nadj][GAT::ntype+1][ncent+1][GAT::ndpsi+1][GAT::nptt][GAT::nxi] {};

    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int icent = 0; icent < ncent+1; icent++)
            {
                for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                {
                    for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    {
                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {
                            ostringstream oss;
                            oss<<"dphi";
                            if (iadj)
                                oss<<"_adj";
                            oss<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<icent<<"_dpsi"<<idpsi<<"_xi"<<ixi;

                            input->GetObject(oss.str().c_str(),h1dphi[iadj][itype][icent][idpsi][iptt][ixi]);
                            h1dphi[iadj][itype][icent][idpsi][iptt][ixi]->Scale(1./2.);//dphi folding factor required

                            tgrdphi[iadj][itype][icent][idpsi][iptt][ixi] = new TGraphErrors(h1dphi[iadj][itype][icent][idpsi][iptt][ixi]);

                            // //remove points
                            // {
                            //  {
                            //      tgrdphi[iadj][itype][icent][idpsi][iptt][ixi]->RemovePoint(ipt);
                            //  }
                            // }


                            // {
                            //  {
                            //      tgrdphi[iadj][itype][icent][idpsi][iptt][ixi]->SetPoint(j,x,y);
                            //      tgrdphi[iadj][itype][icent][idpsi][iptt][ixi]->SetPointError(j,width * ex,ey);
                            //      tgrdphi_sys[iadj][itype][icent][idpsi][iptt][ixi]->SetPoint(j,x,y);
                            //      tgrdphi_sys[iadj][itype][icent][idpsi][iptt][ixi]->SetPointError(j,width * ex,ey);    
                            //      j++;
                            //  }
                            // }

                            tgrdphi_sys[iadj][itype][icent][idpsi][iptt][ixi] = new TGraphErrors(h1dphi[iadj][itype][icent][idpsi][iptt][ixi]);

                            //i thought the way in which the tgraphs were being made was messing up the error widths.
                            //turns out the error width mess was arising below during setting of error width. go figure.
                            //anyway, i changed the way tgraphs were intialized many times. keeping this option tho others could work too.
                            {
                                int j = 0;//tgraph point
                                TAxis *xaxis = h1dphi[iadj][itype][icent][idpsi][iptt][ixi]->GetXaxis();
                                for (Int_t i = 0; i < xaxis->GetNbins(); i++)
                                {
                                    double y = h1dphi[iadj][itype][icent][idpsi][iptt][ixi]->GetBinContent(i + 1);
                                    if (y == 0.)
                                    {
                                        tgrdphi[iadj][itype][icent][idpsi][iptt][ixi]->RemovePoint(j);
                                        tgrdphi_sys[iadj][itype][icent][idpsi][iptt][ixi]->RemovePoint(j);
                                    }
                                    else
                                    {
                                        j++;
                                    }
                                }
                            }

                            // {
                            //  tgrdphi[iadj][itype][icent][idpsi][iptt][ixi]->SetPointError(i,.02,tgrdphi[iadj][itype][icent][idpsi][iptt][ixi]->GetErrorY(i));
                            //  tgrdphi_sys[iadj][itype][icent][idpsi][iptt][ixi]->SetPointError(i,.02,tgrdphi_sys[iadj][itype][icent][idpsi][iptt][ixi]->GetErrorY(i));
                            // }

                            // ---------- sys error from bg ----------
                            // idpsi was erroneously set to IN!
                            double fe0 = errinfo[iadj][itype][iptt][icent][idpsi][ixi].flate[0];
                            double fe1 = errinfo[iadj][itype][iptt][icent][idpsi][ixi].flate[1];
                            double ce0 = errinfo[iadj][itype][iptt][icent][idpsi][ixi].cose[0];
                            double ce1 = errinfo[iadj][itype][iptt][icent][idpsi][ixi].cose[1];

                            double fe = sqrt(fe0*fe0 + fe1*fe1);
                            double ce = sqrt(ce0*ce0 + ce1*ce1);

                            TF1* fbge = new TF1(Form("%s_bge",h1dphi[iadj][itype][icent][idpsi][iptt][ixi]->GetName()),"sqrt(pow([0],2.)+pow([1]*cos(x*2.*TMath::Pi()),2.))",0.,1.);
                            fbge->SetParameter(0,fe);
                            fbge->SetParameter(1,ce);

                            double ew = (1.-h1dphi[iadj][itype][icent][idpsi][iptt][ixi]->GetBinCenter(GAT::ndphibin))/3.;
                            
                            //use data from tgrdphi. fill tgraph sys error
                            //tgrdphi_sys[iadj][itype][icent][idpsi][iptt][ixi]

                            // {

                            //  // sys err from bg


                            //  tgrdphi_sys[iadj][itype][icent][idpsi][iptt][ixi]->SetPointError(idp,ew,tote);
                            // }

                            for (int idp = 0; idp < tgrdphi[iadj][itype][icent][idpsi][iptt][ixi]->GetN(); idp++)
                            {
                                double y = tgrdphi[iadj][itype][icent][idpsi][iptt][ixi]->GetY()[idp];
                                double x = tgrdphi[iadj][itype][icent][idpsi][iptt][ixi]->GetX()[idp];

                                // sys err from bg
                                double bge = fbge->Eval(x);

                                double effe = efferr*y;

                                double smr = gsmearpct[iadj][itype][iptt]->GetY()[ixi]/100.*y;

                                double tote = sqrt(bge*bge + effe*effe + smr*smr);

                                tgrdphi_sys[iadj][itype][icent][idpsi][iptt][ixi]->SetPointError(idp,ew,tote);
                            }

                        }
                    }
                }
            }
        }
    }

    //lets create the I(IN - OUT) / ALL graphs

    TH1F* h1v2j[nadj][GAT::ntype+1][ncent+1][GAT::nptt][GAT::nxi] {};
    TGraphErrors* tgrv2j[nadj][GAT::ntype+1][ncent+1][GAT::nptt][GAT::nxi] {};
    double v2jsyserr[nadj][GAT::ntype+1][ncent+1][GAT::nptt][GAT::nxi] {};

    //----------------------------------------------------------------------
    //v2jet vs xi
    TGraphErrors* tgrixiv2j[nadj][GAT::ntype+1][GAT::nptt][ncent+1] = {};
    TGraphErrors* tgrixiv2j_syserr[nadj][GAT::ntype+1][GAT::nptt][ncent+1] = {};

    //v2jet vs ptt
    TGraphErrors* tgripttv2j[nadj][GAT::ntype+1][ncent+1][GAT::nxi] = {};
    TGraphErrors* tgripttv2j_syserr[nadj][GAT::ntype+1][ncent+1][GAT::nxi] = {};

    //x positions. includes offset. we are copying how integAwaySide.C does things
    double xpos_xi[GAT::nptt][GAT::nxi];
    double xpos_ptt[GAT::nptt][GAT::nxi];

    double sex_xi = 0.;//sys error x for xi bins
    double sex_ptt = 0.;//sys error x for ptt bins

    {
        double xioff[GAT::nptt] = {-.12,-.04,.04,.12};
        double ptoff[GAT::nxi] = {-.09,.0,.09};

        double xibw = GAT::limxi[1]-GAT::limxi[0];
        double pttbw = GAT::limptt[1]-GAT::limptt[0];

        sex_xi = .035*xibw;
        sex_ptt = .04*pttbw;

        for (int iptt = 0; iptt < GAT::nptt; iptt++)
        {
            for (int ixi = 0; ixi < GAT::nxi; ixi++)
            {
                double ximid = (GAT::limxi[ixi+1]+GAT::limxi[ixi])/2.;
                double ptmid = (GAT::limptt[iptt+1]+GAT::limptt[iptt])/2.;

                xpos_xi[iptt][ixi] = ximid + xioff[iptt]*xibw;
                xpos_ptt[iptt][ixi] = ptmid + ptoff[ixi]*pttbw;
            }
        }
    }
    //----------------------------------------------------------------------

    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {           
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;//doesnt exist

            for (int icent = 0; icent < ncent+1; icent++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        //creating tgraphs to hold v2j results
                        if (!tgripttv2j[iadj][itype][icent][ixi])
                        {
                            tgripttv2j[iadj][itype][icent][ixi] = new TGraphErrors();
                            tgripttv2j[iadj][itype][icent][ixi]->SetName(Form("ipttv2jfromfit_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi));
                            tgripttv2j[iadj][itype][icent][ixi]->SetTitle(";p_{Tt};v_{2,away-side}");

                            tgripttv2j_syserr[iadj][itype][icent][ixi] = new TGraphErrors();
                            tgripttv2j_syserr[iadj][itype][icent][ixi]->SetName(Form("ipttv2jfromfit_syserr_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi));
                            tgripttv2j_syserr[iadj][itype][icent][ixi]->SetTitle(";p_{Tt};v_{2,away-side}");

                        }
                        if (!tgrixiv2j[iadj][itype][iptt][icent])
                        {
                            tgrixiv2j[iadj][itype][iptt][icent] = new TGraphErrors();
                            tgrixiv2j[iadj][itype][iptt][icent]->SetName(Form("ixiv2jfromfit_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent));
                            tgrixiv2j[iadj][itype][iptt][icent]->SetTitle(";#xi;v_{2,away-side}");

                            tgrixiv2j_syserr[iadj][itype][iptt][icent] = new TGraphErrors();
                            tgrixiv2j_syserr[iadj][itype][iptt][icent]->SetName(Form("ixiv2jfromfit_syserr_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent));               
                            tgrixiv2j_syserr[iadj][itype][iptt][icent]->SetTitle(";#xi;v_{2,away-side}");
                        }
                        
                        ostringstream oss;
                        oss<<"v2j";
                        if (iadj)
                            oss<<"_adj";
                        oss<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<icent<<"_xi"<<ixi;

                        //do the arithmetic on the histos
                        
                        h1v2j[iadj][itype][icent][iptt][ixi] = (TH1F*) h1dphi[iadj][itype][icent][IN][iptt][ixi]->Clone(oss.str().c_str());
                        h1v2j[iadj][itype][icent][iptt][ixi]->Add(h1dphi[iadj][itype][icent][OUT][iptt][ixi],-1.);

                        h1v2j[iadj][itype][icent][iptt][ixi]->Divide(h1dphi[iadj][itype][icent][ALL][iptt][ixi]);

                        //before scaling calculate the error
                        {
                            //using last bin to determine sys err since v2j is amplitude of cosine function.
                            //also approximating cos(2Deltaphi) as 1 here

                            int nbinsx = h1v2j[iadj][itype][icent][iptt][ixi]->GetNbinsX();
                          
                            double R = h1v2j[iadj][itype][icent][iptt][ixi]->GetBinContent(nbinsx);//ratio
                            double D = h1dphi[iadj][itype][icent][ALL][iptt][ixi]->GetBinContent(nbinsx);//denominator

                            double Ntin = errinfo[iadj][itype][iptt][icent][IN][ixi].ntrg;//Ntrg in plane
                            double Ntout = errinfo[iadj][itype][iptt][icent][OUT][ixi].ntrg;//Ntrg out of plane
                            double Nt = errinfo[iadj][itype][iptt][icent][ALL][ixi].ntrg;//Ntrg total

                            double dRda1 = -((1./Ntin) + (1./Ntout))/D;//dR/d#alpha_1
                            double dRdb1 = dRda1;//dR/d#beta_1; these are equal when cos2Dphi = 1

                            double dRda0 = ((1./Ntout) - (1./Ntin))/D + R/(D*Nt);//dR/d#alpha_0
                            double dRdb0 = dRda0;//dR/d#beta_0; these are equal when cos2Dphi = 1

                            double ea0 = errinfo[iadj][itype][iptt][icent][IN][ixi].flate[0] * Ntin;
                            double ea1 = errinfo[iadj][itype][iptt][icent][IN][ixi].flate[1] * Ntin;
                            double eb0 = errinfo[iadj][itype][iptt][icent][IN][ixi].cose[0] * Ntin;
                            double eb1 = errinfo[iadj][itype][iptt][icent][IN][ixi].cose[1] * Ntin;

                            v2jsyserr[iadj][itype][icent][iptt][ixi] =
                                sqrt( pow(dRda0*ea0,2.) + pow(dRda1*ea1,2.) + pow(dRdb0*eb0,2.) + pow(dRdb1*eb1,2.) );
                        }

                        double sig = tgrexpsig->GetY()[icent];
                        sig *= tgrksig[iadj][itype][iptt]->GetY()[icent];
                        
                        double scale = PI/(16.*sig);
                        h1v2j[iadj][itype][icent][iptt][ixi]->Scale(scale);
                        v2jsyserr[iadj][itype][icent][iptt][ixi] *= scale;

                        //==================== NO! DO NOT DO THIS SCALING. ====================
                        // ========== Scaling not necessary as this is a ratio! ==========
                        // h1v2j[iadj][itype][icent][iptt][ixi]->Scale(1./(h1v2j[iadj][itype][icent][iptt][ixi]->GetBinWidth(1)*PI*2.));
                        // v2jsyserr[iadj][itype][icent][iptt][ixi] /= h1v2j[iadj][itype][icent][iptt][ixi]->GetBinWidth(1)*PI*2.;

                        //create tgraph with result

                        tgrv2j[iadj][itype][icent][iptt][ixi] = new TGraphErrors(h1v2j[iadj][itype][icent][iptt][ixi]);

                        for (int ipoint = 0, ibin = 1; ibin <= h1v2j[iadj][itype][icent][iptt][ixi]->GetNbinsX(); ibin++)
                        {
                            bool setzero = false;
                            for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                            {
                                if (h1dphi[iadj][itype][icent][idpsi][iptt][ixi]->GetBinContent(ibin) == 0.)
                                    setzero = true;
                            }
                            if (setzero)
                                tgrv2j[iadj][itype][icent][iptt][ixi]->RemovePoint(ipoint);
                            else
                                ipoint++;
                        }

                    }
                }
            }
        }
    }

    //in and out same plot
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int icent = 0; icent <= ncent; icent++)
            {
                TCanvas* bigcan[2];// one for 5-7 & 7-9 and the other 9-12 & 12-15
                TPad* bigpad[2][2][3];// 2 pt columns, 3 xi rows
                TLatex* bigtxtmom[2][2][3];
                TLatex* bigtxtcent[2][2][3];
                TMultiGraph* bigmgr[2][2][3];
                double biglim[2][3][2];// x index not needed
                TLegend* bigleg[2][2][3];

                // set up canvas
                double xpos = 0.;
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {

                    int ibc = iptt/2;
                    if (iptt%2==0)
                    {
                        xpos = 0;
                        bigcan[ibc] = new TCanvas(Form("bigcan%i_inout_%s_adj%i_cent%i",ibc,GAT::spec[itype],iadj,icent),"",3000,3500);
                    }
                    int ipx = iptt%2;//pad x pos


                    double left = .02;
                    if (ipx == 0) left = .15;

                    double right = .02;
                    if (ipx == 1) right = .05;


                    double wd = .93/1.76;
                    if (ipx == 1) wd *= .83/.93;

                    double ypos = 0.;               
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        int ipy = ixi;//pad y pos

                        //set pad margins and create
                        double top = .01;
                        if (ipy == 2) top = .15;

                        double bot = .01;
                        if (ipy == 0) bot = .15;


                        double ht = .84/2.80;
                        if (ipy != 1) ht *= .98/.84;

                        bigpad[ibc][ipx][ipy] = new TPad(Form("bigpad%i%i%i",ibc,ipx,ipy),"",xpos,ypos,xpos+wd,ypos+ht);

                        bigmgr[ibc][ipx][ipy] = new TMultiGraph();//make mgr while were at it

                        // initialize y range limits (unecessarily repeated)
                        biglim[ibc][ipy][0] = numeric_limits<double>::max();
                        biglim[ibc][ipy][1] = numeric_limits<double>::min();

                        bigpad[ibc][ipx][ipy]->SetTopMargin(top);
                        bigpad[ibc][ipx][ipy]->SetBottomMargin(bot);
                        bigpad[ibc][ipx][ipy]->SetRightMargin(right);
                        bigpad[ibc][ipx][ipy]->SetLeftMargin(left);

                        // doesnt just "SetTicks" the labels appear too :(
                        // bigpad[ibc][ipx][ipy]->SetTicky(2);
                            
                        bigpad[ibc][ipx][ipy]->Draw();


                        double txtx = .5;
                        if (ipx == 0 && ipy == 2) txtx = .35;
                        
                        bigtxtmom[ibc][ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.075*(1.-top-bot),Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                        bigtxtmom[ibc][ipx][ipy]->SetTextAlign(22);
                        bigtxtmom[ibc][ipx][ipy]->SetNDC();
                        bigtxtmom[ibc][ipx][ipy]->SetTextFont(63);
                        bigtxtmom[ibc][ipx][ipy]->SetTextSize(60);

                        if (icent != ncent)
                            bigtxtcent[ibc][ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.15*(1.-top-bot),Form("Au+Au %i-%i%% 200GeV",limcent[icent],limcent[icent+1]));
                        else
                            bigtxtcent[ibc][ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.15*(1.-top-bot),Form("Au+Au %i-%i%% 200GeV",limcent[0],limcent[ncent]));
                        bigtxtcent[ibc][ipx][ipy]->SetTextAlign(22);
                        bigtxtcent[ibc][ipx][ipy]->SetNDC();
                        bigtxtcent[ibc][ipx][ipy]->SetTextFont(63);
                        bigtxtcent[ibc][ipx][ipy]->SetTextSize(60);

                        bigleg[ibc][ipx][ipy] = new TLegend(1.-right-.33*(1.-right-left),1.-top-.2*(1.-top-bot),1.-right,1.-top);
                        bigleg[ibc][ipx][ipy]->SetTextFont(43);
                        bigleg[ibc][ipx][ipy]->SetTextSize(60);

                        ypos += ht;
                    }
                    xpos += wd;
                }

                
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    int ibc = iptt/2;
                    int ipx = iptt%2;//pad x pos
                        
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        int ipy = ixi;//pad y pos
                        
                        TCanvas can ("can","",2500,2000);//2x3
                        TMultiGraph* mgr = new TMultiGraph;
                        TLegend leg (0.65,0.75,0.9,0.9);
                        leg.SetTextSize(0.04);
                        
                        for (int idpsi = 0; idpsi < GAT::ndpsi; idpsi++)
                        {
                            TGraphErrors* tgr = tgrdphi[iadj][itype][icent][idpsi][iptt][ixi];
                            tgr->SetMarkerColor(TColor::GetColorPalette(440/3*idpsi));
                            tgr->SetLineColor(TColor::GetColorPalette(440/3*idpsi));
                            tgr->SetMarkerStyle(markstyle[idpsi]);

                            TGraphErrors* tgr_sys = tgrdphi_sys[iadj][itype][icent][idpsi][iptt][ixi];
                            tgr_sys->SetFillColorAlpha(TColor::GetColorPalette(440/3*idpsi),.5);

                            mgr->Add(tgr_sys,"2");
                            mgr->Add(tgr,"p");
                            TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                            tgrclone->SetMarkerSize(0);                 
                            mgr->Add(tgrclone,"p");
                            leg.AddEntry(tgr,Form("%s",strplane[idpsi].c_str()),"p");

                            bigmgr[ibc][ipx][ipy]->Add(tgr_sys,"2");
                            bigmgr[ibc][ipx][ipy]->Add(tgr,"p");
                            bigmgr[ibc][ipx][ipy]->Add(tgrclone,"p");
                            bigleg[ibc][ipx][ipy]->AddEntry(tgr,Form("%s",strplane[idpsi].c_str()),"p");
                            
                        }

                        mgr->SetTitle(";#Delta#phi/#pi;1/N_{ t}  dN_{ h }/d#Delta#phi");
                        mgr->Draw("A");
                        leg.Draw();

                        TLatex txtmom (0.4,.85,Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                        txtmom.SetTextAlign(22);
                        txtmom.SetNDC();
                        txtmom.SetTextSize(0.04);
                        txtmom.Draw();

                        TLatex* txtcent;
                        if (icent != ncent)
                        {
                            txtcent = new TLatex(0.4,.8,Form("Au+Au %i-%i%% 200GeV",limcent[icent],limcent[icent+1]));
                        }
                        else
                        {
                            txtcent = new TLatex(0.4,.8,Form("Au+Au %i-%i%% 200GeV",limcent[0],limcent[ncent]));
                            
                        }
                        txtcent->SetTextAlign(22);
                        txtcent->SetNDC();
                        txtcent->SetTextSize(0.04);
                        txtcent->Draw();
                        
                        double ymin = min(0.,mgr->GetHistogram()->GetMinimum());
                        double ymax = max(0.,mgr->GetHistogram()->GetMaximum());
                        double ydiff = ymax-ymin;
                        double bigymin = ymin;
                        double bigymax = ymax;
                        
                        ymin = ymin - (ydiff)*.1;
                        ymax = ymax + (ydiff)*.3;

                        biglim[ibc][ipy][0] = min(bigymin,biglim[ibc][ipy][0]);
                        biglim[ibc][ipy][1] = max(bigymax,biglim[ibc][ipy][1]);
    
                        mgr->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);
                        mgr->GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                        TLine* line = new TLine(0.,0.,1.,0.);
                        line->SetLineWidth(1);
                        line->SetLineStyle(2);
                        line->Draw();

                        mgr->Draw();

                        can.Update();
                        can.SaveAs(Form("dphiplots/inout_%s_adj%i_cent%i_ptt%i_xi%i.pdf",GAT::spec[itype],iadj,icent,iptt,ixi));

                        // draw on big pads too!
                        bigpad[ibc][ipx][ipy]->cd();
                        bigmgr[ibc][ipx][ipy]->SetTitle(";#Delta#phi (1/#pi);1/N_{ t}  dN_{ h }/d#Delta#phi");                  
                        bigmgr[ibc][ipx][ipy]->Draw("A");

                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleFont(43);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelFont(43);
                        if (ipy == 0)
                        {
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleOffset(2.75);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(70);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(60);
                        }
                        else
                        {
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(0);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(0);
                        }

                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleFont(43);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelFont(43);
                        if (ipx == 0)
                        {
                            
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleOffset(3.5);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(70);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(60);
                        }
                        else
                        {
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(0);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(0);
                        }
                        
                        line->Draw();
                        bigtxtmom[ibc][ipx][ipy]->Draw();
                        // bigtxtcent[ibc][ipx][ipy]->Draw();

                        if (ipx == 0 && ipy == 2)
                        {
                            bigleg[ibc][ipx][ipy]->Draw();
                            bigtxtcent[ibc][ipx][ipy]->Draw();
                        }
                    }//ixi
                }//iptt

                for (int ibc = 0; ibc < 2; ibc++)
                {
                    for (int ipx = 0; ipx < 2; ipx++)
                    {
                        for (int ipy = 0; ipy < 3; ipy++)
                        {
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetRangeUser(biglim[ibc][ipy][0],biglim[ibc][ipy][1]);
                            bigpad[ibc][ipx][ipy]->Modified();
                            bigpad[ibc][ipx][ipy]->Update();
                        }
                    }
                    
                    bigcan[ibc]->SaveAs(Form("dphiplots/%s.pdf",bigcan[ibc]->GetName()));
                }
            }
        }
    }

    //(IN - OUT) / ALL
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int icent = 0; icent <= ncent; icent++)
            {
                TCanvas* bigcan[2];// one for 5-7 & 7-9 and the other 9-12 & 12-15
                TPad* bigpad[2][2][3];// 2 pt columns, 3 xi rows
                TLatex* bigtxtmom[2][2][3];
                TLatex* bigtxtcent[2][2][3];
                TMultiGraph* bigmgr[2][2][3];
                double biglim[2][3][2];// x index not needed
                TLegend* bigleg[2][2][3];

                // set up canvas
                double xpos = 0.;
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {

                    int ibc = iptt/2;
                    if (iptt%2==0)
                    {
                        xpos = 0;
                        bigcan[ibc] = new TCanvas(Form("bigcan%i_v2j_%s_adj%i_cent%i",ibc,GAT::spec[itype],iadj,icent),"",3000,3500);
                    }
                    int ipx = iptt%2;//pad x pos

                    double left = .02;
                    if (ipx == 0) left = .15;

                    double right = .02;
                    if (ipx == 1) right = .05;

                    double wd = .93/1.76;
                    if (ipx == 1) wd *= .83/.93;

                    double ypos = 0.;               
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        int ipy = ixi;//pad y pos

                        //set pad margins and create
                        double top = .01;
                        if (ipy == 2) top = .15;

                        double bot = .01;
                        if (ipy == 0) bot = .15;

                        double ht = .84/2.80;
                        if (ipy != 1) ht *= .98/.84;

                        bigpad[ibc][ipx][ipy] = new TPad(Form("bigpad%i%i%i",ibc,ipx,ipy),"",xpos,ypos,xpos+wd,ypos+ht);

                        bigmgr[ibc][ipx][ipy] = new TMultiGraph();//make mgr while were at it

                        // initialize y range limits (unecessarily repeated)
                        biglim[ibc][ipy][0] = numeric_limits<double>::max();
                        biglim[ibc][ipy][1] = numeric_limits<double>::min();

                        bigpad[ibc][ipx][ipy]->SetTopMargin(top);
                        bigpad[ibc][ipx][ipy]->SetBottomMargin(bot);
                        bigpad[ibc][ipx][ipy]->SetRightMargin(right);
                        bigpad[ibc][ipx][ipy]->SetLeftMargin(left);
                            
                        bigpad[ibc][ipx][ipy]->Draw();


                        double txtx = .5;
                        if (ipx == 0 && ipy == 2) txtx = .40;
                        
                        bigtxtmom[ibc][ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.075*(1.-top-bot),Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                        bigtxtmom[ibc][ipx][ipy]->SetTextAlign(22);
                        bigtxtmom[ibc][ipx][ipy]->SetNDC();
                        bigtxtmom[ibc][ipx][ipy]->SetTextFont(63);
                        bigtxtmom[ibc][ipx][ipy]->SetTextSize(60);

                        if (icent != ncent)
                            bigtxtcent[ibc][ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.15*(1.-top-bot),Form("Au+Au %i-%i%% 200GeV",limcent[icent],limcent[icent+1]));
                        else
                            bigtxtcent[ibc][ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.15*(1.-top-bot),Form("Au+Au %i-%i%% 200GeV",limcent[0],limcent[ncent]));
                        bigtxtcent[ibc][ipx][ipy]->SetTextAlign(22);
                        bigtxtcent[ibc][ipx][ipy]->SetNDC();
                        bigtxtcent[ibc][ipx][ipy]->SetTextFont(63);
                        bigtxtcent[ibc][ipx][ipy]->SetTextSize(60);

                        bigleg[ibc][ipx][ipy] = new TLegend(1.-right-.225*(1.-right-left),1.-top-.1*(1.-top-bot),1.-right,1.-top);
                        bigleg[ibc][ipx][ipy]->SetTextFont(43);
                        bigleg[ibc][ipx][ipy]->SetTextSize(60);

                        ypos += ht;
                    }
                    xpos += wd;
                }

                
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    int ibc = iptt/2;
                    int ipx = iptt%2;//pad x pos
                        
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        int ipy = ixi;//pad y pos
                        
                        TCanvas can ("can","",2500,2000);//2x3
                        TMultiGraph* mgr = new TMultiGraph;
                        TLegend leg (0.75,0.75,0.9,0.9);
                        leg.SetTextSize(0.04);

                        TGraphErrors* tgr = tgrv2j[iadj][itype][icent][iptt][ixi];
                        tgr->SetMarkerColor(TColor::GetColorPalette(440/3*2));
                        tgr->SetLineColor(TColor::GetColorPalette(440/3*2));
                        tgr->SetMarkerStyle(markstyle[0]);


                        mgr->Add(tgr,"p0");
                        TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                        tgrclone->SetMarkerSize(0);                     
                        mgr->Add(tgrclone,"p0");
                        if (icent != ncent)
                            leg.AddEntry(tgr,Form("%i-%i%%",limcent[icent],limcent[icent+1]),"p");
                        else
                            leg.AddEntry(tgr,Form("%i-%i%%",limcent[0],limcent[ncent]),"p");

                        // bigmgr[ibc][ipx][ipy]->Add(tgr_sys,"2");
                        bigmgr[ibc][ipx][ipy]->Add(tgr,"p0");
                        bigmgr[ibc][ipx][ipy]->Add(tgrclone,"p0");
                        if (icent != ncent)
                            bigleg[ibc][ipx][ipy]->AddEntry(tgr,Form("%i-%i%%",limcent[icent],limcent[icent+1]),"p");
                        else
                            bigleg[ibc][ipx][ipy]->AddEntry(tgr,Form("%i-%i%%",limcent[0],limcent[ncent]),"p");
                        
                        mgr->SetTitle(";#Delta#phi/#pi;v_{2} cos(2#Delta#phi)/d#Delta#phi");
                        mgr->Draw("A");
                        leg.Draw();

                        TLatex txtmom (0.4,.85,Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                        txtmom.SetTextAlign(22);
                        txtmom.SetNDC();
                        txtmom.SetTextSize(0.04);
                        txtmom.Draw();

                        // {
                        // }
                        // {
                        // }

                        TLatex txtcent (0.5,.8,"Au+Au 200GeV");
                        txtcent.SetTextAlign(22);
                        txtcent.SetNDC();
                        txtcent.SetTextSize(0.04);
                        txtcent.Draw();
                        
                        double ymin = -1.45;
                        double ymax = 1.45;
                        if (itype != GAT::PI0)
                        {
                            ymin = -5.95;
                            ymax = 5.95;
                        }
                        double ydiff = ymax-ymin;

                        double bigymin = ymin;
                        double bigymax = ymax;
                        
                        ymin = ymin - (ydiff)*.1;
                        ymax = ymax + (ydiff)*.3;

                        biglim[ibc][ipy][0] = min(bigymin,biglim[ibc][ipy][0]);
                        biglim[ibc][ipy][1] = max(bigymax,biglim[ibc][ipy][1]);
    
                        mgr->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);
                        mgr->GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                        TLine* line = new TLine(0.,0.,1.,0.);
                        line->SetLineWidth(1);
                        line->SetLineStyle(2);
                        line->Draw();

                        mgr->Draw();

                        can.Update();
                        can.SaveAs(Form("dphiplots/v2j_%s_adj%i_cent%i_ptt%i_xi%i.pdf",GAT::spec[itype],iadj,icent,iptt,ixi));

                        // draw on big pads too!
                        bigpad[ibc][ipx][ipy]->cd();
                        bigmgr[ibc][ipx][ipy]->SetTitle(";#Delta#phi (1/#pi);v_{2} cos(2#Delta#phi)/d#Delta#phi");
                        bigmgr[ibc][ipx][ipy]->Draw("A");

                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleFont(43);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelFont(43);
                        if (ipy == 0)
                        {
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleOffset(2.75);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(70);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(60);
                        }
                        else
                        {
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(0);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(0);
                        }

                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleFont(43);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelFont(43);
                        if (ipx == 0)
                        {
                            
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleOffset(3.5);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(70);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(60);
                        }
                        else
                        {
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(0);
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(0);
                        }
                        
                        line->Draw();
                        bigtxtmom[ibc][ipx][ipy]->Draw();
                        // bigtxtcent[ibc][ipx][ipy]->Draw();

                        if (ipx == 0 && ipy == 2)
                        {
                            bigleg[ibc][ipx][ipy]->Draw();
                            // bigtxtcent[ibc][ipx][ipy]->Draw();
                        }
                    }//ixi
                }//iptt

                for (int ibc = 0; ibc < 2; ibc++)
                {
                    for (int ipx = 0; ipx < 2; ipx++)
                    {
                        for (int ipy = 0; ipy < 3; ipy++)
                        {
                            bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetRangeUser(biglim[ibc][ipy][0],biglim[ibc][ipy][1]);
                            bigpad[ibc][ipx][ipy]->Modified();
                            bigpad[ibc][ipx][ipy]->Update();
                        }
                    }
                    
                    bigcan[ibc]->SaveAs(Form("dphiplots/%s.pdf",bigcan[ibc]->GetName()));
                }
            }
        }
    }
    //END (IN - OUT) / ALL

    //(IN - OUT) / ALL zoomed
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int icent = 0; icent <= ncent; icent++)
            {

                TCanvas* bigcan = new TCanvas(Form("bigcan_zoom_v2j_%s_adj%i_cent%i",GAT::spec[itype],iadj,icent),"",3000,3500);
                TPad* bigpad[4][3];// 4 pt columns, 3 xi rows
                TLatex* bigtxtmom[4][3];
                TLatex* bigtxtcent[4][3];
                TMultiGraph* bigmgr[4][3];
                double biglim[3][2];// x index not needed
                TLegend* bigleg[4][3];

                // set up canvas
                double xpos = 0.;
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    int ipx = iptt;//pad x pos

                    double left = 0.;
                    if (ipx == 0) left = .25;

                    double right = 0.;
                    if (ipx == 3) right = .15;

                    double wd = .85/2.875;
                    if (ipx == 1 || ipx == 2) wd *= .75;
                    if (ipx == 3) wd *= .75/.85;

                    double ypos = 0.;               
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        int ipy = ixi;//pad y pos

                        //set pad margins and create
                        double top = 0.;
                        if (ipy == 2) top = .15;

                        double bot = 0.;
                        if (ipy == 0) bot = .15;

                        double ht = .85/2.85;
                        if (ipy != 1) ht *= 1./.85;

                        bigpad[ipx][ipy] = new TPad(Form("bigpad%i%i",ipx,ipy),"",xpos,ypos,xpos+wd,ypos+ht);
                        bigmgr[ipx][ipy] = new TMultiGraph();//make mgr while were at it

                        //ok quick point, bigmgr not need anymore.. but we left it bc we are copying code and it doesn't hinder us
                        //although down below we try to redraw bigmgr over our error bands but it is not possible

                        // initialize y range limits (unecessarily repeated)
                        biglim[ipy][0] = numeric_limits<double>::max();
                        biglim[ipy][1] = numeric_limits<double>::min();

                        bigpad[ipx][ipy]->SetTopMargin(top);
                        bigpad[ipx][ipy]->SetBottomMargin(bot);
                        bigpad[ipx][ipy]->SetRightMargin(right);
                        bigpad[ipx][ipy]->SetLeftMargin(left);
                            
                        bigpad[ipx][ipy]->Draw();

                        double txtx = .5;
                        if (ipx == 0 && ipy == 2) txtx = .40;
                        
                        bigtxtmom[ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.075*(1.-top-bot),Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                        bigtxtmom[ipx][ipy]->SetTextAlign(22);
                        bigtxtmom[ipx][ipy]->SetNDC();
                        bigtxtmom[ipx][ipy]->SetTextFont(63);
                        bigtxtmom[ipx][ipy]->SetTextSize(60);

                        bigtxtcent[ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.15*(1.-top-bot),"Au+Au 200GeV");
                        bigtxtcent[ipx][ipy]->SetTextAlign(22);
                        bigtxtcent[ipx][ipy]->SetNDC();
                        bigtxtcent[ipx][ipy]->SetTextFont(63);
                        bigtxtcent[ipx][ipy]->SetTextSize(60);

                        bigleg[ipx][ipy] = new TLegend(1.-right-.225*(1.-right-left),1.-top-.1*(1.-top-bot),1.-right,1.-top);
                        bigleg[ipx][ipy]->SetTextFont(43);
                        bigleg[ipx][ipy]->SetTextSize(60);

                        ypos += ht;
                    }
                    xpos += wd;
                }


                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    int ipx = iptt;//pad x pos
                
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        int ipy = ixi;//pad y pos
                    
                        TCanvas can ("can","",2500,2000);//2x3
                        TMultiGraph* mgr = new TMultiGraph;
                        
                        TGraphErrors* tgr = tgrv2j[iadj][itype][icent][iptt][ixi];
                        tgr->SetMarkerColor(TColor::GetColorPalette(440/3*2));
                        tgr->SetLineColor(TColor::GetColorPalette(440/3*2));
                        tgr->SetMarkerStyle(markstyle[0]);

                        TF1* f1v2j = new TF1(Form("fit_%s",h1v2j[iadj][itype][icent][iptt][ixi]->GetName()),"[0]*cos(x*2.*TMath::Pi())",.5,1.);
                        f1v2j->SetLineWidth(1);;
                        //eventually must fit distribution with sys errors added to stat
                        tgr->Fit(f1v2j,"R0");

                        double v2jerr = f1v2j->GetParError(0);

                        double v2j = f1v2j->Eval(1);

                        tgrixiv2j[iadj][itype][iptt][icent]->SetPoint(ixi,xpos_xi[iptt][ixi],v2j);
                        tgrixiv2j[iadj][itype][iptt][icent]->SetPointError(ixi,0.,v2jerr);

                        tgrixiv2j_syserr[iadj][itype][iptt][icent]->SetPoint(ixi,xpos_xi[iptt][ixi],v2j);
                        tgrixiv2j_syserr[iadj][itype][iptt][icent]->SetPointError(ixi,sex_xi,v2jsyserr[iadj][itype][icent][iptt][ixi]);

                        tgripttv2j[iadj][itype][icent][ixi]->SetPoint(iptt,xpos_ptt[iptt][ixi],v2j);
                        tgripttv2j[iadj][itype][icent][ixi]->SetPointError(iptt,0.,v2jerr);

                        tgripttv2j_syserr[iadj][itype][icent][ixi]->SetPoint(iptt,xpos_ptt[iptt][ixi],v2j);
                        tgripttv2j_syserr[iadj][itype][icent][ixi]->SetPointError(iptt,sex_ptt,v2jsyserr[iadj][itype][icent][iptt][ixi]);

                        
                        mgr->Add(tgr,"p0");
                        TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                        tgrclone->SetMarkerSize(0);                     
                        mgr->Add(tgrclone,"p0");

                        // bigmgr[ipx][ipy]->Add(tgr_sys,"2");
                        bigmgr[ipx][ipy]->Add(tgr,"p0");
                        bigmgr[ipx][ipy]->Add(tgrclone,"p0");

                        mgr->SetTitle(";#Delta#phi/#pi;v_{2} cos(2#Delta#phi)/d#Delta#phi");
                    
                        mgr->Draw("A");

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

                        can.Update();
                        mgr->GetHistogram()->GetXaxis()->SetLimits(0.5,1.);
                        mgr->GetHistogram()->GetXaxis()->SetNdivisions(505);
                    
                        double ymin = INT_MAX;
                        double ymax = INT_MIN;
                        {
                            int n = tgr->GetN();
                            for (int i = n/2; i < n; i++)
                            {
                                ymin = min(ymin,tgr->GetY()[i] - tgr->GetEY()[i]);
                            

                                ymax = max(ymax,tgr->GetY()[i] + tgr->GetEY()[i]);
                            
                            }
                        }
                        ymin = -.7;
                        ymax = .7;
                        if (itype != GAT::PI0)
                        {
                            ymin = -2.9;
                            ymax = 2.9;
                        }

                        double bigymin = ymin;
                        double bigymax = ymax;

                        biglim[ipy][0] = min(bigymin,biglim[ipy][0]);
                        biglim[ipy][1] = max(bigymax,biglim[ipy][1]);

                        double ydiff = ymax-ymin;
                    
                        ymin = ymin - (ydiff)*.1;
                        ymax = ymax + (ydiff)*.3;
    
                        mgr->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);

                        double xlims[2] = {0.501,.999};
                        TLine* line = new TLine(xlims[0],0.,xlims[1],0.);
                        line->SetLineWidth(1);
                        line->SetLineStyle(2);
                        line->Draw();

                        mgr->Draw();

                        can.Update();
                        can.SaveAs(Form("dphiplots/v2j_zoom_%s_adj%i_cent%i_ptt%i_xi%i.pdf",GAT::spec[itype],iadj,icent,iptt,ixi));

                        // draw on big pads too!
                        bigpad[ipx][ipy]->cd();
                        // bigmgr[ipx][ipy]->SetTitle(";#Delta#phi/#pi;1/N_{t}  dN_{ h }/d#Delta#phi");
                        if (ipx == 0 && ipy == 2)
                            bigmgr[ipx][ipy]->SetTitle(";;v_{2} cos(2#Delta#phi)/d#Delta#phi");
                        else if (ipx == 3 && ipy == 0)
                            bigmgr[ipx][ipy]->SetTitle(";#Delta#phi/#pi;");
                        else
                            bigmgr[ipx][ipy]->SetTitle("");

                        //ticks
                        if (ipx != 3)
                            bigpad[ipx][ipy]->SetTicky(2);
                        if (ipy != 2)
                            bigpad[ipx][ipy]->SetTickx(2);

                    
                        bigmgr[ipx][ipy]->Draw("A");

                        // DrawError(f1v2j,v2jerr,.5);
                        // DrawError(f1v2j,v2jsyserr[iadj][itype][icent][iptt][ixi],.5,kRed,.30);

                        DrawError(f1v2j,v2jerr,.5,kOrange,.3);
                        DrawError(f1v2j,v2jsyserr[iadj][itype][icent][iptt][ixi],.5,kRed,.3);


                        //store v2j and errors in tgraphs

                        tgr->Draw("p0 same");
                        tgrclone->Draw("p0 same");

                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLimits(xlims[0],xlims[1]);
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetNdivisions(505);
                    
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleFont(43);
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelFont(43);
                        if (ipy == 0)
                        {
                            bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleOffset(2.75);
                            bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(70);
                            bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(60);
                        }
                        else
                        {
                            bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(0);
                            bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(0);
                        }

                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleFont(43);
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelFont(43);
                        if (ipx == 0)
                        {
                            
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleOffset(6.0);
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelOffset(.01);
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(70);
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(60);
                        }
                        else
                        {
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(0);
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(0);
                        }

                        if (ipy == 2)
                        {                       
                            double rt = bigpad[ipx][ipy]->GetRightMargin();
                            double lt = bigpad[ipx][ipy]->GetLeftMargin();
                            double xpos = ((1.-rt)+lt)/2.;

                            TLatex* txt = new TLatex(xpos,.89,Form("%.1f #leq p_{Tt} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
                        
                            txt->SetNDC();
                            txt->SetTextAlign(22);
                            txt->SetTextFont(43);
                            txt->SetTextSize(60);
                            txt->Draw();
                        }

                        if (ipx == 3)
                        {
                            TLatex* txt;
                            double tp = bigpad[ipx][ipy]->GetTopMargin();
                            double bt = bigpad[ipx][ipy]->GetBottomMargin();
                            double ypos = ((1.-tp)+bt)/2.;

                            txt = new TLatex(.92,ypos,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]));

                            txt->SetNDC();
                            txt->SetTextAlign(22);
                            txt->SetTextAngle(90);
                            txt->SetTextFont(43);
                            txt->SetTextSize(60);
                            txt->Draw();
                        }
                    
                        line->Draw();
                    
                    }//xi
                }//ptt

                for (int ipy = 0; ipy < 3; ipy++)
                {
                    double ydiff = biglim[ipy][1]-biglim[ipy][0];
                    
                    biglim[ipy][0] = biglim[ipy][0] - (ydiff)*.05;
                    biglim[ipy][1] = biglim[ipy][1] + (ydiff)*.05;
                    
                    for (int ipx = 0; ipx < 4; ipx++)
                    {
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetRangeUser(biglim[ipy][0],biglim[ipy][1]);
                        bigpad[ipx][ipy]->Modified();
                        bigpad[ipx][ipy]->Update();
                    }
                }
                    
                bigcan->SaveAs(Form("dphiplots/%s.pdf",bigcan->GetName()));

            }
        }
    }
    //END (IN - OUT) / ALL zoomed


    //central and mid central same plot
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            TCanvas* bigcan[2];// one for 5-7 & 7-9 and the other 9-12 & 12-15
            TPad* bigpad[2][2][3];// 2 pt columns, 3 xi rows
            TLatex* bigtxtmom[2][2][3];
            TLatex* bigtxtcent[2][2][3];
            TMultiGraph* bigmgr[2][2][3];
            double biglim[2][3][2];// x index not needed
            TLegend* bigleg[2][2][3];

            // set up canvas
            double xpos = 0.;
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {

                int ibc = iptt/2;
                if (iptt%2==0)
                {
                    xpos = 0;
                    bigcan[ibc] = new TCanvas(Form("bigcan%i_central_%s_adj%i",ibc,GAT::spec[itype],iadj),"",3000,3500);
                }
                int ipx = iptt%2;//pad x pos

                double left = .02;
                if (ipx == 0) left = .15;

                double right = .02;
                if (ipx == 1) right = .05;

                double wd = .93/1.76;
                if (ipx == 1) wd *= .83/.93;

                double ypos = 0.;                   
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    int ipy = ixi;//pad y pos

                    //set pad margins and create
                    double top = .01;
                    if (ipy == 2) top = .15;

                    double bot = .01;
                    if (ipy == 0) bot = .15;

                    double ht = .84/2.80;
                    if (ipy != 1) ht *= .98/.84;

                    bigpad[ibc][ipx][ipy] = new TPad(Form("bigpad%i%i%i",ibc,ipx,ipy),"",xpos,ypos,xpos+wd,ypos+ht);

                    bigmgr[ibc][ipx][ipy] = new TMultiGraph();//make mgr while were at it

                    // initialize y range limits (unecessarily repeated)
                    biglim[ibc][ipy][0] = numeric_limits<double>::max();
                    biglim[ibc][ipy][1] = numeric_limits<double>::min();

                    bigpad[ibc][ipx][ipy]->SetTopMargin(top);
                    bigpad[ibc][ipx][ipy]->SetBottomMargin(bot);
                    bigpad[ibc][ipx][ipy]->SetRightMargin(right);
                    bigpad[ibc][ipx][ipy]->SetLeftMargin(left);
                            
                    bigpad[ibc][ipx][ipy]->Draw();


                    double txtx = .5;
                    if (ipx == 0 && ipy == 2) txtx = .40;
                        
                    bigtxtmom[ibc][ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.075*(1.-top-bot),Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                    bigtxtmom[ibc][ipx][ipy]->SetTextAlign(22);
                    bigtxtmom[ibc][ipx][ipy]->SetNDC();
                    bigtxtmom[ibc][ipx][ipy]->SetTextFont(63);
                    bigtxtmom[ibc][ipx][ipy]->SetTextSize(60);

                    bigtxtcent[ibc][ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.15*(1.-top-bot),"Au+Au 200GeV");
                    bigtxtcent[ibc][ipx][ipy]->SetTextAlign(22);
                    bigtxtcent[ibc][ipx][ipy]->SetNDC();
                    bigtxtcent[ibc][ipx][ipy]->SetTextFont(63);
                    bigtxtcent[ibc][ipx][ipy]->SetTextSize(60);

                    bigleg[ibc][ipx][ipy] = new TLegend(1.-right-.225*(1.-right-left),1.-top-.2*(1.-top-bot),1.-right,1.-top);
                    bigleg[ibc][ipx][ipy]->SetTextFont(43);
                    bigleg[ibc][ipx][ipy]->SetTextSize(60);

                    ypos += ht;
                }
                xpos += wd;
            }


            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                int ibc = iptt/2;
                int ipx = iptt%2;//pad x pos
                
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    int ipy = ixi;//pad y pos
                    
                    TCanvas can ("can","",2500,2000);//2x3
                    TMultiGraph* mgr = new TMultiGraph;
                    TLegend leg (0.75,0.75,0.9,0.9);
                    leg.SetTextSize(0.04);
                        
                    for (int icent = 0; icent < ncent; icent++)
                    {
                        TGraphErrors* tgr = tgrdphi[iadj][itype][icent][ALL][iptt][ixi];
                        tgr->SetMarkerColor(TColor::GetColorPalette(440/3*icent));
                        tgr->SetLineColor(TColor::GetColorPalette(440/3*icent));
                        tgr->SetMarkerStyle(markstyle[icent]);

                        TGraphErrors* tgr_sys = tgrdphi_sys[iadj][itype][icent][ALL][iptt][ixi];
                        tgr_sys->SetFillColorAlpha(TColor::GetColorPalette(440/3*icent),.5);

                        mgr->Add(tgr_sys,"2");
                        mgr->Add(tgr,"p");
                        TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                        tgrclone->SetMarkerSize(0);                     
                        mgr->Add(tgrclone,"p");
                        leg.AddEntry(tgr,Form("%i-%i%%",limcent[icent],limcent[icent+1]),"p");

                        bigmgr[ibc][ipx][ipy]->Add(tgr_sys,"2");
                        bigmgr[ibc][ipx][ipy]->Add(tgr,"p");
                        bigmgr[ibc][ipx][ipy]->Add(tgrclone,"p");
                        bigleg[ibc][ipx][ipy]->AddEntry(tgr,Form("%i-%i%%",limcent[icent],limcent[icent+1]),"p");

                    }

                    mgr->SetTitle(";#Delta#phi/#pi;1/N_{ t}  dN_{ h }/d#Delta#phi");
                    mgr->Draw("A");
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

                    double ymin = min(0.,mgr->GetHistogram()->GetMinimum());
                    double ymax = max(0.,mgr->GetHistogram()->GetMaximum());
                    double ydiff = ymax-ymin;

                    double bigymin = ymin;
                    double bigymax = ymax;
                    
                    ymin = ymin - (ydiff)*.1;
                    ymax = ymax + (ydiff)*.3;

                    biglim[ibc][ipy][0] = min(bigymin,biglim[ibc][ipy][0]);
                    biglim[ibc][ipy][1] = max(bigymax,biglim[ibc][ipy][1]);
    
                    mgr->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);
                    mgr->GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                    TLine* line = new TLine(0.,0.,1.,0.);
                    line->SetLineWidth(1);
                    line->SetLineStyle(2);
                    line->Draw();

                    mgr->Draw();

                    can.Update();
                    can.SaveAs(Form("dphiplots/central_%s_adj%i_ptt%i_xi%i.pdf",GAT::spec[itype],iadj,iptt,ixi));

                    // draw on big pads too!
                    bigpad[ibc][ipx][ipy]->cd();
                    bigmgr[ibc][ipx][ipy]->SetTitle(";#Delta#phi (1/#pi);1/N_{ t}  dN_{ h }/d#Delta#phi");
                    bigmgr[ibc][ipx][ipy]->Draw("A");

                    bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                    bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleFont(43);
                    bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelFont(43);
                    if (ipy == 0)
                    {
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleOffset(2.75);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(70);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(60);
                    }
                    else
                    {
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(0);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(0);
                    }

                    bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleFont(43);
                    bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelFont(43);
                    if (ipx == 0)
                    {
                            
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleOffset(3.5);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(70);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(60);
                    }
                    else
                    {
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(0);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(0);
                    }
                        
                    line->Draw();
                    bigtxtmom[ibc][ipx][ipy]->Draw();
                    // bigtxtcent[ibc][ipx][ipy]->Draw();

                    if (ipx == 0 && ipy == 2)
                    {
                        bigleg[ibc][ipx][ipy]->Draw();
                        bigtxtcent[ibc][ipx][ipy]->Draw();
                    }
                }//xi
            }//ptt

            for (int ibc = 0; ibc < 2; ibc++)
            {
                for (int ipx = 0; ipx < 2; ipx++)
                {
                    for (int ipy = 0; ipy < 3; ipy++)
                    {
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetRangeUser(biglim[ibc][ipy][0],biglim[ibc][ipy][1]);
                        bigpad[ibc][ipx][ipy]->Modified();
                        bigpad[ibc][ipx][ipy]->Update();
                    }
                }
                    
                bigcan[ibc]->SaveAs(Form("dphiplots/%s.pdf",bigcan[ibc]->GetName()));
            }

        }
    }

    //0-60% only. copied from cent ratio
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            TCanvas* bigcan[2];// one for 5-7 & 7-9 and the other 9-12 & 12-15
            TPad* bigpad[2][2][3];// 2 pt columns, 3 xi rows
            TLatex* bigtxtmom[2][2][3];
            TLatex* bigtxtcent[2][2][3];
            TMultiGraph* bigmgr[2][2][3];
            double biglim[2][3][2];// x index not needed
            TLegend* bigleg[2][2][3];

            // set up canvas
            double xpos = 0.;
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {

                int ibc = iptt/2;
                if (iptt%2==0)
                {
                    xpos = 0;
                    bigcan[ibc] = new TCanvas(Form("bigcan%i_060_%s_adj%i",ibc,GAT::spec[itype],iadj),"",3000,3500);
                }
                int ipx = iptt%2;//pad x pos

                double left = .02;
                if (ipx == 0) left = .15;

                double right = .02;
                if (ipx == 1) right = .05;

                double wd = .93/1.76;
                if (ipx == 1) wd *= .83/.93;

                double ypos = 0.;                   
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    int ipy = ixi;//pad y pos

                    //set pad margins and create
                    double top = .01;
                    if (ipy == 2) top = .15;

                    double bot = .01;
                    if (ipy == 0) bot = .15;

                    double ht = .84/2.80;
                    if (ipy != 1) ht *= .98/.84;

                    bigpad[ibc][ipx][ipy] = new TPad(Form("bigpad%i%i%i",ibc,ipx,ipy),"",xpos,ypos,xpos+wd,ypos+ht);

                    bigmgr[ibc][ipx][ipy] = new TMultiGraph();//make mgr while were at it

                    // initialize y range limits (unecessarily repeated)
                    biglim[ibc][ipy][0] = numeric_limits<double>::max();
                    biglim[ibc][ipy][1] = numeric_limits<double>::min();

                    bigpad[ibc][ipx][ipy]->SetTopMargin(top);
                    bigpad[ibc][ipx][ipy]->SetBottomMargin(bot);
                    bigpad[ibc][ipx][ipy]->SetRightMargin(right);
                    bigpad[ibc][ipx][ipy]->SetLeftMargin(left);
                            
                    bigpad[ibc][ipx][ipy]->Draw();


                    double txtx = .5;
                    if (ipx == 0 && ipy == 2) txtx = .40;
                        
                    bigtxtmom[ibc][ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.075*(1.-top-bot),Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                    bigtxtmom[ibc][ipx][ipy]->SetTextAlign(22);
                    bigtxtmom[ibc][ipx][ipy]->SetNDC();
                    bigtxtmom[ibc][ipx][ipy]->SetTextFont(63);
                    bigtxtmom[ibc][ipx][ipy]->SetTextSize(60);

                    bigtxtcent[ibc][ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.15*(1.-top-bot),"Au+Au 200GeV");
                    bigtxtcent[ibc][ipx][ipy]->SetTextAlign(22);
                    bigtxtcent[ibc][ipx][ipy]->SetNDC();
                    bigtxtcent[ibc][ipx][ipy]->SetTextFont(63);
                    bigtxtcent[ibc][ipx][ipy]->SetTextSize(60);

                    bigleg[ibc][ipx][ipy] = new TLegend(1.-right-.225*(1.-right-left),1.-top-.1*(1.-top-bot),1.-right,1.-top);
                    bigleg[ibc][ipx][ipy]->SetTextFont(43);
                    bigleg[ibc][ipx][ipy]->SetTextSize(60);

                    ypos += ht;
                }
                xpos += wd;
            }


            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                int ibc = iptt/2;
                int ipx = iptt%2;//pad x pos
                
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    int ipy = ixi;//pad y pos
                    
                    TCanvas can ("can","",2500,2000);//2x3
                    TMultiGraph* mgr = new TMultiGraph;
                    TLegend leg (0.75,0.75,0.9,0.9);
                    leg.SetTextSize(0.04);
                        
                    TGraphErrors* tgr = tgrdphi[iadj][itype][ncent][ALL][iptt][ixi];
                    tgr->SetMarkerColor(TColor::GetColorPalette(440/3*ncent));
                    tgr->SetLineColor(TColor::GetColorPalette(440/3*ncent));
                    tgr->SetMarkerStyle(markstyle[0]);

                    TGraphErrors* tgr_sys = tgrdphi_sys[iadj][itype][ncent][ALL][iptt][ixi];
                    tgr_sys->SetFillColorAlpha(TColor::GetColorPalette(440/3*ncent),.5);

                    mgr->Add(tgr_sys,"2");
                    mgr->Add(tgr,"p");
                    TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                    tgrclone->SetMarkerSize(0);                 
                    mgr->Add(tgrclone,"p");
                    leg.AddEntry(tgr,Form("%i-%i%%",limcent[0],limcent[ncent]),"p");

                    bigmgr[ibc][ipx][ipy]->Add(tgr_sys,"2");
                    bigmgr[ibc][ipx][ipy]->Add(tgr,"p");
                    bigmgr[ibc][ipx][ipy]->Add(tgrclone,"p");
                    bigleg[ibc][ipx][ipy]->AddEntry(tgr,Form("%i-%i%%",limcent[0],limcent[ncent]),"p");

                    mgr->SetTitle(";#Delta#phi/#pi;1/N_{ t}  dN_{ h }/d#Delta#phi");
                    mgr->Draw("A");
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

                    double ymin = min(0.,mgr->GetHistogram()->GetMinimum());
                    double ymax = max(0.,mgr->GetHistogram()->GetMaximum());
                    double ydiff = ymax-ymin;

                    double bigymin = ymin;
                    double bigymax = ymax;
                    
                    ymin = ymin - (ydiff)*.1;
                    ymax = ymax + (ydiff)*.3;

                    biglim[ibc][ipy][0] = min(bigymin,biglim[ibc][ipy][0]);
                    biglim[ibc][ipy][1] = max(bigymax,biglim[ibc][ipy][1]);
    
                    mgr->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);
                    mgr->GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                    TLine* line = new TLine(0.,0.,1.,0.);
                    line->SetLineWidth(1);
                    line->SetLineStyle(2);
                    line->Draw();

                    mgr->Draw();

                    can.Update();
                    can.SaveAs(Form("dphiplots/060_%s_adj%i_ptt%i_xi%i.pdf",GAT::spec[itype],iadj,iptt,ixi));

                    // draw on big pads too!
                    bigpad[ibc][ipx][ipy]->cd();
                    bigmgr[ibc][ipx][ipy]->SetTitle(";#Delta#phi (1/#pi);1/N_{ t}  dN_{ h }/d#Delta#phi");
                    bigmgr[ibc][ipx][ipy]->Draw("A");

                    bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                    bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleFont(43);
                    bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelFont(43);
                    if (ipy == 0)
                    {
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleOffset(2.75);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(70);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(60);
                    }
                    else
                    {
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(0);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(0);
                    }

                    bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleFont(43);
                    bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelFont(43);
                    if (ipx == 0)
                    {
                            
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleOffset(3.5);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(70);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(60);
                    }
                    else
                    {
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(0);
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(0);
                    }
                        
                    line->Draw();
                    bigtxtmom[ibc][ipx][ipy]->Draw();
                    // bigtxtcent[ibc][ipx][ipy]->Draw();

                    if (ipx == 0 && ipy == 2)
                    {
                        bigleg[ibc][ipx][ipy]->Draw();
                        bigtxtcent[ibc][ipx][ipy]->Draw();
                    }
                }//xi
            }//ptt

            for (int ibc = 0; ibc < 2; ibc++)
            {
                for (int ipx = 0; ipx < 2; ipx++)
                {
                    for (int ipy = 0; ipy < 3; ipy++)
                    {
                        bigmgr[ibc][ipx][ipy]->GetHistogram()->GetYaxis()->SetRangeUser(biglim[ibc][ipy][0],biglim[ibc][ipy][1]);
                        bigpad[ibc][ipx][ipy]->Modified();
                        bigpad[ibc][ipx][ipy]->Update();
                    }
                }
                    
                bigcan[ibc]->SaveAs(Form("dphiplots/%s.pdf",bigcan[ibc]->GetName()));
            }

        }
    }

    //0-60% only. all ptt one plot. zoomed.
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            TCanvas* bigcan = new TCanvas(Form("bigcan_zoom_060_%s_adj%i",GAT::spec[itype],iadj),"",3000,3500);
            TPad* bigpad[4][3];// 4 pt columns, 3 xi rows
            TLatex* bigtxtmom[4][3];
            TLatex* bigtxtcent[4][3];
            TMultiGraph* bigmgr[4][3];
            double biglim[3][2];// x index not needed
            TLegend* bigleg[4][3];

            
            // set up canvas
            double xpos = 0.;
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {

                int ipx = iptt;//pad x pos


                // //(1-left-right)w


                double left = 0.;
                if (ipx == 0) left = .25;

                double right = 0.;
                if (ipx == 3) right = .15;

                //(1-left-right)w

                double wd = .85/2.875;
                if (ipx == 1 || ipx == 2) wd *= .75;
                if (ipx == 3) wd *= .75/.85;

                double ypos = 0.;                   
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    int ipy = ixi;//pad y pos

                    // //set pad margins and create


                    // //set pad margins and create


                    //set pad margins and create
                    double top = 0.;
                    if (ipy == 2) top = .15;

                    double bot = 0.;
                    if (ipy == 0) bot = .15;


                    double ht = .85/2.85;
                    if (ipy != 1) ht *= 1./.85;


                    bigpad[ipx][ipy] = new TPad(Form("bigpad%i%i",ipx,ipy),"",xpos,ypos,xpos+wd,ypos+ht);
                    bigmgr[ipx][ipy] = new TMultiGraph();//make mgr while were at it

                    // initialize y range limits (unecessarily repeated)
                    biglim[ipy][0] = numeric_limits<double>::max();
                    biglim[ipy][1] = numeric_limits<double>::min();

                    bigpad[ipx][ipy]->SetTopMargin(top);
                    bigpad[ipx][ipy]->SetBottomMargin(bot);
                    bigpad[ipx][ipy]->SetRightMargin(right);
                    bigpad[ipx][ipy]->SetLeftMargin(left);
                            
                    bigpad[ipx][ipy]->Draw();


                    double txtx = .5;
                    if (ipx == 0 && ipy == 2) txtx = .40;
                        
                    bigtxtmom[ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.075*(1.-top-bot),Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                    bigtxtmom[ipx][ipy]->SetTextAlign(22);
                    bigtxtmom[ipx][ipy]->SetNDC();
                    bigtxtmom[ipx][ipy]->SetTextFont(63);
                    bigtxtmom[ipx][ipy]->SetTextSize(60);

                    bigtxtcent[ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.15*(1.-top-bot),"Au+Au 200GeV");
                    bigtxtcent[ipx][ipy]->SetTextAlign(22);
                    bigtxtcent[ipx][ipy]->SetNDC();
                    bigtxtcent[ipx][ipy]->SetTextFont(63);
                    bigtxtcent[ipx][ipy]->SetTextSize(60);

                    bigleg[ipx][ipy] = new TLegend(1.-right-.225*(1.-right-left),1.-top-.1*(1.-top-bot),1.-right,1.-top);
                    bigleg[ipx][ipy]->SetTextFont(43);
                    bigleg[ipx][ipy]->SetTextSize(60);

                    ypos += ht;
                }
                xpos += wd;
            }


            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                int ipx = iptt;//pad x pos
                
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    int ipy = ixi;//pad y pos
                    
                    TCanvas can ("can","",2500,2000);//2x3
                    TMultiGraph* mgr = new TMultiGraph;
                        
                    TGraphErrors* tgr = tgrdphi[iadj][itype][ncent][ALL][iptt][ixi];
                    tgr->SetMarkerColor(TColor::GetColorPalette(440/3*ncent));
                    tgr->SetLineColor(TColor::GetColorPalette(440/3*ncent));
                    tgr->SetMarkerStyle(markstyle[0]);

                    TGraphErrors* tgr_sys = tgrdphi_sys[iadj][itype][ncent][ALL][iptt][ixi];
                    tgr_sys->SetFillColorAlpha(TColor::GetColorPalette(440/3*ncent),.5);

                    mgr->Add(tgr_sys,"2");
                    mgr->Add(tgr,"p");
                    TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                    tgrclone->SetMarkerSize(0);                 
                    mgr->Add(tgrclone,"p");

                    bigmgr[ipx][ipy]->Add(tgr_sys,"2");
                    bigmgr[ipx][ipy]->Add(tgr,"p");
                    bigmgr[ipx][ipy]->Add(tgrclone,"p");

                    mgr->SetTitle(";#Delta#phi/#pi;1/N_{ t}  dN_{ h }/d#Delta#phi");
                    
                    mgr->Draw("A");

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

                    can.Update();
                    mgr->GetHistogram()->GetXaxis()->SetLimits(0.5,1.);
                    mgr->GetHistogram()->GetXaxis()->SetNdivisions(505);
                    
                    //will getminimum use current range?

                    // this method somehow fucks up limits on big graph..
                    // {
                    //  {
                    //  }
                    // }

                    double ymin = INT_MAX;
                    double ymax = INT_MIN;
                    {
                        int n = tgr->GetN();
                        for (int i = n/2; i < n; i++)
                        {
                            ymin = min(ymin,tgr->GetY()[i] - tgr->GetEY()[i]);
                            
                            ymin = min(ymin,tgr_sys->GetY()[i] - tgr_sys->GetEY()[i]);

                            
                            ymax = max(ymax,tgr->GetY()[i] + tgr->GetEY()[i]);
                            
                            ymax = max(ymax,tgr_sys->GetY()[i] + tgr_sys->GetEY()[i]);
                        }
                    }

                    double bigymin = ymin;
                    double bigymax = ymax;

                    biglim[ipy][0] = min(bigymin,biglim[ipy][0]);
                    biglim[ipy][1] = max(bigymax,biglim[ipy][1]);

                    double ydiff = ymax-ymin;
                    
                    ymin = ymin - (ydiff)*.1;
                    ymax = ymax + (ydiff)*.3;
    
                    mgr->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);

                    double xlims[2] = {0.501,.999};
                    TLine* line = new TLine(xlims[0],0.,xlims[1],0.);
                    line->SetLineWidth(1);
                    line->SetLineStyle(2);
                    line->Draw();

                    mgr->Draw();

                    can.Update();
                    can.SaveAs(Form("dphiplots/060_zoom_%s_adj%i_ptt%i_xi%i.pdf",GAT::spec[itype],iadj,iptt,ixi));

                    // draw on big pads too!
                    bigpad[ipx][ipy]->cd();
                    // bigmgr[ipx][ipy]->SetTitle(";#Delta#phi/#pi;1/N_{t}  dN_{ h }/d#Delta#phi");
                    if (ipx == 0 && ipy == 2)
                        bigmgr[ipx][ipy]->SetTitle(";;1/N_{ t}  dN_{ h }/d#Delta#phi");
                    else if (ipx == 3 && ipy == 0)
                        bigmgr[ipx][ipy]->SetTitle(";#Delta#phi/#pi;");
                    else
                        bigmgr[ipx][ipy]->SetTitle("");

                    //ticks
                    if (ipx != 3)
                        bigpad[ipx][ipy]->SetTicky(2);
                    if (ipy != 2)
                        bigpad[ipx][ipy]->SetTickx(2);

                    
                    bigmgr[ipx][ipy]->Draw("A");

                    bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLimits(xlims[0],xlims[1]);
                    bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetNdivisions(505);
                    
                    bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleFont(43);
                    bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelFont(43);
                    if (ipy == 0)
                    {
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleOffset(2.75);
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(70);
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(60);
                    }
                    else
                    {
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(0);
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(0);
                    }

                    bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleFont(43);
                    bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelFont(43);
                    if (ipx == 0)
                    {
                            
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleOffset(6.0);
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelOffset(.01);
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(70);
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(60);
                    }
                    else
                    {
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(0);
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(0);
                    }

                    if (ipy == 2)
                    {                   
                        double rt = bigpad[ipx][ipy]->GetRightMargin();
                        double lt = bigpad[ipx][ipy]->GetLeftMargin();
                        double xpos = ((1.-rt)+lt)/2.;

                        TLatex* txt = new TLatex(xpos,.89,Form("%.1f #leq p_{Tt} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
                        
                        txt->SetNDC();
                        txt->SetTextAlign(22);
                        txt->SetTextFont(43);
                        txt->SetTextSize(60);
                        txt->Draw();
                    }

                    if (ipx == 3)
                    {
                        TLatex* txt;
                        double tp = bigpad[ipx][ipy]->GetTopMargin();
                        double bt = bigpad[ipx][ipy]->GetBottomMargin();
                        double ypos = ((1.-tp)+bt)/2.;

                        txt = new TLatex(.92,ypos,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]));

                        txt->SetNDC();
                        txt->SetTextAlign(22);
                        txt->SetTextAngle(90);
                        txt->SetTextFont(43);
                        txt->SetTextSize(60);
                        txt->Draw();

                    }


                    line->Draw();
                    // bigtxtmom[ipx][ipy]->Draw();
                    // bigtxtcent[ipx][ipy]->Draw();

                    // {
                    //  bigtxtcent[ipx][ipy]->Draw();
                    // }
                }//xi
            }//ptt

            for (int ipy = 0; ipy < 3; ipy++)
            {
                double ydiff = biglim[ipy][1]-biglim[ipy][0];
                    
                biglim[ipy][0] = biglim[ipy][0] - (ydiff)*.05;
                biglim[ipy][1] = biglim[ipy][1] + (ydiff)*.05;
                    
                for (int ipx = 0; ipx < 4; ipx++)
                {
                    bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetRangeUser(biglim[ipy][0],biglim[ipy][1]);
                    bigpad[ipx][ipy]->Modified();
                    bigpad[ipx][ipy]->Update();
                }
            }
                    
            bigcan->SaveAs(Form("dphiplots/%s.pdf",bigcan->GetName()));

        }
    }

    //central and mid central same plot. all ptt one plot. zoomed. copied from 0-60 zoomed
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            TCanvas* bigcan = new TCanvas(Form("bigcan_zoom_central_%s_adj%i",GAT::spec[itype],iadj),"",3000,3500);
            TPad* bigpad[4][3];// 4 pt columns, 3 xi rows
            TLatex* bigtxtmom[4][3];
            TLatex* bigtxtcent[4][3];
            TMultiGraph* bigmgr[4][3];
            double biglim[3][2];// x index not needed
            TLegend* bigleg[4][3];          

            // set up canvas
            double xpos = 0.;
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {

                int ipx = iptt;//pad x pos

                //more widht options shown in preivious plottings

                double left = 0.;
                if (ipx == 0) left = .25;

                double right = 0.;
                if (ipx == 3) right = .15;

                //(1-left-right)w

                double wd = .85/2.875;
                if (ipx == 1 || ipx == 2) wd *= .75;
                if (ipx == 3) wd *= .75/.85;

                double ypos = 0.;                   
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    int ipy = ixi;//pad y pos

                    //more pad height options shown in previous plottings

                    //set pad margins and create
                    double top = 0.;
                    if (ipy == 2) top = .15;

                    double bot = 0.;
                    if (ipy == 0) bot = .15;


                    double ht = .85/2.85;
                    if (ipy != 1) ht *= 1./.85;


                    bigpad[ipx][ipy] = new TPad(Form("bigpad%i%i",ipx,ipy),"",xpos,ypos,xpos+wd,ypos+ht);
                    bigmgr[ipx][ipy] = new TMultiGraph();//make mgr while were at it

                    // initialize y range limits (unecessarily repeated)
                    biglim[ipy][0] = numeric_limits<double>::max();
                    biglim[ipy][1] = numeric_limits<double>::min();

                    bigpad[ipx][ipy]->SetTopMargin(top);
                    bigpad[ipx][ipy]->SetBottomMargin(bot);
                    bigpad[ipx][ipy]->SetRightMargin(right);
                    bigpad[ipx][ipy]->SetLeftMargin(left);
                            
                    bigpad[ipx][ipy]->Draw();


                    double txtx = .5;
                    if (ipx == 0 && ipy == 2) txtx = .40;
                        
                    bigtxtmom[ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.075*(1.-top-bot),Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                    bigtxtmom[ipx][ipy]->SetTextAlign(22);
                    bigtxtmom[ipx][ipy]->SetNDC();
                    bigtxtmom[ipx][ipy]->SetTextFont(63);
                    bigtxtmom[ipx][ipy]->SetTextSize(60);

                    bigtxtcent[ipx][ipy] = new TLatex(txtx*(1.-left-right)+left,1.-top-.15*(1.-top-bot),"Au+Au 200GeV");
                    bigtxtcent[ipx][ipy]->SetTextAlign(22);
                    bigtxtcent[ipx][ipy]->SetNDC();
                    bigtxtcent[ipx][ipy]->SetTextFont(63);
                    bigtxtcent[ipx][ipy]->SetTextSize(60);

                    bigleg[ipx][ipy] = new TLegend(.95-right-.45*(.95-right-left),.98-top-.2*(.98-top-bot),.95-right,.98-top);
                    bigleg[ipx][ipy]->SetBorderSize(0);
                    bigleg[ipx][ipy]->SetTextFont(43);
                    bigleg[ipx][ipy]->SetTextSize(60);

                    ypos += ht;
                }
                xpos += wd;
            }


            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                int ipx = iptt;//pad x pos
                
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    int ipy = ixi;//pad y pos
                    
                    TCanvas can ("can","",2500,2000);//2x3
                    TMultiGraph* mgr = new TMultiGraph;

                    TLegend leg (0.75,0.75,0.9,0.9);
                    leg.SetTextSize(0.04);

                    double ymin = INT_MAX;
                    double ymax = INT_MIN;

                    for (int icent = 0; icent < ncent; icent++)
                    {                   
                        TGraphErrors* tgr = tgrdphi[iadj][itype][icent][ALL][iptt][ixi];
                        tgr->SetMarkerColor(TColor::GetColorPalette(440/3*icent));
                        tgr->SetLineColor(TColor::GetColorPalette(440/3*icent));
                        tgr->SetMarkerStyle(markstyle[0]);

                        TGraphErrors* tgr_sys = tgrdphi_sys[iadj][itype][icent][ALL][iptt][ixi];
                        tgr_sys->SetFillColorAlpha(TColor::GetColorPalette(440/3*icent),.5);                    

                        mgr->Add(tgr_sys,"2");
                        mgr->Add(tgr,"p");
                        TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                        tgrclone->SetMarkerSize(0);                     
                        mgr->Add(tgrclone,"p");

                        leg.AddEntry(tgr,Form("%i-%i%%",limcent[icent],limcent[icent+1]),"p");

                        bigmgr[ipx][ipy]->Add(tgr_sys,"2");
                        bigmgr[ipx][ipy]->Add(tgr,"p");
                        bigmgr[ipx][ipy]->Add(tgrclone,"p");
                        bigleg[ipx][ipy]->AddEntry(tgr,Form("%i-%i%%",limcent[icent],limcent[icent+1]),"p");

                        int n = tgr->GetN();
                        for (int i = n/2; i < n; i++)
                        {
                            ymin = min(ymin,tgr->GetY()[i] - tgr->GetEY()[i]);
                            
                            ymin = min(ymin,tgr_sys->GetY()[i] - tgr_sys->GetEY()[i]);

                            
                            ymax = max(ymax,tgr->GetY()[i] + tgr->GetEY()[i]);
                            
                            ymax = max(ymax,tgr_sys->GetY()[i] + tgr_sys->GetEY()[i]);
                        }
                    }

                    mgr->SetTitle(";#Delta#phi/#pi;1/N_{ t}  dN_{ h }/d#Delta#phi");
                    
                    mgr->Draw("A");
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

                    can.Update();
                    mgr->GetHistogram()->GetXaxis()->SetLimits(0.5,1.);
                    mgr->GetHistogram()->GetXaxis()->SetNdivisions(505);

                    double bigymin = ymin;
                    double bigymax = ymax;

                    biglim[ipy][0] = min(bigymin,biglim[ipy][0]);
                    biglim[ipy][1] = max(bigymax,biglim[ipy][1]);

                    double ydiff = ymax-ymin;
                    
                    ymin = ymin - (ydiff)*.1;
                    ymax = ymax + (ydiff)*.3;
    
                    mgr->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);

                    double xlims[2] = {0.501,.999};
                    TLine* line = new TLine(xlims[0],0.,xlims[1],0.);
                    line->SetLineWidth(1);
                    line->SetLineStyle(2);
                    line->Draw();

                    mgr->Draw();

                    can.Update();
                    can.SaveAs(Form("dphiplots/central_zoom_%s_adj%i_ptt%i_xi%i.pdf",GAT::spec[itype],iadj,iptt,ixi));

                    // draw on big pads too!
                    bigpad[ipx][ipy]->cd();
                    // bigmgr[ipx][ipy]->SetTitle(";#Delta#phi/#pi;1/N_{t}  dN_{ h }/d#Delta#phi");
                    if (ipx == 0 && ipy == 2)
                        bigmgr[ipx][ipy]->SetTitle(";;1/N_{ t}  dN_{ h }/d#Delta#phi");
                    else if (ipx == 3 && ipy == 0)
                        bigmgr[ipx][ipy]->SetTitle(";#Delta#phi/#pi;");
                    else
                        bigmgr[ipx][ipy]->SetTitle("");

                    //ticks
                    if (ipx != 3)
                        bigpad[ipx][ipy]->SetTicky(2);
                    if (ipy != 2)
                        bigpad[ipx][ipy]->SetTickx(2);

                    
                    bigmgr[ipx][ipy]->Draw("A");

                    bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLimits(xlims[0],xlims[1]);
                    bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetNdivisions(505);
                    
                    bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleFont(43);
                    bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelFont(43);
                    if (ipy == 0)
                    {
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleOffset(2.75);
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(70);
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(60);
                    }
                    else
                    {
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(0);
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(0);
                    }

                    bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleFont(43);
                    bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelFont(43);
                    if (ipx == 0)
                    {
                            
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleOffset(6.0);
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelOffset(.01);
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(70);
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(60);
                    }
                    else
                    {
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(0);
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(0);
                    }

                    if (ipy == 2)
                    {                   
                        double rt = bigpad[ipx][ipy]->GetRightMargin();
                        double lt = bigpad[ipx][ipy]->GetLeftMargin();
                        double xpos = ((1.-rt)+lt)/2.;

                        TLatex* txt = new TLatex(xpos,.89,Form("%.1f #leq p_{Tt} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
                        
                        txt->SetNDC();
                        txt->SetTextAlign(22);
                        txt->SetTextFont(43);
                        txt->SetTextSize(60);
                        txt->Draw();
                    }

                    if (ipx == 3)
                    {
                        TLatex* txt;
                        double tp = bigpad[ipx][ipy]->GetTopMargin();
                        double bt = bigpad[ipx][ipy]->GetBottomMargin();
                        double ypos = ((1.-tp)+bt)/2.;

                        txt = new TLatex(.92,ypos,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]));

                        txt->SetNDC();
                        txt->SetTextAlign(22);
                        txt->SetTextAngle(90);
                        txt->SetTextFont(43);
                        txt->SetTextSize(60);
                        txt->Draw();

                    }


                    line->Draw();
                    // bigtxtmom[ipx][ipy]->Draw();
                    // bigtxtcent[ipx][ipy]->Draw();

                    if (ipx == 0 && ipy == 2)
                    {
                        bigleg[ipx][ipy]->Draw();
                        // bigtxtcent[ipx][ipy]->Draw();
                    }
                }//xi
            }//ptt

            for (int ipy = 0; ipy < 3; ipy++)
            {
                double ydiff = biglim[ipy][1]-biglim[ipy][0];
                    
                biglim[ipy][0] = biglim[ipy][0] - (ydiff)*.05;
                biglim[ipy][1] = biglim[ipy][1] + (ydiff)*.05;
                    
                for (int ipx = 0; ipx < 4; ipx++)
                {
                    bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetRangeUser(biglim[ipy][0],biglim[ipy][1]);
                    bigpad[ipx][ipy]->Modified();
                    bigpad[ipx][ipy]->Update();
                }
            }
                    
            bigcan->SaveAs(Form("dphiplots/%s.pdf",bigcan->GetName()));

        }
    }

    //in and out same plot. all ptt one plot. zoomed. copied from central and mid central same plot zoomed
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int icent = 0; icent <= ncent; icent++)
            {
                TCanvas* bigcan = new TCanvas(Form("bigcan_zoom_inout_%s_adj%i_cent%i",GAT::spec[itype],iadj,icent),"",3000,3500);
                TPad* bigpad[4][3];// 4 pt columns, 3 xi rows
                TMultiGraph* bigmgr[4][3];
                double biglim[3][2];// x index not needed
                TLegend* bigleg[4][3];      

                // set up canvas
                double xpos = 0.;
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {

                    int ipx = iptt;//pad x pos

                    //more widht options shown in preivious plottings

                    double left = 0.;
                    if (ipx == 0) left = .25;

                    double right = 0.;
                    if (ipx == 3) right = .15;

                    //(1-left-right)w

                    double wd = .85/2.875;
                    if (ipx == 1 || ipx == 2) wd *= .75;
                    if (ipx == 3) wd *= .75/.85;

                    double ypos = 0.;               
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        int ipy = ixi;//pad y pos

                        //more pad height options shown in previous plottings

                        //set pad margins and create
                        double top = 0.;
                        if (ipy == 2) top = .15;

                        double bot = 0.;
                        if (ipy == 0) bot = .15;


                        double ht = .85/2.85;
                        if (ipy != 1) ht *= 1./.85;


                        bigpad[ipx][ipy] = new TPad(Form("bigpad%i%i",ipx,ipy),"",xpos,ypos,xpos+wd,ypos+ht);

                        bigmgr[ipx][ipy] = new TMultiGraph();//make mgr while were at it

                        // initialize y range limits (unecessarily repeated)
                        biglim[ipy][0] = numeric_limits<double>::max();
                        biglim[ipy][1] = numeric_limits<double>::min();

                        bigpad[ipx][ipy]->SetTopMargin(top);
                        bigpad[ipx][ipy]->SetBottomMargin(bot);
                        bigpad[ipx][ipy]->SetRightMargin(right);
                        bigpad[ipx][ipy]->SetLeftMargin(left);
                            
                        bigpad[ipx][ipy]->Draw();


                        double txtx = .5;
                        if (ipx == 0 && ipy == 2) txtx = .40;              

                        bigleg[ipx][ipy] = new TLegend(.95-right-.65*(.95-right-left),.98-top-.2*(.98-top-bot),.95-right,.98-top);
                        bigleg[ipx][ipy]->SetBorderSize(0);
                        bigleg[ipx][ipy]->SetTextFont(43);
                        bigleg[ipx][ipy]->SetTextSize(60);

                        ypos += ht;
                    }
                    xpos += wd;
                }


                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    int ipx = iptt;//pad x pos
                
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        int ipy = ixi;//pad y pos
                    
                        TCanvas can ("can","",2500,2000);//2x3
                        TMultiGraph* mgr = new TMultiGraph;

                        TLegend leg (0.65,0.75,0.9,0.9);
                        leg.SetTextSize(0.04);

                        double ymin = INT_MAX;
                        double ymax = INT_MIN;

                        for (int idpsi = 0; idpsi < GAT::ndpsi; idpsi++)
                        {
                            TGraphErrors* tgr = tgrdphi[iadj][itype][icent][idpsi][iptt][ixi];
                            tgr->SetMarkerColor(TColor::GetColorPalette(440/3*idpsi));
                            tgr->SetLineColor(TColor::GetColorPalette(440/3*idpsi));
                            tgr->SetMarkerStyle(markstyle[0]);

                            TGraphErrors* tgr_sys = tgrdphi_sys[iadj][itype][icent][idpsi][iptt][ixi];
                            tgr_sys->SetFillColorAlpha(TColor::GetColorPalette(440/3*idpsi),.5);                        

                            mgr->Add(tgr_sys,"2");
                            mgr->Add(tgr,"p");
                            TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                            tgrclone->SetMarkerSize(0);                 
                            mgr->Add(tgrclone,"p");

                            leg.AddEntry(tgr,Form("%s",strplane[idpsi].c_str()),"p");

                            bigmgr[ipx][ipy]->Add(tgr_sys,"2");
                            bigmgr[ipx][ipy]->Add(tgr,"p");
                            bigmgr[ipx][ipy]->Add(tgrclone,"p");
                            bigleg[ipx][ipy]->AddEntry(tgr,Form("%s",strplane[idpsi].c_str()),"p");

                            int n = tgr->GetN();
                            for (int i = n/2; i < n; i++)
                            {
                                ymin = min(ymin,tgr->GetY()[i] - tgr->GetEY()[i]);
                            
                                ymin = min(ymin,tgr_sys->GetY()[i] - tgr_sys->GetEY()[i]);

                            
                                ymax = max(ymax,tgr->GetY()[i] + tgr->GetEY()[i]);
                            
                                ymax = max(ymax,tgr_sys->GetY()[i] + tgr_sys->GetEY()[i]);
                            }
                        }

                        mgr->SetTitle(";#Delta#phi/#pi;1/N_{ t}  dN_{ h }/d#Delta#phi");
                    
                        mgr->Draw("A");
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

                        can.Update();
                        mgr->GetHistogram()->GetXaxis()->SetLimits(0.5,1.);
                        mgr->GetHistogram()->GetXaxis()->SetNdivisions(505);

                        double bigymin = ymin;
                        double bigymax = ymax;

                        biglim[ipy][0] = min(bigymin,biglim[ipy][0]);
                        biglim[ipy][1] = max(bigymax,biglim[ipy][1]);

                        double ydiff = ymax-ymin;
                    
                        ymin = ymin - (ydiff)*.1;
                        ymax = ymax + (ydiff)*.3;
    
                        mgr->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);

                        double xlims[2] = {0.501,.999};
                        TLine* line = new TLine(xlims[0],0.,xlims[1],0.);
                        line->SetLineWidth(1);
                        line->SetLineStyle(2);
                        line->Draw();

                        mgr->Draw();

                        can.Update();
                        can.SaveAs(Form("dphiplots/inout_zoom_%s_adj%i_cent%i_ptt%i_xi%i.pdf",GAT::spec[itype],iadj,icent,iptt,ixi));

                        // draw on big pads too!
                        bigpad[ipx][ipy]->cd();
                        // bigmgr[ipx][ipy]->SetTitle(";#Delta#phi/#pi;1/N_{t}  dN_{ h }/d#Delta#phi");
                        if (ipx == 0 && ipy == 2)
                            bigmgr[ipx][ipy]->SetTitle(";;1/N_{ t}  dN_{ h }/d#Delta#phi");
                        else if (ipx == 3 && ipy == 0)
                            bigmgr[ipx][ipy]->SetTitle(";#Delta#phi/#pi;");
                        else
                            bigmgr[ipx][ipy]->SetTitle("");

                        //ticks
                        if (ipx != 3)
                            bigpad[ipx][ipy]->SetTicky(2);
                        if (ipy != 2)
                            bigpad[ipx][ipy]->SetTickx(2);

                    
                        bigmgr[ipx][ipy]->Draw("A");

                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLimits(xlims[0],xlims[1]);
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetNdivisions(505);
                    
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleFont(43);
                        bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelFont(43);
                        if (ipy == 0)
                        {
                            bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleOffset(2.75);
                            bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(70);
                            bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(60);
                        }
                        else
                        {
                            bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetTitleSize(0);
                            bigmgr[ipx][ipy]->GetHistogram()->GetXaxis()->SetLabelSize(0);
                        }

                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleFont(43);
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelFont(43);
                        if (ipx == 0)
                        {
                            
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleOffset(6.0);
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelOffset(.01);
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(70);
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(60);
                        }
                        else
                        {
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetTitleSize(0);
                            bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetLabelSize(0);
                        }

                        if (ipy == 2)
                        {                       
                            double rt = bigpad[ipx][ipy]->GetRightMargin();
                            double lt = bigpad[ipx][ipy]->GetLeftMargin();
                            double xpos = ((1.-rt)+lt)/2.;

                            TLatex* txt = new TLatex(xpos,.89,Form("%.1f #leq p_{Tt} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
                        
                            txt->SetNDC();
                            txt->SetTextAlign(22);
                            txt->SetTextFont(43);
                            txt->SetTextSize(60);
                            txt->Draw();
                        }

                        if (ipx == 3)
                        {
                            TLatex* txt;
                            double tp = bigpad[ipx][ipy]->GetTopMargin();
                            double bt = bigpad[ipx][ipy]->GetBottomMargin();
                            double ypos = ((1.-tp)+bt)/2.;

                            txt = new TLatex(.92,ypos,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]));

                            txt->SetNDC();
                            txt->SetTextAlign(22);
                            txt->SetTextAngle(90);
                            txt->SetTextFont(43);
                            txt->SetTextSize(60);
                            txt->Draw();

                        }


                        line->Draw();
                        // bigtxtmom[ipx][ipy]->Draw();
                        // bigtxtcent[ipx][ipy]->Draw();

                        if (ipx == 0 && ipy == 2)
                        {
                            bigleg[ipx][ipy]->Draw();
                            // bigtxtcent[ipx][ipy]->Draw();
                        }
                    }//xi
                }//ptt

                for (int ipy = 0; ipy < 3; ipy++)
                {
                    double ydiff = biglim[ipy][1]-biglim[ipy][0];
                    
                    biglim[ipy][0] = biglim[ipy][0] - (ydiff)*.05;
                    biglim[ipy][1] = biglim[ipy][1] + (ydiff)*.05;
                    
                    for (int ipx = 0; ipx < 4; ipx++)
                    {
                        bigmgr[ipx][ipy]->GetHistogram()->GetYaxis()->SetRangeUser(biglim[ipy][0],biglim[ipy][1]);
                        bigpad[ipx][ipy]->Modified();
                        bigpad[ipx][ipy]->Update();
                    }
                }
                    
                bigcan->SaveAs(Form("dphiplots/%s.pdf",bigcan->GetName()));
            }

        }
    }


    //ppg comparison
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            //ppg comparison to direct photons
            if (itype != GAT::DIR) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    for (int icent = 0; icent <= ncent; icent++)
                    {
                        TCanvas can ("can","",2500,2000);//2x3
                        TMultiGraph* mgr = new TMultiGraph;
                        TLegend leg (0.75,0.75,0.9,0.9);
                        leg.SetTextSize(0.04);

                        //MY DATA
                        TGraphErrors* tgr = tgrdphi[iadj][itype][icent][ALL][iptt][ixi];
                        tgr->SetMarkerColor(TColor::GetColorPalette(440/3*0));
                        tgr->SetLineColor(TColor::GetColorPalette(440/3*0));
                        tgr->SetMarkerStyle(markstyle[0]);

                        mgr->Add(tgr,"p");
                        TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                        tgrclone->SetMarkerSize(0);                     
                        mgr->Add(tgrclone,"p");

                        if (icent < ncent)
                            leg.AddEntry(tgr,Form("%i-%i%%",limcent[icent],limcent[icent+1]),"p");
                        else
                            leg.AddEntry(tgr,Form("%i-%i%%",limcent[0],limcent[ncent]),"p");

                        //PPG210
                        TGraphErrors* tgrppg = new TGraphErrors(nbinppg,vx.data(),vyxim[ixi].data(),0,veyxim[ixi].data());
                        tgrppg->SetMarkerColor(TColor::GetColorPalette(440/3*1));
                        tgrppg->SetLineColor(TColor::GetColorPalette(440/3*1));
                        tgrppg->SetMarkerStyle(markstyle[1]);

                        mgr->Add(tgrppg,"p");
                        TGraphErrors* tgrppgclone = (TGraphErrors*) tgrppg->Clone();
                        tgrppgclone->SetMarkerSize(0);                  
                        mgr->Add(tgrppgclone,"p");

                        leg.AddEntry(tgrppg,"PPG210 0-40%","p");
                        
                        mgr->SetTitle(";#Delta#phi/#pi;1/N_{ t}  dN_{ h }/d#Delta#phi");
                        mgr->Draw("A");
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

                        double ymin = min(0.,mgr->GetHistogram()->GetMinimum());
                        double ymax = max(0.,mgr->GetHistogram()->GetMaximum());
                        double ydiff = ymax-ymin;
                        ymin = ymin - (ydiff)*.1;
                        ymax = ymax + (ydiff)*.3;
    
                        mgr->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);
                        mgr->GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                        TLine line(0.,0.,1.,0.);
                        line.SetLineWidth(2);
                        line.SetLineStyle(2);
                        line.Draw();

                        mgr->Draw();

                        can.Update();
                        can.SaveAs(Form("dphiplots/ppg_%s_adj%i_cent%i_ptt%i_xi%i.pdf",GAT::spec[itype],iadj,icent,iptt,ixi));
                    }
                }
            }
        }
    }

    //-------------------- v3 shenanigans --------------------
    //open v3 tgraph input
    const char* inname_v3 = "v3_1020.root";
    TFile* input_v3 = new TFile(inname_v3,"READ");
    if (!input_v3 || input_v3->IsZombie())
    {
        cerr<<"Could not open "<<inname_v3<<endl;
        return;
    }
    const char* name_tgrv3 = "v3_1020";
    TGraph* tgrv3 = (TGraph*) gROOT->FindObject(name_tgrv3);

    double v3t[GAT::nptt] {};
    double v3terr[GAT::nptt] {};
    double v3h[GAT::nptt][GAT::nxi] {};
    double v3herr[GAT::nptt][GAT::nxi] {};

    const double v3ptmax = 7.;//v3 tgraph ends here. use constant instead of extrapolating down
    for (int iptt = 0; iptt < GAT::nptt; iptt++)
    {
        double ptt = (GAT::limptt[iptt] + GAT::limptt[iptt+1])/2.;
        v3t[iptt] = (ptt < v3ptmax) ? tgrv3->Eval(ptt,nullptr,"S") : tgrv3->Eval(v3ptmax,nullptr,"S");
        
        //set error to rough max difference         
        for (int i = 0; i <= 1000; i++)//1000 seems fine..
        {
            double ptttmp = GAT::limptt[iptt] + (GAT::limptt[iptt+1]-GAT::limptt[iptt])*i/1000.;
            double v3tmp = (ptttmp < v3ptmax) ? tgrv3->Eval(ptttmp,nullptr,"S") : tgrv3->Eval(v3ptmax,nullptr,"S");
            v3terr[iptt] = max(v3terr[iptt],fabs(v3t[iptt]-v3tmp));
        }
        v3terr[iptt]*=1.1;//include 10% which is about error from ppg173


        // v3terr[iptt] = max(fabs(v3t[iptt]-v3t0),fabs(v3t[iptt]-v3t1));

        for (int ixi = 0; ixi < GAT::nxi; ixi++)
        {
            //we will use pta limits associated with xi limits of the lower limit ptt
            double pta0 = GAT::limptt[iptt+1]/exp(GAT::limxi[ixi]);
            double pta1 = GAT::limptt[iptt]/exp(GAT::limxi[ixi+1]);
            double pta = (pta0+pta1)/2.;

            v3h[iptt][ixi] = (pta < v3ptmax) ? tgrv3->Eval(pta,nullptr,"S") : tgrv3->Eval(v3ptmax,nullptr,"S");
            
            //set error to rough max difference     
            for (int i = 0; i <= 1000; i++)//1000 seems fine..
            {
                double ptatmp = pta0 + (pta1-pta0)*i/1000.;
                double v3tmp = (ptatmp < v3ptmax) ? tgrv3->Eval(ptatmp,nullptr,"S") : tgrv3->Eval(v3ptmax,nullptr,"S");
                v3herr[iptt][ixi] = max(v3herr[iptt][ixi],fabs(v3h[iptt][ixi]-v3tmp));
            }
            v3herr[iptt][ixi]*=1.1;//include 10% which is about error from ppg173

            
            // v3herr[iptt][ixi] = max(fabs(v3h[iptt][ixi]-v3h0),fabs(v3h[iptt][ixi]-v3h1));
        }
    }

    
    TFile* input_k = new TFile(inname_k,"READ");
    if (!input_k || input_k->IsZombie())
    {
        cerr<<"Could not open "<<inname_k<<endl;
        return;
    }
    
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int icent = 0; icent <= ncent; icent++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    ostringstream sufg;
                    sufg<<"adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt;

                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        ostringstream sufh;
                        sufh<<sufg.str()<<"_xi"<<ixi;

                        TGraphErrors* ag = (TGraphErrors*) gROOT->FindObject(Form("ag_Np_%s",sufg.str().c_str()));
                        TGraphErrors* ah = (TGraphErrors*) gROOT->FindObject(Form("ah_Np_%s",sufh.str().c_str()));
                    
                        TCanvas can ("can","",2500,2000);//2x3
                        can.SetLeftMargin(.15);
                        

                        TGraphErrors* tgr = tgrdphi[iadj][itype][icent][ALL][iptt][ixi];
                        tgr->SetMarkerColor(TColor::GetColorPalette(440/3*0));
                        tgr->SetLineColor(TColor::GetColorPalette(440/3*0));
                        tgr->SetMarkerStyle(markstyle[0]);
                        gStyle->SetEndErrorSize(30);
                        tgr->SetMarkerSize(5);

                        tgr->GetYaxis()->SetTitleSize(.05);
                        tgr->GetYaxis()->SetLabelSize(.045);
                        tgr->GetYaxis()->SetTitleOffset(1.4);
                        tgr->GetXaxis()->SetTitleSize(.05);
                        tgr->GetXaxis()->SetLabelSize(.045);
                        
                        tgr->SetTitle(";#Delta#phi/#pi;1/N_{ t}  dN_{ h }/d#Delta#phi");
                        tgr->Draw("AP");

                        TGraphErrors* tgrclone = (TGraphErrors*) tgr->Clone();
                        tgrclone->SetMarkerSize(0);                     
                        tgrclone->Draw("P");


                        TLatex txtmom (0.525,.85,Form("%.1f #leq p_{T} < %.1f #times %.1f #leq #xi < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                        txtmom.SetTextAlign(22);
                        txtmom.SetNDC();
                        txtmom.SetTextSize(0.045);
                        txtmom.Draw();


                        int c0,c1;
                        if (icent != ncent)
                        {
                            c0 = limcent[icent];
                            c1 = limcent[icent+1];
                        }
                        else
                        {
                            c0 = limcent[0];
                            c1 = limcent[ncent];
                        }
                        TLatex txtcent (0.525,.8,Form("Au+Au %i-%i%% 200GeV",c0,c1));
                        txtcent.SetTextAlign(22);
                        txtcent.SetNDC();
                        txtcent.SetTextSize(0.045);
                        txtcent.Draw();

                        double ymin = min(0.,tgr->GetHistogram()->GetMinimum());
                        double ymax = max(0.,tgr->GetHistogram()->GetMaximum());
                        double ydiff = ymax-ymin;
                        ymin = ymin - (ydiff)*.05;
                        ymax = ymax + (ydiff)*.15;
    
                        tgr->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax);
                        tgr->GetHistogram()->GetXaxis()->SetLimits(0.,1.);

                        TLine line(0.,0.,1.,0.);
                        line.SetLineWidth(2);
                        line.SetLineStyle(2);
                        line.Draw();

                        can.Update();
                        can.SaveAs(Form("dphiplots/single_%s_adj%i_cent%i_ptt%i_xi%i.pdf",GAT::spec[itype],iadj,icent,iptt,ixi));

                        //draw v3 contribution
                        //per trig normalization removes 2piag. just need ah.
                        //sum ah cent bins
                        double ahsum = 0.;
                        int nbin = 0;
                        for (int jcent = 0; jcent < GAT::ncent; jcent++)
                        {
                            if (GAT::limcent[jcent] >= c0 && GAT::limcent[jcent+1] <= c1)
                            {
                                ahsum += ah->GetY()[jcent];
                                nbin++;
                            }
                        }
                        double ahavg = ahsum/nbin;
                        {
                            const Int_t n = 24;
                            Double_t x[n], y[n],ymin[n], ymax[n];
                            for (int i=0;i<n;i++) {
                                x[i] = i/(n-1.);

                                // ymax[i] = ahavg*pow(v3+.01,2.)*cos(3.*x[i]*PI);
                                // ymin[i] = ahavg*pow(v3-.01,2.)*cos(3.*x[i]*PI);
                                // y[i] = ahavg*pow(v3,2.)*cos(3.*x[i]*PI);

                                ymax[i] = ahavg*(v3t[iptt]+v3terr[iptt])*(v3h[iptt][ixi]+v3herr[iptt][ixi])*cos(3.*x[i]*PI);
                                ymin[i] = ahavg*(v3t[iptt]-v3terr[iptt])*(v3h[iptt][ixi]-v3herr[iptt][ixi])*cos(3.*x[i]*PI);
                                y[i] = ahavg*v3t[iptt]*v3h[iptt][ixi]*cos(3.*x[i]*PI);
                            }
                            TGraph *grmin = new TGraph(n,x,ymin);
                            TGraph *grmax = new TGraph(n,x,ymax);
                            TGraph *gr    = new TGraph(n,x,y);
                            TGraph *grshade = new TGraph(2*n);
                            for (int i=0;i<n;i++) {
                                grshade->SetPoint(i,x[i],ymax[i]);
                                grshade->SetPoint(n+i,x[n-i-1],ymin[n-i-1]);
                            }
                            grshade->SetFillStyle(3013);
                            grshade->SetFillColor(16);
                            grshade->Draw("f");

                            grmin->SetLineColor(16);
                            grmin->SetLineColor(TColor::GetColorPalette(440/3*1));
                            grmin->Draw("l");

                            grmax->SetLineColor(16);
                            grmax->SetLineColor(TColor::GetColorPalette(440/3*1));
                            grmax->Draw("l");

                            gr->SetLineColor(TColor::GetColorPalette(440/3*1));
                            gr->SetLineWidth(2);
                            gr->Draw("l");
                        }

                        //redraw points over v3 graphs
                        tgr->GetYaxis()->SetTitleSize(.045);
                        tgr->GetYaxis()->SetTitleOffset(1.1);
                        tgr->GetXaxis()->SetTitleSize(.045);
                        tgr->GetXaxis()->SetTitleOffset(0.9);
                        tgr->Draw("P");
                        tgrclone->Draw("P");

                        can.Update();
                        can.SaveAs(Form("dphiplots/v3_%s_adj%i_cent%i_ptt%i_xi%i.pdf",GAT::spec[itype],iadj,icent,iptt,ixi));
                    }
                }
            }
        }
    }
    cout<<"Finished drawing!"<<endl;

    cout<<"Writing some output to "<<outname<<endl;

    //output file
    TFile* output = new TFile(outname,"RECREATE");

    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int icent = 0; icent < ncent+1; icent++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    tgrixiv2j[iadj][itype][iptt][icent]->Write();
                    tgrixiv2j_syserr[iadj][itype][iptt][icent]->Write();
                }

                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    tgripttv2j[iadj][itype][icent][ixi]->Write();
                    tgripttv2j_syserr[iadj][itype][icent][ixi]->Write();
                }
            }
        }
    }

    output->Close();
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
        cout<<"adj"<<info.iadj<<"_"<<GAT::spec[info.itype]<<"_ptt"<<info.iptt<<"_cent"<<info.icent<<"_dpsi"<<info.idpsi<<"_xi"<<info.ixi<<": fe0 = "<<info.flate[0]<<"  fe1 = "<<info.flate[1]<<" ce0 = "<<info.cose[0]<<"  ce1 = "<<info.cose[1]<<endl;
    }

    /* some notes on error. how to get from stored values to desired error
       In-Plane:
       info.flate[0] = kf0err;// = sqrt(kf0errmsum[icent])/(2.*ntrg)
       info.flate[1] = kf1err;// = sqrt(kf1errmsum[icent])/ntrg
       info.cose[0] = kc0err;// = sqrt(kc0errmsum[icent])/(2.*ntrg)
       info.cose[1] = kc1err;// = sqrt(kc1errmsum[icent])/ntrg

       fe_in = sqrt(kf0err^2 + kf1err^2);
       //check: sqrt(kf0errsum/4 + kf1errsum)/ntrg : CORRECT!
       ce_in = sqrt(kc0errsum/4 + kc1errsum)/ntrg;
       ce_in = sqrt(kc0err^2 + kc1err^2);
       
       All-Plane
       fe_all = sqrt(kf0errsum)/ntrg;
       info.flat[0] = kf0err = sqrt(kf0errsum)/2.;//note factor 2 above doesnt apply
       fe_all = kf0err;
       
       OK. For all cases just do sqrt(kf0e^2 + kf1e^2) and 
       sqrt(kc0e^2 + kc1e^2) 

     */
}

double GetError(TGraphErrors* tgr, TF1* f1, double xmin)
{
    double num = 0., den = 0.;

    for (int i = 0; i < tgr->GetN(); i++)
    {
        double x = tgr->GetX()[i];

        if (x < xmin) continue;
        
        double y = tgr->GetY()[i];
        double ey = tgr->GetEY()[i];

        double fval = f1->Eval(x);

        num += pow((y-fval)/ey,2.);
        den += pow(1./ey,2.);
    }

    return sqrt(num/den);
}

void DrawError(TF1* f1, double err, double xmin, int kColor = kRed, double alpha = .15, int ipar = 0)
{
    double par = f1->GetParameter(ipar);
    
    TF1* f1min = (TF1*) f1->Clone("f1min");
    f1min->SetParameter(ipar,par-err);

    TF1* f1max = (TF1*) f1->Clone("f1max");
    f1max->SetParameter(ipar,par+err);
        
    const Int_t n = 24;
    Double_t x[n], ymin[n], ymax[n];
    for (int i=0;i<n;i++) {
        x[i] = xmin + i*(1.-xmin)/(n-1.);

        ymax[i] = f1max->Eval(x[i]);
        ymin[i] = f1min->Eval(x[i]);

    }
    TGraph *grmin = new TGraph(n,x,ymin);
    TGraph *grmax = new TGraph(n,x,ymax);
    TGraph *grshade = new TGraph(2*n);
    for (int i=0;i<n;i++) {
        grshade->SetPoint(i,x[i],ymax[i]);
        grshade->SetPoint(n+i,x[n-i-1],ymin[n-i-1]);
    }
    

    grshade->SetFillStyle(1001);
    grshade->SetFillColorAlpha(kColor,alpha);
        
    grshade->Draw("f");

    grmin->SetLineColor(kColor);
    grmin->SetLineWidth(1);
    grmin->Draw("l");

    grmax->SetLineColor(kColor);
    grmax->SetLineWidth(1);
    grmax->Draw("l");

    f1->Draw("same");
}

