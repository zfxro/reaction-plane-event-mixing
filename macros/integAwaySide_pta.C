enum RPO {IN,OUT,ALL};

const float PI = TMath::Pi();
const int nadj = 2;
const int nmix = 2;
const int ntype = 3;

const int ncent = 8;
const int nsumcent = 2;

const int markstyle[] = {20,24,21,25,22,26,23,32,34,47};

namespace GAT {
    #include "gamAnaTree.h"
}

//TGraphs:
//central to midcentral ratio as function of pta
TGraphErrors* tgriptayield[nadj][ntype][GAT::nptt][nsumcent][GAT::ndpsi+1] = {};
TGraphErrors* tgripttyield[nadj][ntype][GAT::npta][nsumcent][GAT::ndpsi+1] = {};
TGraphErrors* tgriptacentrat[nadj][ntype][GAT::nptt] = {};
//v2jet vs pta
TGraphErrors* tgriptav2j[nadj][ntype][GAT::nptt][nsumcent] = {};
TGraphErrors* tgricentv2j[nadj][ntype][GAT::nptt][GAT::npta] = {};

double yield[nadj][ntype][GAT::nptt][ncent][GAT::ndpsi+1][GAT::npta] = {};
double yieldsum[nadj][ntype][GAT::nptt][nsumcent][GAT::ndpsi+1][GAT::npta] = {};

double yielderr[nadj][ntype][GAT::nptt][ncent][GAT::ndpsi+1][GAT::npta] = {};
double yielderrsum[nadj][ntype][GAT::nptt][nsumcent][GAT::ndpsi+1][GAT::npta] = {};

double v2j[nadj][ntype][GAT::nptt][ncent][GAT::npta] = {};
double v2jsum[nadj][ntype][GAT::nptt][nsumcent][GAT::npta] = {};

double v2jerr[nadj][ntype][GAT::nptt][ncent][GAT::npta] = {};
double v2jerrsum[nadj][ntype][GAT::nptt][nsumcent][GAT::npta] = {};

TGraphErrors* tgrsig = nullptr;
TGraphErrors* tgrexpsig = nullptr;
TGraphErrors* tgrksig[nadj][ntype][GAT::nptt] = {};

TFile* input = nullptr;
TFile* inputsig = nullptr;

int ReadFile(const char* name, TFile*& file);
void GetYield(int iadj, int itype, int iptt, int idpsi, int ipta, int icent, bool issum = false);

void integAwaySide_pta(
    const char* inname = "weighted_bg_subtracted_pta.root",
    const char* innamesig = "kComponents/sigmak.root",
    const char* outname = "integawaysideyield_pta.root"
    )
{
    TFile* output = new TFile(outname,"RECREATE");
    
    if (ReadFile(inname,input)<0) return;
    input->GetObject("npartsigma",tgrsig);

    if (ReadFile(innamesig,inputsig)<0) return;
    inputsig->GetObject("expsig",tgrexpsig);

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        ostringstream osssuf;
        osssuf<<"dphi_";
        if (iadj)
            osssuf<<"adj_";
        
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype != GAT::PI0) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {

                inputsig->GetObject(Form("ksigma_Np_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt),tgrksig[iadj][itype][iptt]);
                
                for (int ipta = 0; ipta < GAT::npta; ipta++)
                {                   
                    for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                    {                   
                        for (int icent = 0; icent < ncent; icent++)
                        {
                            GetYield(iadj,itype,iptt,idpsi,ipta,icent);

                            if (idpsi!=ALL) continue;

                            double iny = yield[iadj][itype][iptt][icent][IN][ipta];
                            double inye = yielderr[iadj][itype][iptt][icent][IN][ipta];

                            double outy = yield[iadj][itype][iptt][icent][OUT][ipta];
                            double outye = yielderr[iadj][itype][iptt][icent][OUT][ipta];

                            double diff = iny - outy;
                            double error = sqrt(inye*inye + outye*outye);

                            double sigma = tgrsig->GetY()[icent];

                            diff/=sigma;
                            error/=sigma;

                            double ally = yield[iadj][itype][iptt][icent][ALL][ipta];
                            double allye = yielderr[iadj][itype][iptt][icent][ALL][ipta];
                            
                            diff = diff/ally;
                            error = error/ally;
                            error = sqrt( pow(error,2.) + pow(allye*diff/ally,2.) );

                            v2j[iadj][itype][iptt][icent][ipta] = diff;
                            v2jerr[iadj][itype][iptt][icent][ipta] = error;
                        }

                        for (int icent = 0; icent < nsumcent; icent++)
                        {
                            GetYield(iadj,itype,iptt,idpsi,ipta,icent,true);

                            if (idpsi!=ALL) continue;

                            double iny = yieldsum[iadj][itype][iptt][icent][IN][ipta];
                            double inye = yielderrsum[iadj][itype][iptt][icent][IN][ipta];

                            double outy = yieldsum[iadj][itype][iptt][icent][OUT][ipta];
                            double outye = yielderrsum[iadj][itype][iptt][icent][OUT][ipta];

                            double diff = iny - outy;
                            double error = sqrt(inye*inye + outye*outye);

                            double sigma = tgrexpsig->GetY()[icent];

                            sigma *= tgrksig[iadj][itype][iptt]->GetY()[icent];

                            diff/=sigma;
                            error/=sigma;

                            double ally = yieldsum[iadj][itype][iptt][icent][ALL][ipta];
                            double allye = yielderrsum[iadj][itype][iptt][icent][ALL][ipta];
                            
                            diff = diff/ally;
                            error = error/ally;
                            error = sqrt( pow(error,2.) + pow(allye*diff/ally,2.) );

                            v2jsum[iadj][itype][iptt][icent][ipta] = diff;
                            v2jerrsum[iadj][itype][iptt][icent][ipta] = error;
                        }//icent < nsumcent
                    }//idpsi
                }//ipta
            }//iptt
        }//itype
    }//iadj

    //Create TGraphs of interest
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype != GAT::PI0) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                for (int ipta = 0; ipta < GAT::npta; ipta++)
                {
                    double centyield = yieldsum[iadj][itype][iptt][0][ALL][ipta];
                    double midcentyield = yieldsum[iadj][itype][iptt][1][ALL][ipta];
                    double ratio = centyield/midcentyield;

                    double e0 = yielderrsum[iadj][itype][iptt][0][ALL][ipta];
                    double e1 = yielderrsum[iadj][itype][iptt][1][ALL][ipta];
                    e0 /= midcentyield;
                    e1 = -e1 * ratio / midcentyield;

                    double error = sqrt(e0*e0 + e1*e1);

                    if (!tgriptacentrat[iadj][itype][iptt])
                    {
                        tgriptacentrat[iadj][itype][iptt] = new TGraphErrors();
                        tgriptacentrat[iadj][itype][iptt]->SetName(Form("iptacentrat_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt));
                        tgriptacentrat[iadj][itype][iptt]->SetTitle(";p_{Ta};Y_{0-20%}/Y_{20-60%}");
                    }
                    
                    tgriptacentrat[iadj][itype][iptt]->SetPoint(ipta,(GAT::limpta[ipta+1]+GAT::limpta[ipta])/2.,ratio);
                    tgriptacentrat[iadj][itype][iptt]->SetPointError(ipta,0,error);
                                                               
                    for (int icent = 0; icent < nsumcent; icent++)
                    {
                        for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                        {
                            if (!tgriptayield[iadj][itype][iptt][icent][idpsi])
                            {
                                tgriptayield[iadj][itype][iptt][icent][idpsi] = new TGraphErrors();
                                tgriptayield[iadj][itype][iptt][icent][idpsi]->SetName(Form("iptayield_adj%i_%s_ptt%i_cent%i_dpsi%i",iadj,GAT::spec[itype],iptt,icent,idpsi));
                                tgriptayield[iadj][itype][iptt][icent][idpsi]->SetTitle(";p_{Ta};1/N_{#gamma} dN/dp_{Ta}");
                            }
                            tgriptayield[iadj][itype][iptt][icent][idpsi]->SetPoint(ipta,(GAT::limpta[ipta+1]+GAT::limpta[ipta])/2.,yieldsum[iadj][itype][iptt][icent][idpsi][ipta]);
                            tgriptayield[iadj][itype][iptt][icent][idpsi]->SetPointError(ipta,0,yielderrsum[iadj][itype][iptt][icent][idpsi][ipta]);

                            if (!tgripttyield[iadj][itype][ipta][icent][idpsi])
                            {
                                tgripttyield[iadj][itype][ipta][icent][idpsi] = new TGraphErrors();
                                tgripttyield[iadj][itype][ipta][icent][idpsi]->SetName(Form("ipttyield_adj%i_%s_pta%i_cent%i_dpsi%i",iadj,GAT::spec[itype],ipta,icent,idpsi));
                                tgripttyield[iadj][itype][ipta][icent][idpsi]->SetTitle(";p_{Tt};1/N_{#gamma} dN/dp_{Tt}");
                            }
                            tgripttyield[iadj][itype][ipta][icent][idpsi]->SetPoint(iptt,(GAT::limptt[iptt+1]+GAT::limptt[iptt])/2.,yieldsum[iadj][itype][iptt][icent][idpsi][ipta]);
                            tgripttyield[iadj][itype][ipta][icent][idpsi]->SetPointError(iptt,0,yielderrsum[iadj][itype][iptt][icent][idpsi][ipta]);
                        }

                        if (!tgriptav2j[iadj][itype][iptt][icent])
                        {
                            tgriptav2j[iadj][itype][iptt][icent] = new TGraphErrors();
                            tgriptav2j[iadj][itype][iptt][icent]->SetName(Form("iptav2j_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent));
                            tgriptav2j[iadj][itype][iptt][icent]->SetTitle(";p_{Ta};#pi/(4#sqrt{3}) (1/N_{#gamma} dN_{h}/dp_{Ta})_{in-out}/#kappa_{#sigma}<#sigma>");
                        }                       
                        tgriptav2j[iadj][itype][iptt][icent]->SetPoint(ipta,(GAT::limpta[ipta+1]+GAT::limpta[ipta])/2.,v2jsum[iadj][itype][iptt][icent][ipta]*(PI/(4.*sqrt(3.))));
                        tgriptav2j[iadj][itype][iptt][icent]->SetPointError(ipta,0,v2jerrsum[iadj][itype][iptt][icent][ipta]*(PI/(4.*sqrt(3.))));
                        
                        if (!tgricentv2j[iadj][itype][iptt][ipta])
                        {
                            tgricentv2j[iadj][itype][iptt][ipta] = new TGraphErrors();
                            tgricentv2j[iadj][itype][iptt][ipta]->SetName(Form("icentv2j_adj%i_%s_ptt%i_pta%i",iadj,GAT::spec[itype],iptt,ipta));
                            tgricentv2j[iadj][itype][iptt][ipta]->SetTitle(";icent;(1/N_{#gamma} dN_{h}/dp_{Ta})_{in-out}/#kappa_{#sigma}<#sigma>");
                        }
                        tgricentv2j[iadj][itype][iptt][ipta]->SetPoint(icent,icent,v2jsum[iadj][itype][iptt][icent][ipta]);
                        tgricentv2j[iadj][itype][iptt][ipta]->SetPointError(icent,0,v2jerrsum[iadj][itype][iptt][icent][ipta]);
                    }//icent
                }//ipta
            }//iptt
        }//itype
    }//iadj

    //write tgraphs
    output->cd();
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype != GAT::PI0) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                tgriptacentrat[iadj][itype][iptt]->Write();
                for (int ipta = 0; ipta < GAT::npta; ipta++)
                {
                    tgricentv2j[iadj][itype][iptt][ipta]->Write();
                }
                for (int icent = 0; icent < nsumcent; icent++)
                {
                    tgriptav2j[iadj][itype][iptt][icent]->Write();
                    for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                    {
                        tgriptayield[iadj][itype][iptt][icent][idpsi]->Write();
                        if (iptt==0)
                        {
                            for (int ipta = 0; ipta < GAT::npta; ipta++)
                            {
                                tgripttyield[iadj][itype][ipta][icent][idpsi]->Write();
                            }
                        }
                    }
                }
            }
        }
    }

    cout<<"Output written to "<<output->GetName()<<endl;
    output->Close();
    cout<<"Done!"<<endl;
    return;
}

int ReadFile(const char* name, TFile*& file)
{
    file = new TFile(name,"READ");
    if (!file)
    {
        cerr<<"Could not open "<<name<<endl;
        return -1;
    }
    if (file->IsZombie())
    {
        cerr<<"File "<<name<<" is zombie!"<<endl;
        return -2;
    }
    cout<<"Reading file "<<name<<endl;
    return 0;
}

void GetYield(int iadj, int itype, int iptt, int idpsi, int ipta, int icent, bool issum)
{
            
    ostringstream oss;
    oss<<"dphi";
    if (iadj)
        oss<<"_adj";
    oss<<"_"<<GAT::spec[itype];
    if (!issum)
        oss<<"_mix0";
    oss<<"_ptt"<<iptt<<"_";
    if (issum)
        oss<<"sum";
    oss<<"cent"<<icent<<"_dpsi"<<idpsi<<"_pta"<<ipta;

    TH1F* h1dphi = nullptr;
    input->GetObject(oss.str().c_str(),h1dphi);
    if (!h1dphi)
    {
        cerr<<"Could not find "<<oss.str()<<endl;
    }
    else
    {
        cout<<"Found "<<oss.str()<<endl;
    }

    const int bin0 = 9;
    const double expec = .67;

    double binlowedge = h1dphi->GetBinLowEdge(bin0);
    if (fabs(binlowedge - expec) > .01)
    {
        cerr<<"ERROR: binlowedge = "<<binlowedge<<endl;
        return;
    }

    int nbin = h1dphi->GetNbinsX();
    
    double error = -9999.;
    double y = h1dphi->IntegralAndError(bin0,nbin,error,"width");
    //width option multiplies by bin width.  the pi has been dropped and must be included manually...
    y*=PI;
    error*=PI;

    cout<<"Scaling yield and error by 1/ptabinwidth"<<endl;
    y/= (GAT::limpta[ipta+1]-GAT::limpta[ipta]);
    error /= (GAT::limpta[ipta+1]-GAT::limpta[ipta]);

    if (issum)
    {
        yieldsum[iadj][itype][iptt][icent][idpsi][ipta] = y;
        yielderrsum[iadj][itype][iptt][icent][idpsi][ipta] = error;
    }
    else
    {
        yield[iadj][itype][iptt][icent][idpsi][ipta] = y;
        yielderr[iadj][itype][iptt][icent][idpsi][ipta] = error;
    }
}
