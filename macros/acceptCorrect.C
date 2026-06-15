
// read gamtree output with rp mixing
// also read gametree output without rp mixing
// read expnh values from hadronaccept directory
// calculate acceptance
// apply to histos of gamtree with rp mixing
// save to new file

const int verbosity = 1;

const int ncent_merge = 2;
const int limcent_merge[ncent_merge+1] = {0,20,60};

const float PI = TMath::Pi();

const int nadj = 2;
const int nmix = 2;
const int ntype = 4;

const double detah = 0.7;

int iptamax_acc = 2;

TFile* gamtree_rp;
TFile* gamtree_norp;
TFile* output;

int GetNhExpXi();
int GetNhExpPta();
int GetNoRpHistos();
int GetRpHistos();
void CalcAcc(const char* name_out_raw);

void GetKComp(const char* name, bool);
void CalcCorrect();
void CalcCorrHelper(int iadj, int itype, int iptt, int idpsi, int ixi, bool ismerge);
int ApplyCorrect(bool docor, bool dopsicor);

void ApplyAcc();

namespace GAT {
     #include "gamAnaTree.h"
}

double expnhpta[GAT::ncent][GAT::npta] = {};
TH1F* h1accpta[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1][GAT::npta];
TH1F* h1dphipta[nadj][ntype][nmix][GAT::nptt][GAT::ncent][GAT::ndpsi+1][GAT::npta];

double expnhxi[nadj][ntype][GAT::ncent][GAT::ndpsi+1][GAT::nptt][GAT::nxi] = {};
TH1F* h1acc[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1][GAT::nxi];
TH1F* h1dphi[nadj][ntype][nmix][GAT::nptt][GAT::ncent][GAT::ndpsi+1][GAT::nxi];
TH2F* h2centptt[nadj][ntype][GAT::ndpsi+1];

TGraphErrors* tgrng[nadj][ntype][GAT::nptt][GAT::ndpsi+1];
TGraphErrors* tgreff[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];

vector<double> corr[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];
vector<double> corre[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];
TGraphErrors* tgrcorr[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];

vector<double> corr_merge[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];
vector<double> corre_merge[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];
TGraphErrors* tgrcorr_merge[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];

vector<double> psimod;
vector<double> psimode;
TGraphErrors* tgrpsimod;

TH1F* hnpart[GAT::ncent];
vector<double> vnpart;
vector<double> vnparte;
int GetNpartDists();

TH1F* hnpart_merge[ncent_merge];
vector<double> vnpart_merge;
vector<double> vnparte_merge;
int GetNpartDistsMerge();

const int npartmin = 2;
const int npartmax = 390;
TGraphErrors* MakeProd(TGraphErrors* tgra, TGraphErrors* tgrae, TGraphErrors* tgrb, TGraphErrors* tgrbe, string name);
void DrawKComps(TGraphErrors* gra, TGraphErrors* grae, string laba, TGraphErrors* grb, TGraphErrors* grbe, string labb, TGraphErrors* grp, string labp, TH1F* hnpart[GAT::ncent], string kname = "", TGraphErrors* grk = nullptr, TGraphErrors* gra0 = nullptr, TGraphErrors* grae0 = nullptr, string laba0 = "", TGraphErrors* grb0 = nullptr, TGraphErrors* grbe0 = nullptr, string labb0 = "");
void CalcExp(TF1* fa, TGraphErrors* grae, TH1F* hnp, double& expa, double& expae);
    

void acceptCorrect(
    bool dopsimod = true,
    const char* name_out = "gamtree_13466_acceptcorrect.root",
    const char* name_kcomp = "kComponents/kComponents.root",
    const char* name_gamtree_rp = "rootfiles/gamtree_13466_psimix.root",
    const char* name_gamtree_norp = "rootfiles/gamtree_13466_nopsimix.root",
    const char* name_out_raw = "gamtree_13466_rawnh.root"
    )
{
    cout<<"Running acceptCorrect.C"<<endl;
    cout<<"dopsimod set to "<<dopsimod<<endl;
    
    cout<<"ATTENTION: nh exp values should be updated before analysis finalized. see calc_acc_binmode.C and run hadronaccept/wrk/xi/drive_calc_acc_binmode.sh"<<endl;
    
    cout<<"ATTENTION: Acceptance calculated for ipta bin "<<iptamax_acc<<" ("<<GAT::limpta[iptamax_acc]<<" < pta < "<<GAT::limpta[iptamax_acc+1]<<") will be used for all higher bins as well"<<endl;
    if (GAT::limpta[iptamax_acc] != GAT::ptamixmax0 || GAT::limpta[iptamax_acc+1] != GAT::ptamixmax1)
    {
        cerr<<"iptamax_acc not set properly"<<endl;
    }
    cout<<"These bin limits agree with the ptamixmax variables in gamAnaTree.h.  Proceed."<<endl;
    
    
    if (GetNpartDists()<0) return;

    if (GetNpartDistsMerge()<0) return;

    gamtree_rp = new TFile(name_gamtree_rp,"READ");
    if (!gamtree_rp || gamtree_rp->IsZombie())
    {
        cerr<<"Could not open "<<name_gamtree_rp<<endl;
        return;
    }

    gamtree_norp = new TFile(name_gamtree_norp,"READ");
    if (!gamtree_norp || gamtree_norp->IsZombie())
    {
        cerr<<"Could not open "<<name_gamtree_norp<<endl;
        return;
    }

    output = new TFile(name_out,"RECREATE");
    
    if (GetNhExpXi() < 0) return;

    GetNhExpPta();

    gamtree_norp->cd();
    cout<<"Changed directory to "<<gamtree_norp->GetName()<<"  gDirectory->GetFile()->GetName() = "<<gDirectory->GetFile()->GetName()<<endl;
    if (GetNoRpHistos() < 0) return;

    gamtree_rp->cd();
    cout<<"Changed directory to "<<gamtree_rp->GetName()<<"  gDirectory->GetFile()->GetName() = "<<gDirectory->GetFile()->GetName()<<endl;
    if (GetRpHistos() < 0) return;
    
    CalcAcc(name_out_raw);

    GetKComp(name_kcomp,dopsimod);
    CalcCorrect();
    if (ApplyCorrect(false,dopsimod) < 0)
    {
        cout<<"Issue with ApplyCorrect(). Exiting."<<endl;
        return;
    }
    
    ApplyAcc();

    cout<<"Writing to "<<output->GetName()<<endl;
    output->Write();
    output->Close();
    cout<<"Done"<<endl;
}

int GetNhExpXi()
{
    cout<<"Getting NhExpXi"<<endl;

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            for (int icent = 0; icent < GAT::ncent; icent++)
            {
                for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                {
                    for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    {       
                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {
                            expnhxi[iadj][itype][icent][idpsi][iptt][ixi] = -9999.;
                        }
                    }
                }
            }
        }
    }

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int icent = 0; icent < GAT::ncent; icent++)
            {
                for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                {
                    ostringstream ossdir;
                    ossdir << "adj" << iadj;
                    if (itype == 0)
                        ossdir<<"_dir";
                    else if (itype == 1)
                        ossdir<<"_dec";
                    else if (itype == 2)
                        ossdir<<"_pi0";
                    else if (itype == 3)
                        ossdir<<"_inc";
                    ossdir << "_cent" << icent;
                    ossdir << "_dpsi" << idpsi;

                    ostringstream ossin;
                    ossin<<"hadronaccept/wrk/xi/"<<ossdir.str()<<"/integrated_xi_"<<ossdir.str()<<".txt";
                    
                    ifstream ifs(ossin.str(),ifstream::in);
                    if (!ifs.is_open())
                    {
                        cerr<<"Cannot open "<<ossin.str()<<endl;
                        return -1;
                    }

                    if (verbosity > 1)
                        cout<<"Reading "<<ossin.str()<<endl;
                    
                    string line;
                    while (getline(ifs,line))
                    {
                        if (!isdigit(line[0])) continue;                        
                        istringstream iss(line.data());

                        int iptt, ixi;
                        float tmpnh;
                        
                        iss >> iptt >> ixi >> tmpnh;
                        if (verbosity > 1)
                            cout << " " << iptt << " " << ixi << " " << tmpnh << endl;

                        expnhxi[iadj][itype][icent][idpsi][iptt][ixi] = tmpnh*detah*2.*PI;
                    }
                }
            }
        }
    }
    cout<<"Got NhExpXi"<<endl;
    return 1;
}

int GetNhExpPta()
{
    cout<<"Getting NhExpPta"<<endl;

    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        for (int ipta = 0; ipta < GAT::npta; ipta++)
        {           
            expnhpta[icent][ipta] = -9999.;
        }       
    }
            
    for (int icent = 0; icent < GAT::ncent; icent++)
    {

        ostringstream ossin;
        ossin<<"hadronaccept/wrk/pta/integrated_pta_cent"<<icent<<".txt";
                    
        ifstream ifs(ossin.str(),ifstream::in);
        if (!ifs.is_open())
        {
            cerr<<"Cannot open "<<ossin.str()<<endl;
            return -1;
        }

        if (verbosity > 1)
            cout<<"Reading "<<ossin.str()<<endl;
                    
        string line;
        while (getline(ifs,line))
        {
            if (!isdigit(line[0])) continue;                    
            istringstream iss(line.data());

            int ipta;
            float tmpnh;
                        
            iss >> ipta >> tmpnh;
            if (verbosity > 1)
                cout << " " << ipta << " " << tmpnh << endl;

            expnhpta[icent][ipta] = tmpnh*detah*2.*PI;
        }
    }
    cout<<"Got NhExpPta"<<endl;
    return 1;
}
                    
int GetNoRpHistos()
{
    cout<<"Getting No RP Histos"<<endl;
    
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {       
                for (int icent = 0; icent < GAT::ncent; icent++)
                {
                    for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                    {
                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {
                            if (idpsi == GAT::ndpsi)
                            {
                                h1acc[iadj][itype][iptt][icent][idpsi][ixi] = (TH1F*) h1acc[iadj][itype][iptt][icent][0][ixi]->Clone(Form("acc_adj%i_%s_ptt%i_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],iptt,icent,idpsi,ixi));
                                h1acc[iadj][itype][iptt][icent][idpsi][ixi]->Add(h1acc[iadj][itype][iptt][icent][1][ixi]);
                                h1acc[iadj][itype][iptt][icent][idpsi][ixi]->SetDirectory(output);
                                continue;
                            }
                            
                            ostringstream oss;
                            oss<<"dphiacc";
                            if (iadj)
                                oss<<"_adj";
                            oss<<"_"<<GAT::spec[itype];
                            oss<<"_ptt"<<iptt;
                            oss<<"_cent"<<icent;
                            oss<<"_dpsi"<<idpsi;
                            oss<<"_xi"<<ixi;

                            TH1F* h1tmp = (TH1F*) gROOT->FindObject(oss.str().c_str());
                            if (!h1tmp)
                            {
                                cerr<<"Could not find "<<oss.str()<<endl;
                                return -1;
                            }

                            h1acc[iadj][itype][iptt][icent][idpsi][ixi] = (TH1F*) h1tmp->Clone(Form("acc_adj%i_%s_ptt%i_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],iptt,icent,idpsi,ixi));
                            h1acc[iadj][itype][iptt][icent][idpsi][ixi]->SetDirectory(output);
                        }

                        for (int ipta = 0; ipta < GAT::npta; ipta++)
                        {
                            if (idpsi == GAT::ndpsi)
                            {
                                h1accpta[iadj][itype][iptt][icent][idpsi][ipta] = (TH1F*) h1accpta[iadj][itype][iptt][icent][0][ipta]->Clone(Form("acc_adj%i_%s_ptt%i_cent%i_dpsi%i_pta%i",iadj,GAT::spec[itype],iptt,icent,idpsi,ipta));
                                h1accpta[iadj][itype][iptt][icent][idpsi][ipta]->Add(h1accpta[iadj][itype][iptt][icent][1][ipta]);
                                h1accpta[iadj][itype][iptt][icent][idpsi][ipta]->SetDirectory(output);
                                continue;
                            }
                            
                            ostringstream oss;
                            oss<<"dphi";
                            if (iadj)
                                oss<<"_adj";
                            oss<<"_"<<GAT::spec[itype];
                            oss<<"_mix"<<1;
                            oss<<"_ptt"<<iptt;
                            oss<<"_cent"<<icent;
                            oss<<"_dpsi"<<idpsi;
                            oss<<"_pta"<<ipta;

                            TH1F* h1tmp = (TH1F*) gROOT->FindObject(oss.str().c_str());
                            if (!h1tmp)
                            {
                                cerr<<"Could not find "<<oss.str()<<endl;
                                return -1;
                            }

                            h1accpta[iadj][itype][iptt][icent][idpsi][ipta] = (TH1F*) h1tmp->Clone(Form("acc_adj%i_%s_ptt%i_cent%i_dpsi%i_pta%i",iadj,GAT::spec[itype],iptt,icent,idpsi,ipta));
                            h1accpta[iadj][itype][iptt][icent][idpsi][ipta]->SetDirectory(output);
                        }

                    }
                }
            }
        }
    }
    cout<<"Got No RP Histos"<<endl;
    return 1;
}


int GetRpHistos()
{
    cout<<"Getting RP Histos"<<endl;

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
            {
                ostringstream oss;
                oss<<"cent_ptt";
                if (iadj)
                    oss<<"_adj";
                oss<<"_"<<GAT::spec[itype];
                oss<<"_dpsi"<<idpsi;

                if (idpsi == GAT::ndpsi)
                {
                    h2centptt[iadj][itype][idpsi] = (TH2F*) h2centptt[iadj][itype][0]->Clone(oss.str().c_str());

                    h2centptt[iadj][itype][idpsi]->Add(h2centptt[iadj][itype][1]);
                    h2centptt[iadj][itype][idpsi]->SetDirectory(output);
                    continue;
                }       

                h2centptt[iadj][itype][idpsi] = (TH2F*) gROOT->FindObject(oss.str().c_str());
                if (!h2centptt[iadj][itype][idpsi])
                {
                    cerr<<"Could not find "<<oss.str()<<endl;
                    return -1;
                }
                h2centptt[iadj][itype][idpsi]->SetDirectory(output);
            }
        }
    }
                

    
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int imix = 0; imix < nmix; imix++)
            {

                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {           
                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {
                        for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                        {
                            for (int ixi = 0; ixi < GAT::nxi; ixi++)
                            {
                                ostringstream oss;
                                oss<<"dphi";
                                if (iadj)
                                    oss<<"_adj";
                                oss<<"_"<<GAT::spec[itype];
                                oss<<"_mix"<<imix;
                                oss<<"_ptt"<<iptt;
                                oss<<"_cent"<<icent;
                                oss<<"_dpsi"<<idpsi;
                                oss<<"_xi"<<ixi;
                            
                                if (idpsi == GAT::ndpsi)
                                {
                                    h1dphi[iadj][itype][imix][iptt][icent][idpsi][ixi] = (TH1F*) h1dphi[iadj][itype][imix][iptt][icent][0][ixi]->Clone(oss.str().c_str());

                                    h1dphi[iadj][itype][imix][iptt][icent][idpsi][ixi]->Add(h1dphi[iadj][itype][imix][iptt][icent][1][ixi]);
                                    h1dphi[iadj][itype][imix][iptt][icent][idpsi][ixi]->SetDirectory(output);
                                    continue;
                                }                           

                                h1dphi[iadj][itype][imix][iptt][icent][idpsi][ixi] = (TH1F*) gROOT->FindObject(oss.str().c_str());
                                if (!h1dphi[iadj][itype][imix][iptt][icent][idpsi][ixi])
                                {
                                    cerr<<"Could not find "<<oss.str()<<endl;
                                    return -1;
                                }
                                h1dphi[iadj][itype][imix][iptt][icent][idpsi][ixi]->SetDirectory(output);

                            }

                            for (int ipta = 0; ipta < GAT::npta; ipta++)
                            {
                                ostringstream oss;
                                oss<<"dphi";
                                if (iadj)
                                    oss<<"_adj";
                                oss<<"_"<<GAT::spec[itype];
                                oss<<"_mix"<<imix;
                                oss<<"_ptt"<<iptt;
                                oss<<"_cent"<<icent;
                                oss<<"_dpsi"<<idpsi;
                                oss<<"_pta"<<ipta;
                            
                                if (idpsi == GAT::ndpsi)
                                {
                                    h1dphipta[iadj][itype][imix][iptt][icent][idpsi][ipta] = (TH1F*) h1dphipta[iadj][itype][imix][iptt][icent][0][ipta]->Clone(oss.str().c_str());

                                    h1dphipta[iadj][itype][imix][iptt][icent][idpsi][ipta]->Add(h1dphipta[iadj][itype][imix][iptt][icent][1][ipta]);
                                    h1dphipta[iadj][itype][imix][iptt][icent][idpsi][ipta]->SetDirectory(output);
                                    continue;
                                }                           

                                h1dphipta[iadj][itype][imix][iptt][icent][idpsi][ipta] = (TH1F*) gROOT->FindObject(oss.str().c_str());
                                if (!h1dphipta[iadj][itype][imix][iptt][icent][idpsi][ipta])
                                {
                                    cerr<<"Could not find "<<oss.str()<<endl;
                                    return -1;
                                }
                                h1dphipta[iadj][itype][imix][iptt][icent][idpsi][ipta]->SetDirectory(output);

                            }
                        }
                    }
                }
            }
        }
    }
    cout<<"Got RP Histos"<<endl;
    return 1;
}


void CalcAcc(const char* name_out_raw)
{
    cout<<"Calculating Acceptance"<<endl;

    TFile* output_rawnh = new TFile(name_out_raw,"RECREATE");
    cout<<"Writing raw nh distributions to "<<name_out_raw<<endl;
    
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
            {
                TH1F* h1ptt_tmp = (TH1F*) h2centptt[iadj][itype][idpsi]->ProjectionY();

                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    int pttbinlo = h1ptt_tmp->FindBin(GAT::limptt[iptt]);
                    int pttbinhi = h1ptt_tmp->FindBin(GAT::limptt[iptt+1])-1;
                    if (pttbinhi > h1ptt_tmp->GetNbinsX())
                    {
                        cerr<<"THIS SHOULD NEVER HAPPEN"<<endl;
                        pttbinhi = h1ptt_tmp->GetNbinsX();
                    }

                    TH1F* h1cent = (TH1F*) h2centptt[iadj][itype][idpsi]->ProjectionX(Form("cent_adj%i_%s_ptt%i_dpsi%i",iadj,GAT::spec[itype],iptt,idpsi),pttbinlo,pttbinhi);
                    h1cent->SetDirectory(output);
                
                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {
                        double ntrg = 0.;
                        for (int jcent = GAT::limcent[icent+1]; jcent > GAT::limcent[icent]; jcent--)
                        {
                            int tmpbin = h1cent->FindBin(jcent);
                            ntrg += h1cent->GetBinContent(tmpbin);
                        }
                            
                        
                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {
                            h1acc[iadj][itype][iptt][icent][idpsi][ixi]->Scale(1./ntrg);

                            cout<<"Writing "<<h1acc[iadj][itype][iptt][icent][idpsi][ixi]->GetName()<<" to "<<gDirectory->GetFile()->GetName()<<endl;
                            h1acc[iadj][itype][iptt][icent][idpsi][ixi]->Write();
                            
                            double scale = GAT::ndphibin/expnhxi[iadj][itype][icent][idpsi][iptt][ixi];
                            h1acc[iadj][itype][iptt][icent][idpsi][ixi]->Scale(scale);
                        }

                        for (int ipta = 0; ipta < GAT::npta; ipta++)
                        {
                            h1accpta[iadj][itype][iptt][icent][idpsi][ipta]->Scale(1./ntrg);

                            cout<<"Writing "<<h1accpta[iadj][itype][iptt][icent][idpsi][ipta]->GetName()<<" to "<<gDirectory->GetFile()->GetName()<<endl;
                            h1accpta[iadj][itype][iptt][icent][idpsi][ipta]->Write();
                            
                            double scale = GAT::ndphibin/expnhpta[icent][ipta];
                            h1accpta[iadj][itype][iptt][icent][idpsi][ipta]->Scale(scale);
                        }

                    }
                }
                
                delete h1ptt_tmp;
            }
        }
    }

    cout<<"Closing "<<output_rawnh->GetName()<<endl;
    output_rawnh->Close();

    cout<<"Acceptance Calculated"<<endl;
}

void ApplyAcc()
{
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int imix = 0; imix < nmix; imix++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {           
                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {
                        for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                        {
                            for (int ixi = 0; ixi < GAT::nxi; ixi++)
                            {
                                h1dphi[iadj][itype][imix][iptt][icent][idpsi][ixi]->Divide(h1acc[iadj][itype][iptt][icent][idpsi][ixi]);
                            }

                            for (int ipta = 0; ipta < GAT::npta; ipta++)
                            {
                                TH1F* h1accptr = h1accpta[iadj][itype][iptt][icent][idpsi][ipta];
                                if (ipta > iptamax_acc)
                                    h1accptr = h1accpta[iadj][itype][iptt][icent][idpsi][iptamax_acc];
                                h1dphipta[iadj][itype][imix][iptt][icent][idpsi][ipta]->Divide(h1accptr);
                            }
                        }
                    }
                }
            }
        }
    }
}
                                

int GetNpartDists()
{
    
    const char* centdistsname = "data/CentDists_9bins.root";
    TFile* fcentdists = new TFile(centdistsname,"READ");
    if (!fcentdists || fcentdists->IsZombie())
    {
        cerr<<"Could not open "<<centdistsname<<endl;
        return -1;
    }
    cout<<"Reading Npart distributions from "<<centdistsname<<endl;
    
    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        hnpart[icent] = (TH1F*) gROOT->FindObject(Form("AuAu_Npart_%i_%i",GAT::limcent[icent],GAT::limcent[icent+1]));
        if (!hnpart[icent])
        {
            cerr<<"Could not find AuAu_Npart_"<<GAT::limcent[icent]<<"_"<<GAT::limcent[icent+1]<<endl;
            return -1;
        }
        
        hnpart[icent]->SetDirectory(0);
        vnpart.push_back(hnpart[icent]->GetMean());
        vnparte.push_back(hnpart[icent]->GetRMS());
    }
    
    fcentdists->Close();
    return 0;
}

int GetNpartDistsMerge()
{
    const char* centdistsname = "data/CentDists_2bins.root";
    TFile* fcentdists = new TFile(centdistsname,"READ");
    if (!fcentdists || fcentdists->IsZombie())
    {
        cerr<<"Could not open "<<centdistsname<<endl;
        return -1;
    }
    cout<<"Reading Npart distributions (merged) from "<<centdistsname<<endl;
    
    for (int icent = 0; icent < ncent_merge; icent++)
    {
        hnpart_merge[icent] = (TH1F*) gROOT->FindObject(Form("AuAu_Npart_%i_%i",limcent_merge[icent],limcent_merge[icent+1]));
        if (!hnpart_merge[icent])
        {
            cerr<<"Could not find AuAu_Npart_"<<limcent_merge[icent]<<"_"<<limcent_merge[icent+1]<<endl;
            return -1;
        }
        hnpart_merge[icent]->SetDirectory(0);
        vnpart_merge.push_back(hnpart_merge[icent]->GetMean());
        vnparte_merge.push_back(hnpart_merge[icent]->GetRMS());
    }
    fcentdists->Close();
    return 0;
}


void GetKComp(const char* name, bool dopsimod)
{
    TFile* input = new TFile(name,"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Could not open "<<name<<endl;
        return;
    }

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                {                   
                    tgrng[iadj][itype][iptt][idpsi] = (TGraphErrors*) gROOT->FindObject(Form("ng_Np_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt));
                    tgrng[iadj][itype][iptt][idpsi]->SetName(Form("ng_Np_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt));
                    if (verbosity > 1)
                    {
                        cout<<"Printing contents of "<<tgrng[iadj][itype][iptt][idpsi]->GetName()<<endl;
                        tgrng[iadj][itype][iptt][idpsi]->Print();
                    }
                }
            }
        }
    }

    if (dopsimod)
    {
        TGraphErrors* tgrbhs = (TGraphErrors*) gROOT->FindObject("bhs_Np_adj0_pi0_ptt0_xi2");
        TGraphErrors* tgrah = (TGraphErrors*) gROOT->FindObject("ah_Np_adj0_pi0_ptt0_xi2");

        for (int icent = 0; icent < GAT::ncent; icent++)
        {
            double bhs = tgrbhs->GetY()[icent];
            double bhse = tgrbhs->GetEY()[icent];

            double ah = tgrah->GetY()[icent];
            double ahe = tgrah->GetEY()[icent];
                        
            psimod.push_back(bhs/ah);

            double e0 = bhse/ah;
            double e1 = - (bhs * ahe) / (ah*ah);
            psimode.push_back(sqrt( e0*e0 + e1*e1 ));
        }
    }

    input->Close();
    delete input;
}


void CalcCorrect()
{
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {       
                for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                {                       
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        vector<double> vacc;
                        vector<double> verr;

                        for (int icent = 0; icent < GAT::ncent; icent++)
                        {
                            double nbin = h1acc[iadj][itype][iptt][icent][idpsi][ixi]->GetNbinsX();

                            double error = -9999.;
                            double integral = h1acc[iadj][itype][iptt][icent][idpsi][ixi]->IntegralAndError(1,nbin,error);

                            double avgacc = integral/nbin;
                            error = error/nbin;

                            vacc.push_back(avgacc);
                            verr.push_back(error);

                        }

                        tgreff[iadj][itype][iptt][idpsi][ixi] = new TGraphErrors(GAT::ncent,&vnpart[0],&vacc[0],&vnparte[0],&verr[0]);
                        tgreff[iadj][itype][iptt][idpsi][ixi]->SetName(Form("eff_adj%i_%s_ptt%i_dpsi%i_xi%i",iadj,GAT::spec[itype],iptt,idpsi,ixi));
                        tgreff[iadj][itype][iptt][idpsi][ixi]->SetTitle(Form("eff_adj%i_%s_ptt%i_dpsi%i_xi%i;N_{part};#epsilon",iadj,GAT::spec[itype],iptt,idpsi,ixi));

                        output->cd();
                        cout<<"Writing "<<tgreff[iadj][itype][iptt][idpsi][ixi]->GetName()<<" to "<<gDirectory->GetFile()->GetName()<<endl;
                        tgreff[iadj][itype][iptt][idpsi][ixi]->Write();

                        CalcCorrHelper(iadj,itype,iptt,idpsi,ixi, true);
                        cout<<"Writing "<<tgrcorr_merge[iadj][itype][iptt][idpsi][ixi]->GetName()<<" to "<<gDirectory->GetFile()->GetName()<<endl;
                        tgrcorr_merge[iadj][itype][iptt][idpsi][ixi]->Write();
                        
                        CalcCorrHelper(iadj,itype,iptt,idpsi,ixi, false);
                        cout<<"Writing "<<tgrcorr[iadj][itype][iptt][idpsi][ixi]->GetName()<<" to "<<gDirectory->GetFile()->GetName()<<endl;
                        tgrcorr[iadj][itype][iptt][idpsi][ixi]->Write();
                        

                    }
                }
            }
        }
    }
}

void CalcCorrHelper(int iadj, int itype, int iptt, int idpsi, int ixi, bool ismerge)
{
    TGraphErrors* tgra = tgrng[iadj][itype][iptt][idpsi];
    TGraphErrors* tgrb = tgreff[iadj][itype][iptt][idpsi][ixi];
    
    bool usespline = true;

    TH1F** hnp = hnpart;
    int ncent = GAT::ncent;
    if (ismerge)
    {
        hnp = hnpart_merge;
        ncent = ncent_merge;
    }

    for (int icent = 0; icent < ncent; icent++)
    {
        int nbinx = hnp[icent]->GetNbinsX();
        
        double k = 0.;
        double ke = 0.;

        double expa = 0.;
        double expb = 0.;
        
        for (int ibin = 1; ibin <= nbinx; ibin++)
        {
            double w = hnp[icent]->GetBinContent(ibin);
            if (w <= 0.0) continue;
            
            double we = hnp[icent]->GetBinError(ibin);
            
            double npart = hnp[icent]->GetBinCenter(ibin);
            
            double a, ae, b, be;
            
            if (usespline)
            {
                a = tgra->Eval(npart);
                if (a < 0.)
                {
                    ostringstream oss;
                    oss<<tgra->GetName();
                    string str = oss.str();
                    if (str[0] == 'a')
                    {
                        cout<<tgra->GetName()<<"->Eval("<<npart<<") evaluated to "<<a<<". Setting to 0..."<<endl;
                        cout<<hnp[icent]->GetName()<<": hnp["<<icent<<"]->GetBinContent("<<ibin<<") = "<<w<<endl;
                        a = 0;
                    }
                }

                ae = 0.;
                
                b = tgrb->Eval(npart);
                if (b < 0.)
                {
                    ostringstream oss;
                    oss<<tgrb->GetName();
                    string str = oss.str();
                    if (str[0] == 'a')
                    {
                        cout<<tgrb->GetName()<<"->Eval("<<npart<<") evaluated to "<<b<<". Setting to 0..."<<endl;
                        cout<<hnp[icent]->GetName()<<": hnp["<<icent<<"]->GetBinContent("<<ibin<<") = "<<w<<endl;
                        b = 0;
                    }
                }

                be = 0.;
            }
            else
            {
                cerr<<"THIS IS NO ELSE!"<<endl;
                return;
            }       
            
            k += w*a*b;
            
            ke += pow(w*ae*b,2.) + pow(w*a*be,2.);
            
            expa += w*a;
            
            expb += w*b;

        }

        k /= (expa*expb);
        if (ismerge)
            corr_merge[iadj][itype][iptt][idpsi][ixi].push_back(k);
        else
            corr[iadj][itype][iptt][idpsi][ixi].push_back(k);
                            
        ke = sqrt(ke);
        ke /= (expa*expb);
        if (ismerge)
            corre_merge[iadj][itype][iptt][idpsi][ixi].push_back(ke);
        else
            corre[iadj][itype][iptt][idpsi][ixi].push_back(ke);

    }

    if (ismerge)
    {
        tgrcorr_merge[iadj][itype][iptt][idpsi][ixi] = new TGraphErrors(ncent,&vnpart_merge[0],&corr_merge[iadj][itype][iptt][idpsi][ixi][0],&vnparte_merge[0],&corre_merge[iadj][itype][iptt][idpsi][ixi][0]);
        tgrcorr_merge[iadj][itype][iptt][idpsi][ixi]->SetName(Form("corr_merge_adj%i_%s_ptt%i_dpsi%i_xi%i",iadj,GAT::spec[itype],iptt,idpsi,ixi));
        tgrcorr_merge[iadj][itype][iptt][idpsi][ixi]->SetTitle(Form("corr_merge_adj%i_%s_ptt%i_dpsi%i_xi%i;N_{part};#epsilon_{cor}",iadj,GAT::spec[itype],iptt,idpsi,ixi));
    }
    else
    {
        tgrcorr[iadj][itype][iptt][idpsi][ixi] = new TGraphErrors(ncent,&vnpart[0],&corr[iadj][itype][iptt][idpsi][ixi][0],&vnparte[0],&corre[iadj][itype][iptt][idpsi][ixi][0]);
        tgrcorr[iadj][itype][iptt][idpsi][ixi]->SetName(Form("corr_adj%i_%s_ptt%i_dpsi%i_xi%i",iadj,GAT::spec[itype],iptt,idpsi,ixi));
        tgrcorr[iadj][itype][iptt][idpsi][ixi]->SetTitle(Form("corr_adj%i_%s_ptt%i_dpsi%i_xi%i;N_{part};#epsilon_{cor}",iadj,GAT::spec[itype],iptt,idpsi,ixi));
    }

}

int ApplyCorrect(bool docor, bool dopsicor)
{
    cout<<"Applying Correction."<<endl;
    if (docor)
    {
        cout<<"WARNING! Doing corr correction which has been deprecated."<<endl;
        return -1;
    }
    if (!dopsicor)
    {
        cout<<"Not doing psi modulation to acceptance.  Skipping ApplyCorrect()."<<endl;
        return 1;
    }

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int imix = 0; imix < nmix; imix++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {           
                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {
                        for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                        {
                            for (int ixi = 0; ixi < GAT::nxi; ixi++)
                            {
                                if (docor)
                                {
                                    h1acc[iadj][itype][iptt][icent][idpsi][ixi]->Scale(corr[iadj][itype][iptt][idpsi][ixi][icent]);
                                }
                                if (dopsicor && idpsi!=GAT::ndpsi)
                                {
                                    double mod = psimod[icent];
                                    mod*=(2./PI);
                                    if (idpsi==1) mod*=-1.;

                                    double npart = vnpart[icent];
                                    
                                    TF1* f1 = new TF1(Form("fpsicor_adj%i_%s_ptt%i_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],iptt,icent,idpsi,ixi),[=](double* x, double* par){return tgreff[iadj][itype][iptt][idpsi][ixi]->Eval(npart*(1+mod*cos(2.*PI*x[0])))/tgreff[iadj][itype][iptt][idpsi][ixi]->Eval(npart);},0.,1.0,0);
                                    cout<<"Writing "<<f1->GetName()<<" to "<<gDirectory->GetFile()->GetName()<<endl;
                                    f1->Write();

                                    h1acc[iadj][itype][iptt][icent][idpsi][ixi]->Multiply(f1);
                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}

TGraphErrors* MakeProd(TGraphErrors* tgra, TGraphErrors* tgrae, TGraphErrors* tgrb, TGraphErrors* tgrbe, string name)
{
    vector<double> vprod;
    vector<double> vprode;
    vector<double> vnpart;
    for (int npart = npartmin; npart <= npartmax; npart++)
    {
        double a, ae, b, be;
        
        a = tgra->Eval(npart);
        if (a < 0.)
        {
            ostringstream oss;
            oss<<tgra->GetName();
            string str = oss.str();
            if (str[0] == 'a' || str[0] == 'n')
            {
                cout<<tgra->GetName()<<"->Eval("<<npart<<") evaluated to "<<a<<". Setting to 0."<<endl;
                a = 0;
            }
        }

        if (tgrae)
        {
            ae = tgrae->Eval(npart);
            if (ae < 0.) ae = 0.;
        }
        else ae = 0.;
                
        b = tgrb->Eval(npart);
        if (b < 0.)
        {
            ostringstream oss;
            oss<<tgrb->GetName();
            string str = oss.str();
            if (str[0] == 'a' || str[0] == 'n')
            {
                cout<<tgrb->GetName()<<"->Eval("<<npart<<") evaluated to "<<b<<". Setting to 0."<<endl;
                b = 0;
            }
        }

        if (tgrbe)
        {
            be = tgrbe->Eval(npart);
            if (be < 0.) be = 0.;
        }
        else be = 0.;
        
        double prod = a*b;
        double prode = sqrt(pow(a*be,2.) + pow(ae*b,2.));

        vprod.push_back(prod);
        vprode.push_back(prode);
        vnpart.push_back(npart);
    }

    TGraphErrors* tgr = new TGraphErrors(vprod.size(),&vnpart[0],&vprod[0],0,&vprode[0]);
    tgr->SetName(name.c_str());

    return tgr;
}

void DrawKComps(TGraphErrors* gra, TGraphErrors* grae, string laba, TGraphErrors* grb, TGraphErrors* grbe, string labb, TGraphErrors* grp, string labp, TH1F* hnpart[GAT::ncent], string kname = "", TGraphErrors* grk = nullptr, TGraphErrors* gra0 = nullptr, TGraphErrors* grae0 = nullptr, string laba0 = "", TGraphErrors* grb0 = nullptr, TGraphErrors* grbe0 = nullptr, string labb0 = "")
{
    bool flag = labb0.length() != 0;

    gra->SetMarkerColor(kViolet);
    gra->SetMarkerStyle(21);

    TF1* fa = new TF1(Form("%s_eval",gra->GetName()),[=](double* x, double* par){return gra->Eval(x[0]);},npartmin,npartmax,0);
    fa->SetLineColor(kViolet);

    TF1* fa0 = nullptr;
    if (flag)
    {
        gra0->SetMarkerColor(kBlue);
        gra0->SetMarkerStyle(21);

        fa0 = new TF1(Form("%s_eval",gra0->GetName()),[=](double* x, double* par){return gra0->Eval(x[0]);},npartmin,npartmax,0);
        fa0->SetLineColor(kBlue);
    }

    grb->SetMarkerColor(kRed);
    grb->SetMarkerStyle(20);

    TF1* fb = new TF1("fb",[=](double* x, double* par){return grb->Eval(x[0]);},npartmin,npartmax,0);
    fb->SetLineColor(kRed);

    TF1* fb0 = nullptr; 
    if (flag)
    {
        grb0->SetMarkerColor(kOrange);
        grb0->SetMarkerStyle(20);

        fb0 = new TF1(Form("%s_eval",grb0->GetName()),[=](double* x, double* par){return grb0->Eval(x[0]);},npartmin,npartmax,0);
        fb0->SetLineColor(kOrange);
    }


    vector<double> vmean, vrms, vexpa, vexpae, vexpb, vexpbe, vprod, vprode;    

    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        TCanvas* can = new TCanvas("can","",1200,1350);
        can->cd();
        TPad *pad1 = new TPad("pad1","",0,0.5,1,1);
        pad1->SetBottomMargin(0.01);
        pad1->Draw();
        pad1->cd();

        double xmin, xmax, ymin, ymax;

        double mean = hnpart[icent]->GetMean();
        vmean.push_back(mean);
        
        double rms = hnpart[icent]->GetRMS();
        vrms.push_back(rms);
        
        xmin = mean-5.*rms;
        if (xmin < npartmin) xmin = npartmin;
        xmax = mean+5.*rms;
        if (xmax > npartmax) xmax = npartmax;
        hnpart[icent]->GetXaxis()->SetRangeUser(xmin,xmax);

        ymin = gra->Eval(xmin);
        ymax = gra->Eval(xmin);
        for (int ix = xmin; ix < xmax; ++ix)
        {
            double eval = gra->Eval(ix);
            if (eval < ymin) ymin = eval;
            if (eval > ymax) ymax = eval;
        }

        ymax = ymax + fabs(ymax)*.1;
        ymin = ymin - fabs(ymin)*.1;
        ymin = ymin - (ymax - ymin)*.1;
        
        TH1F *hr = can->DrawFrame(xmin,ymin,xmax,ymax);
        hr->SetTitle(Form("%s_cent%i",grp->GetName(),icent));

        gra->Draw("LP");

        fa->Draw("same");

        if (flag)
        {
            gra0->Draw("LP same");

            fa0->Draw("same");
        }

        double expa = -9999.;
        double expae = -9999.;

        if (flag)
            CalcExp(fa0,grae0,hnpart[icent],expa,expae);
        else
            CalcExp(fa,grae,hnpart[icent],expa,expae);

        vexpa.push_back(expa);
        vexpae.push_back(expae);
        
        TGraphErrors* grexpa = new TGraphErrors(1,&mean,&expa,&rms,&expae);
        grexpa->SetMarkerColor(kCyan);
        grexpa->SetMarkerStyle(21);
        grexpa->Draw("P same");



        can->cd();
        TPad *overlay1 = new TPad("overlay1","",0,0.5,1,1);
        overlay1->SetBottomMargin(0.01);
        overlay1->SetFillStyle(0);
        overlay1->SetFillColor(0);
        overlay1->SetFrameFillStyle(0);
        overlay1->Draw("FA");
        overlay1->cd();

        ymin = grb->Eval(xmin);
        ymax = grb->Eval(xmin);
        for (int ix = xmin; ix < xmax; ++ix)
        {
            double eval = grb->Eval(ix);
            if (eval < ymin) ymin = eval;
            if (eval > ymax) ymax = eval;
        }

        ymax = ymax + fabs(ymax)*.1;
        ymin = ymin - fabs(ymin)*.1;
        ymin = ymin - (ymax - ymin)*.1;
        
        TH1F *hframe = overlay1->DrawFrame(xmin,ymin,xmax,ymax);
        hframe->GetXaxis()->SetLabelOffset(99);
        hframe->GetYaxis()->SetTickLength(0);
        hframe->GetYaxis()->SetLabelOffset(99);
        grb->Draw("LP");

        fb->Draw("same");

        if (flag)
        {
            grb0->Draw("LP same");

            fb0->Draw("same");
        }

        double expb = -9999.;
        double expbe = -9999.;       

        if (flag)
            CalcExp(fb0,grbe0,hnpart[icent],expb,expbe);
        else
            CalcExp(fb,grbe,hnpart[icent],expb,expbe);

        vexpb.push_back(expb);
        vexpbe.push_back(expbe);
                
        TGraphErrors* grexpb = new TGraphErrors(1,&mean,&expb,&rms,&expbe);
        grexpb->SetMarkerColor(kYellow);
        grexpb->SetMarkerStyle(20);
        grexpb->Draw("P same");

        TGaxis *axis1 = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
        axis1->SetLineColor(kRed);
        axis1->SetLabelColor(kRed);
        axis1->Draw();

        TLegend* leg1 = new TLegend(0.8,0.1,0.85,0.3);
        leg1->AddEntry(gra,laba.c_str(),"p");
        if (flag)
        {
            leg1->AddEntry(gra0,laba0.c_str(),"p");
            leg1->AddEntry(grexpa,Form("<%s>",laba0.c_str()),"p");
        }           
        else
            leg1->AddEntry(grexpa,Form("<%s>",laba.c_str()),"p");
        
        leg1->AddEntry(grb,labb.c_str(),"p");
        if (flag)
        {
            leg1->AddEntry(grb0,labb0.c_str(),"p");
            leg1->AddEntry(grexpb,Form("<%s>",labb0.c_str()),"p");
        }
        else
            leg1->AddEntry(grexpb,Form("<%s>",labb.c_str()),"p");
        leg1->Draw();

        can->cd();

        TPad *pad2 = new TPad("pad2","",0,0,1,0.5);
        pad2->SetTopMargin(0.05);
        pad2->Draw();
        pad2->cd();

        gStyle->SetOptStat(0);
        hnpart[icent]->SetTitle(";N_{part}");
        hnpart[icent]->Draw();
        hnpart[icent]->GetYaxis()->SetTickLength(0);
        hnpart[icent]->GetYaxis()->SetLabelOffset(99);
        ymax = hnpart[icent]->GetBinContent(hnpart[icent]->GetMaximumBin()) * 2;
        hnpart[icent]->SetMaximum(ymax);
        hnpart[icent]->SetMinimum(0.);
        pad2->Update();

        can->cd();
        TPad *overlay2 = new TPad("overlay2","",0,0,1,0.5);
        overlay2->SetTopMargin(0.05);
        overlay2->SetFillStyle(0);
        overlay2->SetFillColor(0);
        overlay2->SetFrameFillStyle(0);
        overlay2->Draw("FA");
        overlay2->cd();

        ymin = grp->Eval(xmin);
        ymax = grp->Eval(xmin);
        for (int ix = xmin; ix < xmax; ++ix)
        {
            double eval = grp->Eval(ix);
            if (eval < ymin) ymin = eval;
            if (eval > ymax) ymax = eval;
        }
    
        ymin = ymin - fabs(ymin)*.1;
        ymax = ymax + fabs(ymax)*.1;
        ymin = ymin - (ymax - ymin)*.1;
        TH1F *hframe2 = overlay2->DrawFrame(xmin,ymin,xmax,ymax);
        hframe2->GetXaxis()->SetLabelOffset(99);
        hframe2->GetYaxis()->SetTickLength(0);
        hframe2->GetYaxis()->SetLabelOffset(99);
        grp->Draw("L");

        double prod = expa*expb;
        vprod.push_back(prod);
        
        double prode = sqrt(pow(expa*expbe,2.) + pow(expae*expb,2.));
        vprode.push_back(prode);
        
        TGraphErrors* grprod = new TGraphErrors(1,&mean,&prod,&rms,&prode);
        grprod->SetMarkerStyle(21);
        grprod->SetMarkerColor(kMagenta);
        grprod->Draw("LP same");


        TGaxis *axis2 = new TGaxis(xmin,ymin,xmin,ymax, ymin,ymax,510);
        axis2->Draw();
        overlay2->Update();

        TLegend* leg2 = new TLegend(0.77,0.25,0.88,0.4);
        leg2->AddEntry(grp,labp.c_str(),"le");
        if (flag)
            leg2->AddEntry(grprod,Form("<%s><%s>",laba0.c_str(),labb0.c_str()),"p");
        else
            leg2->AddEntry(grprod,Form("<%s><%s>",laba.c_str(),labb.c_str()),"p");

        if (kname.length() != 0)
        {
            ostringstream oss;
            oss<<kname<<" = "<<grk->GetY()[icent];
            leg2->SetHeader(oss.str().c_str(),"C");
        }
        leg2->Draw();

        can->SaveAs(Form("kdraw/%s_cent%i.eps",grp->GetName(),icent));
        delete can;
    }
    
    TCanvas* can = new TCanvas("can","",1200,1350);
    TPad *pad1 = new TPad("pad1","",0,0.65,1,1);
    pad1->SetBottomMargin(0.01);
    pad1->Draw();
    pad1->cd();

    double xmin = npartmin;
    double xmax = npartmax;

    double ymin = gra->Eval(xmin);
    double ymax = gra->Eval(xmin);
    for (int ix = xmin; ix < xmax; ++ix)
    {
        double eval = gra->Eval(ix);
        if (eval < ymin) ymin = eval;
        if (eval > ymax) ymax = eval;
    }

    ymax = ymax + fabs(ymax)*.1;
    ymin = ymin - fabs(ymin)*.1;
    ymin = ymin - (ymax - ymin)*.1;

    TH1F *frame1 = can->DrawFrame(xmin,ymin,xmax,ymax);
    frame1->SetTitle(grp->GetName());

    gra->Draw("LP");

    fa->Draw("same");

    if (flag)
    {
        gra0->Draw("LP same");

        fa0->Draw("same");
    }

    TGraphErrors* grexpaall = new TGraphErrors(GAT::ncent,&vmean[0],&vexpa[0],&vrms[0],&vexpae[0]);
    grexpaall->SetMarkerColor(kCyan);
    grexpaall->SetMarkerStyle(21);
    grexpaall->Draw("P same");


    can->cd();
    TPad *overlay1 = new TPad("overlay1","",0,0.65,1,1);
    overlay1->SetBottomMargin(0.01);
    overlay1->SetFillStyle(0);
    overlay1->SetFillColor(0);
    overlay1->SetFrameFillStyle(0);
    overlay1->Draw("FA");
    overlay1->cd();

    
    ymin = grb->Eval(xmin);
    ymax = grb->Eval(xmin);
    for (int ix = xmin; ix < xmax; ++ix)
    {
        double eval = grb->Eval(ix);
        if (eval < ymin) ymin = eval;
        if (eval > ymax) ymax = eval;
    }

    ymax = ymax + fabs(ymax)*.1;
    ymin = ymin - fabs(ymin)*.1;
    ymin = ymin - (ymax - ymin)*.1;

    TH1F *frame2 = overlay1->DrawFrame(xmin,ymin,xmax,ymax);
    frame2->GetXaxis()->SetLabelOffset(99);
    frame2->GetYaxis()->SetTickLength(0);
    frame2->GetYaxis()->SetLabelOffset(99);
    grb->Draw("LP");

    fb->Draw("same");

    if (flag)
    {
        grb0->Draw("LP same");

        fb0->Draw("same");
    }

    TGraphErrors* grexpball = new TGraphErrors(GAT::ncent,&vmean[0],&vexpb[0],&vrms[0],&vexpbe[0]);
    grexpball->SetMarkerColor(kYellow);
    grexpball->SetMarkerStyle(20);
    grexpball->Draw("P same");

    TGaxis *axis1 = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
    axis1->SetLineColor(kRed);
    axis1->SetLabelColor(kRed);
    axis1->Draw();
    
    TLegend* leg1 = new TLegend(0.8,0.1,0.85,0.3);
    leg1->AddEntry(gra,laba.c_str(),"p");
    if (flag)
    {
        leg1->AddEntry(gra0,laba0.c_str(),"p");
        leg1->AddEntry(grexpaall,Form("<%s>",laba0.c_str()),"p");
    }       
    else
        leg1->AddEntry(grexpaall,Form("<%s>",laba.c_str()),"p");
        
    leg1->AddEntry(grb,labb.c_str(),"p");
    if (flag)
    {
        leg1->AddEntry(grb0,labb0.c_str(),"p");
        leg1->AddEntry(grexpball,Form("<%s>",labb0.c_str()),"p");
    }
    else
        leg1->AddEntry(grexpball,Form("<%s>",labb.c_str()),"p");
    leg1->Draw();

    can->cd();

    TPad *pad2 = new TPad("pad2","",0,.3,1,0.65);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.);
    pad2->Draw();
    pad2->cd();

    gStyle->SetOptStat(0);
    hnpart[GAT::ncent-1]->SetTitle("");
    hnpart[GAT::ncent-1]->Draw();
    hnpart[GAT::ncent-1]->GetYaxis()->SetTickLength(0);
    hnpart[GAT::ncent-1]->GetYaxis()->SetLabelOffset(99);
    ymax = hnpart[GAT::ncent-1]->GetBinContent(hnpart[GAT::ncent-1]->GetMaximumBin()) * 2;
    hnpart[GAT::ncent-1]->SetMaximum(ymax);
    hnpart[GAT::ncent-1]->SetMinimum(0.);
    hnpart[GAT::ncent-1]->GetXaxis()->SetRangeUser(xmin,xmax);
    pad2->Update();

    for (int icent = 0; icent < GAT::ncent-1; icent++)
    {   
        hnpart[icent]->Draw("same");
    }
    pad2->Update();
    
    
    can->cd();
    TPad *overlay2 = new TPad("overlay2","",0,.3,1,0.65);
    overlay2->SetTopMargin(0.05);
    overlay2->SetBottomMargin(0.);
    overlay2->SetFillStyle(0);
    overlay2->SetFillColor(0);
    overlay2->SetFrameFillStyle(0);
    overlay2->Draw("FA");
    overlay2->cd();

    ymin = grp->Eval(xmin);
    ymax = grp->Eval(xmin);
    for (int ix = xmin; ix < xmax; ++ix)
    {
        double eval = grp->Eval(ix);
        if (eval < ymin) ymin = eval;
        if (eval > ymax) ymax = eval;
    }
    
    ymin = ymin - fabs(ymin)*.1;
    ymax = ymax + fabs(ymax)*.1;
    ymin = ymin - (ymax - ymin)*.1;

    TH1F *frame3 = overlay2->DrawFrame(xmin,ymin,xmax,ymax);
    frame3->GetXaxis()->SetLabelOffset(99);
    frame3->GetYaxis()->SetTickLength(0);
    frame3->GetYaxis()->SetLabelOffset(99);
    grp->Draw("L");

    TGraphErrors* grprodall = new TGraphErrors(GAT::ncent,&vmean[0],&vprod[0],&vrms[0],&vprode[0]);
    grprodall->SetMarkerStyle(21);
    grprodall->SetMarkerColor(kMagenta);
    grprodall->Draw("P same");

    TGaxis *axis2 = new TGaxis(xmin,ymin,xmin,ymax, ymin,ymax,510);
    axis2->Draw();

    TLegend* leg2 = new TLegend(0.77,0.25,0.88,0.4);
    leg2->AddEntry(grp,labp.c_str(),"le");
    if (flag)
        leg2->AddEntry(grprodall,Form("<%s><%s>",laba0.c_str(),labb0.c_str()),"p");
    else
        leg2->AddEntry(grprodall,Form("<%s><%s>",laba.c_str(),labb.c_str()),"p");

    leg2->Draw();

    can->cd();
    TPad* pad3 = new TPad("pad3","",0,0,1,.3);
    pad3->SetTopMargin(0);
    pad3->Draw();
    pad3->cd();

    if (flag)
        grk->SetTitle(Form(";N_{part};#frac{<%s%s>}{<%s><%s>}",laba.c_str(),labb.c_str(),laba0.c_str(),labb0.c_str()));
    else
        grk->SetTitle(Form(";N_{part};#frac{<%s%s>}{<%s><%s>}",laba.c_str(),labb.c_str(),laba.c_str(),labb.c_str()));
    grk->SetMarkerStyle(kFullCircle);
    grk->SetMarkerColor(kBlue);
    grk->SetLineColor(kBlue);
    grk->GetXaxis()->SetRangeUser(xmin,xmax);
    grk->Draw("ALP");
    pad3->Update();
    ymin = grk->GetHistogram()->GetMinimum();
    ymax = grk->GetHistogram()->GetMaximum();
    ymin = ymin - (ymax-ymin)*.1;
    if (ymin > 1) ymin = 1;
    grk->GetYaxis()->SetRangeUser(ymin,ymax);
    TLine* tl = new TLine(xmin,1.0,xmax,1.0);
    tl->SetLineStyle(9);
    tl->Draw();
    pad3->Update();

    can->SaveAs(Form("kdraw/%s.eps",grp->GetName()));
}

void CalcExp(TF1* fa, TGraphErrors* grae, TH1F* hnp, double& expa, double& expae)
{
    int nbinx = hnp->GetNbinsX();

    expa = 0.;
    expae = 0.;
        
    for (int ibin = 1; ibin <= nbinx; ibin++)
    {
        double w = hnp->GetBinContent(ibin);
        if (w <= 0.0) continue;
            
        double we = hnp->GetBinError(ibin);
            
        double npart = hnp->GetBinCenter(ibin);     
            
        double a = fa->Eval(npart);
        if (a < 0.)
        {
            ostringstream oss;
            oss<<fa->GetName();
            string str = oss.str();
            if (str[0] == 'a')
            {
                cout<<fa->GetName()<<"->Eval("<<npart<<") evaluated to "<<a<<". Setting to 0."<<endl;
                continue;
            }
        }

        double ae = 0.;
        if (grae)
        {
            double ae = grae->Eval(npart);
            if (ae < 0.) ae = 0.;                           
        }
        
        expa += w*a;

        expae += w*ae;

    }

    return;
}
