
// Merges centrality bins (0-20% and 20-60%) for the acceptance correction.
// Reads trigger efficiency (tgreff), hadron ratio (tgrbhs/tgrah), and correction (tgrcorr) TGraphs
// from acceptCorrect.C output and kComponents_Merge.C, then combines them.

const int ncent_merge = 2;
const int limcent_merge[ncent_merge+1] = {0,20,60};

const float PI = TMath::Pi();

const int nadj = 2;
const int nmix = 2;
const int ntype = 3;

const double detah = 0.7;

int iptamax_acc = 2;

TFile* gamtree_rp;
TFile* gamtree_norp;
TFile* output;

void GetNhExpXi();
void GetNhExpPta();
int GetNoRpHistos();
int GetRpHistos();
void CalcAcc(const char* name_out_raw);

void GetKComp(const char* name, bool);
int ReadAcceptCorrect(const char* name);
int ApplyCorrect(bool docor, bool dopsicor);

void ApplyAcc();
void CentMerge();

namespace GAT {
     #include "gamAnaTree.h"
}

double expnhpta[ncent_merge][GAT::npta] = {};
TH1F* h1accpta[nadj][ntype][GAT::nptt][ncent_merge][GAT::ndpsi+1][GAT::npta];
TH1F* h1dphipta[nadj][ntype][nmix][GAT::nptt][ncent_merge][GAT::ndpsi+1][GAT::npta];

double expnhxi[nadj][ntype][ncent_merge][GAT::ndpsi+1][GAT::nptt][GAT::nxi] = {};
TH1F* h1acc[nadj][ntype][GAT::nptt][ncent_merge][GAT::ndpsi+1][GAT::nxi];
TH1F* h1dphi[nadj][ntype][nmix][GAT::nptt][ncent_merge][GAT::ndpsi+1][GAT::nxi];
TH2F* h2centptt[nadj][ntype][GAT::ndpsi+1];

TGraphErrors* tgrng[nadj][ntype][GAT::nptt][GAT::ndpsi+1];
TGraphErrors* tgreff[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];

vector<double> corr_merge[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];
vector<double> corre_merge[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];
TGraphErrors* tgrcorr_merge[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];

vector<double> psimod;
vector<double> psimode;
TGraphErrors* tgrpsimod;

TH1F* hnpart_merge[ncent_merge];
vector<double> vnpart_merge;
vector<double> vnparte_merge;
int GetNpartDistsMerge();

TH1F* hnpart[GAT::ncent];
vector<double> vnpart;
vector<double> vnparte;
int GetNpartDists();

const int npartmin = 2;
const int npartmax = 390;
TGraphErrors* MakeProd(TGraphErrors* tgra, TGraphErrors* tgrae, TGraphErrors* tgrb, TGraphErrors* tgrbe, string name);
void CalcExp(TF1* fa, TGraphErrors* grae, TH1F* hnpart, double& expa, double& expae);
    

void acceptCorrect_Merge(
    bool dopsimod = true,
    const char* name_out = "gamtree_13466_acceptcorrect_merge.root",
    const char* name_kcomp = "kComponents/kComponents_Merge.root",
    const char* name_gamtree_acc = "data/gamtree_13466_acceptcorrect.root",
    const char* name_gamtree_rp = "rootfiles/gamtree_13466_psimix_asr_merge.root",
    const char* name_gamtree_norp = "rootfiles/gamtree_13466_nopsimix_asr_merge.root",
    const char* name_out_raw = "gamtree_13466_rawnh_merge.root"
    )
{

    if (GetNpartDistsMerge()<0) return;

    gamtree_rp = new TFile(name_gamtree_rp,"READ");
    if (!gamtree_rp || gamtree_rp->IsZombie())
    {
        cerr<<"Could not open "<<name_gamtree_rp<<endl;
        return;
    }
    cout<<"Opened "<<name_gamtree_rp<<" successfully."<<endl;

    gamtree_norp = new TFile(name_gamtree_norp,"READ");
    if (!gamtree_norp || gamtree_norp->IsZombie())
    {
        cerr<<"Could not open "<<name_gamtree_norp<<endl;
        return;
    }
    cout<<"Opened "<<name_gamtree_norp<<" successfully."<<endl;
    
    output = new TFile(name_out,"RECREATE");
    cout<<"Output will be written to "<<name_out<<endl;
    
    GetNhExpXi();

    GetNhExpPta();

    gamtree_norp->cd();
    cout<<"Changed directory to "<<gamtree_norp->GetName()<<"  gDirectory->GetFile()->GetName() = "<<gDirectory->GetFile()->GetName()<<endl;
    if (GetNoRpHistos() < 0) return;

    gamtree_rp->cd();
    cout<<"Changed directory to "<<gamtree_rp->GetName()<<"  gDirectory->GetFile()->GetName() = "<<gDirectory->GetFile()->GetName()<<endl;
    if (GetRpHistos() < 0) return;
    
    CalcAcc(name_out_raw);

    GetKComp(name_kcomp,dopsimod);
    ReadAcceptCorrect(name_gamtree_acc);

    output->cd();
    cout<<"Changed directory to "<<output->GetName()<<"  gDirectory->GetFile()->GetName() = "<<gDirectory->GetFile()->GetName()<<endl;
    if (ApplyCorrect(true,dopsimod) < 0)
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

void GetNhExpXi()
{
    cout<<"Getting NhExpXi"<<endl;

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            for (int icent = 0; icent < ncent_merge; icent++)
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

            for (int icent = 0; icent < ncent_merge; icent++)
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
                    ossdir << "_cent" << icent;
                    ossdir << "_dpsi" << idpsi;

                    ostringstream ossin;
                    ossin<<"hadronaccept/wrk/xi_merge/"<<ossdir.str()<<"/integrated_xi_"<<ossdir.str()<<".txt";

                    ifstream ifs(ossin.str(),ifstream::in);
                    if (!ifs.is_open())
                    {
                        cerr<<"Cannot open "<<ossin.str()<<endl;
                        return;
                    }

                    cout<<"Reading "<<ossin.str()<<endl;
                    
                    string line;
                    while (getline(ifs,line))
                    {
                        if (!isdigit(line[0])) continue;                        
                        istringstream iss(line.data());

                        int iptt, ixi;
                        float tmpnh;
                        
                        iss >> iptt >> ixi >> tmpnh;
                        cout << " " << iptt << " " << ixi << " " << tmpnh << endl;

                        expnhxi[iadj][itype][icent][idpsi][iptt][ixi] = tmpnh*detah*2.*PI;
                    }
                }
            }
        }
    }
    cout<<"Got NhExpXi"<<endl;
}

void GetNhExpPta()
{
    cout<<"Getting NhExpPta"<<endl;

    for (int icent = 0; icent < ncent_merge; icent++)
    {
        for (int ipta = 0; ipta < GAT::npta; ipta++)
        {           
            expnhpta[icent][ipta] = -9999.;
        }       
    }
            
    for (int icent = 0; icent < ncent_merge; icent++)
    {

        ostringstream ossin;
        ossin<<"hadronaccept/wrk/pta_merge/integrated_pta_cent"<<icent<<".txt";
                    
        ifstream ifs(ossin.str(),ifstream::in);
        if (!ifs.is_open())
        {
            cerr<<"Cannot open "<<ossin.str()<<endl;
            return;
        }

        cout<<"Reading "<<ossin.str()<<endl;
                    
        string line;
        while (getline(ifs,line))
        {
            if (!isdigit(line[0])) continue;                    
            istringstream iss(line.data());

            int ipta;
            float tmpnh;
                        
            iss >> ipta >> tmpnh;
            cout << " " << ipta << " " << tmpnh << endl;

            expnhpta[icent][ipta] = tmpnh*detah*2.*PI;
        }
    }
    cout<<"Got NhExpPta"<<endl;
}
                    
int GetNoRpHistos()
{
    cout<<"Getting No RP Histos"<<endl;
    
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {       
                for (int icent = 0; icent < ncent_merge; icent++)
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

            for (int imix = 0; imix < nmix; imix++)
            {

                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {           
                    for (int icent = 0; icent < ncent_merge; icent++)
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
                
                    for (int icent = 0; icent < ncent_merge; icent++)
                    {
                        double ntrg = 0.;
                        for (int jcent = limcent_merge[icent+1]; jcent > limcent_merge[icent]; jcent--)
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

            for (int imix = 0; imix < nmix; imix++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {           
                    for (int icent = 0; icent < ncent_merge; icent++)
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
                                    h1accptr = h1accpta[iadj][itype][iptt][icent][idpsi][3];
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
            
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                {                   
                    tgrng[iadj][itype][iptt][idpsi] = (TGraphErrors*) gROOT->FindObject(Form("ng_Np_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt));
                    tgrng[iadj][itype][iptt][idpsi]->SetName(Form("ng_Np_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt));
                    cout<<"Printing contents of "<<tgrng[iadj][itype][iptt][idpsi]->GetName()<<endl;
                    tgrng[iadj][itype][iptt][idpsi]->Print();
                }
            }
        }
    }

    if (dopsimod)
    {
        TGraphErrors* tgrbhs = (TGraphErrors*) gROOT->FindObject("bhs_Np_adj0_pi0_ptt0_xi2");
        TGraphErrors* tgrah = (TGraphErrors*) gROOT->FindObject("ah_Np_adj0_pi0_ptt0_xi2");

        for (int icent = 0; icent < ncent_merge; icent++)
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

int ReadAcceptCorrect(const char* name)
{
    cout<<"Retrieving tgrcorr_merge and tgreff from "<<name<<endl;

    TFile* input = new TFile(name,"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Could not open "<<name<<endl;
        return -1;
    }
    
    cout<<"Opened "<<name<<endl;

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                {
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        tgrcorr_merge[iadj][itype][iptt][idpsi][ixi] = (TGraphErrors*) gROOT->FindObject(Form("corr_merge_adj%i_%s_ptt%i_dpsi%i_xi%i",iadj,GAT::spec[itype],iptt,idpsi,ixi));
                        
                        cout<<"Printing contents of "<<tgrcorr_merge[iadj][itype][iptt][idpsi][ixi]->GetName()<<endl;
                        tgrcorr_merge[iadj][itype][iptt][idpsi][ixi]->Print();

                        tgreff[iadj][itype][iptt][idpsi][ixi] = (TGraphErrors*) gROOT->FindObject(Form("eff_adj%i_%s_ptt%i_dpsi%i_xi%i",iadj,GAT::spec[itype],iptt,idpsi,ixi));
                        cout<<"Printing contents of "<<tgreff[iadj][itype][iptt][idpsi][ixi]->GetName()<<endl;
                        tgreff[iadj][itype][iptt][idpsi][ixi]->Print();

                    }
                }
            }
        }
    }
    return 0;
}

int ApplyCorrect(bool docor, bool dopsicor)
{
    cout<<"Applying Correction."<<endl;

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;

            for (int imix = 0; imix < nmix; imix++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {           
                    for (int icent = 0; icent < ncent_merge; icent++)
                    {
                        for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                        {
                            for (int ixi = 0; ixi < GAT::nxi; ixi++)
                            {
                                if (docor)
                                {
                                    h1acc[iadj][itype][iptt][icent][idpsi][ixi]->Scale(tgrcorr_merge[iadj][itype][iptt][idpsi][ixi]->GetY()[icent]);
                                }
                                if (dopsicor && idpsi!=GAT::ndpsi)
                                {
                                    double mod = psimod[icent];
                                    mod*=(2./PI);
                                    if (idpsi==1) mod*=-1.;

                                    double npart = vnpart_merge[icent];

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

void CalcExp(TF1* fa, TGraphErrors* grae, TH1F* hnpart, double& expa, double& expae)
{
    int nbinx = hnpart->GetNbinsX();

    expa = 0.;
    expae = 0.;
        
    for (int ibin = 1; ibin <= nbinx; ibin++)
    {
        double w = hnpart->GetBinContent(ibin);
        if (w <= 0.0) continue;
            
        double we = hnpart->GetBinError(ibin);
            
        double npart = hnpart->GetBinCenter(ibin);          
            
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
