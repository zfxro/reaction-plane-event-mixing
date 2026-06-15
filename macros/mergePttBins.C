//Run this before mergeXiBins.
//dphi histograms are added together; ntrg_pttcentdpsi must be added together.

const static int nxi = 6;
const float limxi[nxi+1] = {0.0,0.4,0.8,1.2,1.6,2.0,2.4};

const static int npta = 5;
const float limpta[npta+1] = {0.5,0.75,1.25,2.,3.,5.};

const static int ncent = 3;
const int limcent[ncent+1] = {0,20,60,100};

const static int nptt = 4;
const float limptt[nptt+1] = {5.,7.,9.,12.,15.};

const static int ndpsi = 2;

const static int ntype = 3;
const char* spec[ntype+1] = {"dir","dec","pi0","inc"};
enum TRIG
{
    DIR=0,
    DEC=1,
    PI0=2,
    INC=3
};

void mergePttBins(
    string inname = "gamtree_13466_ncent3_psimix.root",
    int irebin = 2
    )
{
    int dnptt = -9999;
    if (irebin == 1)
    {
        dnptt = 4;
    }
    else if (irebin == 2)
    {
        dnptt = 2;
    }
    else
    {
        cerr<<"Invalid value, "<<irebin<<", for irebin. Must be 1 or 2."<<endl;
    }

    //--------------- Root File Names ---------------
    string outname = inname;

    int pos = outname.find("13466") + 5;

    ostringstream tag;
    tag<<"_nptt"<<irebin;
    outname.insert(pos,tag.str());

    cout<<"Input = "<<inname<<endl;
    cout<<"Output = "<<outname<<endl;

    TFile* output = new TFile(outname.c_str(),"RECREATE");
    TFile* input = new TFile(inname.c_str(),"READ");
    //--------------------------------------------------

    TH1F* h1ntrg[ntype][nptt/dnptt][ncent+1][ndpsi+1];
    TH1F* h1ntrg_adj[ntype][nptt/dnptt][ncent+1][ndpsi+1];

    for (int itype = 0; itype < ntype; itype++)
    {
        for (int jptt = 0; jptt < nptt/dnptt; jptt++)
        {
            for (int icent = 0; icent < ncent+1; icent++)
            {
                for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
                {
                    h1ntrg[itype][jptt][icent][idpsi] = NULL;
                    h1ntrg_adj[itype][jptt][icent][idpsi] = NULL;
                }
            }
        }
    }

    for (int itype = 0; itype < ntype; itype++)
    {
        for (int icent = 0; icent < ncent+1; icent++)
        {
            for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
            {
                for (int iptt = 0; iptt < nptt; iptt++)
                {
                    int jptt = iptt/dnptt;

                    TH1F* h1tmp = (TH1F*) gROOT->FindObject(Form("ntrg_%s_pttb%i_centb%i_dpsib%i",spec[itype],iptt,icent,idpsi));

                    if(!h1ntrg[itype][jptt][icent][idpsi])
                    {
                        h1ntrg[itype][jptt][icent][idpsi] = (TH1F*) h1tmp->Clone(Form("ntrg_new_%s_pttb%i_centb%i_dpsib%i",spec[itype],jptt,icent,idpsi));
                        h1ntrg[itype][jptt][icent][idpsi]->SetDirectory(output);
                    }
                    else
                    {
                        h1ntrg[itype][jptt][icent][idpsi]->Add(h1tmp);
                    }

                    if (itype == PI0) continue;

                    TH1F* h1tmp_adj = (TH1F*) gROOT->FindObject(Form("ntrg_adj_%s_pttb%i_centb%i_dpsib%i",spec[itype],iptt,icent,idpsi));

                    if(!h1ntrg_adj[itype][jptt][icent][idpsi])
                    {
                        h1ntrg_adj[itype][jptt][icent][idpsi] = (TH1F*) h1tmp_adj->Clone(Form("ntrg_adj_new_%s_pttb%i_centb%i_dpsib%i",spec[itype],jptt,icent,idpsi));
                        h1ntrg_adj[itype][jptt][icent][idpsi]->SetDirectory(output);
                    }
                    else
                    {
                        h1ntrg_adj[itype][jptt][icent][idpsi]->Add(h1tmp_adj);
                    }

                }
            }
        }
    }

    //===========================================================================

    vector<vector<vector<vector<vector<vector<TH1F*> > > > > > h1dphi_xi(ntype+1, vector<vector<vector<vector<vector<TH1F*> > > > >(3, vector<vector<vector<vector<TH1F*> > > >(nptt/dnptt, vector<vector<vector<TH1F*> > >(ncent+1, vector<vector<TH1F*> >(ndpsi+1, vector<TH1F*>(nxi))))));

    vector<vector<vector<vector<vector<vector<TH1F*> > > > > > h1dphi_xi_adj(ntype+1, vector<vector<vector<vector<vector<TH1F*> > > > >(3, vector<vector<vector<vector<TH1F*> > > >(nptt/dnptt, vector<vector<vector<TH1F*> > >(ncent+1, vector<vector<TH1F*> >(ndpsi+1, vector<TH1F*>(nxi))))));

    vector<vector<vector<vector<vector<vector<TH1F*> > > > > > h1dphi(ntype+1, vector<vector<vector<vector<vector<TH1F*> > > > >(3, vector<vector<vector<vector<TH1F*> > > >(nptt/dnptt, vector<vector<vector<TH1F*> > >(ncent+1, vector<vector<TH1F*> >(ndpsi+1, vector<TH1F*>(nxi))))));

    vector<vector<vector<vector<vector<vector<TH1F*> > > > > > h1dphi_adj(ntype+1, vector<vector<vector<vector<vector<TH1F*> > > > >(3, vector<vector<vector<vector<TH1F*> > > >(nptt/dnptt, vector<vector<vector<TH1F*> > >(ncent+1, vector<vector<TH1F*> >(ndpsi+1, vector<TH1F*>(nxi))))));

    for (int itype = 0; itype < ntype; itype++)
    {
        for (int imix = 0; imix < 2; imix++)
        {
            for (int icent = 0; icent <= ncent; icent++)
            {
                int locent;
                int hicent;
                if (icent == ncent)
                {
                    locent = limcent[0];
                    hicent = limcent[ncent];
                }
                else
                {
                    locent = limcent[icent];
                    hicent = limcent[icent+1];
                }

                for (int idpsi = 0; idpsi <= ndpsi; idpsi++)
                {
                    for (int iptt = 0; iptt < nptt; iptt++)
                    {
                        int jptt = iptt/dnptt;

                        for (int ixi = 0; ixi < nxi; ixi++)
                        {

                            TH1F* h1tmp = (TH1F*) gROOT->FindObject(Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],imix,iptt,icent,idpsi,ixi));

                            if(!h1dphi_xi[itype][imix][jptt][icent][idpsi][ixi])
                            {
                                h1dphi_xi[itype][imix][jptt][icent][idpsi][ixi] = (TH1F*) h1tmp->Clone(Form("dphi_new_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],imix,jptt,icent,idpsi,ixi));
                                h1dphi_xi[itype][imix][jptt][icent][idpsi][ixi]->SetDirectory(output);
                            }
                            else
                            {
                                h1dphi_xi[itype][imix][jptt][icent][idpsi][ixi]->Add(h1tmp);
                            }

                            if (itype != PI0)
                            {
                                TH1F* h1tmp_adj = (TH1F*) gROOT->FindObject(Form("dphi_adj_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],imix,iptt,icent,idpsi,ixi));

                                if(!h1dphi_xi_adj[itype][imix][jptt][icent][idpsi][ixi])
                                {
                                    h1dphi_xi_adj[itype][imix][jptt][icent][idpsi][ixi] = (TH1F*) h1tmp_adj->Clone(Form("dphi_adj_new_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],imix,jptt,icent,idpsi,ixi));
                                    h1dphi_xi_adj[itype][imix][jptt][icent][idpsi][ixi]->SetDirectory(output);
                                }
                                else
                                {
                                    h1dphi_xi_adj[itype][imix][jptt][icent][idpsi][ixi]->Add(h1tmp_adj);
                                }
                            }
                        }

                        for (int ipta = 0; ipta < npta; ipta++)
                        {

                            TH1F* h1tmp = (TH1F*) gROOT->FindObject(Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",spec[itype],imix,iptt,icent,idpsi,ipta));

                            if(!h1dphi[itype][imix][jptt][icent][idpsi][ipta])
                            {
                                h1dphi[itype][imix][jptt][icent][idpsi][ipta] = (TH1F*) h1tmp->Clone(Form("dphi_new_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",spec[itype],imix,jptt,icent,idpsi,ipta));
                                h1dphi[itype][imix][jptt][icent][idpsi][ipta]->SetDirectory(output);
                            }
                            else
                            {
                                h1dphi[itype][imix][jptt][icent][idpsi][ipta]->Add(h1tmp);
                            }

                            if (itype != PI0)
                            {
                                TH1F* h1tmp_adj = (TH1F*) gROOT->FindObject(Form("dphi_adj_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",spec[itype],imix,iptt,icent,idpsi,ipta));

                                if(!h1dphi_adj[itype][imix][jptt][icent][idpsi][ipta])
                                {
                                    h1dphi_adj[itype][imix][jptt][icent][idpsi][ipta] = (TH1F*) h1tmp_adj->Clone(Form("dphi_adj_new_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",spec[itype],imix,jptt,icent,idpsi,ipta));
                                    h1dphi_adj[itype][imix][jptt][icent][idpsi][ipta]->SetDirectory(output);
                                }
                                else
                                {
                                    h1dphi_adj[itype][imix][jptt][icent][idpsi][ipta]->Add(h1tmp_adj);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //==================================================
    // also for pta histos binned in xi
    //==================================================

    TH1F* h1pta[ntype][nptt/dnptt][ncent+1][ndpsi+1][nxi];
    TH1F* h1pta_adj[ntype][nptt/dnptt][ncent+1][ndpsi+1][nxi];

    for (int itype = 0; itype < ntype; itype++)
    {
        for (int jptt = 0; jptt < nptt/dnptt; jptt++)
        {
            for (int icent = 0; icent < ncent+1; icent++)
            {
                for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
                {
                    for (int ixi = 0; ixi < nxi; ixi++)
                    {
                        h1pta[itype][jptt][icent][idpsi][ixi] = NULL;
                        h1pta_adj[itype][jptt][icent][idpsi][ixi] = NULL;
                    }
                }
            }
        }
    }

    for (int itype = 0; itype < ntype; itype++)
    {
        for (int icent = 0; icent < ncent+1; icent++)
        {
            for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
            {
                for (int ixi = 0; ixi < nxi; ixi++)
                {
                    for (int iptt = 0; iptt < nptt; iptt++)
                    {
                        int jptt = iptt/dnptt;

                        TH1F* h1tmp = (TH1F*) gROOT->FindObject(Form("pta_%s_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],iptt,icent,idpsi,ixi));

                        if(!h1pta[itype][jptt][icent][idpsi][ixi])
                        {
                            h1pta[itype][jptt][icent][idpsi][ixi] = (TH1F*) h1tmp->Clone(Form("pta_new_%s_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],jptt,icent,idpsi,ixi));
                            h1pta[itype][jptt][icent][idpsi][ixi]->SetDirectory(output);
                        }
                        else
                        {
                            h1pta[itype][jptt][icent][idpsi][ixi]->Add(h1tmp);
                        }

                        if (itype == PI0) continue;

                        TH1F* h1tmp_adj = (TH1F*) gROOT->FindObject(Form("pta_adj_%s_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],iptt,icent,idpsi,ixi));

                        if(!h1pta_adj[itype][jptt][icent][idpsi][ixi])
                        {
                            h1pta_adj[itype][jptt][icent][idpsi][ixi] = (TH1F*) h1tmp_adj->Clone(Form("pta_adj_new_%s_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],jptt,icent,idpsi,ixi));
                            h1pta_adj[itype][jptt][icent][idpsi][ixi]->SetDirectory(output);
                        }
                        else
                        {
                            h1pta_adj[itype][jptt][icent][idpsi][ixi]->Add(h1tmp_adj);
                        }

                    }
                }
            }
        }
    }

    //===========================================================================

    input->Close();

    cout<<"List of histograms to be saved:"<<endl;

    for (int itype = 0; itype < ntype; itype++)
    {
        for (int jptt = 0; jptt < nptt/dnptt; jptt++)
        {
            for (int icent = 0; icent < ncent+1; icent++)
            {
                for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
                {
                    h1ntrg[itype][jptt][icent][idpsi]->SetName(Form("ntrg_%s_pttb%i_centb%i_dpsib%i",spec[itype],jptt,icent,idpsi));
                    cout<<h1ntrg[itype][jptt][icent][idpsi]->GetName()<<endl;

                    if (itype != PI0)
                    {
                        h1ntrg_adj[itype][jptt][icent][idpsi]->SetName(Form("ntrg_adj_%s_pttb%i_centb%i_dpsib%i",spec[itype],jptt,icent,idpsi));
                        cout<<h1ntrg_adj[itype][jptt][icent][idpsi]->GetName()<<endl;
                    }

                    for (int ixi = 0; ixi < nxi; ixi++)
                    {
                        h1pta[itype][jptt][icent][idpsi][ixi]->SetName(Form("pta_%s_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],jptt,icent,idpsi,ixi));
                        cout<<h1pta[itype][jptt][icent][idpsi][ixi]->GetName()<<endl;

                        if (itype != PI0)
                        {
                            h1pta_adj[itype][jptt][icent][idpsi][ixi]->SetName(Form("pta_adj_%s_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],jptt,icent,idpsi,ixi));
                            cout<<h1pta_adj[itype][jptt][icent][idpsi][ixi]->GetName()<<endl;
                        }
                    }

                    for (int imix = 0; imix < 2; imix++)
                    {

                        for (int ixi = 0; ixi < nxi; ixi++)
                        {
                            h1dphi_xi[itype][imix][jptt][icent][idpsi][ixi]->SetName(Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],imix,jptt,icent,idpsi,ixi));
                            cout<<h1dphi_xi[itype][imix][jptt][icent][idpsi][ixi]->GetName()<<endl;

                            if (itype != PI0)
                            {
                                h1dphi_xi_adj[itype][imix][jptt][icent][idpsi][ixi]->SetName(Form("dphi_adj_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],imix,jptt,icent,idpsi,ixi));
                                cout<<h1dphi_xi_adj[itype][imix][jptt][icent][idpsi][ixi]->GetName()<<endl;
                            }
                        }

                        for (int ipta = 0; ipta < npta; ipta++)
                        {
                            h1dphi[itype][imix][jptt][icent][idpsi][ipta]->SetName(Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",spec[itype],imix,jptt,icent,idpsi,ipta));
                            cout<<h1dphi[itype][imix][jptt][icent][idpsi][ipta]->GetName()<<endl;

                            if (itype != PI0)
                            {
                                h1dphi_adj[itype][imix][jptt][icent][idpsi][ipta]->SetName(Form("dphi_adj_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",spec[itype],imix,jptt,icent,idpsi,ipta));
                                cout<<h1dphi_adj[itype][imix][jptt][icent][idpsi][ipta]->GetName()<<endl;
                            }
                        }

                    }
                }
            }
        }
    }

    cout<<"Writing "<<output->GetName()<<endl;
    output->Write();
    cout<<"Closing "<<output->GetName()<<endl;
    output->Close();

    return;
}
