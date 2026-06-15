//Make sure binning in gamAnaTree.h matches input root file
//dphi histograms are added together; ntrg_pttcentdpsi must be added together.

namespace GAT {
    #include "gamAnaTree.h"
}

const static int nxi = 3;
const float limxi[nxi+1] = {0.0,0.8,1.6,2.4};

const int nxi_old = 6;
const float limxi_old[nxi_old+1] = {0.0,0.4,0.8,1.2,1.6,2.0,2.4};

void mergeXiBins(
    string inname = "gamtree_13466_nptt2_ncent3_psimix.root"
    )
{

    if (nxi_old%nxi!=0)
    {
        cerr<<"nxi_old "<<nxi_old<<" is not multiple of nxi "<<nxi<<"\nExiting!"<<endl;
        return;
    }

    //--------------- Root File Names ---------------
    string outname = inname;

    int pos = outname.find("13466") + 5;

    ostringstream tag;
    tag<<"_nxi"<<nxi;
    outname.insert(pos,tag.str());

    cout<<"Input = "<<inname<<endl;
    cout<<"Output = "<<outname<<endl;

    TFile* output = new TFile(outname.c_str(),"RECREATE");
    TFile* input = new TFile(inname.c_str(),"READ");
    //--------------------------------------------------

    //===========================================================================

    vector<vector<vector<vector<vector<vector<TH1F*> > > > > > h1dphi_xi(GAT::ntype+1, vector<vector<vector<vector<vector<TH1F*> > > > >(3, vector<vector<vector<vector<TH1F*> > > >(GAT::nptt, vector<vector<vector<TH1F*> > >(GAT::ncent+1, vector<vector<TH1F*> >(GAT::ndpsi+1, vector<TH1F*>(nxi))))));

    vector<vector<vector<vector<vector<vector<TH1F*> > > > > > h1dphi_xi_adj(GAT::ntype+1, vector<vector<vector<vector<vector<TH1F*> > > > >(3, vector<vector<vector<vector<TH1F*> > > >(GAT::nptt, vector<vector<vector<TH1F*> > >(GAT::ncent+1, vector<vector<TH1F*> >(GAT::ndpsi+1, vector<TH1F*>(nxi))))));

    for (int itype = 0; itype < GAT::ntype; itype++)
    {
        for (int imix = 0; imix < 2; imix++)
        {
            for (int icent = 0; icent <= GAT::ncent; icent++)
            {
                int locent;
                int hicent;
                if (icent == GAT::ncent)
                {
                    locent = GAT::limcent[0];
                    hicent = GAT::limcent[GAT::ncent];
                }
                else
                {
                    locent = GAT::limcent[icent];
                    hicent = GAT::limcent[icent+1];
                }

                for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                {
                    for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    {
                        for (int ixi = 0; ixi < nxi_old; ixi++)
                        {
                            int jxi = ixi/(nxi_old/nxi);

                            TH1F* h1tmp = (TH1F*) gROOT->FindObject(Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],imix,iptt,icent,idpsi,ixi));

                            if(!h1dphi_xi[itype][imix][iptt][icent][idpsi][jxi])
                            {
                                h1dphi_xi[itype][imix][iptt][icent][idpsi][jxi] = (TH1F*) h1tmp->Clone(Form("dphi_new_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],imix,iptt,icent,idpsi,jxi));

                                h1dphi_xi[itype][imix][iptt][icent][idpsi][jxi]->SetTitle(Form("#splitline{dphi_%s_mix%i}{Ptt: [%.2f,%.2f)  Centrality: [%i,%i)  d#Psi_{2}: %s  Xi: [%.2f,%.2f)};#Delta#phi;((1/N_{trig}) (dN_{ch}/d#Delta#phi)",GAT::spec[itype],imix,GAT::limptt[iptt],GAT::limptt[iptt+1],locent,hicent,(GAT::dpsistr[idpsi]),limxi[jxi],limxi[jxi+1]));

                                h1dphi_xi[itype][imix][iptt][icent][idpsi][jxi]->SetDirectory(output);
                            }
                            else
                            {
                                h1dphi_xi[itype][imix][iptt][icent][idpsi][jxi]->Add(h1tmp);
                            }

                            if (itype != GAT::PI0)
                            {
                                TH1F* h1tmp_adj = (TH1F*) gROOT->FindObject(Form("dphi_adj_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],imix,iptt,icent,idpsi,ixi));

                                if(!h1dphi_xi_adj[itype][imix][iptt][icent][idpsi][jxi])
                                {
                                    h1dphi_xi_adj[itype][imix][iptt][icent][idpsi][jxi] = (TH1F*) h1tmp_adj->Clone(Form("dphi_adj_new_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],imix,iptt,icent,idpsi,jxi));

                                    h1dphi_xi_adj[itype][imix][iptt][icent][idpsi][jxi]->SetTitle(Form("#splitline{dphi_adj_%s_mix%i}{Ptt: [%.2f,%.2f)  Centrality: [%i,%i)  d#Psi_{2}: %s  Xi: [%.2f,%.2f)};#Delta#phi;((1/N_{trig}) (dN_{ch}/d#Delta#phi)",GAT::spec[itype],imix,GAT::limptt[iptt],GAT::limptt[iptt+1],locent,hicent,(GAT::dpsistr[idpsi]),limxi[jxi],limxi[jxi+1]));

                                    h1dphi_xi_adj[itype][imix][iptt][icent][idpsi][jxi]->SetDirectory(output);
                                }
                                else
                                {
                                    h1dphi_xi_adj[itype][imix][iptt][icent][idpsi][jxi]->Add(h1tmp_adj);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //===========================================================================
    // pta histos binned in xi
    //===========================================================================

    TH1F* h1pta[GAT::ntype][GAT::nptt][GAT::ncent+1][GAT::ndpsi+1][nxi];
    TH1F* h1pta_adj[GAT::ntype][GAT::nptt][GAT::ncent+1][GAT::ndpsi+1][nxi];

    for (int itype = 0; itype < GAT::ntype; itype++)
    {
        for (int icent = 0; icent <= GAT::ncent; icent++)
        {
            for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    for (int ixi = 0; ixi < nxi; ixi++)
                    {
                        h1pta[itype][iptt][icent][idpsi][ixi] = NULL;
                        h1pta_adj[itype][iptt][icent][idpsi][ixi] = NULL;
                    }
                }
            }
        }
    }


    for (int itype = 0; itype < GAT::ntype; itype++)
    {
        for (int icent = 0; icent <= GAT::ncent; icent++)
        {
            int locent;
            int hicent;
            if (icent == GAT::ncent)
            {
                locent = GAT::limcent[0];
                hicent = GAT::limcent[GAT::ncent];
            }
            else
            {
                locent = GAT::limcent[icent];
                hicent = GAT::limcent[icent+1];
            }

            for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    for (int ixi = 0; ixi < nxi_old; ixi++)
                    {
                        int jxi = ixi/(nxi_old/nxi);

                        TH1F* h1tmp = (TH1F*) gROOT->FindObject(Form("pta_%s_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],iptt,icent,idpsi,ixi));

                        if(!h1pta[itype][iptt][icent][idpsi][jxi])
                        {
                            h1pta[itype][iptt][icent][idpsi][jxi] = (TH1F*) h1tmp->Clone(Form("pta_new_%s_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],iptt,icent,idpsi,jxi));
                            h1pta[itype][iptt][icent][idpsi][jxi]->SetTitle(Form("#splitline{pta_%s}{Ptt: [%.2f,%.2f)  Centrality: [%i,%i)  d#Psi_{2,BBC}: %s  Xi: [%.2f,%.2f)};p_{T,a}",GAT::spec[itype],GAT::limptt[iptt],GAT::limptt[iptt+1],locent,hicent,(GAT::dpsistr[idpsi]),limxi[jxi],limxi[jxi+1]));
                            h1pta[itype][iptt][icent][idpsi][jxi]->SetDirectory(output);
                        }
                        else
                        {
                            h1pta[itype][iptt][icent][idpsi][jxi]->Add(h1tmp);
                        }

                        if (itype != GAT::PI0)
                        {
                            TH1F* h1tmp_adj = (TH1F*) gROOT->FindObject(Form("pta_adj_%s_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],iptt,icent,idpsi,ixi));

                            if(!h1pta_adj[itype][iptt][icent][idpsi][jxi])
                            {
                                h1pta_adj[itype][iptt][icent][idpsi][jxi] = (TH1F*) h1tmp->Clone(Form("pta_adj_new_%s_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],iptt,icent,idpsi,jxi));
                                h1pta_adj[itype][iptt][icent][idpsi][jxi]->SetTitle(Form("#splitline{pta_%s}{Ptt: [%.2f,%.2f)  Centrality: [%i,%i)  d#Psi_{2,BBC}: %s  Xi: [%.2f,%.2f)};p_{T,a}",GAT::spec[itype],GAT::limptt[iptt],GAT::limptt[iptt+1],locent,hicent,(GAT::dpsistr[idpsi]),limxi[jxi],limxi[jxi+1]));
                                h1pta_adj[itype][iptt][icent][idpsi][jxi]->SetDirectory(output);
                            }
                            else
                            {
                                h1pta_adj[itype][iptt][icent][idpsi][jxi]->Add(h1tmp);
                            }
                        }
                    }
                }
            }
        }
    }


    //===========================================================================

    //get all other nec histos
    for (int itype = 0; itype < GAT::ntype; itype++)
    {
        for (int icent = 0; icent < GAT::ncent+1; icent++)
        {
            for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    TH1F* h1tmp = (TH1F*) gROOT->FindObject(Form("ntrg_%s_pttb%i_centb%i_dpsib%i",GAT::spec[itype],iptt,icent,idpsi));
                    h1tmp->SetDirectory(output);

                    if (itype != GAT::PI0)
                    {
                        TH1F* h1tmp_adj = (TH1F*) gROOT->FindObject(Form("ntrg_adj_%s_pttb%i_centb%i_dpsib%i",GAT::spec[itype],iptt,icent,idpsi));
                        h1tmp_adj->SetDirectory(output);
                    }

                    for (int imix = 0; imix < 2; imix++)
                    {
                        for (int ipta = 0; ipta < GAT::npta; ipta++)
                        {
                            TH1F* h1ptatmp = (TH1F*) gROOT->FindObject(Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",GAT::spec[itype],imix,iptt,icent,idpsi,ipta));
                            h1ptatmp->SetDirectory(output);

                            if (itype != GAT::PI0)
                            {
                                TH1F* h1ptatmp_adj = (TH1F*) gROOT->FindObject(Form("dphi_adj_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",GAT::spec[itype],imix,iptt,icent,idpsi,ipta));
                                h1ptatmp_adj->SetDirectory(output);
                            }
                        }
                    }
                }
            }
        }
    }

    //===========================================================================

    input->Close();

    cout<<"List of histograms to be saved:"<<endl;

    for (int itype = 0; itype < GAT::ntype; itype++)
    {
        for (int iptt = 0; iptt < GAT::nptt; iptt++)
        {
            for (int icent = 0; icent < GAT::ncent+1; icent++)
            {
                for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                {
                    for (int ixi = 0; ixi < nxi; ixi++)
                    {
                        h1pta[itype][iptt][icent][idpsi][ixi]->SetName(Form("pta_%s_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],iptt,icent,idpsi,ixi));
                        cout<<h1pta[itype][iptt][icent][idpsi][ixi]->GetName()<<endl;

                        if (itype != GAT::PI0)
                        {
                            h1pta_adj[itype][iptt][icent][idpsi][ixi]->SetName(Form("pta_adj_%s_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],iptt,icent,idpsi,ixi));
                            cout<<h1pta_adj[itype][iptt][icent][idpsi][ixi]->GetName()<<endl;
                        }

                        for (int imix = 0; imix < 2; imix++)
                        {
                            h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi]->SetName(Form("dphi_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],imix,iptt,icent,idpsi,ixi));
                            cout<<h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi]->GetName()<<endl;

                            if (itype != GAT::PI0)
                            {
                                h1dphi_xi_adj[itype][imix][iptt][icent][idpsi][ixi]->SetName(Form("dphi_adj_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",GAT::spec[itype],imix,iptt,icent,idpsi,ixi));
                                cout<<h1dphi_xi_adj[itype][imix][iptt][icent][idpsi][ixi]->GetName()<<endl;
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
