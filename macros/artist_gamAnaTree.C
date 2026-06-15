// if rebinning ptt make sure to change:
//      - input file in GetNhExpXi
//      - gamAnaTree.h nptt
//      - as well of course gamAnaTree input files
// Look for keyword REBIN
//
// if rebinning xi make sure to change:
//      - gamAnaTree.h nxi
//      - as well of course gamAnaTree input files

#include "artist_gamAnaTree.h"

void artist_gamAnaTree(
    const char* inname = "input.root",
    const char* inname_nopsimix = "input_nopsimix.root",
    const char* outname = "artist_out.root",
    bool isadj = false
    )
{
    verbosity = 2;

    TFile* input_nopsimix = new TFile(inname_nopsimix,"READ");
    if (!input_nopsimix || input_nopsimix->IsZombie())
    {
        cout<<"Could not open "<<inname_nopsimix<<endl;
        return;
    }

    TFile *input = TFile::Open(inname);

    if (!input)
    {
        cout<<"Could not open file!"<<endl;
        return;
    }

    getCounters(isadj);

    //==================================================
    //get expected number of hadrons
    //==================================================

    expnh.resize(ncent+1);
    for (int icent = 0; icent <= ncent; icent++)
    {
        expnh[icent].resize(npta);
    }

    cout<<"Skipping any acceptance correction. expnh will be set to 1 for all centrality and hadron pt."<<endl;
    for (int icent = 0; icent <= ncent; icent++)
    {
        for (int ipta = 0; ipta < npta; ipta++)
        {
            expnh[icent][ipta] = 1.0;
        }
    }

    TH1F* h1expnh[ncent+1];
    for (int icent = 0; icent <= ncent; icent++)
    {
        h1expnh[icent] = new TH1F(Form("expnh%i",icent),Form("expnh%i",icent),npta,0,npta);
        for (int ipta = 0; ipta < npta; ipta++)
        {
            h1expnh[icent]->SetBinContent(ipta+1,expnh[icent][ipta]);
        }
        h1expnh[icent]->Draw();
    }

    expnhxi.resize(ntype+1);
    for (int itype = 0; itype <= ntype; itype++)
    {
        expnhxi[itype].resize(ncent+1);
        for (int icent = 0; icent <= ncent; icent++)
        {
            expnhxi[itype][icent].resize(ndpsi+1);
            for (int idpsi = 0; idpsi <= ndpsi; idpsi++)
            {
                expnhxi[itype][icent][idpsi].resize(nptt);
                for (int iptt = 0; iptt < nptt; iptt++)
                {
                    expnhxi[itype][icent][idpsi][iptt].resize(nxi);

                    for (int ixi = 0; ixi < nxi; ixi++)
                    {
                        expnhxi[itype][icent][idpsi][iptt][ixi] = 0.0;
                    }
                }
            }

        }
    }

    GetNhExpXi(isadj);

    //==================================================
    //obtained expected number of hadrons
    //==================================================


    //============================================================
    // Retrieve/Add dPhi histograms
    //============================================================
    cout<<"Retrieving/adding dphi histograms..."<<endl;
    dPhiHistos(input,input_nopsimix,isadj);
    dPhiHistosXi(input,input_nopsimix,isadj);
    cout<<"Retrieved/added dphi histograms successfully."<<endl;
    //============================================================

    //============================================================
    //Make inclusive histos
    //============================================================

    for (int iptt = 0; iptt < nptt; iptt++)
    {
        for (int icent = 0; icent < ncent+1; icent++)
        {
            for (int idpsi = 0; idpsi <= ndpsi; idpsi++)
            {
                for (int ipta = 0; ipta < npta; ipta++)
                {

                    h1dphi_nopsimix[ntype][iptt][icent][idpsi][ipta]->Add(h1dphi_nopsimix[0][iptt][icent][idpsi][ipta]);
                    h1dphi_nopsimix[ntype][iptt][icent][idpsi][ipta]->Add(h1dphi_nopsimix[1][iptt][icent][idpsi][ipta]);

                    for (int imix = 0; imix < 2; imix++)
                    {
                        h1dphi[ntype][imix][iptt][icent][idpsi][ipta]->Add(h1dphi[0][imix][iptt][icent][idpsi][ipta]);
                        h1dphi[ntype][imix][iptt][icent][idpsi][ipta]->Add(h1dphi[1][imix][iptt][icent][idpsi][ipta]);

                    }
                }

                for (int ixi = 0; ixi < nxi; ixi++)
                {

                    h1dphi_xi_nopsimix[ntype][iptt][icent][idpsi][ixi]->Add(h1dphi_xi_nopsimix[0][iptt][icent][idpsi][ixi]);
                    h1dphi_xi_nopsimix[ntype][iptt][icent][idpsi][ixi]->Add(h1dphi_xi_nopsimix[1][iptt][icent][idpsi][ixi]);

                    if (useksic)
                    {
                        h1dphi_xi_nopsimix_ksic[ntype][iptt][icent][idpsi][ixi]->Add(h1dphi_xi_nopsimix_ksic[0][iptt][icent][idpsi][ixi]);
                        h1dphi_xi_nopsimix_ksic[ntype][iptt][icent][idpsi][ixi]->Add(h1dphi_xi_nopsimix_ksic[1][iptt][icent][idpsi][ixi]);
                    }

                    for (int imix = 0; imix < 2; imix++)
                    {
                        h1dphi_xi[ntype][imix][iptt][icent][idpsi][ixi]->Add(h1dphi_xi[0][imix][iptt][icent][idpsi][ixi]);
                        h1dphi_xi[ntype][imix][iptt][icent][idpsi][ixi]->Add(h1dphi_xi[1][imix][iptt][icent][idpsi][ixi]);

                    }
                }
            }
        }
    }
    //============================================================


    //Calculate remaining flow background
    GetRemBg();

    const float dphibw = PI/(ndphibin/rebinint);

    //============================================================
    // Normalize dphi histos (after inclusive dists made)
    //============================================================

    cout<<"Normalizing dphi histos"<<endl;

    for (int iptt = 0; iptt < nptt; iptt++)
    {
        for (int icent = 0; icent < ncent+1; icent++)
        {
            for (int idpsi = 0; idpsi <= ndpsi; idpsi++)
            {
                for (int ipta = 0; ipta < npta; ipta++)
                {
                    float tmpscl = 1./(ntrg_pttcentdpsi[0][iptt][icent][idpsi] + ntrg_pttcentdpsi[1][iptt][icent][idpsi]);
                    for (int imix = 0; imix < 2; imix++)
                    {
                        h1dphi[ntype][imix][iptt][icent][idpsi][ipta]->Scale(tmpscl);
                        h1dphi[ntype][imix][iptt][icent][idpsi][ipta]->Scale(1./dphibw);
                    }

                    h1dphi_nopsimix[ntype][iptt][icent][idpsi][ipta]->Scale(tmpscl);
                    h1dphi_nopsimix[ntype][iptt][icent][idpsi][ipta]->Scale(1./dphibw);
                    h1dphi_nopsimix[ntype][iptt][icent][idpsi][ipta]->Scale(1./expnh[icent][ipta]);

                    h1dphi_bg[ntype][iptt][icent][idpsi][ipta]->Scale(tmpscl);
                    h1dphi_bg[ntype][iptt][icent][idpsi][ipta]->Scale(1./dphibw);

                    for (int itype = 0; itype < ntype; itype++)
                    {
                        if (skippi0 && itype==PI0) continue;

                        tmpscl = 1./ntrg_pttcentdpsi[itype][iptt][icent][idpsi];
                        for (int imix = 0; imix < 2; imix++)
                        {
                            h1dphi[itype][imix][iptt][icent][idpsi][ipta]->Scale(tmpscl);
                            h1dphi[itype][imix][iptt][icent][idpsi][ipta]->Scale(1./dphibw);
                        }

                        h1dphi_nopsimix[itype][iptt][icent][idpsi][ipta]->Scale(tmpscl);
                        h1dphi_nopsimix[itype][iptt][icent][idpsi][ipta]->Scale(1./dphibw);
                        h1dphi_nopsimix[itype][iptt][icent][idpsi][ipta]->Scale(1./expnh[icent][ipta]);


                        h1dphi_bg[itype][iptt][icent][idpsi][ipta]->Scale(tmpscl);
                        h1dphi_bg[itype][iptt][icent][idpsi][ipta]->Scale(1./dphibw);
                    }
                }
            }
        }
    }

    //=============== Also for xi binned ===============
    for (int iptt = 0; iptt < nptt; iptt++)
    {
        for (int icent = 0; icent < ncent+1; icent++)
        {
            for (int idpsi = 0; idpsi <= ndpsi; idpsi++)
            {
                for (int ixi = 0; ixi < nxi; ixi++)
                {
                    float tmpscl = 1./(ntrg_pttcentdpsi[0][iptt][icent][idpsi] + ntrg_pttcentdpsi[1][iptt][icent][idpsi]);
                    for (int imix = 0; imix < 2; imix++)
                    {
                        h1dphi_xi[ntype][imix][iptt][icent][idpsi][ixi]->Scale(tmpscl);
                        h1dphi_xi[ntype][imix][iptt][icent][idpsi][ixi]->Scale(1./dphibw);
                    }

                    h1dphi_xi_nopsimix[ntype][iptt][icent][idpsi][ixi]->Scale(tmpscl);
                    h1dphi_xi_nopsimix[ntype][iptt][icent][idpsi][ixi]->Scale(1./dphibw);
                    h1dphi_xi_nopsimix[ntype][iptt][icent][idpsi][ixi]->Scale(1./expnhxi[ntype][icent][idpsi][iptt][ixi]);

                    h1dphi_xi_bg[ntype][iptt][icent][idpsi][ixi]->Scale(tmpscl);
                    h1dphi_xi_bg[ntype][iptt][icent][idpsi][ixi]->Scale(1./dphibw);

                    for (int itype = 0; itype < ntype; itype++)
                    {
                        if (skippi0 && itype==PI0) continue;

                        tmpscl = 1./ntrg_pttcentdpsi[itype][iptt][icent][idpsi];
                        for (int imix = 0; imix < 2; imix++)
                        {
                            h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi]->Scale(tmpscl);
                            h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi]->Scale(1./dphibw);
                        }

                        h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi]->Scale(tmpscl);
                        h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi]->Scale(1./dphibw);
                        h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi]->Scale(1./expnhxi[itype][icent][idpsi][iptt][ixi]);


                        h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi]->Scale(tmpscl);
                        h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi]->Scale(1./dphibw);
                    }
                }
            }
        }
    }


    //============================================================

    MixSubNorm();

    SubRemBg();

    TFile* output = new TFile(outname,"RECREATE");

    for (int itype = 0; itype < ntype+1; itype++)
    {
        if (skippi0 && itype==PI0) continue;

        for (int icent = 0; icent < ncent+1; icent++)
        {
            for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
            {
                for (int iptt = 0; iptt < nptt; iptt++)
                {

                    for (int ipta = 0; ipta < npta; ipta++)
                    {
                        for (int imixt = 0; imixt < NMIXTYPE; imixt++)
                        {
                            h1dphi[itype][imixt][iptt][icent][idpsi][ipta]->SetDirectory(output);
                        }
                    }

                    for (int ixi = 0; ixi < nxi; ixi++)
                    {
                        h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi]->SetDirectory(output);
                        h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi]->SetDirectory(output);
                        for (int imixt = 0; imixt < NMIXTYPE; imixt++)
                        {
                            h1dphi_xi[itype][imixt][iptt][icent][idpsi][ixi]->SetDirectory(output);

                        }

                    }

                }
            }
        }
    }

    output->Write();
    output->Close();

    cout<<"All done!"<<endl;
    return;
}

void GetNhExpXi(int isadj)
{
    cout<<"Getting NhExpXi"<<endl;

    for (int itype = 0; itype < ntype; itype++)
    {
        for (int icent = 0; icent <= ncent; icent++)
        {
            for (int idpsi = 0; idpsi <= ndpsi; idpsi++)
            {
                ostringstream ossdir;
                if (itype == PI0)
                    ossdir << "adj0_pi0";
                else
                {
                    ossdir << "adj" << isadj;
                    if (itype == DIR)
                        ossdir<<"_dir";
                    else if (itype == DEC)
                        ossdir<<"_dec";
                }

                ossdir << "_cent" << icent;
                ossdir << "_dpsi" << idpsi;

                ostringstream ossin;
                ossin<<"./expnh_xi/"<<ossdir.str()<<"/expnh_xi_"<<ossdir.str()<<".txt";

                ifstream ifs(ossin.str(),ifstream::in);
                if (!ifs.is_open())
                {
                    cerr<<"Cannot open "<<ossin.str()<<endl;
                    return;
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

                    expnhxi[itype][icent][idpsi][iptt][ixi] = tmpnh*detah*2.;
                }
            }
        }
    }

    for (int icent = 0; icent <= ncent; icent++)
    {
        for (int idpsi = 0; idpsi <= ndpsi; idpsi++)
        {
            for (int iptt = 0; iptt < nptt; iptt++)
            {
                for (int ixi = 0; ixi < nxi; ixi++)
                {

                    float tmpnh = 0.;
                    float tmpntrg = 0.;
                    for (int itype = 0; itype < 2; itype++)
                    {
                        tmpnh +=  ntrg_pttcentdpsi[itype][iptt][icent][idpsi]*expnhxi[itype][icent][idpsi][iptt][ixi];
                        tmpntrg += ntrg_pttcentdpsi[itype][iptt][icent][idpsi];
                    }

                    if (verbosity > 1)
                    {
                        cout<<"icent: "<<icent<<"  idpsi: "<<idpsi<<"  iptt: "<<iptt<<"  ixi: "<<ixi<<endl;
                        cout<<"tmpnh: "<<tmpnh<<"  tmpntrg: "<<tmpntrg<<endl;
                    }

                    tmpnh/=tmpntrg;

                    if (verbosity > 1)
                    {
                        cout<<"expnhxi["<<ntype<<"]["<<icent<<"]["<<idpsi<<"]["<<iptt<<"]["<<ixi<<"] = tmpnh/tmpntrg = "<<tmpnh<<endl;
                    }

                    expnhxi[ntype][icent][idpsi][iptt][ixi] = tmpnh;
                }
            }
        }
    }

    cout<<"Got NhExpXi"<<endl;
}

void dPhiHistos(TFile* input, TFile* input_nopsimix, bool isadj)
{
    for (int itype = 0; itype < ntype+1; itype++)
    {
        if (skippi0 && itype==PI0) continue;

        string prefix = "dphi";
        if (isadj && itype!=PI0)
            prefix.append("_adj");

        for (int iptt = 0; iptt < nptt; iptt++)
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
                    for (int ipta = 0; ipta < npta; ipta++)
                    {
                        if (itype == ntype)
                        {
                            h1dphi_nopsimix[itype][iptt][icent][idpsi][ipta] = new TH1F (
                                Form("%s_nopsimix_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",prefix.c_str(),spec[itype],MIX,iptt,icent,idpsi,ipta),"",ndphibin,mindphi,maxdphi);
                            h1dphi_nopsimix[itype][iptt][icent][idpsi][ipta]->Sumw2();
                        }
                        else
                        {
                            input_nopsimix->cd();

                            TH1F* h1tmp_nopsimix = (TH1F*) gROOT->FindObject(Form("%s_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",prefix.c_str(),spec[itype],MIX,iptt,icent,idpsi,ipta));
                            if (h1tmp_nopsimix == NULL)
                            {
                                cout<<__LINE__;
                                printf(": NoPsiMix %s_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i not found!\n",prefix.c_str(),spec[itype],MIX,iptt,icent,idpsi,ipta);
                                return;
                            }
                            h1dphi_nopsimix[itype][iptt][icent][idpsi][ipta] = (TH1F*) h1tmp_nopsimix->Clone(Form("%s_nopsimix_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",prefix.c_str(),spec[itype],MIX,iptt,icent,idpsi,ipta));

                            input->cd();
                        }

                        h1dphi_nopsimix[itype][iptt][icent][idpsi][ipta]->SetTitle(Form("NoPsiMix! %s p_{tt}:[%.2f,%.2f) Cent:[%i,%i) d#Psi_{2,BBC}:%s  p_{ta}:[%.2f,%.2f);#Delta#phi (#pi rad);#frac{1}{N_{trg}}#frac{dN_{pair}}{d(#Delta#phi)}",spec[itype],limptt[iptt],limptt[iptt+1],locent,hicent,(dpsistr[idpsi]),limpta[ipta],limpta[ipta+1]));
                        h1dphi_nopsimix[itype][iptt][icent][idpsi][ipta]->Rebin(rebinint);


                        for (int imix = 0; imix < NMIXTYPE; imix++)
                        {
                            if (imix == MSN || itype == INC)
                            {
                                h1dphi[itype][imix][iptt][icent][idpsi][ipta] = new TH1F (
                                    Form("%s_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",prefix.c_str(),spec[itype],imix,iptt,icent,idpsi,ipta),
                                    Form("p_{tt}:[%.2f,%.2f) Cent:[%i,%i) d#Psi_{2,BBC}:%s  p_{ta}:[%.2f,%.2f);#Delta#phi (#pi rad);#frac{1}{N_{trg}}#frac{dN_{pair}}{d(#Delta#phi)}",limptt[iptt],limptt[iptt+1],locent,hicent,(dpsistr[idpsi]),limpta[ipta],limpta[ipta+1]),
                                    ndphibin/rebinint,mindphi,maxdphi);
                                h1dphi[itype][imix][iptt][icent][idpsi][ipta]->Sumw2();

                                continue;
                            }

                            h1dphi[itype][imix][iptt][icent][idpsi][ipta] = (TH1F*) gROOT->FindObject(Form("%s_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i",prefix.c_str(),spec[itype],imix,iptt,icent,idpsi,ipta));
                            if (h1dphi[itype][imix][iptt][icent][idpsi][ipta] == NULL)
                            {
                                cout<<__LINE__;
                                printf(": %s_%s_mix%i_ptt%i_cent%i_dpsi%i_pta%i not found!\n",prefix.c_str(),spec[itype],imix,iptt,icent,idpsi,ipta);
                                return;
                            }
                            h1dphi[itype][imix][iptt][icent][idpsi][ipta]->SetTitle(Form("p_{tt}:[%.2f,%.2f) Cent:[%i,%i) d#Psi_{2,BBC}:%s  p_{ta}:[%.2f,%.2f);#Delta#phi (#pi rad);#frac{1}{N_{trg}}#frac{dN_{pair}}{d(#Delta#phi)}",limptt[iptt],limptt[iptt+1],locent,hicent,(dpsistr[idpsi]),limpta[ipta],limpta[ipta+1]));
                            h1dphi[itype][imix][iptt][icent][idpsi][ipta]->Rebin(rebinint);
                        }
                    }
                }
            }
        }
    }
}

void dPhiHistosXi(TFile* input, TFile* input_nopsimix, bool isadj)
{
    for (int itype = 0; itype < ntype+1; itype++)
    {
        if (skippi0 && itype==PI0) continue;

        string prefix = "dphi";
        string prefix_ksic = "dphi_ksic";
        if (isadj && itype!=PI0)
        {
            prefix.append("_adj");
            prefix_ksic.append("_adj");
        }

        for (int iptt = 0; iptt < nptt; iptt++)
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
                    for (int ixi = 0; ixi < nxi; ixi++)
                    {

                        if (itype == ntype)
                        {
                            h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi] = new TH1F (
                                Form("%s_nopsimix_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",prefix.c_str(),spec[itype],MIX,iptt,icent,idpsi,ixi),"",ndphibin,mindphi,maxdphi);
                            h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi]->Sumw2();

                            if (useksic)
                            {
                                h1dphi_xi_nopsimix_ksic[itype][iptt][icent][idpsi][ixi] = new TH1F (
                                    Form("%s_nopsimix_ksic_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",prefix.c_str(),spec[itype],MIX,iptt,icent,idpsi,ixi),"",ndphibin,mindphi,maxdphi);
                                h1dphi_xi_nopsimix_ksic[itype][iptt][icent][idpsi][ixi]->Sumw2();
                            }

                        }
                        else
                        {
                            input_nopsimix->cd();

                            TH1F* h1tmp_nopsimix = (TH1F*) gROOT->FindObject(Form("%s_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",prefix.c_str(),spec[itype],MIX,iptt,icent,idpsi,ixi));
                            if (h1tmp_nopsimix == NULL)
                            {
                                cout<<__LINE__;
                                printf(": NoPsiMix %s_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i not found!\n",prefix.c_str(),spec[itype],MIX,iptt,icent,idpsi,ixi);
                                return;
                            }
                            h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi] = (TH1F*) h1tmp_nopsimix->Clone(Form("%s_nopsimix_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",prefix.c_str(),spec[itype],MIX,iptt,icent,idpsi,ixi));

                            if (useksic)
                            {
                                TH1F* h1tmp_nopsimix_ksic = (TH1F*) gROOT->FindObject(Form("%s_%s_ptt%i_cent%i_dpsi%i_xi%i",prefix_ksic.c_str(),spec[itype],iptt,icent,idpsi,ixi));
                                if (h1tmp_nopsimix_ksic == NULL)
                                {
                                    cout<<__LINE__;
                                    printf(": Nopsimix_Ksic %s_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i not found!\n",prefix.c_str(),spec[itype],MIX,iptt,icent,idpsi,ixi);
                                    return;
                                }
                                h1dphi_xi_nopsimix_ksic[itype][iptt][icent][idpsi][ixi] = (TH1F*) h1tmp_nopsimix_ksic->Clone(Form("%s_nopsimix_ksic_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",prefix.c_str(),spec[itype],MIX,iptt,icent,idpsi,ixi));
                            }

                            input->cd();
                        }

                        h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi]->SetTitle(Form("NoPsiMix! %s p_{tt}:[%.2f,%.2f) Cent:[%i,%i) d#Psi_{2,BBC}:%s  #xi:[%.2f,%.2f);#Delta#phi (#pi rad);#frac{1}{N_{trg}}#frac{dN_{pair}}{d(#Delta#phi)}",spec[itype],limptt[iptt],limptt[iptt+1],locent,hicent,(dpsistr[idpsi]),limxi[ixi],limxi[ixi+1]));
                        h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi]->Rebin(rebinint);

                        if (useksic)
                            h1dphi_xi_nopsimix_ksic[itype][iptt][icent][idpsi][ixi]->Rebin(rebinint);


                        for (int imix = 0; imix < NMIXTYPE; imix++)
                        {
                            if (imix == MSN || itype == INC)
                            {
                                h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi] = new TH1F (
                                    Form("%s_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",prefix.c_str(),spec[itype],imix,iptt,icent,idpsi,ixi),
                                    Form("p_{tt}:[%.2f,%.2f) Cent:[%i,%i) d#Psi_{2,BBC}:%s  #xi:[%.2f,%.2f);#Delta#phi (#pi rad);#frac{1}{N_{trg}}#frac{dN_{pair}}{d(#Delta#phi)}",limptt[iptt],limptt[iptt+1],locent,hicent,(dpsistr[idpsi]),limxi[ixi],limxi[ixi+1]),
                                    ndphibin/rebinint,mindphi,maxdphi);
                                h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi]->Sumw2();

                                continue;
                            }

                            if (useksic && imix == MIX)
                            {
                                h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi] = (TH1F*) gROOT->FindObject(Form("%s_%s_ptt%i_cent%i_dpsi%i_xi%i",prefix_ksic.c_str(),spec[itype],iptt,icent,idpsi,ixi));
                                if (!h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi])
                                {
                                    cout<<__LINE__;
                                    printf(": %s_%s_ptt%i_cent%i_dpsi%i_xi%i not found!",prefix_ksic.c_str(),spec[itype],iptt,icent,idpsi,ixi);
                                    return;
                                }
                            }
                            else
                            {
                                h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi] = (TH1F*) gROOT->FindObject(Form("%s_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i",prefix.c_str(),spec[itype],imix,iptt,icent,idpsi,ixi));
                                if (h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi] == NULL)
                                {
                                    cout<<__LINE__;
                                    printf(": %s_%s_mix%i_ptt%i_cent%i_dpsi%i_xi%i not found!\n",prefix.c_str(),spec[itype],imix,iptt,icent,idpsi,ixi);
                                    return;
                                }
                            }
                            h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi]->SetTitle(Form("p_{tt}:[%.2f,%.2f) Cent:[%i,%i) d#Psi_{2,BBC}:%s  #xi:[%.2f,%.2f);#Delta#phi (#pi rad);#frac{1}{N_{trg}}#frac{dN_{pair}}{d(#Delta#phi)}",limptt[iptt],limptt[iptt+1],locent,hicent,(dpsistr[idpsi]),limxi[ixi],limxi[ixi+1]));
                            h1dphi_xi[itype][imix][iptt][icent][idpsi][ixi]->Rebin(rebinint);
                        }
                    }
                }
            }
        }
}


void MixSubNorm()
{
    cout<<"Subtracting mixed events and correcting for hadron acceptance"<<endl;

    for (int itype = 0; itype < ntype+1; itype++)
    {
        if (skippi0 && itype==PI0) continue;
        for (int iptt = 0; iptt < nptt; iptt++)
        {

            for (int icent = 0; icent < ncent+1; icent++)
            {
                for (int idpsi = 0; idpsi <= ndpsi; idpsi++)
                {
                    for (int ipta = 0; ipta < npta; ipta++)
                    {
                        h1dphi[itype][MSN][iptt][icent][idpsi][ipta]->Add(h1dphi[itype][RAW][iptt][icent][idpsi][ipta]);

                        h1dphi[itype][MSN][iptt][icent][idpsi][ipta]->Add(h1dphi[itype][MIX][iptt][icent][idpsi][ipta],-1.);

                        h1dphi[itype][MSN][iptt][icent][idpsi][ipta]->Divide(h1dphi_nopsimix[itype][iptt][icent][idpsi][ipta]);

                        h1dphi[itype][MIX][iptt][icent][idpsi][ipta]->Divide(h1dphi_nopsimix[itype][iptt][icent][idpsi][ipta]);

                        h1dphi[itype][RAW][iptt][icent][idpsi][ipta]->Divide(h1dphi_nopsimix[itype][iptt][icent][idpsi][ipta]);

                    }

                    for (int ixi = 0; ixi < nxi; ixi++)
                    {
                        h1dphi_xi[itype][MSN][iptt][icent][idpsi][ixi]->Add(h1dphi_xi[itype][RAW][iptt][icent][idpsi][ixi]);

                        h1dphi_xi[itype][MSN][iptt][icent][idpsi][ixi]->Add(h1dphi_xi[itype][MIX][iptt][icent][idpsi][ixi],-1.);

                        h1dphi_xi[itype][MSN][iptt][icent][idpsi][ixi]->Divide(h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi]);

                        h1dphi_xi[itype][MIX][iptt][icent][idpsi][ixi]->Divide(h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi]);

                        h1dphi_xi[itype][RAW][iptt][icent][idpsi][ixi]->Divide(h1dphi_xi_nopsimix[itype][iptt][icent][idpsi][ixi]);

                    }

                }
            }
        }
    }
    return;

}

int getipad(int ic, int ip)
{
    if (ic == ncent) ic = -1;
    return (ip+1) + ((ncent-ic)-1)*(npta);
}

void getCounters(bool isadj)
{
    //============================================================
    // Get histo counter for triggers in pt and cent bins
    // Fill our counter array
    //============================================================

    memset(ntrg_pttcentdpsi,0,sizeof(ntrg_pttcentdpsi));

    for (int itype = 0; itype < ntype; itype++)
    {
        if (skippi0 && itype==PI0) continue;

        string prefix = "ntrg";
        if (isadj && itype!=PI0)
            prefix.append("_adj");

        for (int iptt = 0; iptt < nptt; iptt++)
        {
            for (int icent = 0; icent < ncent; icent++)
            {
                for (int idpsi = 0; idpsi < ndpsi; idpsi++)
                {
                    TH1F* h1tmp = (TH1F*) gROOT->FindObject(Form("%s_%s_pttb%i_centb%i_dpsib%i",prefix.c_str(),spec[itype],iptt,icent,idpsi));

                    if (h1tmp == NULL)
                    {
                        cout<<__LINE__;
                        printf(": Could not find %s_%s_pttb%i_centb%i_dpsib%i",prefix.c_str(),spec[itype],iptt,icent,idpsi);
                        return;
                    }

                    float tmpbc = h1tmp->GetBinContent(1);
                    ntrg_pttcentdpsi[itype][iptt][icent][idpsi] = tmpbc;
                    ntrg_pttcentdpsi[itype][iptt][ncent][idpsi] += tmpbc;
                    ntrg_pttcentdpsi[itype][iptt][ncent][ndpsi] += tmpbc;
                    ntrg_pttcentdpsi[itype][iptt][icent][ndpsi] += tmpbc;

                }
            }
        }
    }


    for (int itype = 0; itype < ntype; itype++)
    {
        if (skippi0 && itype==PI0) continue;
        for (int iptt = 0; iptt < nptt; iptt++)
        {
            for (int icent = 0; icent < ncent+1; icent++)
            {
                for (int idpsi = 0; idpsi <= ndpsi; idpsi++)
                {
                    cout<<"ntrg_pttcentdpsi["<<itype<<"]["<<iptt<<"]["<<icent<<"]["<<idpsi<<"] = "<<ntrg_pttcentdpsi[itype][iptt][icent][idpsi]<<endl;
                }
            }
        }
    }

    //============================================================
    // Get histo counter for n mixed events in pt and cent bins
    // Fill our counter array
    //============================================================

    //============================================================
}

void GetRemBg()
{
    cout<<"Getting the remaining background"<<endl;

    for (int itype = 0; itype < ntype+1; itype++)
    {

        if (skippi0 && itype==PI0) continue;

        for (int iptt = 0; iptt < nptt; iptt++)
        {

            for (int icent = 0; icent < ncent+1; icent++)
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
                    for (int ipta = 0; ipta < npta; ipta++)
                    {

                        h1dphi_bg[itype][iptt][icent][idpsi][ipta] = (TH1F*) h1dphi[itype][MIX][iptt][icent][ndpsi][ipta]->Clone(Form("dphi_bg_%s_ptt%i_cent%i_dpsi%i_pta%i",spec[itype],iptt,icent,idpsi,ipta));

                        h1dphi_bg[itype][iptt][icent][idpsi][ipta]->SetTitle(Form("BG p_{tt}:[%.2f,%.2f) Cent:[%i,%i) d#Psi_{2,BBC}:%s  p_{ta}:[%.2f,%.2f);#Delta#phi (#pi rad);#frac{1}{N_{trg}}#frac{dN_{pair}}{d(#Delta#phi)}",limptt[iptt],limptt[iptt+1],locent,hicent,(dpsistr[idpsi]),limpta[ipta],limpta[ipta+1]));

                        h1dphi_bg[itype][iptt][icent][idpsi][ipta]->Add(h1dphi_nopsimix[itype][iptt][icent][ndpsi][ipta],-1.);

                    }
                    for (int ixi = 0; ixi < nxi; ixi++)
                    {

                        h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi] = (TH1F*) h1dphi_xi[itype][MIX][iptt][icent][ndpsi][ixi]->Clone(Form("dphi_xi_bg_%s_ptt%i_cent%i_dpsi%i_xi%i",spec[itype],iptt,icent,idpsi,ixi));

                        h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi]->SetTitle(Form("BG p_{tt}:[%.2f,%.2f) Cent:[%i,%i) d#Psi_{2,BBC}:%s  #xi:[%.2f,%.2f);#Delta#phi (#pi rad);#frac{1}{N_{trg}}#frac{dN_{pair}}{d(#Delta#phi)}",limptt[iptt],limptt[iptt+1],locent,hicent,(dpsistr[idpsi]),limxi[ixi],limxi[ixi+1]));

                        if (useksic)
                            h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi]->Add(h1dphi_xi_nopsimix_ksic[itype][iptt][icent][ndpsi][ixi],-1.);
                        else
                            h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi]->Add(h1dphi_xi_nopsimix[itype][iptt][icent][ndpsi][ixi],-1.);

                    }
                }
            }
        }
    }
    cout<<"Successfully got remaining background"<<endl;
}

void SubRemBg()
{
    cout<<"Subtracting the remaining background"<<endl;

    for (int itype = 0; itype < ntype+1; itype++)
    {
        if (skippi0 && itype==PI0) continue;
        for (int iptt = 0; iptt < nptt; iptt++)
        {
            for (int icent = 0; icent < ncent+1; icent++)
            {
                for (int idpsi = 0; idpsi <= ndpsi; idpsi++)
                {
                    for (int ipta = 0; ipta < npta; ipta++)
                    {
                        h1dphi_bg[itype][iptt][icent][idpsi][ipta]->Divide(h1dphi_nopsimix[itype][iptt][icent][ndpsi][ipta]);

                        if (idpsi != ndpsi)
                            h1dphi_bg[itype][iptt][icent][idpsi][ipta]->Scale(0.5);

                        h1dphi_bg[itype][iptt][icent][idpsi][ipta]->Scale(1./pow(reso[icent],2) - 1.);

                        TF1* fit = new TF1(Form("fit_%i%i%i%i%i",itype,iptt,ipta,icent,idpsi),"[0]*cos(x*2*TMath::Pi())",0.,1.0);
                        cout<<"\nFitting "<<h1dphi_bg[itype][iptt][icent][idpsi][ipta]->GetName()<<" with function of form "<<fit->GetExpFormula()<<endl;
                        h1dphi_bg[itype][iptt][icent][idpsi][ipta]->Fit(fit,"M");

                        h1dphi[itype][MSN][iptt][icent][idpsi][ipta]->Add(fit,-1.);

                    }

                    for (int ixi = 0; ixi < nxi; ixi++)
                    {
                        h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi]->Divide(h1dphi_xi_nopsimix[itype][iptt][icent][ndpsi][ixi]);

                        if (idpsi != ndpsi)
                            h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi]->Scale(0.5);

                        h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi]->Scale(1./pow(reso[icent],2) - 1.);

                        TF1* fit = new TF1(Form("fit_xi_%i%i%i%i%i",itype,iptt,ixi,icent,idpsi),"[0]*cos(x*2*TMath::Pi())",0.,1.0);
                        cout<<"\nFitting "<<h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi]->GetName()<<" with function of form "<<fit->GetExpFormula()<<endl;
                        h1dphi_xi_bg[itype][iptt][icent][idpsi][ixi]->Fit(fit,"M");

                        h1dphi_xi[itype][MSN][iptt][icent][idpsi][ixi]->Add(fit,-1.);

                    }
                }
            }
        }
    }
}
