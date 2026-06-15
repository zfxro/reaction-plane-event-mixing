void incphotonyieldrebin()
{

    //vector for pt limits
    std::vector<float> v_lim;
    float lim_tmp =.75;
    v_lim.push_back(lim_tmp);
    float max = 35.;

    while (lim_tmp < max)
    {
        float bininc = .25;

        if (lim_tmp >= 2) bininc = .5;
        if (lim_tmp >= 5) bininc = 1.;
        if (lim_tmp >= 10) bininc = 3.;
        if (lim_tmp >= 15) bininc = 5.;

        lim_tmp+=bininc;
        v_lim.push_back(lim_tmp);
    }

    for (std::vector<float>::iterator it = v_lim.begin(); it != v_lim.end(); it++)
    {
        printf("i = %li; lim = %f;\n",it - v_lim.begin(), *it);
    }

    std::vector<float> v_x;
    for (std::vector<float>::iterator it = v_lim.begin(); it != v_lim.end() - 1; it++)
    {
        float tmpx = *(it+1) + *it;
        tmpx = tmpx/2.;
        printf("i = %li; x = %f; limits = [%f,%f);\n",it - v_lim.begin(),tmpx,*it,*(it+1));
        v_x.push_back(tmpx);
    }

    const int nc = 4;
    int centlim[nc+1] = {0,20,40,60,100};
    std::vector<float> v_bc[nc];
    TGraph* tgr[nc];
    TCanvas* can[nc];

    TFile* file = TFile::Open("IncGamSpectra.root");

    for (int ic = 0; ic < nc; ic++)
    {
        TH1F* he = (TH1F*) gROOT->FindObject(Form("e%i",ic));
        if (!he)
        {
            printf("Could not find histogram e%i. Exiting.\n",ic);
            return;
        }
        int nb = he->GetNbinsX();

        for (std::vector<float>::iterator it = v_lim.begin(); it != v_lim.end() - 1; it++)
        {
            float bc = 0.;
            for (int ib = 0; ib < nb; ib++)
            {
                float center = he->GetBinCenter(ib+1);
                if (center >= *it && center < *(it+1))
                    bc+=he->GetBinContent(ib+1);
            }
            float binw = *(it+1) - *it;
            bc = bc/binw;
            v_bc[ic].push_back(bc);
        }

        tgr[ic] = new TGraph(v_x.size(),&(v_x[0]),&(v_bc[ic][0]));
        tgr[ic]->Print();

        can[ic] = new TCanvas(Form("can%i",ic),Form("can%i",ic),1000,600);
        tgr[ic]->SetTitle(Form("centrality %i",ic));
        tgr[ic]->SetMarkerStyle(20);
        tgr[ic]->SetMarkerSize(1);
        tgr[ic]->Draw("ap");

        delete he;

    }

    //print to txt
    for (int ic = 0; ic < nc; ic++)
    {
        ofstream* ofs = new ofstream(Form("incgam_%i_%i.txt",centlim[ic],centlim[ic+1]));
        for (std::vector<float>::iterator it = v_x.begin(); it != v_x.end(); it++)
        {
            *ofs<<*it<<" "<<v_bc[ic][it - v_x.begin()]<<endl;
        }
    }

}
