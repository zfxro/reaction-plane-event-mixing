#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TTree.h>

using namespace std;

void save_nevent(const char* inlist = "filelist.txt", const char* outname = "EventCount.root")
{
    cout<<"Running save_nevent."<<endl;
    
    TFile* output = new TFile(outname,"RECREATE");
    if (!output)
    {
        cout<<"Could not open "<<outname<<endl;
        return;
    }
    if (output->IsZombie()) 
    {
        cout<<"Output "<<outname<<" is zombie."<<endl;
        return;
    }
    cout<<"Output = "<<outname<<endl;

    ifstream ifs (inlist);
    string inname;
    int nfile = 0;
    while(getline(ifs,inname))
    {
        if (nfile%100==0)
        {
            cout<<"Reading file "<<nfile<<" "<<inname<<endl;
        }
        
        int dotpos = inname.find(".");
        string runno = inname.substr(dotpos-6,6);

        TFile* input = new TFile(inname.c_str(),"READ");
        if (!input)
        {
            cout<<"Could not open "<<inname<<endl;
            return;
        }
        if (input->IsZombie()) 
        {
            cout<<"Input "<<inname<<" is zombie."<<endl;
            return;
        }

        TH2F* h2bz0cent_tmp = (TH2F*)gROOT->FindObject("bz0cent");
        h2bz0cent = (TH2F*) h2bz0cent_tmp->Clone(Form("bz0cent_%s",runno));
        output->cd();
        h2bz0cent->Write();

        delete h2bz0cent_tmp;
        delete h2bz0cent;

        input->Close();
        delete input;

        nfile++;
    }

    output->Write();
    output->Close();
    delete output;
}
