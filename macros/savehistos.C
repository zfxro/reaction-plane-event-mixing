#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TTree.h>

using namespace std;

void savehistos(const char* filelist = "inlist.txt", const char* outdir = "/phenix/scratch/zrowan/")
{
    ifstream ifs (filelist);
    string filename;
    while(getline(ifs,filename))
    {
        int dotpos = filename.find(".");
        string runno = filename.substr(dotpos-6,6);

        ostringstream oss;
        oss<<outdir<<"histos"<<"_"<<runno<<".root";
        cout<<oss.str()<<endl;

        TFile* file = new TFile(filename.c_str(),"READ");
        if (file->IsZombie()) 
        {
            cout<<"File is zombie."<<endl;
            continue;
        }
        TIter next(file->GetListOfKeys());
        TKey *key;

        TFile* outfile = new TFile(oss.str().c_str(),"RECREATE");

        while ((key = (TKey*)next())) 
        {
            TClass *cl = gROOT->GetClass(key->GetClassName());
            if (cl->InheritsFrom("TH1"))
            {
                TH1 *h1 = (TH1*)key->ReadObj();
                h1->SetDirectory(0);
                outfile->cd();
                h1->Write();
            }
            if (cl->InheritsFrom("TH2"))
            {
                TH2 *h2 = (TH2*)key->ReadObj();
                h2->SetDirectory(0);
                outfile->cd();
                h2->Write();
            }
        }
    }
}
