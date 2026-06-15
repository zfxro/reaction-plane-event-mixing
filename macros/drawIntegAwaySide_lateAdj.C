//draw results from lateAdj.C

// //----- from lateAdj.C -----
// //yields
// //v2j

const float PI = TMath::Pi();
const int nadj = 2;
const int ntype = 3;
const int nsumcent = 2;

const int limsumcent[nsumcent+1] = {0,20,60};

//-------------------- yields --------------------
TGraphErrors* tgrixiyield[nadj][ntype][GAT::nptt][nsumcent+1][GAT::ndpsi+1] = {};
TGraphErrors* tgrixiyield_syserr[nadj][ntype][GAT::nptt][nsumcent+1][GAT::ndpsi+1] = {};

TGraphErrors* tgripttyield[nadj][ntype][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};
TGraphErrors* tgripttyield_syserr[nadj][ntype][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};

