#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "iostream"
#include "fstream"
#include "stdlib.h"

using namespace std;


int macro()
{
    TString samples[6] = {"pp","pomp","pompom","phop","phopom","phopho"};
    map<TString,TFile*> files;
    map<TString,TTree*> trees;
    map<TString,TH1F*>  heta;
    map<TString,TH1F*>  hetab1;
    map<TString,TH1F*>  hptb1;
    map<TString,TH1F*>  hpt2bb;
    map<TString,TH1F*>  hetab2;
    map<TString,TH1F*>  hnbackward;

    for(int i=0;i<6;i++) heta[samples[i]]=new TH1F("eta_"+samples[i],"",100,-10,10);
    for(int i=0;i<6;i++) hnbackward[samples[i]]=new TH1F("nbackward_"+samples[i],"",100,-0.5,100.5);
    for(int i=0;i<6;i++) hetab1[samples[i]]=new TH1F("b1eta_"+samples[i],"",100,-10,10);
    for(int i=0;i<6;i++) hptb1[samples[i]]=new TH1F("b1pt_"+samples[i],"",100,0,100);
    for(int i=0;i<6;i++) hpt2bb[samples[i]]=new TH1F("bb2pt_"+samples[i],"",100,0,100);
    for(int i=0;i<6;i++) hetab2[samples[i]]=new TH1F("b2eta_"+samples[i],"",100,-10,10);



    Float_t cs[6] = {
        359746153.83010089, 
        26272939.090299036,
        1506425.5262538597,
        11501.607068787487,
        1745.5448611727602,
        0.10936557517629096};

    Int_t ngen;
    Float_t px[5000],py[5000],pz[5000],e[5000];
    Int_t id[5000], ist[5000];

    for(int i=0;i<6;i++) {
        cout << samples[i] << endl;
        TFile *file = new TFile("/data/ntuples/FPMC/data"+samples[i]+"_hq.root");
        TTree *t    =(TTree*)file->Get("h777");
        //trees[samples[i]]=(TTree*)files[samples[i]]->Get("h777");
        t->SetBranchAddress("ngen",&ngen);
        t->SetBranchAddress("id",id);
        t->SetBranchAddress("ist",ist);
        t->SetBranchAddress("px",px);
        t->SetBranchAddress("py",py);
        t->SetBranchAddress("pz",pz);
        t->SetBranchAddress("e",e);

        TLorentzVector *b1 = new TLorentzVector(); //bottom
        TLorentzVector *b2 = new TLorentzVector(); //bbar
        TLorentzVector *tmp = new TLorentzVector(); //bbar
        Int_t lhcb_nobackward=0;

        for(int ie=0;ie<t->GetEntries();ie++) {
            t->GetEntry(ie);
            //if(ie%100000==0) cout << ie << " ngen = " << ngen << endl;
            b1->SetPxPyPzE(0,0,0,0);
            b2->SetPxPyPzE(0,0,0,0);
            int backward=0, forward=0;
            for(int ip=0;ip<ngen;ip++) {
                //if(ist[ip]!=158) continue;
                if(ist[ip]==123 || ist[ip]==124) {
                    if(id[ip]==5  && b1->E()<0.1) b1->SetPxPyPzE(px[ip],py[ip],pz[ip],e[ip]);
                    if(id[ip]==-5 && b2->E()<0.1) b2->SetPxPyPzE(px[ip],py[ip],pz[ip],e[ip]);
                }
                if(ist[ip]!=1 || sqrt(px[ip]*px[ip]+py[ip]*py[ip])<.01) continue;
                tmp->SetPxPyPzE(px[ip],py[ip],pz[ip],e[ip]);
                if(tmp->Eta()>-4.5&&tmp->Eta()<-1.5) backward++;
                if(tmp->Eta()<4.5&&tmp->Eta()>1.5) forward++;
                if(tmp->Eta()>-7.&&tmp->Eta()<-5) backward++;
                if(tmp->Eta()<7.&&tmp->Eta()>5) forward++;

            }
            if(b1->E()<0.1 || b2->E()<0.1) {
                //cout << "not found b's. ie = " << ie << endl; 
                //    for(int ip=0;ip<ngen;ip++) {
                //        cout << ip << " " << id[ip] << " " << ist[ip] << endl;
                //    }
                continue;
            }
            //if(b1->Pt()<5) continue;
            hetab1[samples[i]]->Fill(b1->Eta(),cs[i]/t->GetEntries());
            hptb1[samples[i]]->Fill(b1->Pt(),cs[i]/t->GetEntries());
            hpt2bb[samples[i]]->Fill(pow((*b1+*b2).Pt(),2.),cs[i]/t->GetEntries());
            hetab2[samples[i]]->Fill(b2->Eta(),cs[i]/t->GetEntries());
            if(b1->Eta()>1.5&&b1->Eta()<4.5&&b2->Eta()>1.5&&b2->Eta()<4.5) {
              if(backward==0) lhcb_nobackward++;
              hnbackward[samples[i]]->Fill(backward);
            }
            if(samples[i]=="pomp" || samples[i]=="phopom" || samples[i]=="phop") {
                hetab1[samples[i]]->Fill(-b1->Eta(),cs[i]/t->GetEntries());
                hptb1[samples[i]]->Fill(b1->Pt(),cs[i]/t->GetEntries());
                hpt2bb[samples[i]]->Fill(pow((*b1+*b2).Pt(),2.),cs[i]/t->GetEntries());
                hetab2[samples[i]]->Fill(-b2->Eta(),cs[i]/t->GetEntries());
                if((b1->Eta()<-1.5&&b1->Eta()>-4.5)&&(b2->Eta()<-1.5&&b2->Eta()>-4.5)) {
                    if(forward==0) lhcb_nobackward++;
                    hnbackward[samples[i]]->Fill(forward);
                }
            }
        }
        cout << samples[i] << " cs[pb] = " << cs[i]*hetab1[samples[i]]->GetEntries()/t->GetEntries() << endl;
        cout << samples[i] << " cs[pb] in LHCb with backward gap = " << cs[i]*lhcb_nobackward/t->GetEntries() << endl;
    }

    TCanvas *c2 = new TCanvas("c2","lego");
    c2->cd();
    for(int i=0;i<6;i++) {
        if(i==0)hptb1[samples[i]]->Draw();
        else hptb1[samples[i]]->Draw("same");
    }

    TCanvas *c3 = new TCanvas("c3","");
    c3->cd();
    cout << "finished" << endl;
    TFile *fsave = new TFile("plots.root","recreate");
    for(int i=0;i<6;i++) hnbackward[samples[i]]->Write();
    for(int i=0;i<6;i++) hptb1[samples[i]]->Write();
    for(int i=0;i<6;i++) hpt2bb[samples[i]]->Write();
    fsave->Close();

    return 1;
}

/*
   /data/ntuples/FPMC/dataphopho_hq.log: Cross section[pb]=  0.10936557517629096       NEVENTS=      500000
   /data/ntuples/FPMC/dataphop_hq.log: Cross section[pb]=   11501.607068787487       NEVENTS=      500000
   /data/ntuples/FPMC/dataphopom_hq.log: Cross section[pb]=   1745.5448611727602       NEVENTS=      500000
   /data/ntuples/FPMC/datapomp_hq.log: Cross section[pb]=   26272939.090299036       NEVENTS=      500000
   /data/ntuples/FPMC/datapompom_hq.log: Cross section[pb]=   1506425.5262538597       NEVENTS=      500000
   /data/ntuples/FPMC/datapp_hq.log: Cross section[pb]=   359746153.83010089       NEVENTS=      300000
   */
