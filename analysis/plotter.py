#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

import ROOT

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPadRightMargin(0.10);


############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--type',action="store",type="string",dest="type",default="gg")
parser.add_option('--ptbin',action="store",type="int",dest="ptbin",default=500)

(options, args) = parser.parse_args()


############################################################

def makeCanvas(hists, names, canname):
    
    bin = options.ptbin;
    directory = "figs_bin"+str(bin);
    colors = [2,4,1,6,7];
    
    max = -999.;
    for hist in hists:
        if max < hist.GetMaximum(): max = hist.GetMaximum();

    leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    for i in range(len(names)):
        hists[i].SetLineColor(colors[i]);
        leg.AddEntry(hists[i], names[i], "l");


    can = ROOT.TCanvas("can"+canname,"can"+canname,1000,800);
    hists[0].SetMaximum( 1.2*max );
    hists[0].SetMinimum( 0 );    
    hists[0].Draw();    
    for i in range(1,len(hists)):
        hists[i].Draw("sames");
    leg.Draw();
    can.SaveAs(directory+"/"+canname+".eps");
    can.SaveAs(directory+"/"+canname+".png");

    for hist in hists:
        hist.Scale(1./hist.Integral());
        
def makeSingleHistCanvas( h_2D ):
    
    directory = "figs_bin"+str(bin);

    cant = ROOT.TCanvas("cant","cant",800,800);
    h_2D.Draw("colz");
    cant.SaveAs(directory+"/"+h_2D.GetName()+".png");
    cant.SaveAs(directory+"/"+h_2D.GetName()+".eps");

    
############################################################

if __name__ == '__main__':

    type = options.type;
    bin = options.ptbin;
    binp1 = options.ptbin + 100;
        
    directory = "figs_bin"+str(bin);
    if not os.path.exists(directory):
        os.makedirs(directory)

    types = ["gg","qq"];

    files = [];
    files.append("data/boost2013-pythia81-tune4c-lhc8-gg-pt0"+str(bin)+"-0"+str(binp1)+".root");
    files.append("data/boost2013-pythia81-tune4c-lhc8-qq-pt0"+str(bin)+"-0"+str(binp1)+".root");    

    h_pt1 = [];
    h_pt2 = [];
    h_eta1 = [];
    h_eta2 = [];
    h_mass1 = [];
    h_mass2 = [];
    h_tau1_b1 = [];
    h_tau1_b2 = [];
    h_c1_b0 = [];
    h_c1_b1 = [];
    h_c1_b2 = [];
    
    #observables to correlate
    obs2Corr = ["jmass","jtau1_b1","jtau1_b2","jc1_b0","jc1_b1","jc1_b2"];
    obs2Corr_name = ["mass","#tau_{1}^{#beta=1}","#tau_{1}^{#beta=2}","C_{1}^{#beta=0}","C_{1}^{#beta=1}","C_{1}^{#beta=2}"]
    obs2Corr_nb = [40,40,40,40,40,40];
    obs2Corr_lo = [0,0,0,0.2,0,0];
    obs2Corr_hi = [200,0.5,0.5,0.5,0.3,0.3];
            
    h_2D = [];
    for a in range(len(types)):
        h_2D_type = []
        for i in range(len(obs2Corr)):
            for j in range(i+1,len(obs2Corr)):
                hname = "h2d_"+obs2Corr[i]+"_"+obs2Corr[j]+"_"+types[a];
                htitle = ";"+obs2Corr_name[i]+";"+obs2Corr_name[j];
                print "hname = ",hname
                h_2D_type.append( ROOT.TH2F(hname,htitle,obs2Corr_nb[i],obs2Corr_lo[i],obs2Corr_hi[i],obs2Corr_nb[j],obs2Corr_lo[j],obs2Corr_hi[j]) );
        h_2D.append( h_2D_type );
        
        
    for i in range(len(types)):
        h_pt1.append( ROOT.TH1F("h_pt1_"+types[i],";jet 1 pT (GeV); N",50,0,2000) );
        h_pt2.append( ROOT.TH1F("h_pt2_"+types[i],";jet 2 pT (GeV); N",50,0,2000) );
        h_eta1.append( ROOT.TH1F("h_eta1_"+types[i],";jet 1 eta; N",20,-5,5) );
        h_eta2.append( ROOT.TH1F("h_eta2_"+types[i],";jet 2 eta; N",20,-5,5) );
        h_mass1.append( ROOT.TH1F("h_mass1_"+types[i],";jet 1 mass; N",40,0,200) );
        h_mass2.append( ROOT.TH1F("h_mass2_"+types[i],";jet 2 mass; N",40,0,200) );
    
        h_tau1_b1.append( ROOT.TH1F("h_tau1_b1"+types[i],";tau1 (beta=1); N",40,0,0.5) );
        h_tau1_b2.append( ROOT.TH1F("h_tau1_b2"+types[i],";tau1 (beta=2); N",40,0,0.5) );        
        h_c1_b0.append( ROOT.TH1F("h_c1_b0"+types[i],";c1 (beta=0); N",40,0.2,0.5) );
        h_c1_b1.append( ROOT.TH1F("h_c1_b1"+types[i],";c1 (beta=1); N",40,0,0.3) );        
        h_c1_b2.append( ROOT.TH1F("h_c1_b2"+types[i],";c1 (beta=2); N",40,0,0.3) );        
            
    for a in range(len(types)):
    
        f = ROOT.TFile(files[a]);
        t = f.Get("t");
    
        for i in range(t.GetEntries()):

            t.GetEntry(i);
            h_pt1[a].Fill( t.jpt[0] );
            h_pt2[a].Fill( t.jpt[1] );
            h_eta1[a].Fill( t.jeta[0] );
            h_eta2[a].Fill( t.jeta[1] );
            h_mass1[a].Fill( t.jmass[0] );
            h_mass2[a].Fill( t.jmass[1] );
            
            h_tau1_b1[a].Fill( t.jtau1_b1[0] );
            h_tau1_b2[a].Fill( t.jtau1_b2[0] );            
            
            h_c1_b0[a].Fill( t.jc1_b0[0] );                        
            h_c1_b1[a].Fill( t.jc1_b1[0] );                        
            h_c1_b2[a].Fill( t.jc1_b2[0] );   
            
            #2D plots
            obsctr = 0;
            for j in range(len(obs2Corr)):
                for k in range(j+1,len(obs2Corr)):
                    h_2D[a][obsctr].Fill( getattr(t,obs2Corr[j])[0], getattr(t,obs2Corr[k])[0] );
                    obsctr += 1;
            
        del f;
        del t;
                
    # plot
    makeCanvas(h_pt1, types, "jpt1");
    makeCanvas(h_pt2, types, "jpt2");
    makeCanvas(h_eta1, types, "jeta1");
    makeCanvas(h_eta2, types, "jeta2");
    makeCanvas(h_mass1, types, "jmass1");
    makeCanvas(h_mass2, types, "jmass2");
                
    makeCanvas(h_tau1_b1, types, "h_tau1_b1");
    makeCanvas(h_tau1_b2, types, "h_tau1_b2");                    
    makeCanvas(h_c1_b0, types, "h_c1_b0");                    
    makeCanvas(h_c1_b1, types, "h_c1_b1");                    
    makeCanvas(h_c1_b2, types, "h_c1_b2");                            
     
    for a in range(len(types)):
        #2D plots
        obsctr = 0;
        for j in range(len(obs2Corr)):
            for k in range(j+1,len(obs2Corr)):
                makeSingleHistCanvas(h_2D[a][obsctr]);        
                obsctr += 1;
            

        