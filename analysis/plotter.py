#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time

import ROOT

from TMVAhelper import *

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPadRightMargin(0.20);
ROOT.gStyle.SetPalette(1);
ROOT.gStyle.SetPaintTextFormat("1.1f");

#ROOT.gROOT.ProcessLine(".L PerformanceFinal.C");

#from ROOT import myutils
ROOT.gROOT.ProcessLine(".L myutils.C++");
ROOT.gSystem.Load("myUtils_C.so");
#ROOT.gSystem.Load("PerformanceFinal_C.so");
from ROOT import myutils


############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--type',action="store",type="string",dest="type",default="gg")
parser.add_option('--ptbin',action="store",type="int",dest="ptbin",default=500)
parser.add_option('--rPar',action="store",type="int",dest="rPar",default=500)
parser.add_option('--sigEff',action="store",type="float",dest="sigEff",default=500)

parser.add_option('--do2D',action="store_true",dest="do2D",default=False,help='making correlations, do MVA')
parser.add_option('--train',action="store_true",dest="train",default=False,help='do training')
parser.add_option('--train3',action="store_true",dest="train3",default=False,help='do training')
parser.add_option('--trainAllOnly',action="store_true",dest="trainAllOnly",default=False,help='do training')
parser.add_option('--tmvaname',action="store",type="string",dest="tmvaname",default="qg")
parser.add_option('--tmvamethod',action="store",type="string",dest="tmvamethod",default="BDT")

(options, args) = parser.parse_args()

############################################################

# global vars
type = options.type;
bin = options.ptbin;
binp1 = options.ptbin + 100;

bin_cut = bin;
binp1_cut = binp1;
# if bin == 1000: binp1_cut = binp1 + 100;
# else: binp1_cut = binp1 + 50;

cone = "ak%02i" % (options.rPar);

directory = "plots/figs011915_"+options.tmvaname+"_bin"+str(bin)+"_"+cone;
if not os.path.exists(directory): os.makedirs(directory)
if not os.path.exists(directory+"/rocs"): os.makedirs(directory+"/rocs");
if not os.path.exists(directory+"/oneD"): os.makedirs(directory+"/oneD");
if not os.path.exists(directory+"/twoD"): os.makedirs(directory+"/twoD");

globalBannerText = "%s, pT = %i-%i GeV, AK%i" % (options.tmvaname,bin_cut,binp1_cut,options.rPar);
gbanner = ROOT.TLatex(0.17,0.96,globalBannerText);
gbanner.SetNDC();
gbanner.SetTextSize(0.030);
gbanner.SetTextFont(42);
############################################################

def makeCanvas(hists, names, canname, normalize=False):
    
    #bin = options.ptbin;
    #directory = "figs_bin"+str(bin);
    colors = [2,4,1,6,7];
    
    leg = ROOT.TLegend(0.65,0.77,0.85,0.92);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.05);    
    leg.SetTextFont(42);
    for i in range(len(names)):
        hists[i].SetLineColor(colors[i]);
        hists[i].SetLineWidth(2);
        if normalize: hists[i].Scale(1./hists[i].Integral());
        leg.AddEntry(hists[i], names[i], "l");

    banner = ROOT.TLatex(0.77,0.96,"BOOST13WG");
    banner.SetNDC();
    banner.SetTextSize(0.030);  
    banner.SetTextFont(52);      

    max = -999.;
    for hist in hists:
        if max < hist.GetMaximum(): max = hist.GetMaximum();

    can = ROOT.TCanvas("can"+canname,"can"+canname,800,800);
    can.SetRightMargin(0.05);
    hists[0].SetMaximum( 1.4*max );
    if "mass" in canname: hists[0].SetMaximum( 0.7 );
    hists[0].SetMinimum( 0 );    
    hists[0].Draw();    
    for i in range(1,len(hists)):
        hists[i].Draw("sames");
    leg.Draw();
    banner.Draw();
    gbanner.Draw();
    can.SaveAs(directory+"/oneD/"+canname+".eps");
    can.SaveAs(directory+"/oneD/"+canname+".pdf");
    can.SaveAs(directory+"/oneD/"+canname+".png");

    for hist in hists:
        hist.Scale(1./hist.Integral());

def makeCanvas2D(hist0, hist1, types):
    
    hist1.SetLineColor(2);
    
    leg = ROOT.TLegend(0.30,0.8,0.5,0.9);
    leg.SetBorderSize(1);
    leg.SetFillStyle(1001);
    leg.SetTextFont(42);
    leg.SetTextSize(0.035);
    leg.AddEntry(hist0, types[0], "");

    leg1 = ROOT.TLegend(0.30,0.8,0.5,0.9);
    leg1.SetBorderSize(1);
    leg1.SetFillStyle(1001);
    leg1.SetTextSize(0.035);
    leg1.AddEntry(hist1, types[1], "");

    hist0.Scale(1./hist0.Integral());
    hist1.Scale(1./hist1.Integral());
    max = 0;
    if hist0.GetMaximum() < hist1.GetMaximum(): max = hist1.GetMaximum();
    if hist0.GetMaximum() > hist1.GetMaximum(): max = hist0.GetMaximum();

    can = ROOT.TCanvas("canOS"+hist0.GetName(),"canOS"+hist0.GetName(),1600,800);
    #hist0.SetMaximum( max*1.1 );
    #hist1.SetMaximum( max*1.1 );
    can.Divide(2,1);
    can.cd(1);
    hist0.Draw("COLZ");
    leg.Draw();
    can.cd(2);
    hist1.Draw("COLZ");
    leg1.Draw();
    can.SaveAs(directory+"/twoD/"+hist0.GetName()+"_onSame.eps");
    can.SaveAs(directory+"/twoD/"+hist0.GetName()+"_onSame.pdf");    
    can.SaveAs(directory+"/twoD/"+hist0.GetName()+"_onSame.png");

def makeSingleHistCanvas( h_2D ):
    
    #directory = "figs_bin"+str(bin);

    cant = ROOT.TCanvas("cant","cant",800,800);
    h_2D.Draw("colz");
    cant.SaveAs(directory+"/twoD/"+h_2D.GetName()+".png");
    cant.SaveAs(directory+"/twoD/"+h_2D.GetName()+".eps");
    cant.SaveAs(directory+"/twoD/"+h_2D.GetName()+".pdf");    

def makeROCFromHisto(hsig,hbkg,LtoR):

    nbins = hsig.GetNbinsX();
    binsize = hsig.GetBinWidth(1);
    lowedge = hsig.GetBinLowEdge(1);

    #print "lowedge: ",lowedge

    hsigIntegral = hsig.Integral();
    hbkgIntegral = hbkg.Integral();

    xval = array('d', [])
    yval = array('d', [])
    ctr = 0;
    effBkgPrev = -9999;
    for i in range(1,nbins+1):

        effBkg = 0;
        effSig = 0;

        if LtoR: effBkg = hbkg.Integral( i, nbins )/hbkgIntegral;
        else: effBkg = hbkg.Integral( 1, i )/hbkgIntegral;

        if LtoR: effSig = hsig.Integral( i, nbins )/hsigIntegral;
        else: effSig = hsig.Integral( 1, i )/hsigIntegral;

        #if not effBkg == 0.: print "cut: ",(lowedge+(i-1)*binsize),"effBkg: ", 1/effBkg, ", effSig: ", effSig;

        xval.append( effSig );
        yval.append( effBkg );

        #effBkgPrev = effBkg;
        ctr = ctr + 1;

    #print nbins, "and ", ctr
    tg = ROOT.TGraph( nbins, xval, yval );
    tg.SetName( "tg"+hsig.GetName() );
    return tg;

def SetLineColorStyleWidth(tg, color, style, width):
    tg.SetLineColor( color );
    tg.SetLineStyle( style );
    tg.SetLineWidth( width );
    
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if __name__ == '__main__':

    ROOT.TMVA.Tools.Instance();
    tmvamethod = options.tmvamethod;

    types = [];
    types_names = [];
    if options.tmvaname == "qg": types = ["gg","qq"];
    if options.tmvaname == "Wg": types = ["WW","gg"];
    if options.tmvaname == "Wq": types = ["WW","qq"];
    if options.tmvaname == "qg": types_names = ["gluon jet","quark jet"];
    if options.tmvaname == "Wg": types_names = ["W jet","gluon jet"];
    if options.tmvaname == "Wq": types_names = ["W jet","quark jet"];
    
    files = [];
    if options.tmvaname == "Wq" or options.tmvaname == "Wg": files.append("dataProcessedFinal/boost2013-jhu452-lhc8-wlwl-pt%04i-%04i_%s.root" % (bin,binp1,cone));
    if options.tmvaname == "qg" or options.tmvaname == "Wg": files.append("dataProcessedFinal/boost2013-pythia81-tune4c-lhc8-gg-pt%04i-%04i_%s.root" % (bin,binp1,cone));
    if options.tmvaname == "Wq" or options.tmvaname == "qg": files.append("dataProcessed/boost2013-pythia81-tune4c-lhc8-qq-pt%04i-%04i_%s-100k.root" % (bin,binp1,cone));

    trainingFiles = [];
    if options.tmvaname == "Wq" or options.tmvaname == "Wg": trainingFiles.append("dataProcessedFinal/boost2013-jhu452-lhc8-wlwl-pt%04i-%04i_xcheck_%s.root" % (bin,binp1,cone));
    if options.tmvaname == "qg" or options.tmvaname == "Wg": trainingFiles.append("dataProcessedFinal/boost2013-pythia81-tune4c-lhc8-gg-pt%04i-%04i_xcheck_%s.root" % (bin,binp1,cone));
    if options.tmvaname == "Wq" or options.tmvaname == "qg": trainingFiles.append("dataProcessed/boost2013-pythia81-tune4c-lhc8-qq-pt%04i-%04i_xcheck_%s-100k.root" % (bin,binp1,cone));
    # if options.tmvaname == "Wq" or options.tmvaname == "Wg": trainingFiles.append("data/boost2013-jhu452-lhc8-wlwl-pt%04i-%04i_%s.root" % (bin,binp1,cone));
    # if options.tmvaname == "qg" or options.tmvaname == "Wg": trainingFiles.append("data/boost2013-pythia81-tune4c-lhc8-gg-pt%04i-%04i_%s.root" % (bin,binp1,cone));
    # if options.tmvaname == "Wq" or options.tmvaname == "qg": trainingFiles.append("data/boost2013-pythia81-tune4c-lhc8-qq-pt%04i-%04i_%s.root" % (bin,binp1,cone));

    # print options.tmvaname
    # print "Number of testing  files = ", len(files);
    # print "Number of training files = ", len(trainingFiles);

    bdt_nb = 1000;
    obs_nb = 50;
    mass_obs_nb = 50;

    print files

    #######-----------------------------------------------
    ####### setup
    #######-----------------------------------------------

    ## =====================================================================================
    ## 1D setup
    ## =====================================================================================
    
    h_pt1 = [];
    h_pt2 = [];
    h_eta1 = [];
    h_eta2 = [];
    h_mass1 = [];
    h_mass2 = [];
    h_tau1_b1 = [];
    h_tau1_b2 = [];
    h_tau21_b1 = [];
    h_tau21_b2 = [];
    h_c1_b0 = [];
    h_c1_b1 = [];
    h_c1_b2 = [];
    h_c2_b1 = [];
    h_c2_b2 = [];
    h_qjetVol = [];
    h_mass_trim = [];
    h_mass_mmdt = [];
    h_mass_prun = [];
    h_mass_sdb2 = [];
    h_mass_sdm1 = [];
    h_multiplicity = [];
    
    h_bdtall = [];
    
    for i in range(len(types)):
        
        h_pt1.append( ROOT.TH1F("h_jpt1_"+types[i],";pT_{j1} (GeV); fraction of events",50,0,2000) );
        h_pt2.append( ROOT.TH1F("h_jpt2_"+types[i],";pT_{j2} (GeV); fraction of events",50,0,2000) );
        h_eta1.append( ROOT.TH1F("h_jeta1_"+types[i],";#eta_{j1}; fraction of events",20,-5,5) );
        h_eta2.append( ROOT.TH1F("h_jeta2_"+types[i],";#eta_{j2}; fraction of events",20,-5,5) );
        h_mass1.append( ROOT.TH1F("h_jmass1_"+types[i],";m (GeV); fraction of events",mass_obs_nb,0,200) );
        h_mass2.append( ROOT.TH1F("h_jmass2_"+types[i],";m (GeV); fraction of events",obs_nb,0,200) );
        
        h_tau1_b1.append( ROOT.TH1F("h_jtau1_b1"+types[i],";#tau_{1}^{#beta=1}; fraction of events",obs_nb,0,0.5) );
        h_tau1_b2.append( ROOT.TH1F("h_jtau1_b2"+types[i],";#tau_{1}^{#beta=2}; fraction of events",obs_nb,0,0.5) );
        h_tau21_b1.append( ROOT.TH1F("h_jtau21_b1"+types[i],";#tau_{21}^{#beta=1}; fraction of events",obs_nb,0,1.0) );
        h_tau21_b2.append( ROOT.TH1F("h_jtau21_b2"+types[i],";#tau_{21}^{#beta=2}; fraction of events",obs_nb,0,1.0) );
        h_c1_b0.append( ROOT.TH1F("h_jc1_b0"+types[i],";C_{1}^{#beta=0}; fraction of events",obs_nb,0.2,0.5) );
        h_c1_b1.append( ROOT.TH1F("h_jc1_b1"+types[i],";C_{1}^{#beta=1}; fraction of events",obs_nb,0,0.3) );
        h_c1_b2.append( ROOT.TH1F("h_jc1_b2"+types[i],";C_{1}^{#beta=2}; fraction of events",obs_nb,0,0.3) );
        h_c2_b1.append( ROOT.TH1F("h_jc2_b1"+types[i],";C_{2}^{#beta=1}; fraction of events",obs_nb,0,1.0) );
        h_c2_b2.append( ROOT.TH1F("h_jc2_b2"+types[i],";C_{2}^{#beta=2}; fraction of events",obs_nb,0,1.0) );
        
        h_qjetVol.append( ROOT.TH1F("h_j_qjetVol"+types[i],";#Gamma_{Qjet}; fraction of events",obs_nb,0,1) );
        h_mass_trim.append( ROOT.TH1F("h_j_mass_trim"+types[i],";m_{trim} (GeV); fraction of events",mass_obs_nb,0,240) );
        h_mass_mmdt.append( ROOT.TH1F("h_j_mass_mmdt"+types[i],";m_{mmdt} (GeV); fraction of events",mass_obs_nb,0,240) );
        h_mass_prun.append( ROOT.TH1F("h_j_mass_prun"+types[i],";m_{prun} (GeV); fraction of events",mass_obs_nb,0,240) );
        h_mass_sdb2.append( ROOT.TH1F("h_j_mass_sdb2"+types[i],";m_{sd}^{#beta=2} (GeV); fraction of events",mass_obs_nb,0,240) );
        h_mass_sdm1.append( ROOT.TH1F("h_mass_sdm1"+types[i],";m_{sd}^{#beta=-1} (GeV); fraction of events",mass_obs_nb,0,240) );
        h_multiplicity.append( ROOT.TH1F("h_multiplicity"+types[i],";n_{constits}; fraction of events",50,0,150) );

        h_bdtall.append( ROOT.TH1F("h_bdtall"+types[i],";bdt all; N",bdt_nb,-1,1) );

    ## =====================================================================================
    ## 2D setup
    ## =====================================================================================

    if options.do2D:
        obs2Corr_all = [];
        obs2Corr = [];
        obs2Cor_name = [];
        obs2Corr_nb = [];
        obs2Corr_lo = [];
        obs2Corr_hi = [];
        
        if options.tmvaname == "qg":
            #observables to correlate -- quark vs. gluon
            obs2Corr_all = ["jmass","jtau1_b1","jtau1_b2","jc1_b0","jc1_b1","jc1_b2","j_multiplicity","j_qjetVol"];
            obs2Corr = ["jmass","jtau1_b1","jtau1_b2","jc1_b0","jc1_b1","jc1_b2","j_multiplicity","j_qjetVol"];
            obs2Corr_name = ["m","#tau_{1}^{#beta=1}","#tau_{1}^{#beta=2}",
                             "C_{1}^{#beta=0}","C_{1}^{#beta=1}","C_{1}^{#beta=2}",
                             "n_{constits}","#Gamma_{Qjet}"]
            #obs2Corr_nb = [40,40,40,40,40,40,40,40,40,40];
            obs2Corr_nb = [40]*8;
            obs2Corr_lo = [  0,  0,  0,0.2,  0,  0,  0,  0];
            obs2Corr_hi = [200,0.5,0.5,0.5,0.3,0.3,150,  1];


            # obs2Corr_all = ["j_multiplicity","j_qjetVol","jc1_b0","jtau1_b1","jc1_b1","jtau1_b2","jc1_b2","jmass"];
            # obs2Corr = ["j_multiplicity","j_qjetVol","jc1_b0","jtau1_b1","jc1_b1","jtau1_b2","jc1_b2","jmass"];
            # obs2Corr_name = ["n_{constits}","#Gamma_{Qjet}","C_{1}^{#beta=0}","#tau_{1}^{#beta=1}","C_{1}^{#beta=1}","#tau_{1}^{#beta=2}","C_{1}^{#beta=2}","m"]
            # obs2Corr_nb = [40]*8;
            # obs2Corr_lo = [  0,  0,0.2,  0,  0,  0,  0,  0];
            # obs2Corr_hi = [150,  1,0.5,0.5,0.3,0.5,0.3,200];

        if options.tmvaname == "Wq" or options.tmvaname == "Wg":
            
            ### FULL SET
            # #observables to correlate -- W vs. quark/gluon
            #obs2Corr_all = ["jmass","jtau21_b1","jtau21_b2",
                        # "jc2_b1","jc2_b2",
                        # "j_qjetVol",
                        # "j_mass_trim","j_mass_mmdt","j_mass_prun","j_mass_sdb2","j_mass_sdm1"];

            obs2Corr_all = ["jmass","jtau21_b1",
                        "jc2_b1",
                        "j_qjetVol",
                        "j_mass_trim","j_mass_mmdt","j_mass_prun","j_mass_sdb2"];                        
            # # obs2Corr_name = ["m","#tau_{21}^{#beta=1}","#tau_{21}^{#beta=2}",
            #                  "C_{2}^{#beta=1}","C_{2}^{#beta=2}",
            #                  "#Gamma_{Qjet}",
            #                  "m_{trim}","m_{mmdt}","m_{prun}","m_{sd}^{#beta=2}","m_{sd}^{#beta=-1}"]
            # obs2Corr_nb = [40]*11;
            # obs2Corr_lo = [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0];
            # obs2Corr_hi = [200,1.0,1.0,0.6,0.6,  1,200,200,200,200,200];

            ### REDUCED SET
            obs2Corr = ["jmass","jtau21_b1",
                        "jc2_b1",
                        "j_qjetVol",
                        "j_mass_trim","j_mass_mmdt","j_mass_prun","j_mass_sdb2"];
            obs2Corr_name = ["m","#tau_{21}^{#beta=1}",
                             "C_{2}^{#beta=1}",
                             "#Gamma_{Qjet}",
                             "m_{trim}","m_{mmdt}","m_{prun}","m_{sd}^{#beta=2}"]
            obs2Corr_nb = [40]*8;
            obs2Corr_lo = [  0,  0,  0,  0,  0,  0,  0,  0];
            obs2Corr_hi = [200,1.0,0.6,  1,200,200,200,200];

            varFor3Vars = ["j_mass_sdb2","jc2_b1"]
            varFor3Vars_name = ["m_{sd}^{#beta=2}","C_{2}^{#beta=1}"]

            ### just 2 vars
            # obs2Corr = ["jmass","j_mass_mmdt","j_mass_sdb2"];
            # obs2Corr_name = ["m","m_{mmdt}","m_{sd}^{#beta=2}"];
            # obs2Corr_nb = [40]*3;
            # obs2Corr_lo = [  0,  0,  0];
            # obs2Corr_hi = [200,200,200];


        OFILE = ROOT.TFile(directory+"/PCAtrees.root","RECREATE");
        otrees = [];
        for i in range(len(types)): otrees.append( ROOT.TTree("tree_"+types[i],"tree_"+types[i]) );
        aOvars = [];
        for i in range(len(types)): 
            tmp_aOvars = [];
            for j in range(len(obs2Corr)+1): tmp_aOvars.append( array( 'f', [ 0. ] ) );
            aOvars.append(tmp_aOvars);            
        for i in range(len(types)): 
            for j in range(len(obs2Corr)): otrees[i].Branch("bdt_"+obs2Corr[j], aOvars[i][j] , "bdt_"+obs2Corr[j]+"/F");
            otrees[i].Branch("bdt_all", aOvars[i][len(obs2Corr)] , "bdt_all/F");


        h_bdts = [];
        h_bdts_vars = [];
        h_2D = [];
        for a in range(len(types)):
            h_2D_type = []
            h_bdts_type = [];
            for i in range(len(obs2Corr)):
                for j in range(i,len(obs2Corr)):
                    hname = "h2d_"+obs2Corr[i]+"_"+obs2Corr[j]+"_"+types[a];
                    htitle = ";"+obs2Corr_name[i]+";"+obs2Corr_name[j];
                    #print "hname = ",hname
                    h_2D_type.append( ROOT.TH2F(hname,htitle,obs2Corr_nb[i],obs2Corr_lo[i],obs2Corr_hi[i],obs2Corr_nb[j],obs2Corr_lo[j],obs2Corr_hi[j]) );
                    
                    # bdt stuff
                    bdttitle = ";"+obs2Corr_name[i]+"+"+obs2Corr_name[j]+"; N";
                    h_bdts_type.append( ROOT.TH1F("h_bdt_"+hname,bdttitle,bdt_nb,-1,1) )
                    if a == 0: 
                        if i == j: h_bdts_vars.append( [obs2Corr[i]] );
                        else: h_bdts_vars.append( [obs2Corr[i],obs2Corr[j]] );
    
            if options.tmvaname == "Wq" or options.tmvaname == "Wg":                      
                for i in range(len(obs2Corr)):
                    if not (obs2Corr[i] in varFor3Vars):
                        hname = "h2d_"+obs2Corr[i]+"_"+obs2Corr[j]+"_"+types[a];
                        bdttitle = ";"+varFor3Vars_name[0]+"+"+varFor3Vars_name[1]+"+"+obs2Corr_name[i]+"; N";
                        h_bdts_type.append( ROOT.TH1F("h_bdt_"+hname,bdttitle,bdt_nb,-1,1) )
                        h_bdts_vars.append( [varFor3Vars[0],varFor3Vars[1],obs2Corr[i]] );

            h_2D.append( h_2D_type );
            h_bdts.append( h_bdts_type );
                    
        #######-----------------------------------------------
        #######-----------------------------------------------
        ## TMVA stuff
        ## classifier
        ## inputs, list of variables, list of ntuples
        bdtClasses = [];
        bdtctr = 0;
        for i in range(len(obs2Corr)):
            for j in range(i,len(obs2Corr)):
                if i == j: vars_tmp = [obs2Corr[i]+"[0]"];
                else: vars_tmp = [obs2Corr[i]+"[0]",obs2Corr[j]+"[0]"];

                tmvaname = "tmva_"+options.tmvaname+"_"+obs2Corr[i]+"_"+obs2Corr[j]+"_"+cone+"_pt%04i" % (binp1);
                #"weights/TMVAClassification_tmva_Wg_j_mass_mmdt_j_mass_sdb2_1000_ak08_BDTG.weights.xml"
                tmva_tmp = TMVAhelper("tmva_"+options.tmvaname+"_"+obs2Corr[i]+"_"+obs2Corr[j]+"_"+str(bin)+"_"+cone,vars_tmp,trainingFiles);
                bdtClasses.append(tmva_tmp);
                if options.train: 
                    bdtClasses[bdtctr].train(bin_cut,binp1_cut);
                bdtClasses[bdtctr].read(tmvamethod);            
                bdtctr += 1;

        if options.tmvaname == "Wq" or options.tmvaname == "Wg":                      
            for i in range(len(obs2Corr)):
                vars_tmp = [varFor3Vars[0]+"[0]",varFor3Vars[1]+"[0]"]
                if not (obs2Corr[i] in varFor3Vars):
                    vars_tmp.append(obs2Corr[i]+"[0]")

                    tmvaname = "tmva_"+options.tmvaname+"_"+varFor3Vars[0]+"_"+varFor3Vars[1]+"_"+obs2Corr[i]+"_"+cone+"_pt%04i" % (binp1);
                    #"weights/TMVAClassification_tmva_Wg_j_mass_mmdt_j_mass_sdb2_1000_ak08_BDTG.weights.xml"
                    tmva_tmp = TMVAhelper( "tmva_"+options.tmvaname+"_"+varFor3Vars[0]+"_"+varFor3Vars[1]+"_"+obs2Corr[i]+"_"+str(bin)+"_"+cone,vars_tmp,trainingFiles);
                    bdtClasses.append(tmva_tmp);
                    if options.train or options.train3: 
                        bdtClasses[bdtctr].train(bin_cut,binp1_cut);
                    bdtClasses[bdtctr].read(tmvamethod);            
                    bdtctr += 1;
                    print "bdtctr = ", bdtctr

        tmva_all = TMVAhelper("tmva_"+options.tmvaname+"_all"+"_"+cone+"_pt%04i" % (bin),obs2Corr_all,trainingFiles);
        if options.train or options.trainAllOnly: tmva_all.train(bin_cut,binp1_cut);
        tmva_all.read(tmvamethod);
            
    ## =====================================================================================
    ## =====================================================================================
    
    #######-----------------------------------------------
    ####### loop
    #######-----------------------------------------------
    
    for a in range(len(types)):
    
        f = ROOT.TFile(files[a]);
        t = f.Get("t");
    
        for i in range(t.GetEntries()):

            if i % 1000 == 0: print "file #",a,", event #",i
            if i > 10000000: break;
            
            t.GetEntry(i);
            h_pt1[a].Fill( t.jpt[0] );
            h_pt2[a].Fill( t.jpt[1] );
            h_eta1[a].Fill( t.jeta[0] );
            h_eta2[a].Fill( t.jeta[1] );
            h_mass1[a].Fill( t.jmass[0] );
            h_mass2[a].Fill( t.jmass[1] );

            if t.jpt[0] < bin_cut or t.jpt[0] > binp1_cut: continue;
            
            h_tau1_b1[a].Fill( t.jtau1_b1[0] );
            h_tau1_b2[a].Fill( t.jtau1_b2[0] );            
            h_tau21_b1[a].Fill( t.jtau21_b1[0] );
            h_tau21_b2[a].Fill( t.jtau21_b2[0] );
            
            h_c1_b0[a].Fill( t.jc1_b0[0] );                        
            h_c1_b1[a].Fill( t.jc1_b1[0] );                        
            h_c1_b2[a].Fill( t.jc1_b2[0] ); 
            h_c2_b1[a].Fill( t.jc2_b1[0] );
            h_c2_b2[a].Fill( t.jc2_b2[0] );
            
            h_qjetVol[a].Fill( t.j_qjetVol[0] );
            h_mass_trim[a].Fill( t.j_mass_trim[0] );
            h_mass_mmdt[a].Fill( t.j_mass_mmdt[0] );
            h_mass_prun[a].Fill( t.j_mass_prun[0] );
            h_mass_sdb2[a].Fill( t.j_mass_sdb2[0] );
            h_mass_sdm1[a].Fill( t.j_mass_sdm1[0] );
            h_multiplicity[a].Fill( t.j_multiplicity[0] );
            
            ## =====================================================================================
            ## 2D stuff
            if options.do2D:
                #2D plots
                obsctr = 0;
                for j in range(len(obs2Corr)):
                    for k in range(j,len(obs2Corr)):
                        h_2D[a][obsctr].Fill( getattr(t,obs2Corr[j])[0], getattr(t,obs2Corr[k])[0] );

                        tmpList = [];
                        for l in range(len(h_bdts_vars[obsctr])): 
                            tmpList.append( getattr(t,h_bdts_vars[obsctr][l])[0] );
                        bdtOutput = bdtClasses[obsctr].evaluate( tmpList, tmvamethod );
                        h_bdts[a][obsctr].Fill(bdtOutput);
                        
                        if j == k: 
                            aOvars[a][j][0] = bdtOutput;

                        obsctr += 1;

                if options.tmvaname == "Wq" or options.tmvaname == "Wg":                      
                    for j in range(len(obs2Corr)):
                        if not (obs2Corr[j] in varFor3Vars):
                            tmpList = [];
                            for l in range(len(h_bdts_vars[obsctr])): 
                                tmpList.append( getattr(t,h_bdts_vars[obsctr][l])[0] );
                            #print "obsctr = ", obsctr,", total length = ", len(bdtClasses)
                            bdtOutput = bdtClasses[obsctr].evaluate( tmpList, tmvamethod );
                            h_bdts[a][obsctr].Fill(bdtOutput);
                            
                            obsctr += 1;
            
                tmpList = [];
                for j in range(len(obs2Corr_all)):
                    tmpList.append( getattr(t,obs2Corr_all[j])[0] );
                h_bdtall[a].Fill( tmva_all.evaluate( tmpList, tmvamethod ) );
                aOvars[a][len(obs2Corr)][0] = tmva_all.evaluate( tmpList, tmvamethod );

            ## =====================================================================================
            otrees[a].Fill();

        del f;
        del t;
    
    # OFILE.cd();
    # for a in range(len(otrees)): otrees[a].Write();
    # OFILE.Close();

    #######-----------------------------------------------
    ####### plot
    #######-----------------------------------------------

    ## =====================================================================================
    ## 1D stuff

    makeCanvas(h_pt1, types_names, "jpt1", True);
    makeCanvas(h_pt2, types_names, "jpt2", True);
    makeCanvas(h_eta1, types_names, "jeta1", True);
    makeCanvas(h_eta2, types_names, "jeta2", True);
    makeCanvas(h_mass1, types_names, "jmass1", True);
    makeCanvas(h_mass2, types_names, "jmass2", True);
                
    makeCanvas(h_tau1_b1, types_names, "h_tau1_b1", True);
    makeCanvas(h_tau1_b2, types_names, "h_tau1_b2", True);
    makeCanvas(h_tau21_b1, types_names, "h_tau21_b1", True);
    makeCanvas(h_tau21_b2, types_names, "h_tau21_b2", True);
    makeCanvas(h_c1_b0, types_names, "h_c1_b0", True);
    makeCanvas(h_c1_b1, types_names, "h_c1_b1", True);
    makeCanvas(h_c1_b2, types_names, "h_c1_b2", True);
    makeCanvas(h_c2_b1, types_names, "h_c2_b1", True);
    makeCanvas(h_c2_b2, types_names, "h_c2_b2", True);

    makeCanvas(h_qjetVol, types_names, "h_qjetVol", True);
    makeCanvas(h_mass_trim, types_names, "h_mass_trim", True);
    makeCanvas(h_mass_mmdt, types_names, "h_mass_mmdt", True);
    makeCanvas(h_mass_prun, types_names, "h_mass_prun", True);
    makeCanvas(h_mass_sdb2, types_names, "h_mass_sdb2", True);
    makeCanvas(h_mass_sdm1, types_names, "h_mass_sdm1", True);
    makeCanvas(h_multiplicity, types_names, "h_multiplicity", True);

    makeCanvas(h_bdtall, types_names, "h_bdtall", True);

    ## =====================================================================================
    ## 2D stuff
    if options.do2D:
        
        #plot bdts
        obsctr = 0;
        for j in range(len(obs2Corr)):
            for k in range(j,len(obs2Corr)):
                tmpList = [];
                for a in range(len(types)): tmpList.append( h_bdts[a][obsctr] );
                makeCanvas( tmpList, types_names, h_bdts[0][obsctr].GetName(), True );
                obsctr += 1;

        for a in range(len(types)):
            #2D plots
            obsctr = 0;
            for j in range(len(obs2Corr)):
                for k in range(j+1,len(obs2Corr)):
                    makeSingleHistCanvas(h_2D[a][obsctr]);        
                    obsctr += 1;
        obsctr = 0;
        for j in range(len(obs2Corr)):
            for k in range(j+1,len(obs2Corr)):
                makeCanvas2D(h_2D[0][obsctr],h_2D[1][obsctr],types_names);
                obsctr += 1;

        # for Ben's utility
        theutils = myutils()
        ngp = 25;

        # single ROCs
        roc_mass    = theutils.ROC( h_mass1[0], h_mass1[1], ngp );
        SetLineColorStyleWidth( roc_mass, 1, 1, 6 );

        roc_tau1_b1 = makeROCFromHisto( h_tau1_b1[0], h_tau1_b1[1], 1 );
        SetLineColorStyleWidth( roc_tau1_b1, 2, 9, 6 );
        roc_tau1_b2 = makeROCFromHisto( h_tau1_b2[0], h_tau1_b2[1], 1 );
        SetLineColorStyleWidth( roc_tau1_b2, 4, 9, 6 );
        roc_tau21_b1 = makeROCFromHisto( h_tau21_b1[0], h_tau21_b1[1], 0 );
        SetLineColorStyleWidth( roc_tau21_b1, 2, 9, 6 );
        roc_tau21_b2 = makeROCFromHisto( h_tau21_b2[0], h_tau21_b2[1], 0  );
        SetLineColorStyleWidth( roc_tau21_b2, 4, 9, 6 );
        
        roc_c1_b0   = makeROCFromHisto( h_c1_b0[0], h_c1_b0[1], 1 );
        SetLineColorStyleWidth( roc_c1_b0, 6, 9, 6 );
        roc_c1_b1   = makeROCFromHisto( h_c1_b1[0], h_c1_b1[1], 1 );
        SetLineColorStyleWidth( roc_c1_b1, 7, 9, 6 );
        roc_c1_b2   = makeROCFromHisto( h_c1_b2[0], h_c1_b2[1], 1 );
        SetLineColorStyleWidth( roc_c1_b2, 3, 9, 6 );
        roc_c2_b1   = makeROCFromHisto( h_c2_b1[0], h_c2_b1[1], 0 );
        SetLineColorStyleWidth( roc_c2_b1, 7, 9, 6 );
        roc_c2_b2   = makeROCFromHisto( h_c2_b2[0], h_c2_b2[1], 0 );
        SetLineColorStyleWidth( roc_c2_b2, 3, 9, 6 );

        roc_qjetVol = makeROCFromHisto( h_qjetVol[0], h_qjetVol[1], 0 );
        SetLineColorStyleWidth( roc_qjetVol, 3, 1, 6 );

        roc_mass_trim = theutils.ROC( h_mass_trim[0], h_mass_trim[1], ngp );
        SetLineColorStyleWidth( roc_mass_trim, 2, 1, 6 );
        roc_mass_mmdt = theutils.ROC( h_mass_mmdt[0], h_mass_mmdt[1], ngp );
        SetLineColorStyleWidth( roc_mass_mmdt, 46, 1, 6 );
        roc_mass_prun = theutils.ROC( h_mass_prun[0], h_mass_prun[1], ngp );
        SetLineColorStyleWidth( roc_mass_mmdt, 6, 1, 6 );
        roc_mass_sdb2 = theutils.ROC( h_mass_sdb2[0], h_mass_sdb2[1], ngp );
        SetLineColorStyleWidth( roc_mass_sdb2, 7, 1, 6 );
        roc_mass_sdm1 = theutils.ROC( h_mass_sdm1[0], h_mass_sdm1[1], ngp );
        SetLineColorStyleWidth( roc_mass_sdm1, 4, 1, 6 );

        roc_multiplicity = makeROCFromHisto( h_multiplicity[0], h_multiplicity[1], 1 );
        SetLineColorStyleWidth( roc_mass_sdm1, 1, 9, 6 );

        roc_bdtall = makeROCFromHisto( h_bdtall[0], h_bdtall[1], 1 );
        SetLineColorStyleWidth( roc_bdtall, 28, 1, 6 );


        #bdt styles...
        bdtcolors = [1,2,3,4,6,7]
        bdtstyles = [1,3,5,7,9]
        bdtwidths = [2,4];

        bdtlines = [];
        for j in range(len(bdtstyles)):
            for k in range(len(bdtcolors)):
                for m in range(len(bdtwidths)):
                    bdtlines.append( [bdtcolors[k], bdtstyles[j], bdtwidths[m]] );
        roc_bdts = [];
        obsctr = 0;
        for j in range(len(obs2Corr)):
            for k in range(j,len(obs2Corr)):
                roc_bdts.append( makeROCFromHisto( h_bdts[0][obsctr], h_bdts[1][obsctr], 1 ) );
                roc_bdts[obsctr].SetLineColor( bdtlines[obsctr][0] );
                roc_bdts[obsctr].SetLineStyle( bdtlines[obsctr][1] );
                roc_bdts[obsctr].SetLineWidth( bdtlines[obsctr][2] );
                obsctr += 1;
        
        if options.tmvaname == "Wq" or options.tmvaname == "Wg":                      
            for j in range(len(obs2Corr)):
                if not (obs2Corr[j] in varFor3Vars):
                    roc_bdts.append( makeROCFromHisto( h_bdts[0][obsctr], h_bdts[1][obsctr], 1 ) );
                    roc_bdts[obsctr].SetLineColor( bdtlines[obsctr][0] );
                    roc_bdts[obsctr].SetLineStyle( bdtlines[obsctr][1] );
                    roc_bdts[obsctr].SetLineWidth( bdtlines[obsctr][2] );
                    obsctr += 1;

        ## combine all ROCs into one array...
        #### this is deprecated??
        allROCs = [];
        allROCsNames = [];
        # if options.tmvaname == "Wq" or options.tmvaname == "Wg" or options.tmvaname == "qg":
        #     allROCs.append( roc_mass );
        #     allROCsNames.append( "m" );
        # if options.tmvaname == "qg":
        #     allROCs.append( roc_tau1_b1 );
        #     allROCsNames.append( "#tau_{1}^{#beta=1}" );
        #     allROCs.append( roc_tau1_b2 );  
        #     allROCsNames.append( "#tau_{1}^{#beta=2}" );
        # if options.tmvaname == "Wq" or options.tmvaname == "Wg":
        #     allROCs.append( roc_tau21_b1 );
        #     allROCsNames.append( "#tau_{21}^{#beta=1}" );
        #     # allROCs.append( roc_tau21_b2 );
        #     # allROCsNames.append( "#tau_{21}^{#beta=2}" );
        # if options.tmvaname == "qg":
        #     allROCs.append( roc_c1_b0 );
        #     allROCsNames.append( "C_{1}^{#beta=0}" );
        #     allROCs.append( roc_c1_b1 );  
        #     allROCsNames.append( "C_{1}^{#beta=1}" );
        #     allROCs.append( roc_c1_b2 );  
        #     allROCsNames.append( "C_{1}^{#beta=2}" );
        # if options.tmvaname == "Wq" or options.tmvaname == "Wg":
        #     allROCs.append( roc_c2_b1 );
        #     allROCsNames.append( "C_{2}^{#beta=1}" );
        #     # allROCs.append( roc_c2_b2 );
        #     # allROCsNames.append( "C_{2}^{#beta=2}" );
        # if options.tmvaname == "Wq" or options.tmvaname == "Wg" or options.tmvaname == "qg":
        #     allROCs.append( roc_qjetVol );
        #     allROCsNames.append( "#Gamma_{Qjet}" );
        # if options.tmvaname == "Wq" or options.tmvaname == "Wg":
        #     allROCs.append( roc_mass_trim );
        #     allROCsNames.append( "m_{trim}" );
        #     allROCs.append( roc_mass_mmdt );
        #     allROCsNames.append( "m_{mmdt}" );
        #     allROCs.append( roc_mass_prun );  
        #     allROCsNames.append( "m_{prun}" );
        #     allROCs.append( roc_mass_sdb2 );  
        #     allROCsNames.append( "m_{sd}^{#beta=2}" );
        #     # allROCs.append( roc_mass_sdm1 );
        #     # allROCsNames.append( "m_{sd}^{#beta=-1}" );
        # if options.tmvaname == "qg":
        #     allROCs.append( roc_multiplicity );
        #     allROCsNames.append( "n_{constits}" );

        if options.tmvaname == "qg": vallo = 1e-3;
        if options.tmvaname == "Wg": vallo = 1e-4;
        if options.tmvaname == "Wq": vallo = 1e-4;

        obsctr = 0;
        for j in range(len(obs2Corr)):
            for k in range(j,len(obs2Corr)):
                allROCs.append( roc_bdts[obsctr] );  
                if j == k: allROCsNames.append( obs2Corr_name[j] );
                else: allROCsNames.append( obs2Corr_name[j]+"+"+obs2Corr_name[k] );
                obsctr += 1;

        if options.tmvaname == "Wq" or options.tmvaname == "Wg":                      
            for j in range(len(obs2Corr)):
                if not (obs2Corr[j] in varFor3Vars):
                    allROCs.append( roc_bdts[obsctr] );  
                    allROCsNames.append( varFor3Vars_name[0] + "-" + varFor3Vars_name[1] + "-" + obs2Corr_name[j] );
                    obsctr += 1;

        allROCs.append( roc_bdtall );
        allROCsNames.append( "allvars" );

        ## --------------------
        ## start drawing
        leg = ROOT.TLegend( 0.0, 0.0, 1.0, 1.0 );
        leg.SetBorderSize( 0 );
        leg.SetFillStyle( 0 );
        #leg.SetTextSize( 0.02 );
        leg.SetNColumns(4);
        for i in range(len(allROCs)):
            leg.AddEntry( allROCs[i], allROCsNames[i], "l" );

        canRoc = ROOT.TCanvas("canRoc","canRoc",2500,1250);    
        pad1 = ROOT.TPad("pad1","pad1",0.00,0.00,0.5,1.00)
        pad2 = ROOT.TPad("pad2","pad2",0.5,0.00,1.00,1.00)
        pad1.SetFillColor(0)
        pad2.SetFillColor(0)
        pad1.Draw()
        pad2.Draw()
        pad1.SetRightMargin(0.05);
        #pad2.SetTopMargin(0.05);
        
        pad2.cd();
        leg.Draw();
        pad1.cd();
        hrl = canRoc.DrawFrame(0,vallo,1.0,1);
        hrl.GetXaxis().SetTitle("#varepsilon_{sig}");
        hrl.GetYaxis().SetTitle("#varepsilon_{bkg}");
        for i in range(len(allROCs)):
            allROCs[i].Draw();

        canRoc.cd();
        gbanner.Draw();
        canRoc.SaveAs(directory+"/rocs/Rocs_1D_lin.eps");
        canRoc.SaveAs(directory+"/rocs/Rocs_1D_lin.png");
        canRoc.SaveAs(directory+"/rocs/Rocs_1D_lin.pdf");
        pad1.cd();
        ROOT.gPad.SetLogy();
        canRoc.SaveAs(directory+"/rocs/Rocs_1D.eps");
        canRoc.SaveAs(directory+"/rocs/Rocs_1D.png");
        canRoc.SaveAs(directory+"/rocs/Rocs_1D.pdf");

        ## ------------
        
        for a in range(len(obs2Corr)):

            wgbanner = ROOT.TLatex(0.65,0.96,"BOOST13WG");
            wgbanner.SetNDC();
            wgbanner.SetTextSize(0.030);  
            wgbanner.SetTextFont(52);              
            
            legSingle = ROOT.TLegend( 0.0, 0.1, 1.0, 1.0 );
            legSingle.SetBorderSize( 0 );
            legSingle.SetFillStyle( 0 );
            legSingle.SetTextSize( 0.12 );
            legSingle.SetTextFont( 42 );
            legSingle.SetNColumns(1);
            for i in range(len(allROCs)):
                if obs2Corr[a] in allROCs[i].GetName(): legSingle.AddEntry( allROCs[i], allROCsNames[i], "l" );

            canRocSingle = ROOT.TCanvas("canRocSingle","canRocSingle",2000,1300); 
            canRocSingle.SetRightMargin(0.05);   
            pad1s = ROOT.TPad("pad1","pad1",0.00,0.00,0.80,1.00)
            pad2s = ROOT.TPad("pad2","pad2",0.80,0.00,1.00,1.00)
            pad1s.SetFillColor(0)
            pad2s.SetFillColor(0)
            pad1s.Draw()
            pad2s.Draw()
            pad1s.SetRightMargin(0.05);
            pad1s.cd();
            hrlSingle = canRocSingle.DrawFrame(0,vallo,1.0,1);
            hrlSingle.GetXaxis().SetTitle("#varepsilon_{sig}");
            hrlSingle.GetYaxis().SetTitle("#varepsilon_{bkg}");
            for i in range(len(allROCs)):
                print allROCs[i].GetName()
                if obs2Corr[a] in allROCs[i].GetName(): 
                    allROCs[i].Draw();
                    print "match  = ",obs2Corr[a],"and",allROCs[i].GetName()
            ROOT.gPad.SetLogy();

            pad2s.cd();
            legSingle.Draw();
            canRocSingle.cd();
            gbanner.Draw();
            wgbanner.Draw();
            canRocSingle.SaveAs(directory+"/rocs/Rocs_1D_"+obs2Corr[a]+".eps");
            canRocSingle.SaveAs(directory+"/rocs/Rocs_1D_"+obs2Corr[a]+".png");
            canRocSingle.SaveAs(directory+"/rocs/Rocs_1D_"+obs2Corr[a]+".pdf");

            del legSingle
            del hrlSingle
            del canRocSingle
            
        ## ------------    
        bdtallIndex = len(allROCs)-1;
        leg1 = ROOT.TLegend( 0.0, 0.1, 1.0, 1.0 );
        leg1.SetBorderSize( 0 );
        leg1.SetFillStyle( 0 );
        leg1.SetTextSize( 0.15 );
        leg1.SetTextFont( 42 );
        leg1.SetNColumns(1);
        obsctr = 0;
        for j in range(len(obs2Corr)):
            for k in range(j,len(obs2Corr)):
                if j == k: leg1.AddEntry( allROCs[obsctr], allROCsNames[obsctr], "l" );
                obsctr += 1;
        leg1.AddEntry( allROCs[bdtallIndex], allROCsNames[bdtallIndex], "l" );


        wgbanner2 = ROOT.TLatex(0.65,0.96,"BOOST13WG");
        wgbanner2.SetNDC();
        wgbanner2.SetTextSize(0.030);  
        wgbanner2.SetTextFont(52);     

        canRoc1 = ROOT.TCanvas("canRoc1","canRoc1",2000,1300);    
        pad1a = ROOT.TPad("pad1","pad1",0.00,0.00,0.80,1.00)
        pad2a = ROOT.TPad("pad2","pad2",0.80,0.00,1.00,1.00)
        pad1a.SetFillColor(0)
        pad2a.SetFillColor(0)
        pad1a.Draw()
        pad2a.Draw()
        pad1a.SetRightMargin(0.05);
        pad2a.cd();
        leg1.Draw();
        pad1a.cd();
        hrl1 = canRoc1.DrawFrame(0,vallo,1.0,1);
        hrl1.GetXaxis().SetTitle("#varepsilon_{sig}");
        hrl1.GetYaxis().SetTitle("#varepsilon_{bkg}");
        obsctr = 0;
        for j in range(len(obs2Corr)):
            for k in range(j,len(obs2Corr)):
                if j == k: allROCs[obsctr].Draw("l");
                obsctr += 1;
        allROCs[bdtallIndex].Draw("l");
        canRoc1.cd();
        gbanner.Draw();
        wgbanner2.Draw();
        canRoc1.SaveAs(directory+"/rocs/Rocs_1D_single_lin.eps");
        canRoc1.SaveAs(directory+"/rocs/Rocs_1D_single_lin.png");
        canRoc1.SaveAs(directory+"/rocs/Rocs_1D_single_lin.pdf");
        pad1a.cd();
        ROOT.gPad.SetLogy();
        canRoc1.SaveAs(directory+"/rocs/Rocs_1D_single.eps");
        canRoc1.SaveAs(directory+"/rocs/Rocs_1D_single.png");
        canRoc1.SaveAs(directory+"/rocs/Rocs_1D_single.pdf");


        ## save pca tree
        OFILE.cd();
        for a in range(len(otrees)): otrees[a].Write();
        OFILE.Close();


        # ################# store all graphs?
        # graphHolderFN = "graphHolder_pt%04i-%04i_%s.root" % (bin,binp1,cone);
        # fileOfGraphs = ROOT.TFile(graphHolderFN,"RECREATE");
        # fileOfGraphs.cd();
        # for i in range(len(allROCs)): allROCs[i].Write();
        # fileOfGraphs.Close();

        ################# conglomerate
        # plot
        xlen = len(obs2Corr);
        if options.tmvaname == "Wq" or options.tmvaname == "Wg": xlen += 1;

        h2_effBkg = ROOT.TH2F("h2_eff",";;",xlen,0,xlen,len(obs2Corr),0,len(obs2Corr));
        # h2_effBkg.GetZaxis().SetTitle("#varepsilon_{bkg}^{(m only)}/#varepsilon_{bkg}");
        for i in range(len(obs2Corr)): h2_effBkg.GetXaxis().SetBinLabel(i+1,obs2Corr_name[i]);
        for i in range(len(obs2Corr)): h2_effBkg.GetYaxis().SetBinLabel(i+1,obs2Corr_name[i]);
        if options.tmvaname == "Wq" or options.tmvaname == "Wg": h2_effBkg.GetXaxis().SetBinLabel(len(obs2Corr)+1,varFor3Vars_name[0]+"+"+varFor3Vars_name[1]);

        print "obs2Corr_name = ",obs2Corr_name
        allVarsEffBkg = 0;
        theSigEff = options.sigEff;

        theNORMALIZATION = 1.;
        # for i in range(len(allROCs)):
        #     if "+" in allROCsNames[i]: continue;
        #     elif allROCsNames[i] == "allvars": continue;
        #     else: 
        #         singleindex = obs2Corr_name.index(allROCsNames[i]);                
        #         if singleindex == 0: theNORMALIZATION = allROCs[i].Eval(0.5);

        vars3EffBkg = [];
        for i in range(len(allROCs)):
            print allROCsNames[i], " = ", allROCs[i].Eval(0.5);
            # it's a 2 combo
            if "+" in allROCsNames[i]:
                splitted = allROCsNames[i].split("+");
                #print splitted;
                index1 = obs2Corr_name.index(splitted[0]);
                index2 = obs2Corr_name.index(splitted[1]);
                ##print "index1, index2 = ", index1, index2
                h2_effBkg.SetBinContent(index1+1,index2+1, theNORMALIZATION/allROCs[i].Eval(theSigEff));
            # it's all vars
            elif allROCsNames[i] == "allvars":
                allVarsEffBkg = theNORMALIZATION/allROCs[i].Eval(theSigEff);
                #print "allvars"
            elif "-" in allROCsNames[i]:
                splitted = allROCsNames[i].split("-");
                #vars3EffBkg.append( theNORMALIZATION/allROCs[i].Eval(theSigEff) );
                index3 = obs2Corr_name.index(splitted[2]);
                h2_effBkg.SetBinContent(len(obs2Corr)+1,index3+1, theNORMALIZATION/allROCs[i].Eval(theSigEff));
                #print "allvars"
            # its' single
            else:
                singleindex = obs2Corr_name.index(allROCsNames[i]);
                #print "singleindex = ", singleindex
                h2_effBkg.SetBinContent(singleindex+1,singleindex+1,theNORMALIZATION/allROCs[i].Eval(theSigEff));

        banner = ROOT.TLatex(0.45,0.25,"fixed #varepsilon_{sig} = %1.2f"%(theSigEff));
        banner.SetNDC();
        banner.SetTextSize(0.035);
        banner_all = ROOT.TLatex(0.35,0.2,"all variables, 1/#varepsilon_{bkg} = %1.1f"%(allVarsEffBkg));
        banner_all.SetNDC();
        banner_all.SetTextSize(0.035);

        #banner_z = ROOT.TLatex(0.95,0.9,"#varepsilon_{bkg}^{ (m only)}/#varepsilon_{bkg}");
        banner_z = ROOT.TLatex(0.95,0.9,"bkg. rejection (1/#varepsilon_{bkg})");
        banner_z.SetNDC();
        banner_z.SetTextSize(0.055);  
        banner_z.SetTextAngle(270);

        wgbanner2 = ROOT.TLatex(0.65,0.96,"BOOST13WG");
        wgbanner2.SetNDC();
        wgbanner2.SetTextSize(0.030);  
        wgbanner2.SetTextFont(52);   

        canEffBkg = ROOT.TCanvas("canEffBkg","canEffBkg",1000,800);
        canEffBkg.SetLeftMargin(0.13);
        h2_effBkg.SetMinimum(1.);
        if options.tmvaname == "qg": h2_effBkg.SetMaximum(20.);
        if options.tmvaname == "Wg": h2_effBkg.SetMaximum(1000.);
        if options.tmvaname == "Wq": h2_effBkg.SetMaximum(1000.);

        h2_effBkg.Draw("colz text");
        banner.Draw();
        banner_all.Draw();
        banner_z.Draw();
        gbanner.Draw();
        wgbanner2.Draw();
        ROOT.gPad.SetLogz();
        canEffBkg.SaveAs(directory+"/rocs/effBkg2D.eps");
        canEffBkg.SaveAs(directory+"/rocs/effBkg2D.png");
        canEffBkg.SaveAs(directory+"/rocs/effBkg2D.pdf");

        print "vars3EffBkg = ",vars3EffBkg


