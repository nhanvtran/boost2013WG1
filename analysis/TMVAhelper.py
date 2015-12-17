#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time

import ROOT

from ROOT import *
from array import array

class TMVAhelper:
    
    # -------------------------------------
    def __init__(self,named,inputVars,files):

        self._named = named;
        self._inputVars = inputVars;
        self._files = files;
        self._nInputVars = len(self._inputVars);

        #print "NUMBER OF FILES = ",len(self._files)

        ## classifier
        ## inputs, list of variables, list of ntuples
        fout = ROOT.TFile(self._named+".root","RECREATE")        
        self._factory = ROOT.TMVA.Factory("TMVAClassification_"+self._named, fout,
                                    ":".join(["!V",
                                              "!Silent",
                                              "Color",
                                              "DrawProgressBar",
                                              "Transformations=I",
                                              "AnalysisType=Classification"]
                                             ));        
        
        self._reader = ROOT.TMVA.Reader()
        self._varRefs = [];
        self._ptRef = array('f',[0]);

    # -------------------------------------    
    def train(self,ptlo,pthi):
        
        print "training..."
        
        # define inputs
        for i in range(self._nInputVars):
            self._factory.AddVariable(self._inputVars[i],"F");
        self._factory.AddSpectator("jpt[0]","F");

        fs = ROOT.TFile(self._files[0]);
        fb = ROOT.TFile(self._files[1]);
        ts = fs.Get("t");
        tb = fb.Get("t");
        self._factory.AddSignalTree( ts );
        self._factory.AddBackgroundTree( tb );

        cutstring = "(jpt[0] > "+str(ptlo)+") && (jpt[0] < "+str(pthi)+")";
        cuts = ROOT.TCut("");
        self._factory.PrepareTrainingAndTestTree(cuts,   # signal events
                                           ":".join([
                                                     "nTrain_Signal=0",
                                                     "nTrain_Background=0",
                                                     "SplitMode=Random",
                                                     "NormMode=NumEvents",
                                                     "!V"
                                                     ]));

        #method = self._factory.BookMethod(ROOT.TMVA.Types.kCuts, "CutsGA","!H:!V:FitMethod=GA:EffSel:Steps=40:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95");
        #method = self._factory.BookMethod(ROOT.TMVA.Types.kCuts, "CutsMC","!H:!V:FitMethod=MC:EffSel:");

        # method = self._factory.BookMethod(ROOT.TMVA.Types.kCuts, "CutsSA","!H:!V:FitMethod=SA:EffSel:KernelTemp=IncAdaptive:Eps=1e-10:UseDefaultScale");
        # method = self._factory.BookMethod(ROOT.TMVA.Types.kPDEFoam, "PDEFoam","!H:!V::VarTransform=I,N:CreateMVAPdfs:IgnoreNegWeightsInTraining:SigBgSeparate=F:TailCut=0.001"
        #                                                          ":VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T");
        # method = self._factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT",
        #                             ":".join([
        #                                       "!H",
        #                                       "!V",
        #                                       "NTrees=500",
        #                                       "MinNodeSize=5%",
        #                                       "MaxDepth=5",
        #                                       "BoostType=AdaBoost",
        #                                       "AdaBoostBeta=0.5",
        #                                       "SeparationType=GiniIndex",
        #                                       "nCuts=200",
        #                                       "PruneMethod=CostComplexity",
        #                                       "PruneStrength=5"
        #                                       ]));

        ###method = self._factory.BookMethod(ROOT.TMVA.Types.kBDT,"BDTG","!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=F:nCuts=25000:MaxDepth=3:SeparationType=GiniIndex");
        ###if self._inputVars <= 2: method = self._factory.BookMethod(ROOT.TMVA.Types.kBDT,"BDTG","!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=F:nCuts=25000:MaxDepth=3:SeparationType=GiniIndex");
        ###else:                    method = self._factory.BookMethod(ROOT.TMVA.Types.kBDT,"BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=F:nCuts=25000:MaxDepth=5:SeparationType=GiniIndex");
        method = self._factory.BookMethod(ROOT.TMVA.Types.kBDT,"BDTG","!H:!V:NTrees=400:BoostType=Grad:Shrinkage=0.1:UseBaggedGrad=F:nCuts=2000:NNodesMax=10000:MaxDepth=5:UseYesNoLeaf=F:nEventsMin=200");
        #method = self._factory.BookMethod(ROOT.TMVA.Types.kBDT,"BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.1:UseBaggedGrad=F:nCuts=10000:MaxDepth=3:UseYesNoLeaf=F:nEventsMin=200")
        ##method = self._factory.BookMethod(ROOT.TMVA.Types.kBDT,"BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad=F:nCuts=25000:MaxDepth=5:SeparationType=GiniIndex");

        self._factory.TrainAllMethods()
        self._factory.TestAllMethods()
        self._factory.EvaluateAllMethods()
        
        
    # -------------------------------------    
    def read(self,method="BDT"):
    
        print "reading"
        for i in range(self._nInputVars): self._varRefs.append( array('f',[0]) );
        
        for i in range(self._nInputVars):
            self._reader.AddVariable(self._inputVars[i],self._varRefs[i]);        
        self._reader.AddSpectator("jpt[0]",self._ptRef);

        self._reader.BookMVA(method,"weights/TMVAClassification_"+self._named+"_"+method+".weights.xml")        
        #self._reader.BookMVA("BDT","weights/TMVAClassification_BDT.weights.xml")        
        
    # -------------------------------------            
    def evaluate(self,val,method="BDT"):
            
        for i in range(self._nInputVars):
            self._varRefs[i][0] = val[i];  
        
        bdtOutput = self._reader.EvaluateMVA(method)
        return bdtOutput;      
        