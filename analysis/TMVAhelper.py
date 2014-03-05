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

        ## classifier
        ## inputs, list of variables, list of ntuples
        fout = ROOT.TFile("tmva_"+self._named+".root","RECREATE")        
        self._factory = ROOT.TMVA.Factory("TMVAClassification_"+self._named, fout,
                                    ":".join(["!V",
                                              "!Silent",
                                              "Color",
                                              "DrawProgressBar",
                                              "Transformations=I;D;P;G,D",
                                              "AnalysisType=Classification"]
                                             ));        
        
        self._reader = ROOT.TMVA.Reader()
        self._varRefs = [];

    # -------------------------------------    
    def train(self):
        
        print "training..."
        
        # define inputs
        for i in range(self._nInputVars):
            self._factory.AddVariable(self._inputVars[i],"F");

        fs = ROOT.TFile(self._files[0]);
        fb = ROOT.TFile(self._files[1]);
        ts = fs.Get("t");
        tb = fb.Get("t");
        self._factory.AddSignalTree( ts );
        self._factory.AddBackgroundTree( tb );

        cuts = ROOT.TCut("");
        self._factory.PrepareTrainingAndTestTree(cuts,   # signal events
                                           ":".join([
                                                     "nTrain_Signal=0",
                                                     "nTrain_Background=0",
                                                     "SplitMode=Random",
                                                     "NormMode=NumEvents",
                                                     "!V"
                                                     ]));

        method = self._factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT",
                                    ":".join([
                                              "!H",
                                              "!V",
                                              "NTrees=850",
                                              "nEventsMin=150",
                                              "MaxDepth=3",
                                              "BoostType=AdaBoost",
                                              "AdaBoostBeta=0.5",
                                              "SeparationType=GiniIndex",
                                              "nCuts=20",
                                              "PruneMethod=NoPruning",
                                              ]))

        self._factory.TrainAllMethods()
        self._factory.TestAllMethods()
        self._factory.EvaluateAllMethods()
        
        
    # -------------------------------------    
    def read(self):
    
        print "reading"
        for i in range(self._nInputVars): self._varRefs.append( array('f',[0]) );
        
        for i in range(self._nInputVars):
            self._reader.AddVariable(self._inputVars[i],self._varRefs[i]);        

        self._reader.BookMVA("BDT","weights/TMVAClassification_"+self._named+"_BDT.weights.xml")        
        #self._reader.BookMVA("BDT","weights/TMVAClassification_BDT.weights.xml")        
        
    # -------------------------------------            
    def evaluate(self,val):
            
        for i in range(self._nInputVars):
            self._varRefs[i][0] = val[i];  
        
        bdtOutput = self._reader.EvaluateMVA("BDT")
        return bdtOutput;      
        