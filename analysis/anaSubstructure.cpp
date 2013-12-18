// for the paper with Bho


#include "LHEF.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include <cmath>
#include "TTree.h"

#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"

#include "LHEF.h"

//#ifdef __MAKECINT__
//#pragma link C++ class vector<float>+;
//#endif

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

/*
 
 TO COMPILE:
 
 export ROOTSYS=~/Desktop/root
 export PATH=$ROOTSYS/bin:$PATH
 export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
 export DYLD_LIBRARY_PATH=$ROOTSYS/lib:$DYLD_LIBRARY_PATH 
 
 c++ -o anaSubstructure `root-config --glibs --cflags` `/Users/ntran/Research/CMS/PhysicsAnalysis/boost2013/fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -lvectorDict -lEnergyCorrelator anaSubstructure.cpp
 
 TO RUN:     
 
 ./analysis01 ../madevent_3/Events/pietro_test14_unweighted_events.lhe ../Z2j/run_01_unweighted_events.lhe 
 
 */

//! ========================================================================================

////////////////////-----------------------------------------------
// Global variables
int njets;
std::vector<float> jpt;
std::vector<float> jeta;
std::vector<float> jmass;
std::vector<float> jtau1_b1;
std::vector<float> jtau2_b1;
std::vector<float> jtau1_b2;
std::vector<float> jtau2_b2;
std::vector<float> jc1_b0;
std::vector<float> jc1_b1;
std::vector<float> jc1_b2;

void analyzeEvent(std::vector < fastjet::PseudoJet > particles);
////////////////////-----------------------------------------------


int main (int argc, char **argv) {
    
    std::string type = argv[1];   // type "gg" or "qq"
    int bin = atoi(argv[2]);
    int binp1 = bin+100;          // pt bin
    int max = atoi(argv[3]);      // events to run over
    
    char inName[192];
    sprintf( inName, "data/pythia81-tune4c-lhc8-%s-pt0%i-0%i.lhe", type.c_str(), bin, bin+100 );
    std::cout << "fname = " << inName << std::endl;
    std::ifstream ifsbkg (inName) ;
    LHEF::Reader reader(ifsbkg) ;

    char outName[192];
    sprintf( outName, "data/boost2013-pythia81-tune4c-lhc8-%s-pt0%i-0%i.root", type.c_str(), bin, bin+100 );    
    TFile *f = TFile::Open(outName,"RECREATE");
    TTree *t = new TTree("t","Tree with vectors");
    t->Branch("njets"      , &njets      );
    t->Branch("jpt"        , &jpt        );
    t->Branch("jeta"       , &jeta       );
    t->Branch("jmass"      , &jmass      );
    t->Branch("jtau1_b1"   , &jtau1_b1   );
    t->Branch("jtau2_b1"   , &jtau2_b1   );
    t->Branch("jtau1_b2"   , &jtau1_b2   );
    t->Branch("jtau2_b2"   , &jtau2_b2   );
    t->Branch("jc1_b0"      , &jc1_b0      );
    t->Branch("jc1_b1"      , &jc1_b1      );
    t->Branch("jc1_b2"      , &jc1_b2      );

            
    int evtCtr = 0;
    std::vector < fastjet::PseudoJet > particles;

    // loop over events
    while ( reader.readEvent () ) {
        ++evtCtr;
        if (evtCtr % 100 == 0) std::cout << "event " << evtCtr << "\n";
        if (evtCtr > max) break;
        
        // per event
        particles.clear();
        jpt.clear();
        jeta.clear();
        jmass.clear();
        jtau1_b1.clear();
        jtau2_b1.clear();
        jtau1_b2.clear();
        jtau2_b2.clear();
        jc1_b0.clear();
        jc1_b1.clear();
        jc1_b2.clear();
                        
        for (unsigned int i = 0 ; i < reader.hepeup.IDUP.size(); ++i){

            if (reader.hepeup.ISTUP.at(i) == 1){
                float px = reader.hepeup.PUP.at(i).at(0);
                float py = reader.hepeup.PUP.at(i).at(1);
                float pz = reader.hepeup.PUP.at(i).at(2);
                float e  = reader.hepeup.PUP.at(i).at(3);                                    
                particles.push_back( fastjet::PseudoJet( px, py, pz, e ) );
            }   
        }
        analyzeEvent( particles );
        t->Fill();
        
    }
    
    std::cout << "finish loop" << std::endl;
    
    f->cd();
    t->Write();
    f->Close();
    
    delete f;
    return 0 ;
}

void analyzeEvent(std::vector < fastjet::PseudoJet > particles){
    
        // recluster on the fly....
    double rParam = 0.8;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    
    
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    
    fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(particles, jetDef, fjAreaDefinition);
    std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(25.0));
    
    fastjet::Pruner pruner1( fastjet::cambridge_algorithm, 0.1, 0.5 );
    fastjet::Filter trimmer1( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03)) );
        
    // n-subjettiness    
    double beta = 1; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
    double R0 = rParam; // Characteristic jet radius for normalization              
    double Rcut = rParam; // maximum R particles can be from axis to be included in jet   
    // beta = 1                   
    fastjet::contrib::Nsubjettiness nSub1KT_b1(1, fastjet::contrib::Njettiness::onepass_kt_axes, 1, R0, Rcut);
    fastjet::contrib::Nsubjettiness nSub2KT_b1(2, fastjet::contrib::Njettiness::onepass_kt_axes, 1, R0, Rcut);
    // beta = 2
    fastjet::contrib::Nsubjettiness nSub1KT_b2(1, fastjet::contrib::Njettiness::onepass_kt_axes, 2, R0, Rcut);
    fastjet::contrib::Nsubjettiness nSub2KT_b2(2, fastjet::contrib::Njettiness::onepass_kt_axes, 2, R0, Rcut);
    
    // ECF
    fastjet::contrib::EnergyCorrelatorDoubleRatio C2_b0(1,0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio C2_b1(1,1,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio C2_b2(1,2,fastjet::contrib::EnergyCorrelator::pt_R);
                    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    // FILL IN THE TREE
    njets = out_jets.size();
    for (unsigned int i = 0; i < out_jets.size(); i++){
        jpt.push_back( out_jets[i].pt() );
        jeta.push_back( out_jets[i].eta() );
        jmass.push_back( out_jets[i].m() );                

        jtau1_b1.push_back( nSub1KT_b1(out_jets.at(i)) );        
        jtau2_b1.push_back( nSub2KT_b1(out_jets.at(i)) );        
        jtau1_b2.push_back( nSub1KT_b2(out_jets.at(i)) );        
        jtau2_b2.push_back( nSub2KT_b2(out_jets.at(i)) );  
                
        jc1_b0.push_back( C2_b0(out_jets.at(i)) );            
        jc1_b1.push_back( C2_b1(out_jets.at(i)) );            
        jc1_b2.push_back( C2_b2(out_jets.at(i)) );                    
                
        //std::cout << " out_jets[i].pt() = " << out_jets[i].pt() << std::endl;           
    }
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    thisClustering->delete_self_when_unused();
    
}


