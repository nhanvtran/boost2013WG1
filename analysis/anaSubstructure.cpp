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
//#include "fastjet/contrib/ModifiedMassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "LHEF.h"
#include "QjetsPlugin.h"
#include "Qjets.h"

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
int evtCtr;

int njets;
std::vector<float> jpt;
std::vector<float> jeta;
std::vector<float> jmass;
std::vector<float> jtau1_b1;
std::vector<float> jtau2_b1;
std::vector<float> jtau1_b2;
std::vector<float> jtau2_b2;
std::vector<float> jtau21_b2;
std::vector<float> jtau21_b1;
std::vector<float> jc1_b0;
std::vector<float> jc1_b1;
std::vector<float> jc1_b2;
std::vector<float> jc2_b1;
std::vector<float> jc2_b2;
std::vector<float> j_qjetVol;
std::vector<float> j_mass_trim;
std::vector<float> j_mass_mmdt;
std::vector<float> j_mass_prun;
std::vector<float> j_mass_sdb2;
std::vector<float> j_mass_sdm1;
std::vector<float> j_multiplicity;

void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal);

float FindRMS( std::vector< float > qjetmasses );
float FindMean( std::vector< float > qjetmasses );
double getQjetVolatility(std::vector < fastjet::PseudoJet > constits, int QJetsN = 25, int seed = 12345);

////////////////////-----------------------------------------------

int main (int argc, char **argv) {
    
    std::string type = argv[1];   // type "gg" or "qq"
    int bin = atoi(argv[2]);
    int binp1 = bin+100;          // pt bin
    int min = atoi(argv[3]);      // events to run over
    int max = atoi(argv[4]);      // events to run over
    std::string tag = argv[5];
    float rVal = atof(argv[6]);
    
    char inName[192];
    sprintf( inName, "data/%s-pt%04i-%04i.lhe", type.c_str(), bin, bin+100 );
    std::cout << "fname = " << inName << std::endl;
    std::ifstream ifsbkg (inName) ;
    LHEF::Reader reader(ifsbkg) ;

    char outName[192];
    int rInt = (int) (rVal*10.);
    sprintf( outName, "data/boost2013-%s-pt%04i-%04i%s_ak%02i.root", type.c_str(), bin, bin+100, tag.c_str(), rInt );
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
    t->Branch("jtau21_b1"   , &jtau21_b1   );
    t->Branch("jtau21_b2"   , &jtau21_b2   );
    t->Branch("jc1_b0"      , &jc1_b0      );
    t->Branch("jc1_b1"      , &jc1_b1      );
    t->Branch("jc1_b2"      , &jc1_b2      );
    t->Branch("jc2_b1"      , &jc2_b1      );
    t->Branch("jc2_b2"      , &jc2_b2      );
    t->Branch("j_qjetVol"   , &j_qjetVol      );
    t->Branch("j_mass_trim"      , &j_mass_trim      );
    t->Branch("j_mass_mmdt"      , &j_mass_mmdt      );
    t->Branch("j_mass_prun"      , &j_mass_prun      );
    t->Branch("j_mass_sdb2"      , &j_mass_sdb2      );
    t->Branch("j_mass_sdm1"      , &j_mass_sdm1      );
    t->Branch("j_multiplicity"   , &j_multiplicity      );
            
    evtCtr = 0;
    std::vector < fastjet::PseudoJet > particles;

    // loop over events
    while ( reader.readEvent () ) {
        ++evtCtr;
        if (evtCtr % 1000 == 0) std::cout << "event " << evtCtr << "\n";
        if (evtCtr < min) continue;
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
        jtau21_b1.clear();
        jtau21_b2.clear();
        jc1_b0.clear();
        jc1_b1.clear();
        jc1_b2.clear();
        jc2_b1.clear();
        jc2_b2.clear();
        j_qjetVol.clear();
        j_mass_trim.clear();
        j_mass_mmdt.clear();
        j_mass_prun.clear();
        j_mass_sdb2.clear();
        j_mass_sdm1.clear();
        j_multiplicity.clear();
        
        for (unsigned int i = 0 ; i < reader.hepeup.IDUP.size(); ++i){

            if (reader.hepeup.ISTUP.at(i) == 1){
                float px = reader.hepeup.PUP.at(i).at(0);
                float py = reader.hepeup.PUP.at(i).at(1);
                float pz = reader.hepeup.PUP.at(i).at(2);
                float e  = reader.hepeup.PUP.at(i).at(3);                                    
                particles.push_back( fastjet::PseudoJet( px, py, pz, e ) );
            }   
        }
        analyzeEvent( particles, rVal );
        t->Fill();
        
    }
    
    std::cout << "finish loop" << std::endl;
    
    f->cd();
    t->Write();
    f->Close();
    
    delete f;
    return 0 ;
}

void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal){
    
        // recluster on the fly....
    double rParam = rVal;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    
    
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    
    fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(particles, jetDef, fjAreaDefinition);
    std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(25.0));
    
    // groomers/taggers
    fastjet::Pruner pruner1( fastjet::cambridge_algorithm, 0.1, 0.5 );
    fastjet::Filter trimmer1( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03)) );
    
    // use just a symmetry cut for the tagger, with no mass-drop requirement
    //double z_cut = 0.10;
    //ModifiedMassDropTagger tagger(z_cut);
    
    double beta_sd = 1.0;
    double zcut_sd = 0.1;
    double mu_sd   = 1.0;
    contrib::SoftDropTagger soft_drop_mmdt(0.0, zcut_sd, mu_sd);
    contrib::SoftDropTagger soft_drop_sdb2(2.0, zcut_sd, mu_sd);
    contrib::SoftDropTagger soft_drop_sdm1(-1.0, zcut_sd, mu_sd);
    
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
    fastjet::contrib::EnergyCorrelatorDoubleRatio C1_b0(1,0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio C1_b1(1,1,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio C1_b2(1,2,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio C2_b1(2,1,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio C2_b2(2,2,fastjet::contrib::EnergyCorrelator::pt_R);
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    // FILL IN THE TREE
    njets = out_jets.size();
    for (unsigned int i = 0; i < out_jets.size(); i++){
    
        jpt.push_back( out_jets[i].pt() );
        jeta.push_back( out_jets[i].eta() );
        jmass.push_back( out_jets[i].m() );                

        // N-subjettiness
        jtau1_b1.push_back( nSub1KT_b1(out_jets.at(i)) );        
        jtau2_b1.push_back( nSub2KT_b1(out_jets.at(i)) );        
        jtau1_b2.push_back( nSub1KT_b2(out_jets.at(i)) );        
        jtau2_b2.push_back( nSub2KT_b2(out_jets.at(i)) );  
        jtau21_b1.push_back( nSub2KT_b1(out_jets.at(i)) / nSub1KT_b1(out_jets.at(i)) );
        jtau21_b2.push_back( nSub2KT_b2(out_jets.at(i)) / nSub1KT_b2(out_jets.at(i)) );
        
        // energy correlator        
        jc1_b0.push_back( C1_b0(out_jets.at(i)) );
        jc1_b1.push_back( C1_b1(out_jets.at(i)) );
        jc1_b2.push_back( C1_b2(out_jets.at(i)) );
        jc2_b1.push_back( C2_b1(out_jets.at(i)) );
        jc2_b2.push_back( C2_b2(out_jets.at(i)) );
        
        // Qjets computation
        int QJetsPreclustering = 999;
        std::vector<fastjet::PseudoJet> constits;        
        if (i == 0){ // only fill qjets for the hardest jet in the event to save time
            unsigned int nqjetconstits = out_jets.at(i).constituents().size();
            if (nqjetconstits < (unsigned int) QJetsPreclustering) constits = out_jets.at(i).constituents();
            else constits = out_jets.at(i).associated_cluster_sequence()->exclusive_subjets_up_to(out_jets.at(i),QJetsPreclustering);        

            j_qjetVol.push_back( getQjetVolatility(constits, 25, evtCtr*25) );            
        }
        else{
            j_qjetVol.push_back( -1. );
        }
        constits.clear();
         
        // groomed masses
//        PseudoJet tagged_jet = tagger(out_jets.at(i));
//        if (tagged_jet != 0){
//            double Rfilt = min(0.3, tagged_jet.structure_of<ModifiedMassDropTagger>().delta_R()*0.5);
//            int    nfilt = 3;
//            Filter filter(Rfilt, SelectorNHardest(nfilt));
//            PseudoJet filtered_jet = filter(tagged_jet);
//            j_mass_mmdt.push_back( filtered_jet.m() );
//        }
//        else{
//            j_mass_mmdt.push_back( -1 );
//        }
        
        
        j_mass_trim.push_back( trimmer1( out_jets.at(i) ).m() );
        j_mass_prun.push_back( pruner1( out_jets.at(i) ).m() );    
        j_mass_mmdt.push_back( soft_drop_mmdt( out_jets.at(i) ).m() );
        j_mass_sdb2.push_back( soft_drop_sdb2( out_jets.at(i) ).m() );
        j_mass_sdm1.push_back( soft_drop_sdm1( out_jets.at(i) ).m() );
        
        j_multiplicity.push_back( (float) out_jets.at(i).constituents().size() );

    }
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    thisClustering->delete_self_when_unused();
    
}

// ----------------------------------------------------------------------------------

// -----------------------------------------
// other helpers
// -----------------------------------------

float FindRMS( std::vector< float > qjetmasses ){
    
    float total = 0.;
    float ctr = 0.;
    for (unsigned int i = 0; i < qjetmasses.size(); i++){
        total = total + qjetmasses[i];
        ctr++;
    }
    float mean =  total/ctr;
    
    float totalsquared = 0.;
    for (unsigned int i = 0; i < qjetmasses.size(); i++){
        totalsquared += (qjetmasses[i] - mean)*(qjetmasses[i] - mean) ;
    }
    float RMS = sqrt( totalsquared/ctr );
    return RMS;
}

float FindMean( std::vector< float > qjetmasses ){
    float total = 0.;
    float ctr = 0.;
    for (unsigned int i = 0; i < qjetmasses.size(); i++){
        total = total + qjetmasses[i];
        ctr++;
    }
    return total/ctr;
}

double getQjetVolatility(std::vector < fastjet::PseudoJet > constits, int QJetsN, int seed){
    
    std::vector< float > qjetmasses;
    
    double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1), truncationFactor(0.01);          
    
    for(unsigned int ii = 0 ; ii < (unsigned int) QJetsN ; ii++){
        QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncationFactor);
        qjet_plugin.SetRandSeed(seed+ii); // new feature in Qjets to set the random seed
        fastjet::JetDefinition qjet_def(&qjet_plugin);
        fastjet::ClusterSequence qjet_seq(constits, qjet_def);
        vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(5.0));
        
        if (inclusive_jets2.size()>0) { qjetmasses.push_back( inclusive_jets2[0].m() ); }
        
    }
    
    // find RMS of a vector
    float qjetsRMS = FindRMS( qjetmasses );
    // find mean of a vector
    float qjetsMean = FindMean( qjetmasses );
    float qjetsVolatility = qjetsRMS/qjetsMean;
    return qjetsVolatility;
}



