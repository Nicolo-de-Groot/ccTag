/*
This macro shows how to compute jet energy scale.
root -l extract.C'("delphes_output.root", "outputtype")'

with outputtype = "gluon", "cc" or "zqq"
The output ROOT file contains the pT(MC)/pT(Reco) distributions for
various pT(Reco) and |eta| bins. The peak value of such distribution is
interpreted as the jet energy correction to be applied for that
given pT(Reco), |eta| bin.

*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"

#else
class ExRootTreeReader;
class ExRootResult;
#endif

//------------------------------------------------------------------------------


class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------


void PrintTree(TClonesArray *branchTree, int start)
{
    if (start <= 0) return;
    GenParticle *particle = (GenParticle *) branchTree->At(start);
    cout << particle->PID << endl;
    if (particle->D2 >= particle->D1)
    {
        for (int ii = particle->D1; ii <= particle->D2; ii++)
        {
            cout << ">";
            PrintTree(branchTree, ii);
        }
    }
    else
    {
        PrintTree(branchTree, particle->D1);
    }
}

// finds the charmonium particle in the Pythia Event record
int FindCC(TClonesArray *branchParticle)
{
    int pType = -1;
    int cc_index = -1;
    GenParticle *particle;
    // first find out which charmonium state is produced initially (Status = 23)
    for (int i = 0; i < branchParticle->GetEntriesFast(); ++i)
    {
        particle = (GenParticle *) branchParticle->At(i);
        if (particle->Status == 23)
        {
            pType = particle->PID;
            break;
        }
    }
    // than find the decay of this state (Status = 2) min PT of 25 GeV
    for (int i = 0; i < branchParticle->GetEntriesFast(); ++i)
    {
        particle = (GenParticle *) branchParticle->At(i);
        if (particle->PID == pType && particle->Status == 2 && particle->PT > 25)
        {
            cc_index = i;
            break;
        }
    }

    if (abs(particle->Eta) > 2.1)
    {
        cc_index = -1;
    }
    return cc_index;
}

void AnalyseEvents(ExRootTreeReader *treeReader, const char *outputType)
{
    TClonesArray *branchParticle = treeReader->UseBranch("ccParticle");
    TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchTower1 = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchTower2 = treeReader->UseBranch("EFlowNeutralHadron");

    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl; 

    int isSig =  1; 
    int iTgt = 0;
    std::ofstream outf;
    if (strncmp(outputType,"cc",2)==0)
    {
	isSig =1;
        outf.open("cc.csv");
	iTgt = 1;
    }
    else if (strncmp(outputType,"gluon",2) == 0)
    {
        outf.open("gluon.csv");
	isSig =0;
    }
    else 
    {
        outf.open("zqq.csv");
	isSig =-1;
    }
    TRandom3 *myRand = new TRandom3;

    Jet *jet, *genjet;
    GenParticle *particle;
    TObject *object;

    Track *track;
    Tower *tower;

    TLorentzVector jetMomentum, genJetMomentum, trackJet, momentum;

    Float_t deltaR, mass, Eecone[8], Econe[8], Pcone[8];
    Float_t pt, eta, dEta, dPhi;
    Int_t nCharged, nNeutral, jetQ, indexJ;
    Long64_t entry;

    Int_t i, j, iBin;

    // Loop over all events
    for (entry = 0; entry < allEntries; ++entry) {
        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        if (entry % 500 == 0) cout << "Event number: " << entry << endl;

        // find the first charmonium particle
	if (isSig >= 0) {
        	int cc_id = FindCC(branchParticle);
        	if (cc_id < 0) continue; // no J/psi, go to next event
        	particle = (GenParticle *) branchParticle->At(cc_id);
        	genJetMomentum = particle->P4();
	} else {
		genJetMomentum.SetPxPyPzE(0.0, 0.0, 50.0, 50.0);
	}

        // this is simply to avoid warnings from initial state particle
        // having infite rapidity ...

        indexJ = -1;
        bool jetOK;
	TVector3 jetDir;

        // Loop over all reconstructed jets in event
        for (j = 0; j < branchJet->GetEntriesFast(); ++j) {
            jetOK = false;
            int bTag = 0;
            jet = (Jet *) branchJet->At(j);
            jetMomentum = jet->P4();
	    jetDir = jetMomentum.Vect();
	    jetDir.SetMag(1.0);
            if (abs(jet->Eta) > 2.1) continue; // skip endcap jets
            if (jet->PT < 25) continue; // skip low Pt jets
            if (jet->PT > 80) continue; //-- to make the kinematics more similar to J/psi
            if (jet->NCharged < 2) continue; // remove this for calo jets
            deltaR = genJetMomentum.DeltaR(jetMomentum);
            if (isSig == 1) {
                jetOK = (deltaR < 0.2);
            } else if (isSig ==0) {
                jetOK = (jet->Flavor == 21 && deltaR > 1.0);
            } else {
		    jetOK = (jet->Flavor > 0 && jet->Flavor < 6 && deltaR > 1.0);
	    }
            if (!jetOK) continue;
            dEta = jet->DeltaEta;
            dPhi = jet->DeltaPhi;
            nNeutral = jet->NNeutrals;
            nCharged = jet->NCharged;
            jetQ = jet->Charge;
            mass = jet->Mass;

            // now look at jet constututents      

            for (int k = 0; k < 8; k++) Econe[k] = 0.;
            for (int k = 0; k < 8; k++) Eecone[k] = 0.;
            for (int k = 0; k < 8; k++) Pcone[k] = 0.;

            momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            trackJet.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            nCharged = 0;

            bTag = jet->BTag%2;
	    // fix b tag from wrong association, set probability to 0.2%
            if (isSig == 1 ){
                bTag = 0.;
                if (myRand->Rndm() > 0.998) {bTag = 1;}
            }

            Float_t Rem = 0.;
            Float_t Qjet = 0., SumRtPT = 0.;
            Float_t Rtrack = 0, SumPT = 0.;
            Float_t z, theta, SPT=0., LHA=0., MSS=0., WDT=0.;

            for (j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
            {
                object = jet->Constituents.At(j);

                // Check if the constituent is accessible

                if (object == 0) continue;
                if (entry < 5) {
                    cout << object->ClassName() << " found ..." << endl;
                }
                if (object->IsA() == GenParticle::Class()) {
                    particle = (GenParticle *) object;
                    momentum += particle->P4();
                } else if (object->IsA() == Track::Class()) {
                    track = (Track *) object;
                    momentum += track->P4();
                    trackJet += track->P4();
		            deltaR = jetMomentum.DeltaR(track->P4());
                    z = track->PT / jet-> PT;
                    theta = deltaR / 0.4;
                    if (track->PT > 0.4) {
                        nCharged++;
                        Qjet += track->Charge * pow(jetDir.Dot(track->P4().Vect()),0.5);
                        Rtrack += track->PT *deltaR;
                        SumPT += pow(jetDir.Dot(track->P4().Vect()),0.5);
                        SumRtPT += sqrt(track->PT);
                    }
                } else if (object->IsA() == Tower::Class()) {
                    tower = (Tower *) object;
                    momentum += tower->P4();
                    deltaR = jetMomentum.DeltaR(tower->P4());
                    Rem += tower->Eem *deltaR;
                    iBin = floor(deltaR / 0.1);
                    z = tower->ET / jet-> PT;
                    theta = deltaR / 0.4;
                    if (deltaR < 0.4) {
                        for (int ii = iBin; ii < 4; ii++) { Econe[ii] += tower->ET; }
                        for (int ii = iBin; ii < 4; ii++) { 
				            Eecone[ii] += tower->Eem / TMath::CosH(tower->Eta); 
			            }
                    }
                    else if (deltaR < 0.6) {
                        Econe[iBin] += tower->ET;
                        Eecone[iBin] += tower->Eem / TMath::CosH(tower->Eta);
                    }
                }
                LHA += z * sqrt(theta); 
                SPT += z * z;
                WDT += z * theta;
                MSS += z * theta * theta;
            }
            Rem /= (Eecone[3] + 0.00001);
            if (Rem > 1.5) Rem = 1.5;
            Float_t Fcore1 = Econe[0] / (Econe[3] + 0.00001);
            Float_t Fcore2 = Econe[1] / (Econe[3] + 0.00001);
            Float_t Fcore3 = Econe[2] / (Econe[3] + 0.00001);
            Float_t Fem = Eecone[3] / (Econe[3] + 0.00001);

            //
            //-- now loop over tracks and make p-cones
            //
            for (j = 0; j < branchTrack->GetEntriesFast(); ++j) {
                track = (Track *) branchTrack->At(j);
                deltaR = jetMomentum.DeltaR(track->P4());
                if (deltaR < 0.6)
                {
                    iBin = floor(deltaR / 0.1);
                    Pcone[iBin] += track->PT;
                }
            }
            //
            //-- all done, scale cones and fill the historgrams
            //
            for (int k = 0; k < 4; k++) Econe[k] /= Econe[3];
            for (int k = 1; k < 7; k++) Pcone[k] += Pcone[k - 1];
            for (int k = 0; k < 4; k++) Pcone[k] /= trackJet.Perp();

            Float_t Pcore1 = Pcone[0] / (Pcone[3] + 0.00001);
            Float_t Pcore2 = Pcone[1] / (Pcone[3] + 0.00001);
            Qjet /= (SumPT + 0.00001);
            Rtrack /= (SumRtPT + 0.000001);

            // -- output

            outf << trackJet.Mag() << ", " << dEta << ", " << dPhi << ", " << abs(jetQ) << ", ";
            outf << mass << ", " << bTag << ", ";
            outf << LHA << ", " << SPT << ", " << MSS << ", " << WDT << ", ";
            outf << nNeutral << ", " << nCharged << ", " << abs(Qjet) << ", " << Rtrack << ", " << Pcore1 << ", ";
            outf << Pcore2 << ", " << Fem << ", " << Rem << ", " << Fcore1 << ", " << Fcore2 << ", " << iTgt << std::endl;
        }
    }
    outf.close();
}

//------------------------------------------------------------------------------


void extract(const char *inputFile, const char *outputType)
{
    gSystem->Load("libDelphes");

    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);

    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    ExRootResult *result = new ExRootResult();

    AnalyseEvents(treeReader, outputType);

    cout << "** Exiting..." << endl;

    delete result;
    delete treeReader;
    delete chain;
}

//------------------------------------------------------------------------------
