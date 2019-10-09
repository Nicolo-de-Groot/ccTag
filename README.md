# ccTag
Files related to tagging hadronic charmonium decays

Simulation setup

configNoLHE.cmnd lives in the examples Pythia8 directory, edit this to select the process: charmonium + g/gm or quarks

delphes_card_ATLAS.tcl lives in the cards directory

run the simulation with DelphesPythia8 cccc xxxx

Training data extraction

use root

Tensorflow

data.tar.gz contains the .csv files with training and test data

run "python train.py trainfile testfile" to train the network using tensorflow

run "python test.py trainfile testfile" to evaluate network
