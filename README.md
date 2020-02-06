# ccTag
Files related to tagging hadronic charmonium decays

1. Simulation setup

* configNoLHE.cmnd lives in the examples Pythia8 directory, edit this to select the process: charmonium + g/gm or quarks

* delphes_card_ATLAS.tcl lives in the cards directory

* run the simulation with "./DelphesPythia8 cards/delphes_card_ATLAS.tcl examples/Pythia8/configNoLHE.cmnd outputfile.root"

2. Training data extraction from rootfiles

* use "root -l extract.C'("inputfile.root", "outputtype")", with outputtype = "cc", "gluon" or "zqq"

3. Tensorflow

* data.tar.gz contains the .csv files with training and test data

* run "python train.py trainfile testfile" to train the network using tensorflow

* run "python test.py trainfile testfile" to evaluate network
