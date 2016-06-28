# vhHad

This repository contains the 76X versions of code to produce the trees needed to set the limits for HH->bbbb. 
The files rundata_76X_vh.sh , runsig_76X_vh.sh, and runttbar_76X_vh.sh run over the Heppy Ntuples using generalTreeAnalyzer_76X_vh.py
and drawing the name of the files from the .txt files in the TxtFiles/ directory. The options for running, specified in the .sh file, 
include trigger (apply the trigger cut), jets (apply the jet pt and eta cuts), deta (apply the delta eta cut), xsec (cross section), 
syst (systematic being run, if any), and isMC (to distinguish tree information meant to be only in the MC trees). The remaining options
configure which file(s) is (are) being run over.

Once you have these files, you can run them through the secondary analyzer, ttreeAnalyzer_76X_vh.py using runttree_vh.sh. This produces
histograms based on cuts defined in the file, which you can then organize into a .tex file cutflow in cutflowtable_vh.py, which is
run using the python command.
