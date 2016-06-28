import os
import glob
import math    

import ROOT 
#ROOT.gROOT.Macro("rootlogon.C")

import FWCore.ParameterSet.Config as cms

import sys
from DataFormats.FWLite import Events, Handle

from array import array

from optparse import OptionParser
parser = OptionParser()
#output file name
parser.add_option("-o", "--outName", dest="outName",
                  help="output file name")
(options, args) = parser.parse_args()
outputfilename = options.outName

#numberLimit = float(sys.argv[1])
#getting input file
f = ROOT.TFile(sys.argv[1])
#making output file
f1 =  ROOT.TFile(outputfilename, 'recreate')
#print outputfilename
f1.cd()
#getting old tree
treeMine  = f.Get('myTree')
#making new tree
mynewTree = ROOT.TTree('mynewTree', 'mynewTree')

print treeMine
print mynewTree

#getting old branches from old tree
jet1ptB = array('f', [-100.0])
jet2ptB = array('f', [-100.0])
jet1etaB = array('f', [-100.0])
jet2etaB = array('f', [-100.0])
etadiffB = array('f', [-100.0])
dijetmassB = array('f', [-100.0])
jet1tau21B = array('f', [-100.0])
jet2tau21B = array('f', [-100.0])
jet1pmassB = array('f', [-100.0])
jet2pmassB = array('f', [-100.0])
#jet1sdmassB = array('f', [-100.0])
#jet2sdmassB = array('f', [-100.0])
jet1bbtagB = array('f', [-100.0])
jet2bbtagB = array('f', [-100.0])
#jet1mscsvB = array('f', [-100.0])
#jet2mscsvB = array('f', [-100.0])
triggerpassbbB = array('f', [-100.0])
#triggerpassSJB = array('f', [-100.0])
#nHiggsTagsB = array('f', [-100.0])

treeMine.SetBranchAddress("jet1pt", jet1ptB )
treeMine.SetBranchAddress("jet2pt", jet2ptB )
treeMine.SetBranchAddress("jet1eta", jet1etaB )
treeMine.SetBranchAddress("jet2eta", jet2etaB )
treeMine.SetBranchAddress("etadiff", etadiffB )
treeMine.SetBranchAddress("dijetmass", dijetmassB )
treeMine.SetBranchAddress("jet1tau21", jet1tau21B )
treeMine.SetBranchAddress("jet2tau21", jet2tau21B )
treeMine.SetBranchAddress("jet1pmass", jet1pmassB )
treeMine.SetBranchAddress("jet2pmass", jet2pmassB )
#treeMine.SetBranchAddress("jet1sdmass", jet1sdmassB )
#treeMine.SetBranchAddress("jet2sdmass", jet2sdmassB )
treeMine.SetBranchAddress("jet1bbtag", jet1bbtagB )
treeMine.SetBranchAddress("jet2bbtag", jet2bbtagB )
#treeMine.SetBranchAddress("jet1mscsv", jet1mscsvB )
#treeMine.SetBranchAddress("jet2mscsv", jet2mscsvB )
treeMine.SetBranchAddress("triggerpassbb", triggerpassbbB ) 
#treeMine.SetBranchAddress("triggerpasssj", triggerpassSJB)
#treeMine.SetBranchAddress("nHiggsTags", nHiggsTagsB ) 

nevent = treeMine.GetEntries();

#setting histograms we need
c1 = ROOT.TH1F("c1", "After trigger, jet pt eta, deltaeta", 10000, 0, 10000)
c2 = ROOT.TH1F("c2", "After Mass", 10000, 0, 10000)
c3 = ROOT.TH1F("c3", "After Jet tau21 0.6", 10000, 0, 10000)
c4 = ROOT.TH1F("c4", "After Jet tau21 0.45", 10000, 0, 10000)
c5 = ROOT.TH1F("c5", "After Jet tau 21 DDT", 10000, 0, 10000)
c3a = ROOT.TH1F("c3a", "After bbtag 0.3", 10000, 0, 10000)
c3b = ROOT.TH1F("c3b", "After bbtag 0.6", 10000, 0, 10000)
c4a = ROOT.TH1F("c4a", "After bbtag 0.3", 10000, 0, 10000)
c4b = ROOT.TH1F("c4b", "After bbtag 0.6", 10000, 0, 10000)
c5a = ROOT.TH1F("c5a", "After bbtag 0.3", 10000, 0, 10000)
c5b = ROOT.TH1F("c5b", "After bbtag 0.6", 10000, 0, 10000)

for i in range(0, nevent) :
    treeMine.GetEntry(i)

    if triggerpassbbB[0] > 0:
        c1.Fill(dijetmassB[0])
    #jet 1 V jet 2 H
#        if 65 < jet1pmassB[0] < 105 and 105 < jet2pmassB[0] < 135:
#            c2.Fill(dijetmassB[0])
#            jet1tauDDT = jet1tau21B[0] + 0.063*math.log((jet1pmassB[0]**2)/jet1ptB[0])
#            jet2tauDDT = jet2tau21B[0] + 0.063*math.log((jet2pmassB[0]**2)/jet2ptB[0])
#            if jet1tau21B[0] < 0.6:
#                c3.Fill(dijetmassB[0])
#                if jet2bbtagB[0] > 0.3:
#                    c3a.Fill(dijetmassB[0])
#                if jet2bbtagB[0] > 0.6:
#                    c3b.Fill(dijetmassB[0])
#            if jet1tau21B[0] < 0.45:
#                c4.Fill(dijetmassB[0])
#                if jet2bbtagB[0] > 0.3:
#                    c4a.Fill(dijetmassB[0])
#                if jet2bbtagB[0] > 0.6:
#                    c4b.Fill(dijetmassB[0])
#            if jet1tauDDT <= 0.38:
#                c5.Fill(dijetmassB[0])
#                if jet2bbtagB[0] > 0.3:
#                    c5a.Fill(dijetmassB[0])
#                if jet2bbtagB[0] > 0.6:
#                    c5b.Fill(dijetmassB[0])
    #jet 1 H jet 2 V
        if 65 < jet2pmassB[0] < 105 and 105 < jet1pmassB[0] < 135:
            c2.Fill(dijetmassB[0])
            jet1tauDDT = jet1tau21B[0] + 0.063*math.log((jet1pmassB[0]**2)/jet1ptB[0])
            jet2tauDDT = jet2tau21B[0] + 0.063*math.log((jet2pmassB[0]**2)/jet2ptB[0])
            if jet2tau21B[0] < 0.6:
                c3.Fill(dijetmassB[0])
                if jet1bbtagB[0] > 0.3:
                    c3a.Fill(dijetmassB[0])
                if jet1bbtagB[0] > 0.6:
                    c3b.Fill(dijetmassB[0])
            if jet2tau21B[0] < 0.45:
                c4.Fill(dijetmassB[0])
                if jet1bbtagB[0] > 0.3:
                    c4a.Fill(dijetmassB[0])
                if jet1bbtagB[0] > 0.6:
                    c4b.Fill(dijetmassB[0])
            if jet2tauDDT <= 0.38:
                c5.Fill(dijetmassB[0])
                if jet1bbtagB[0] > 0.3:
                    c5a.Fill(dijetmassB[0])
                if jet1bbtagB[0] > 0.6:
                    c5b.Fill(dijetmassB[0])

print "Done with " + outputfilename 

f1.cd()
f1.Write()
f1.Close()

f.Close()


