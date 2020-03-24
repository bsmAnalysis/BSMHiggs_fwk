import ROOT as rt
import CMS_lumi, tdrstyle
import array
import os
import sys

from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double


#set the tdr style
tdrstyle.setTDRStyle()
rt.gStyle.SetOptStat(11111)
#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "" #"13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12

H_ref = 700; 
W_ref = 700; 
W = W_ref
H  = H_ref

# 
# Simple example of macro: plot with CMS name and lumi text
#  (this script does not pretend to work in all configurations)
# iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
# For instance: 
#               iPeriod = 3 means: 7 TeV + 8 TeV
#               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
#               iPeriod = 0 means: free form (uses lumi_sqrtS)
# Initiated by: Gautier Hamel de Monchenault (Saclay)
# Translated in Python by: Joshua Hardenbrook (Princeton)
# Updated by:   Dinko Ferencek (Rutgers)
#

iPeriod = 0

# references for T, B, L, R
T = 0.08*H_ref
B = 0.12*H_ref 
L = 0.08*W_ref
R = 0.12*W_ref

canvas = rt.TCanvas("c2","c2",50,50,W,H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
#canvas.SetTopMargin( T/H )
#canvas.SetBottomMargin( B/H )
canvas.SetTickx(0)
canvas.SetTicky(0)

# configure your card directory below after running ./goodnessfit.sh
directory = "/afs/cern.ch/work/y/yuanc/Analysis/H2a4b/CMSSW_10_2_13/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/JOBS/SB13TeV_SM_Zh_backup/test/out_60_4b"
datafname = "higgsCombineTest.GoodnessOfFit.mH120.root"
toyfname = "higgsCombineTest.GoodnessOfFit.mH120.123456.root"

dataf = rt.TFile.Open(os.path.join(directory, datafname))
datatree = dataf.limit
toyf = rt.TFile.Open(os.path.join(directory, toyfname))
toytree = toyf.limit

toytree.Draw("limit>>toy", "limit>=0","")
#toytree.Draw("limit>>toy({},0,{})".format(100,700), "limit>=0","")
hist = rt.gDirectory.Get("toy")
ymax = hist.GetMaximum()
hist.GetYaxis().SetRangeUser(0, ymax*1.2)

for iev in datatree:
  tobs = iev.limit
if(tobs>=hist.GetXaxis().GetXmax()):
  print("tobs: {}, max on Xaxis:{}".format(tobs, hist.GetXaxis().GetXmax()))
  toytree.Draw("limit>>toy({},0,{})".format(100,int(tobs*1.2)), "limit>=0","")
  hist = rt.gDirectory.Get("toy")
  ymax = hist.GetMaximum()
  hist.GetYaxis().SetRangeUser(0, ymax*1.2)
  #hist.GetXaxis().SetRangeUser(0, 1000)
line = rt.TLine(tobs, 0, tobs, ymax)
line.SetLineColor(rt.kRed)
line.SetLineWidth(3)
line.Draw("same")


#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)


canvas.cd()
canvas.Update()
#canvas.SetLogx(True)
canvas.RedrawAxis()

canvas.SaveAs("GOF.pdf")

canvas.Update()

raw_input("Press Enter to end")
