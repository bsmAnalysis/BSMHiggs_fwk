#!/usr/bin/env python
import os,sys
import glob
import json
import ROOT as r
import getopt
import commands
import subprocess
from array import array
import math
from ROOT import SetOwnership

hists_3b = ['e_A_CR_3b','mu_A_CR_3b']
hists_4b = ['e_A_CR_4b','mu_A_CR_4b']

#startbin = 9 

def drawHist(hist, name):

    c = r.TCanvas(name,name,800,800)
    r.gPad.SetLogy(1)
    hist.SetStats(0)
    hist.SetLineColor(r.kBlack)
    hist.SetMarkerStyle(20)
    hist.SetLineWidth(2)
    hist.Draw("ehist")
    r.gPad.SetLogy(0)

    return c

def applyWeights(hist, name):
    for i in range(1, hist.GetNbinsX()+1):
        sf=1.0
	ht = hist.GetBinCenter(i)
	if ht>=160:
	    if '3b' in name:
	        if ht<380: sf = math.exp(0.11879-0.00080*ht)
		else: sf = -1.831+1.537*pow(10,-2)*ht-2.836*pow(10,-5)*ht**2+1.673*pow(10,-8)*ht**3
	    elif '4b' in name:
	        if ht<380: sf = math.exp(0.14235-0.00084*ht)
		#else: sf = -0.0824+4.805*pow(10,-3)*ht-7.863*pow(10,-6)*ht**2+3.979*pow(10,-9)*ht**3
		#else: sf = 1.237-2.101*pow(10,-3)*ht+3.931*pow(10,-6)*ht**2-2.599*pow(10,-9)*ht**3
		else: sf = 1.412-3.047*pow(10,-3)*ht+5.589*pow(10,-6)*ht**2-3.545*pow(10,-9)*ht**3
	hist.SetBinContent(i, hist.GetBinContent(i)*sf)

        

def weightedAverage(ratio, hNLO, threshold, end=800):
    start = ratio.GetXaxis().FindBin(threshold)
    stop = ratio.GetXaxis().FindBin(end)
    mean_num = 0
    mean_den = 0
    error = 0
    for i in range(start, stop+1):
        weight = hNLO.GetBinContent(i)
	value = ratio.GetBinContent(i)
	err = ratio.GetBinError(i)
	mean_num += weight*value
	mean_den += weight
	error += weight**2*err**2
    if mean_den == 0: # if sum of weights is 0, then set zero to the error and mean
        mean_den = stop+1 - start
	error = 0
	mean_num = 0
    error = math.sqrt(error)/mean_den
    mean = mean_num/mean_den
    for i in range(start,ratio.GetXaxis().GetNbins()+1):
        ratio.SetBinContent(i,mean)
	ratio.SetBinError(i,error)
    return ratio

def ratioPlot(hLO,hNLO,ratio_h,name):

    thred = 800
    if '3b' in name: thred = 400
    elif '4b' in name: thred = 400
#    ratio_h = weightedAverage(ratio_h,hNLO,thred)
    
    # Define the Canvas
    c = r.TCanvas(name,name,800,800)
    
    # Upper plot will be in pad1
    hLO.SetLineColor(r.kBlack)
    hLO.SetMarkerColor(r.kBlack)
    hLO.SetFillStyle(4000)
    hLO.SetLineWidth(2)
    hLO.SetTitle(name)
    hLO.SetMarkerStyle(21)
    hLO.SetMarkerSize(1)
    hNLO.SetLineColor(r.kBlue)
    hNLO.SetMarkerColor(r.kBlue)
    hNLO.SetFillStyle(4000)
    hNLO.SetLineWidth(2)
    hNLO.SetTitle(name)
    hNLO.SetMarkerStyle(21)
    hNLO.SetMarkerSize(1)
    
    ymin = min(hLO.GetMinimum(), hNLO.GetMinimum())
    ymax = 1.1*max(hLO.GetMaximum(), hNLO.GetMaximum())
    if ymin<0: ymin = 1.1*ymin
    else: ymin = 0.9*ymin
    hLO.GetYaxis().SetRangeUser(ymin,ymax)
#    hLO.GetXaxis().SetRangeUser(0,300)
    hLO.GetYaxis().SetTitleOffset(1.6)
    hLO.GetYaxis().SetTitleSize(25)
    hLO.GetYaxis().SetTitleFont(43)
    hLO.GetYaxis().SetLabelFont(43)
    hLO.GetYaxis().SetLabelSize(20)
    hLO.GetXaxis().SetTitleSize(25)
    hLO.GetXaxis().SetTitleFont(43)
    hLO.GetXaxis().SetLabelFont(43)
    hLO.GetXaxis().SetLabelSize(20)
    pad1 = r.TPad("pad1","pad1", 0., 0.3, 1, 1)
    pad1.SetBottomMargin(0) # Upper and lower plot are joined
    pad1.SetGridx()         # Vertical grid
    pad1.SetGridy()         # Vertical grid
    pad1.Draw()
    pad1.cd()		    # pad1 becomes the current pad
    hLO.SetStats(0)	    # No statistics on upper plot
    hLO.Draw("ehist")
    hNLO.Draw("ehistsame")
    leg = r.TLegend(0.75,0.75,0.9,0.85)
    leg.AddEntry(hLO,"Data","lp")
    leg.AddEntry(hNLO,"Bkgd","lp")
    leg.Draw("same")
    SetOwnership( leg, 0 ) # 0 = release (not keep), 1 = keep

    c.cd()		    # Go back to the main canvas before defining pad2
    pad2 = r.TPad("pad2", "pad2", 0., 0, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()

    # Define the ratio plot
    ratio = ratio_h.DrawCopy("ehist_"+name)
#    ratio.GetXaxis().SetRangeUser(0,300)
    ratio.GetYaxis().SetRangeUser(0.6,1.4)
#    ratio = hNLO.Clone(name+"_clone")
    ratio.SetLineColor(r.kBlack)
    ratio.SetLineWidth(2)
#    ratio.SetMinimum(0.8)
#    ratio.SetMaximum(1.35)
#    ratio.Sumw2()
    ratio.SetStats(0);      # No statistics on lower plot
#    ratio.Divide(hLO)
#    ratio.SetMarkerStyle(21)

    ratio.SetTitle("")      # Remove the ratio title
    ratio.GetYaxis().SetTitle("Data/Bkgd")
    ratio.GetYaxis().SetTitleSize(25)
    ratio.GetYaxis().SetTitleFont(43)
    ratio.GetYaxis().SetTitleOffset(1.)
    ratio.GetYaxis().SetLabelFont(43)
    ratio.GetYaxis().SetLabelSize(20)
    ratio.GetXaxis().SetTitle("HT")
    ratio.GetXaxis().SetTitleSize(25)
    ratio.GetXaxis().SetTitleFont(43)
    ratio.GetXaxis().SetTitleOffset(2.5)
    ratio.GetXaxis().SetLabelFont(43)
    ratio.GetXaxis().SetLabelSize(20)

#    fitf_h = r.TF1(name+"_hf", "[0]", thred, 800)
#    ratio.Fit(name+"_hf", "+R")
#    func_h = ratio.GetFunction(name+"_hf")
#    print("Chi square of straight line is: {}".format(func_h.GetChisquare()))
    fitf = r.TF1(name+"_f", "exp([0]+[1]*x)", 160, thred)
    ratio.Fit(name+"_f", "R")
    func = ratio.GetFunction(name+"_f")
    print("Chi square of expo is: {}".format(func.GetChisquare()))
    if '3b' in name:
	fitf_l1 = r.TF1(name+"_f1", "pol3", 380, 800)
	ratio.Fit(name+"_f1", "+R")
	func1 = ratio.GetFunction(name+"_f1")
	print("Chi square of straight line is: {}".format(func1.GetChisquare()))
    elif '4b' in name:
	fitf_l1 = r.TF1(name+"_f1", "pol3", 380, 800)
	ratio.Fit(name+"_f1", "+R")
	func1 = ratio.GetFunction(name+"_f1")
	print("Chi square of straight line is: {}".format(func1.GetChisquare()))
    ratio.Draw("E0")
    
    r.gPad.Update()
    line = r.TLine(r.gPad.GetUxmin(), 1, r.gPad.GetUxmax(), 1)
    line.Draw("same")
    SetOwnership( line, 0 )
    
#    return c,func,func_h
    return c

def computeTopPtSF(channels, inf, outdir, outf):
    bkgd_e = inf.Get(channels[0]+"/total")
    bkgd_mu = inf.Get(channels[1]+"/total")
    data_e = inf.Get(channels[0]+"/data_obs")
    data_mu = inf.Get(channels[1]+"/data_obs")
    if not (bkgd_e and bkgd_mu and data_e and data_mu):
        print("Cannot find total/data histogtams in {} or {}, please check the file ".format(channels[0], channels[1]) + fpath)
        exit(-1)
    bkgd = bkgd_e.Clone("bkdg_tota"+channels[0][-2:])
    bkgd.Add(bkgd_mu)
#    applyWeights(bkgd, channels[0][-2:])
    data = data_e.Clone("data_total"+channels[0][-2:])
    data.Add(data_mu)
    weights = data.Clone("weights"+channels[0][-2:])
    weights.Divide(bkgd)
#    c,f,f1 = ratioPlot(data, bkgd, weights, "topPtSF_"+channels[0][-2:])
    c = ratioPlot(data, bkgd, weights, "topPtSF_"+channels[0][-2:])
    c.SaveAs(os.path.join(outdir,"topPtSF_"+channels[0][-2:]+".pdf"))
    outf.cd()
    c.Write()
#    f.Write()
#    f1.Write()
    weights.Write()


try:
     # retrive command line options
    shortopts  = "s:e:j:d:o:c:l:p:t:g:r:?" #RJ
    opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
    # print help information and exit:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)

CMSSW_BASE=os.environ.get('CMSSW_BASE')
if(len(CMSSW_BASE)==0):
    usage()
    sys.exit(1)

inputdir=''
outdir=''

for o,a in opts:
    if o in('-d'): inputdir = a
    elif o in('-o'): outdir = a

fpath = os.path.join(inputdir, "haa4b_60_13TeV_wh.root")
inf = r.TFile(fpath, "READ")
if not inf.IsOpen():
    print("The file cannot be opened, please check: " + fpath)
    exit(-1)
outf = r.TFile(os.path.join(outdir,"topSF.root"), 'RECREATE')

######## 3b bin ##########
computeTopPtSF(hists_3b, inf, outdir, outf)

######## 4b bin ##########
computeTopPtSF(hists_4b, inf, outdir, outf)


inf.Close()
outf.Close()

raw_input("Press [ENTER] to exit...")
