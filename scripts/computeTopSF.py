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


bkgds = ["otherbkg", "ttbarbba", "ttbarcba", "ttbarlig", "zll", "wlnu", "ddqcd", "qcd"]
#startbin = 9 

low=-1
high=-1
a_3b = 0
b_3b = 0
a_4b = 0
b_4b = 0

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
    global a_3b
    global a_4b
    global b_3b
    global b_4b
    thred_low = 0
    thred_high = 500
    if vh_tag=="zh":
        thred_low = int(low)
	thred_high = int(high)
    for i in range(1, hist.GetNbinsX()+1):
        sf=1.0
	ht = hist.GetBinCenter(i)
	if ht>=thred_low and ht<=thred_high:
	    if '3b' in name:
		#sf = math.exp(0.04182-0.00095*ht) # Wh 2016
		#sf = math.exp(-0.07758-0.00175*ht) # Wh 2017
		sf = math.exp(0.07997+0.00001*ht) # Wh 2018
		
		#sf = math.exp(a_3b+b_3b*ht)
	    elif '4b' in name:
		sf = math.exp(0.00503-0.00036*ht) # Wh 2018
		#sf = math.exp(-0.07055-0.00171*ht) # Wh 2017
		#sf = math.exp(0.02629-0.00098*ht) # Wh 2016
		
		#sf = math.exp(a_4b+b_4b*ht)
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

def ratioOnly(ratio_h, name):
    c = r.TCanvas(name,name,800,600)
    ratio = ratio_h.DrawCopy(name)
    ratio.SetTitle(name.replace("_"," "))
    pad = r.TPad("pad","pad", 0, 0, 1, 1)
    pad.SetBottomMargin(0.1)
    pad.Draw()
    pad.cd()
    fitf = r.TF1(name+"_f", "exp([0]+[1]*x)", 0, 500)
    ratio.Fit(name+"_f", "R")
    ratio.GetYaxis().SetRangeUser(0.4, 1.5)
    ratio.GetXaxis().SetTitleOffset(1.0)
    ratio.Draw("E0")
    return c

def ratioPlot(hLO,hNLO,ratio_h,name):
    global CC
    global a_3b
    global a_4b
    global b_3b
    global b_4b

    thred = 500
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
#    pad1.SetGridx()         # Vertical grid
#    pad1.SetGridy()         # Vertical grid
    pad1.Draw()
    pad1.cd()		    # pad1 becomes the current pad
    hLO.SetStats(0)	    # No statistics on upper plot
    hLO.Draw("ehist")
    hNLO.Draw("ehistsame")
    leg = r.TLegend(0.65,0.75,0.9,0.85)
    leg.AddEntry(hLO,"Data - nonTop","lp")
    leg.AddEntry(hNLO,"Top","lp")
    leg.Draw("same")
    SetOwnership( leg, 0 ) # 0 = release (not keep), 1 = keep

    c.cd()		    # Go back to the main canvas before defining pad2
    pad2 = r.TPad("pad2", "pad2", 0., 0, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.2)
#    pad2.SetGridx()
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()

    # Define the ratio plot
    ratio = ratio_h.DrawCopy("ehist_"+name)
#    ratio.GetXaxis().SetRangeUser(0,300)
    ratio.GetYaxis().SetRangeUser(0.4,1.6)
#    ratio.GetYaxis().SetRangeUser(0.9*ratio.GetMinimum(),2)
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
    ratio.GetYaxis().SetTitle("(Data-nonTop)/Top")
    ratio.GetYaxis().SetTitleSize(20)
    ratio.GetYaxis().SetTitleFont(43)
    ratio.GetYaxis().SetTitleOffset(1.5)
    ratio.GetYaxis().SetLabelFont(43)
    ratio.GetYaxis().SetLabelSize(20)
    ratio.GetXaxis().SetTitle("Pt(W)")
    ratio.GetXaxis().SetTitleSize(25)
    ratio.GetXaxis().SetTitleFont(43)
    ratio.GetXaxis().SetTitleOffset(2.5)
    ratio.GetXaxis().SetLabelFont(43)
    ratio.GetXaxis().SetLabelSize(20)

    c1 = ratioOnly(ratio, name+"_ratio")
#    pad2.cd()
#    fitf = r.TF1(name+"_f", "exp([0]+[1]*x)", int(low), int(high))
#    ratio.Draw("E0")
#    if CC=="False":
#      print("name of function: {}".format(name+"_f"))
#      ratio.Fit(name+"_f", "R")
#      func = ratio.GetFunction(name+"_f")
#      print("Chi square of expo is: {}".format(func.GetChisquare()))
#      if showfull=="True":
#        func.SetRange(0, 500)
#      func.Draw("same")
#      a = func.GetParameter(0)
#      b = func.GetParameter(1)
#      if "3b" in name:
#        a_3b = a
#	b_3b = b
#      elif "4b" in name:
#        a_4b = a
#	b_4b = b

#    r.gPad.Update()
#    line = r.TLine(r.gPad.GetUxmin(), 1, r.gPad.GetUxmax(), 1)
#    line.Draw("same")
#    SetOwnership( line, 0 )
   
#    return c,func
    return c1

def getHists(channels, inf):
    for ch in channels:
	data = inf.Get(ch+"/data_obs")
	if not data:
	    print("Error, cannot find data in the channel: {}".format(ch))
	    exit(-1)
	dataSub = data.Clone("dataSub_h")
	for bkgd in bkgds:
	    if("ttbar" in bkgd): # add in ttbarMC
		h_mc = inf.Get(ch+"/"+bkgd)
		if h_mc:
		    if("ttbarMC" in vars() or "ttbarMC" in globals()): ttbarMC.Add(h_mc)
		    else: ttbarMC = h_mc.Clone("ttbarMC_h")
	    else:   # add in onttbarMC
		h_mc = inf.Get(ch+"/"+bkgd)
		if h_mc:
		    if("nonttbarMC" in vars() or "nonttbarMC" in globals()): nonttbarMC.Add(h_mc)
		    else: nonttbarMC = h_mc.Clone("nonttbarMC_h")
	if not ("ttbarMC" in vars() or "ttbarMC" in globals()):
	    print("Error, cannot find ttbar process in the channel: {}".format(ch))
	    exit(-1)
	if("nonttbarMC" in vars() or "nonttbarMC" in globals()): dataSub.Add(nonttbarMC,-1)
	if("dataSub_ret" in vars() or "dataSub_ret" in globals()): dataSub_ret.Add(dataSub)
	else: dataSub_ret = dataSub.Clone("dataSub_ret_h_"+ch[-2:])
	if("ttbarMC_ret" in vars() or "ttbarMC_ret" in globals()): ttbarMC_ret.Add(ttbarMC)
	else: ttbarMC_ret = ttbarMC.Clone("ttbarMC_ret_h_"+ch[-2:])
	del ttbarMC,dataSub
	if("nonttbarMC" in vars() or "nonttbarMC" in globals()): del nonttbarMC

    return dataSub_ret, ttbarMC_ret

def computeTopPtSF(channels, inf, outdir, outf):
    postfix=""
    if CC=="True": postfix="_CC"
    dataSub_h, ttbarMC_h = getHists(channels, inf)
    if CC=="True":
        applyWeights(ttbarMC_h, "weights"+channels[0][-2:])
    weights = dataSub_h.Clone("weights"+channels[0][-2:]+postfix)
    weights.Divide(ttbarMC_h)
#    c,f = ratioPlot(data, bkgd, weights, "topPtSF_"+channels[0][-2:])
    c = ratioPlot(dataSub_h, ttbarMC_h, weights, "topPtSF_"+channels[0][-2:]+postfix)
    c.SaveAs(os.path.join(outdir,"topPtSF_"+channels[0][-2:]+".pdf")+postfix)
    outf.cd()
    c.Write()
    dataSub_h.Write()
    ttbarMC_h.Write()
    weights.Write()


try:
     # retrive command line options
    shortopts  = "d:o:l:f:u:t:c:?" #RJ
    opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
    # print help information and exit:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
#    usage()
    sys.exit(1)

CMSSW_BASE=os.environ.get('CMSSW_BASE')
if(len(CMSSW_BASE)==0):
    usage()
    sys.exit(1)

inputdir=''
outdir=''
vh_tag=''
CC="False" # do the closure check or not
showfull="True"

for o,a in opts:
    if o in('-d'): inputdir = a
    elif o in('-o'): outdir = a
    elif o in('-t'): vh_tag = a
    elif o in('-l'): low = a
    elif o in('-u'): high = a
    elif o in('-c'): CC = a
    elif o in('-f'): showfull = a

print("vh_tag = {}".format(vh_tag))
print("CC = {}".format(CC))

if vh_tag == "wh":
    hists_3b = ['e_A_CR_3b','mu_A_CR_3b']
    hists_4b = ['e_A_CR_4b','mu_A_CR_4b']
elif vh_tag == "zh":
    hists_3b = ['emu_A_CR_3b']
    hists_4b = ['emu_A_CR_4b']
else:
    print("Error, cannot recognize the tag: {}".format(vh_tag))
    print(vh_tag)
    exit(-1)

fpath = os.path.join(inputdir, "haa4b_60_13TeV_{}.root".format(vh_tag))
inf = r.TFile(fpath, "READ")
if not inf.IsOpen():
    print("The file cannot be opened, please check: " + fpath)
    exit(-1)

outname = "topSF.root"
if CC=="True": outname = "topSF_CC.root"
outf = r.TFile(os.path.join(outdir,outname), 'RECREATE')
#outf_CC = r.TFile(os.path.join(outdir,"topSF_CC.root"), 'RECREATE')
######## 3b bin ##########
computeTopPtSF(hists_3b, inf, outdir, outf)

######## 4b bin ##########
computeTopPtSF(hists_4b, inf, outdir, outf)

#CC="True"
#computeTopPtSF(hists_3b, inf, outdir, outf)
#computeTopPtSF(hists_4b, inf, outdir, outf)

inf.Close()
outf.Close()
#outf_CC.Close()

raw_input("Press [ENTER] to exit...")
