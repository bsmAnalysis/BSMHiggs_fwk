import ROOT as rt
import CMS_lumi, tdrstyle
from array import array
import os
import sys

from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double


# Choose if you need to blind SRs
#iblind=-1
iblind=3

limit_dir="/afs/cern.ch/work/y/yuanc/Analysis/H2a4b/CMSSW_10_2_13/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/cards_SB13TeV_SM_Wh_2017/0060/"
#limit_dir="/afs/cern.ch/work/y/yuanc/Analysis/H2a4b/CMSSW_10_2_13/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/cards_SB13TeV_SM_Wh_2018/0060/"


# Picks up the pre-fit BDT:
dir1="shapes_prefit"
#Picks up the post-fit BDT:
#dir1="shapes_fit_b"

channels = ["mumu_A_CR_3b", "mumu_A_CR_4b", "mumu_A_SR_3b", "mumu_A_SR_4b", "emu_A_CR_3b", "emu_A_CR_4b", "emu_A_SR_3b", "emu_A_SR_4b", "ee_A_CR_3b", "ee_A_CR_4b", "ee_A_SR_3b", "ee_A_SR_4b"]
#channels = ["mumu_A_CR_3b", "mumu_A_CR_4b","emu_A_CR_3b", "emu_A_CR_4b","ee_A_CR_3b", "ee_A_CR_4b"]
#channels = ["mu_A_CR5j_3b", "mu_A_CR5j_4b", "mu_A_CR_3b", "mu_A_CR_4b", "mu_A_SR_3b", "mu_A_SR_4b", "e_A_CR5j_3b", "e_A_CR5j_4b", "e_A_CR_3b", "e_A_CR_4b", "e_A_SR_3b", "e_A_SR_4b"]
#channels = ["mu_A_SR_4b"]
e_mu = ["e", "mu"]
#e_mu = ["e"]
wz = "zh"
#wz = "wh"
CRs = ["mumu_A_CR_3b", "mumu_A_CR_4b","emu_A_CR_3b", "emu_A_CR_4b", "emu_A_SR_3b", "emu_A_SR_4b", "ee_A_CR_3b", "ee_A_CR_4b", "mu_A_CR5j_3b", "mu_A_CR5j_4b", "mu_A_CR_3b", "mu_A_CR_4b", "e_A_CR5j_3b", "e_A_CR5j_4b", "e_A_CR_3b", "e_A_CR_4b"]


#set the tdr style
tdrstyle.setTDRStyle()

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

def printEvtYields(hist, name):
  if not hist: return
  printout = [name+":"]
  for ibin in range(1, hist.GetXaxis().GetNbins()+1):
#    content = hist.GetBinContent(ibin)
#    error = hist.GetBinError(ibin)
    error = Double()
    content = hist.IntegralAndError(ibin, ibin, error)
    yields = "{:.2f}".format(content) +"+-" + "{:.2f}".format(error)
    printout.append(yields)
#  error = Double()
#  content = hist.IntegralAndError(1, hist.GetXaxis().GetNbins(), error)
#  yields = "{:.2f}".format(content) +"+-" + "{:.2f}".format(error)
#  printout.append(yields)
  print("  ".join(printout))

def convertXRange(hist, name, edges):
  out_h = rt.TH1F(name, name, len(edges)-1, array('d', edges))
  for i in range(1, hist.GetNbinsX()+1):
    out_h.SetBinContent(i, hist.GetBinContent(i))
    out_h.SetBinError(i, hist.GetBinError(i))
  return out_h
   

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
L = 0.12*W_ref
R = 0.04*W_ref
for ch in e_mu:
#  file = rt.TFile(limit_dir+"0060/fitDiagnostics_{}.root".format(ch),"READ")
  file = rt.TFile(limit_dir+"fitDiagnostics_{}.root".format(ch),"READ")
#  file = rt.TFile(limit_dir+"fitDiagnostics.root","READ")
  inf = rt.TFile(limit_dir+"haa4b_60_13TeV_{}.root".format(wz),"READ")
  for dir2 in channels:
    blind=iblind
    if dir2 in CRs: blind=-1
    print("blind: " +str(blind))
    dir=dir1+"/"+dir2+"/"
    
    inh = inf.Get(dir2+"/data_obs")
    if not inh:
      print("!!!Warning: cannot find the data_obs in channel {}, skip the channel!".format(dir2))
      continue
    edges = []
    np = inh.GetNbinsX()
    for i in range(1, np+1):
      edges.append(inh.GetXaxis().GetBinLowEdge(i))
    edges.append(inh.GetXaxis().GetBinUpEdge(i))
    print(edges)
    
    canvas = rt.TCanvas("c2","c2",50,50,W,H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(0)
    canvas.SetTicky(0)
    
    h =  rt.TH1F("h","h; BDT; Events",5, 0, 5)
    
    xAxis = h.GetXaxis()
    xAxis.SetNdivisions(6,5,0)
    
    yAxis = h.GetYaxis()
    yAxis.SetNdivisions(6,5,0)
    yAxis.SetTitleOffset(1)
    
    
    bkgd_list = []
    lgname_list = []
    
    data = file.Get(dir+"data")
    if not data: continue
    data.SetMarkerStyle(20)
    data.SetMarkerColor(1)
    
#    total = file.Get(dir+"total_background")
    total = convertXRange(file.Get(dir+"total_background"), "total_background", edges)
    
    sig = file.Get(dir+"wh")
    if sig:
      sig = convertXRange(sig, "wh", edges)
      sig.SetLineColor(rt.kRed+4)
      sig.SetLineWidth(2)
      sig.SetLineStyle(2)
    #sig.Scale(50)
    
    otherbkg = file.Get(dir+"otherbkg")
    if otherbkg:
      otherbkg = convertXRange(otherbkg, "otherbkg", edges)
      otherbkg.SetFillColor(424)
      otherbkg.SetLineColor(1)
      bkgd_list.append(otherbkg)
      lgname_list.append("Other bkgs")
    
    qcd = file.Get(dir+"ddqcd")
    if qcd:
      qcd = convertXRange(qcd, "ddqcd", edges)
      qcd.SetFillColor(634)
      qcd.SetLineColor(1)
      bkgd_list.append(qcd)
      lgname_list.append("DD qcd")
    
    zll = file.Get(dir+"zll")
    if zll:
      zll = convertXRange(zll, "zll", edges)
      zll.SetFillColor(624)
      zll.SetLineColor(1)
      bkgd_list.append(zll)
      lgname_list.append("Z#rightarrow ll")
    
    wlnu = file.Get(dir+"wlnu")
    if wlnu:
      wlnu = convertXRange(wlnu, "wlnu", edges)
      wlnu.SetFillColor(622)
      wlnu.SetLineColor(1)
      bkgd_list.append(wlnu)
      lgname_list.append("W#rightarrow l#nu")
    
    ttbarbba = file.Get(dir+"ttbarbba")
    if ttbarbba:
      ttbarbba = convertXRange(ttbarbba, "ttbarbba", edges)
      ttbarbba.SetFillColor(833)
      ttbarbba.SetLineColor(1)
      bkgd_list.append(ttbarbba)
      lgname_list.append("t#bar{t} + b#bar{b}")
    
    ttbarcba = file.Get(dir+"ttbarcba")
    if ttbarcba:
      ttbarcba = convertXRange(ttbarcba, "ttbarcba", edges)
      ttbarcba.SetFillColor(408)
      ttbarcba.SetLineColor(1)
      bkgd_list.append(ttbarcba)
      lgname_list.append("t#bar{t} + c#bar{c}")
    
    ttbarlig = file.Get(dir+"ttbarlig")
    if ttbarlig:
      ttbarlig = convertXRange(ttbarlig, "ttbarlig", edges)
      ttbarlig.SetFillColor(406)
      ttbarlig.SetLineColor(1)
      bkgd_list.append(ttbarlig)
      lgname_list.append("t#bar{t} + light")
    
    MC = rt.THStack()
    
    #MC.Add(zll)
    #MC.Add(wlnu)
    
    #MC.Add(qcd)
    #MC.Add(otherbkg)
    
    #MC.Add(ttbarbba)
    #MC.Add(ttbarcba)
    #MC.Add(ttbarlig)
    for bkgd in bkgd_list:
      MC.Add(bkgd)
    
    
    t1 = rt.TPad("t1","t1", 0.0, 0.2, 1.0, 1.0)
    t1.SetFillColor(0)
    t1.SetBorderMode(0)
    t1.SetBorderSize(2)
    t1.SetTickx(1)
    t1.SetTicky(1)
    t1.SetLeftMargin(0.10)
    t1.SetRightMargin(0.05)
    t1.SetTopMargin(0.05)
    t1.SetBottomMargin(0.10)
    t1.SetFrameFillStyle(0)
    t1.SetFrameBorderMode(0)
    t1.SetFrameFillStyle(0)
    t1.SetFrameBorderMode(0)
    
    t1.Draw()
    t1.cd()
    
    #MC.Draw("hist")
    hdata = rt.TH1F("hdata", "data bdt", np, array('d',edges))
    
    #xmin=-0.3
    #for i in range(0,hdata.GetXaxis().GetNbins()):
    #    label=str(xmin+i*hdata.GetXaxis().GetBinWidth(i))
    #    hdata.GetXaxis().SetBinLabel(i,label)
        
    #hdata.Rebin(rbin)
    
    htotal=total.Clone()
    
    px, py = Double(), Double()
    pyerr = Double()
    
    nPoints=data.GetN()
    for i in range(0,nPoints):
        data.GetPoint(i,px,py)
        pyerr=data.GetErrorY(i)
        #hdata.Fill(px,py)
        hdata.SetBinContent(i+1,py)
        hdata.SetBinError(i+1,pyerr)
    
    
    
    hdata.GetXaxis().SetLabelOffset(0.007);
    hdata.GetXaxis().SetLabelSize(0.04);
    hdata.GetXaxis().SetTitleOffset(1.2);
    hdata.GetXaxis().SetTitleFont(42);
    hdata.GetXaxis().SetTitleSize(0.04);
    hdata.GetYaxis().SetLabelFont(42);
    hdata.GetYaxis().SetLabelOffset(0.007);
    hdata.GetYaxis().SetLabelSize(0.04);
    hdata.GetYaxis().SetTitleOffset(1.35);
    hdata.GetYaxis().SetTitleFont(42);
    hdata.GetYaxis().SetTitleSize(0.04);
    
    hdata.GetYaxis().SetTitle("Events");
    hdata.GetXaxis().SetTitle("BDT");
    hdata.Draw();
    
    t1.Update();
    
        
    MC.Draw("histsame")
    hdata.Draw("epsame0")
    
    if sig: 
      sig.Draw("histsame")
    
#    hdata.SetMaximum(1.5*max(hdata.GetMaximum(), MC.GetMaximum()))
    hdata.SetMinimum(0.1)
    hmax = max(hdata.GetMaximum(), MC.GetMaximum())
    if hmax > 50000:
	hdata.SetMaximum(1000*hmax)
    elif hmax > 10000:
	hdata.SetMaximum(500*hmax)
    elif hmax > 1000:
        hdata.SetMaximum(100*hmax)
    else:
	hdata.SetMaximum(50*hmax)
    
    if (blind>0):
#        for i in range(hdata.FindBin(blind),hdata.GetNbinsX()+1):
        for i in range(blind,hdata.GetNbinsX()+1):
            hdata.SetBinContent(i,0)
            hdata.SetBinError(i,0)
#        blinding_box = rt.TPave(hdata.GetBinLowEdge(hdata.FindBin(blind)),  hdata.GetMinimum(), hdata.GetXaxis().GetXmax(), hdata.GetMaximum(), 0, "NB" )  
        blinding_box = rt.TPave(hdata.GetBinLowEdge(blind),  hdata.GetMinimum(), hdata.GetXaxis().GetXmax(), hdata.GetMaximum(), 0, "NB" )  
        blinding_box.SetFillColor(15)
        blinding_box.SetFillStyle(3013)
        blinding_box.Draw("same F");
    
    #set the colors and size for the legend
    histLineColor = rt.kOrange+7
    histFillColor = rt.kOrange-2
    markerSize  = 1.0
    
    latex = rt.TLatex()
    n_ = 2
    
    x1_l = 0.95 #0.92
    y1_l = 0.90 #0.60
    
    dx_l = 0.30
    dy_l = 0.18 #0.18
    x0_l = x1_l-dx_l
    y0_l = y1_l-dy_l
    
    legend =  rt.TLegend(0.40,0.74,0.93,0.96, "NDC")
    #legend.SetFillColor( rt.kGray )
    #legend.SetHeader(dir)
    legend.SetNColumns(3)  
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)
    legend.SetLineColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    
    #legend.Draw("same")
    #legend.cd()
    
    ar_l = dy_l/dx_l
    #gap_ = 0.09/ar_l
    gap_ = 1./(n_+1)
    bwx_ = 0.12
    bwy_ = gap_/1.5
    
    x_l = [1.2*bwx_]
    #y_l = [1-(1-0.10)/ar_l]
    y_l = [1-gap_]
    ex_l = [0]
    ey_l = [0.04/ar_l]
    
    ## latex.DrawLatex(xx_+1.*bwx_,yy_,"Data")
    legend.AddEntry(hdata,"Data","LP")
    if sig: legend.AddEntry(sig,"Zh (60)","L")
    print("-"*100)
    print(dir2) 
    for i in range(len(bkgd_list)):
      legend.AddEntry(bkgd_list[i], lgname_list[i], "LF")
      printEvtYields(bkgd_list[i], lgname_list[i])
    printEvtYields(total, "total")
    print("\n")
    printEvtYields(hdata, "data")
    if sig: printEvtYields(sig, "signal")
    
    #legend.AddEntry(wlnu,"W#rightarrow l#nu","LF")
    #legend.AddEntry(zll,"Z#rightarrow ll","LF")
    #legend.AddEntry(ttbarbba,"t#bar{t} + b#bar{b}","LF")
    #legend.AddEntry(ttbarcba,"t#bar{t} + c#bar{c}","LF")
    #legend.AddEntry(ttbarlig,"t#bar{t} + light","LF")
    #legend.AddEntry(otherbkg,"Other bkgs","LF")
    #legend.AddEntry(qcd,"DD qcd","LF")
    
    if(blind>0):
        legend.AddEntry(blinding_box,"Blinded area","F")
    
    legend.Draw("same")
    t1.SetLogy(True)
    
    #draw the lumi text on the canvas
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
    
    
    
    ## ratio plot
    t2 = rt.TPad("t2", "t2",0.0,0.0, 1.0,0.2)
    t2.SetFillColor(0)
    t2.SetBorderMode(0)
    t2.SetBorderSize(2)
    t2.SetGridy()
    t2.SetTickx(1)
    t2.SetTicky(1)
    t2.SetLeftMargin(0.10)
    t2.SetRightMargin(0.05)
    t2.SetTopMargin(0.0)
    t2.SetBottomMargin(0.20)
    t2.SetFrameFillStyle(0)
    t2.SetFrameBorderMode(0)
    t2.SetFrameFillStyle(0)
    t2.SetFrameBorderMode(0)
    t2.Draw()
    t2.cd()
    t2.SetGridy(1)
    t2.SetPad(0,0.0,1.0,0.2)
    
    
    
    hratio=hdata.Clone("myratio")
    hratio.Divide(htotal)
    
    hratio.Draw("same e2")
    #hratio.DrawCopy("histsame"); 
    hratio.SetFillColor(rt.kBlue);
    hratio.SetFillStyle(3018);
    hratio.Draw("same 2");
     
    
    
    #hratio.SetFillColor(rt.kMagenta)
    
    #hratio.Draw("e2")
    #hratio.Draw("hist L same"); # you missed the option L
    #hratio.Draw("same 2 0")
    
    
    hratio.SetMinimum(0.6);
    hratio.SetMaximum(1.4);
    
    hratio.SetLineColor(1)
    hratio.SetFillStyle(3005)
    hratio.SetFillColor(rt.kGray+3)
    hratio.SetMarkerColor(1)
    hratio.SetMarkerStyle(20)
    hratio.GetYaxis().SetLabelSize(0.13)
    hratio.GetYaxis().SetTitleSize(0.13)
    hratio.GetYaxis().SetTitleOffset(1.)
    hratio.GetXaxis().SetLabelSize(0.13)
    hratio.GetXaxis().SetTitleSize(0.13)
    hratio.GetXaxis().SetTitleOffset(1.2)
    hratio.SetMarkerSize(0.5)
    #hratio.GetYaxis().SetTitle("ratio")
    hratio.SetLineWidth(2)
    
    if(blind>0):
#        blinding_box2 = rt.TPave(hdata.GetBinLowEdge(hdata.FindBin(blind)), 0.4, hdata.GetXaxis().GetXmax(), 1.6, 0, "NB" )  
        blinding_box2 = rt.TPave(hdata.GetBinLowEdge(blind), 0.4, hdata.GetXaxis().GetXmax(), 1.6, 0, "NB" )  
        blinding_box2.SetFillColor(15)
        blinding_box2.SetFillStyle(3013)
        blinding_box2.Draw("same F");
    
    
    line = rt.TLine(0.,1,hdata.GetXaxis().GetXmax(),1);
    line.SetLineColor(rt.kBlack)
    line.Draw("same")
    
    t2.Update()
    
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()
    #frame = canvas.GetFrame()
    #frame.Draw()
    
    #canvas.SaveAs(limit_dir+dir1+"_"+dir2+".pdf")
    #canvas.SaveAs(dir2+"_postfit.pdf")
    canvas.SaveAs("{}_{}.pdf".format(dir1,dir2))
    
    canvas.Update()

raw_input("Press Enter to end")
