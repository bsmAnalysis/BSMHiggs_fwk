import ROOT as rt
import CMS_lumi, tdrstyle
from array import array
import os
import sys
import ctypes

import argparse

import math

from ROOT import gROOT, gBenchmark, gRandom, gSystem, gStyle, Double

import ROOT

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()

parser.add_argument( "--do_liny", dest='do_liny', help="Use linear scale for y axis", action='store_true' )
parser.add_argument( "--do_linzoom", dest='do_linzoom', help="Zoom in on the highest bins with a linear vertical scale", action='store_true' )
parser.add_argument( "prodmode", type=str, help="Production mode (wh or zh)")
parser.add_argument( "limit_dir", type=str, help="Input directory for the mass point.  Example = all-batch-output-10x-july24a--new-bdt-binning2--autoMCStats/cards_SB13TeV_SM_Wh_2016_noSoftb/0060/")

args = parser.parse_args()

limit_dir = args.limit_dir

# Choose if you need to blind SRs
iblind=-1
#iblind=4


# Picks up the pre-fit BDT:
#dir1="shapes_prefit"

#Picks up the b post-fit BDT:
dir1="shapes_fit_s"

#Picks up the s+b post-fit BDT:
#dir1="shapes_fit_s"

#wz = "zh"
#wz = "wh"

wz = args.prodmode
if ( not (wz == "wh" or wz == "zh" )):
   print("\n\n *** set prodmode to wh or zh.  ", wz, " is not allowed.\n\n")
   quit()

outdir = "postfit-plots"
if ( "fit_s" in dir1 ): outdir = "postfit_s+b-plots"

if ( wz == "zh" ) :
   outdir += "-zh"
   if ( "2016" in limit_dir ): outdir += "-2016"
   if ( "2017" in limit_dir ): outdir += "-2017"
   if ( "2018" in limit_dir ): outdir += "-2018"
if ( wz == "wh" ) :
   outdir += "-wh"
   if ( "2016" in limit_dir ): outdir += "-2016"
   if ( "2017" in limit_dir ): outdir += "-2017"
   if ( "2018" in limit_dir ): outdir += "-2018"

try:
        os.mkdir( outdir )
except:
        print("\n\n problem making %s" % outdir ) ;


if wz=="wh":
  #channels = ["mu_A_CR5j_3b", "mu_A_CR5j_4b", "mu_A_CR_3b", "mu_A_CR_4b", "mu_A_SR_3b", "mu_A_SR_4b", "e_A_CR5j_3b", "e_A_CR5j_4b", "e_A_CR_3b", "e_A_CR_4b", "e_A_SR_3b", "e_A_SR_4b"]
  channels = ["mu_A_CR_3b", "mu_A_CR_4b", "mu_A_SR_3b", "mu_A_SR_4b", "e_A_CR_3b", "e_A_CR_4b", "e_A_SR_3b", "e_A_SR_4b"]
  e_mu = [""]
elif wz=="zh":
  channels = ["mumu_A_CR_3b", "mumu_A_CR_4b", "mumu_A_SR_3b", "mumu_A_SR_4b", "ee_A_CR_3b", "ee_A_CR_4b", "ee_A_SR_3b", "ee_A_SR_4b"]
  e_mu = ["e", "mu"]
#CRs = ["mumu_A_CR_3b", "mumu_A_CR_4b","emu_A_CR_3b", "emu_A_CR_4b", "emu_A_SR_3b", "emu_A_SR_4b", "ee_A_CR_3b", "ee_A_CR_4b", "mu_A_CR5j_3b", "mu_A_CR5j_4b", "mu_A_CR_3b", "mu_A_CR_4b", "e_A_CR5j_3b", "e_A_CR5j_4b", "e_A_CR_3b", "e_A_CR_4b"]
CRs = ["mumu_A_CR_3b", "mumu_A_CR_4b", "ee_A_CR_3b", "ee_A_CR_4b",  "mu_A_CR_3b", "mu_A_CR_4b", "e_A_CR_3b", "e_A_CR_4b"]

verbose = True

do_liny = False
if args.do_liny:
   do_liny = True
   print("Will use linear axis\n")

do_linzoom = False
if args.do_linzoom:
   do_linzoom = True
   do_liny = True
   print("Will zoom in on senitive bins with linear vertical scale.\n")



#set the tdr style
tdrstyle.setTDRStyle()

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "" #"13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPeriod = 9
if ( "2016" in limit_dir ): iPeriod = 6
if ( "2017" in limit_dir ): iPeriod = 7
if ( "2018" in limit_dir ): iPeriod = 8

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12

H_ref = 700; 
W_ref = 700; 
W = W_ref
H  = H_ref


gStyle.SetErrorX(0.5)  # owen: need to add this to get horizontal bars on data hist


#------------------------------------------------------------

def printEvtYields(hist, name):
  if not hist: return
#  printout = [name+":"]
  printout = ["{0: <19}:".format(name)]
  for ibin in range(1, hist.GetXaxis().GetNbins()+1):
#    content = hist.GetBinContent(ibin)
#    error = hist.GetBinError(ibin)
    error = Double() # owen
    #error = ctypes.c_double()
    content = hist.IntegralAndError(ibin, ibin, error)
    yields = "{0:10.2f}".format(content) +" +- " + "{0:<10.2f}".format(error)
    #yields = "{0:10.2f}".format(content) +" +- " + "{0:<10.2f}".format(error.value)
    #yields = " fixme "
    printout.append(yields)
#  error = Double()
#  content = hist.IntegralAndError(1, hist.GetXaxis().GetNbins(), error)
#  yields = "{:.2f}".format(content) +"+-" + "{:.2f}".format(error)
#  printout.append(yields)
  print("  ".join(printout))


#------------------------------------------------------------

def convertXRange(hist, name, edges):
  out_h = rt.TH1F(name, name, len(edges)-1, array('d', edges))
  for i in range(1, hist.GetNbinsX()+1):
    out_h.SetBinContent(i, hist.GetBinContent(i))
    out_h.SetBinError(i, hist.GetBinError(i))
  return out_h


#------------------------------------------------------------
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

#iPeriod = 0

# references for T, B, L, R
T = 0.08*H_ref
B = 0.12*H_ref 
L = 0.12*W_ref
R = 0.04*W_ref
for ch in e_mu:

  if ( verbose ) : print("\n\n verbose:  ch = ", ch , "\n")

  if wz=="zh":
     file = rt.TFile(limit_dir+"fitDiagnostics_zh_{}.root".format(ch),"READ")
  elif wz=="wh":
     file = rt.TFile(limit_dir+"fitDiagnostics_wh.root","READ")

  if ( verbose ) : print("\n\n verbose: file = ", file.GetName(),"\n" )

  inf = rt.TFile(limit_dir+"haa4b_60_13TeV_{}.root".format(wz),"READ")

  if ( verbose ) : print("\n\n verbose inf = ", inf.GetName(), "\n" )

  for dir2 in channels:

    if ( verbose ) : print("  verbose:  dir2 = ", dir2)

    blind=iblind
    if dir2 in CRs: blind=-1

    print("blind: " +str(blind))

    dir=dir1+"/"+dir2+"/"

    inh = inf.Get(dir2+"/data_obs")

    if ( verbose ) : print("  verbose:  Get arg : ", (dir2+"/data_obs") , " inh = ", inh.GetName(), "  ,  ", inh.GetTitle() )

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
    
    hdata = file.Get(dir+"hdata")


    if not hdata: continue

    if ( verbose ) : print("  verbose:  Get arg: ",(dir+"hdata")," data = ", hdata.GetName(), "  ,  ", hdata.GetTitle() )

    hdata.SetMarkerStyle(20)
    hdata.SetMarkerColor(1)

#    total = file.Get(dir+"total_background")
    tmp_hist = file.Get(dir+"total_background")
    if ( verbose ) : print("  verbose:  Get arg: ", (dir+"total_background"), " tmp_hist = ", tmp_hist.GetName(), "  ,  ", tmp_hist.GetTitle(), "\n" )


    total = convertXRange(file.Get(dir+"total_background"), "total_background", edges)

   #-- owen:  Take the signal from the prefit if this is a fit_b plot (shows example of signal distribution).

    if dir1 == "shapes_fit_b":
       print("\n\n This is shapes_fit_b.  Taking signal from shapes_prefit.\n\n")
       sig = file.Get( "shapes_prefit/" + dir2 + "/wh")
    else:
       sig = file.Get(dir+"wh")


    if sig:
      sig = convertXRange(sig, "wh", edges)
      sig.SetLineColor(rt.kRed+4)
      #sig.SetLineWidth(2)
      sig.SetLineWidth(4)
      sig.SetLineStyle(2)
#      if((dir1 == "shapes_fit_s") and (dir2 not in CRs)): sig.Scale(10)
#      if ( wz == "zh" ): sig.Scale(100)

    otherbkg = file.Get(dir+"otherbkg")
    if otherbkg:
      otherbkg = convertXRange(otherbkg, "otherbkg", edges)
      otherbkg.SetFillColor(852)
      otherbkg.SetLineColor(1)
      bkgd_list.append(otherbkg)
      lgname_list.append("Other bkgs")

    qcd = file.Get(dir+"ddqcd")
    if qcd:
      qcd = convertXRange(qcd, "ddqcd", edges)
      qcd.SetFillColor(634)
      qcd.SetLineColor(1)
      bkgd_list.append(qcd)
      lgname_list.append("QCD")

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

#    hdata = rt.TH1F("hdata", "data bdt", np, array('d',edges))
    
    htotal=total.Clone()
        
    MC.Draw("hist")
#    hdata.Draw("e0psame0")
    
    MC.GetXaxis().SetLabelOffset(0.007);
    MC.GetXaxis().SetLabelSize(0.04);
    MC.GetXaxis().SetTitleOffset(1.2);
    MC.GetXaxis().SetTitleFont(42);
    MC.GetXaxis().SetTitleSize(0.04);
    MC.GetYaxis().SetLabelFont(42);
    MC.GetYaxis().SetLabelOffset(0.007);
    MC.GetYaxis().SetLabelSize(0.04);
    MC.GetYaxis().SetTitleOffset(1.35);
    MC.GetYaxis().SetTitleFont(42);
    MC.GetYaxis().SetTitleSize(0.04);
    MC.GetYaxis().SetMaxDigits(3);

    MC.GetYaxis().SetTitle("Events");
    MC.GetXaxis().SetTitle("BDT");

    t1.Update()          

    if sig: 
      sig.Draw("histsame")
      if ( wz == "wh" ): sig.Scale(100)
      if ( wz == "zh" ):
         if ("3b" in dir2): sig.Scale(100)
         if ("4b" in dir2): sig.Scale(40)

#    hdata.SetMaximum(1.5*max(hdata.GetMaximum(), MC.GetMaximum()))
    MC.SetMinimum(0.1)
    hmax = max(hdata.GetMaximum(), MC.GetMaximum())
    if ( not do_liny ) :
       if hmax > 50000:
           MC.SetMaximum(1000*hmax)
       elif hmax > 10000:
           MC.SetMaximum(500*hmax)
       elif hmax > 1000:
           MC.SetMaximum(100*hmax)
       else:
           MC.SetMaximum(50*hmax) #50*hmax)
    else:
       if ( do_linzoom ):
          zoombin = hdata.GetBinContent(4)
          MC.SetMaximum(3*zoombin)
       else:
          MC.SetMaximum(1.5*hmax)
          if (wz == "wh"):
             if("4b" in dir2): MC.SetMaximum(1.8*hmax) 
          if (wz == "zh"):
             if("4b" in dir2): MC.SetMaximum(2.0*hmax)  

    systBand = rt.TGraphErrors(total.GetNbinsX()); IPoint=0  
    for ibin in range(1,htotal.GetXaxis().GetNbins()+1):  
       systBand.SetPoint(IPoint, htotal.GetBinCenter(ibin), htotal.GetBinContent(ibin) )
       systBand.SetPointError(IPoint, htotal.GetBinWidth(ibin)/2, htotal.GetBinError(ibin) ) 
       IPoint += 1 

    systBand.Set(IPoint);   

    systBand.SetFillStyle(3004); 
    systBand.SetFillColor(rt.kGray+2);  
    systBand.SetMarkerStyle(1);  
    systBand.Draw("same 2 0"); 

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
    
    hdata.Draw("e0psame0")   
    hdata.Draw("e1same")

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
    
    #######legend =  rt.TLegend(0.40,0.74,0.93,0.96, "NDC")
#    legend =  rt.TLegend(0.67,0.96,0.83,0.45, "NDC") 


    legend =  rt.TLegend(0.30,0.70,0.90,0.96, "brNDC")
    #legend.SetFillColor( rt.kGray )
    legend.SetHeader("") #dir)
    legend.SetNColumns(3)  
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
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
    legend.AddEntry(systBand, "Syst + Stat.", "F");  
    #if sig: legend.AddEntry(sig,"Zh (60)","L")
#    if sig: legend.AddEntry(sig,"Zh (60)","L")
    print("-"*100)
    print(dir2) 
    for i in range(len(bkgd_list)):
      legend.AddEntry(bkgd_list[i], lgname_list[i], "F")
      printEvtYields(bkgd_list[i], lgname_list[i])
    printEvtYields(total, "total")
    print("\n")
    printEvtYields(hdata, "data")
    if sig: 
       if ( wz == "zh" ): 
          if ("3b" in dir2): legend.AddEntry(sig,"Signal ZH (60GeV), (x100)","L") 
          if ("4b" in dir2): legend.AddEntry(sig,"Signal ZH (60GeV), (x40)","L")  
       if ( wz == "wh" ): legend.AddEntry(sig,"Signal WH (60GeV), (x100)","L")  
       printEvtYields(sig, "signal")
    
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
    if ( not do_liny ) : t1.SetLogy(True)
    
    #draw the lumi text on the canvas
    CMS_lumi.CMS_lumi(t1, iPeriod, iPos)
    
    t1.Update()

    canvas.cd()
    canvas.Update()
    
    
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
    
    
    hMCtotal = htotal.Clone("mcrelunc")  
#    utils::root::checkSumw2(hMCtotal) 

    # Make ratio with TGraphErrors
    gtotal = rt.TGraphErrors(hMCtotal.GetNbinsX())

    GPoint = 0  
    for i in range(1,hMCtotal.GetXaxis().GetNbins()+1): 
       gtotal.SetPoint(GPoint, hMCtotal.GetBinCenter(i), 1.0)
       if (hMCtotal.GetBinContent(i) != 0): 
          gtotal.SetPointError(GPoint, hMCtotal.GetBinWidth(i)/2, hMCtotal.GetBinError(i)/hMCtotal.GetBinContent(i))
          GPoint += 1  
       else:
          gtotal.SetPointError(GPoint, hMCtotal.GetBinWidth(i)/2, 0) 
          GPoint += 1 ; continue  

       err = hMCtotal.GetBinError(i)/hMCtotal.GetBinContent(i)  
       hMCtotal.SetBinContent(i,1)  
       hMCtotal.SetBinError(i,err)

    gtotal.Set(GPoint) 
    gtotal.SetLineColor(1) 
    gtotal.SetFillStyle(3001) #3004 
    gtotal.SetFillColor(16) #rt.kGray+3)  
    gtotal.SetMarkerColor(1) 
    gtotal.SetMarkerStyle(20)

    hMCtotal.Reset("ICE")
    hMCtotal.SetTitle("") 
    hMCtotal.SetStats(0)
    hMCtotal.Draw()  
    gtotal.Draw("2 0 same")

    yscale = (1.0-0.2)/(0.2)
    hMCtotal.GetYaxis().SetTitle("Data/#Sigma Bkg.") 
    hMCtotal.GetXaxis().SetTitle() # drop the title to gain space
    hMCtotal.SetMinimum(0.2) 
    hMCtotal.SetMaximum(1.8) 

    if "CR" in dir2: 
        hMCtotal.SetMinimum(0.4) 
        hMCtotal.SetMaximum(1.6) 

    hMCtotal.GetXaxis().SetLabelFont(42) 
    hMCtotal.GetXaxis().SetLabelOffset(0.007) 
    hMCtotal.GetXaxis().SetLabelSize(0.04 * yscale) 
    hMCtotal.GetXaxis().SetTitleFont(42) 
    hMCtotal.GetXaxis().SetTitleSize(0.035 * yscale) 
    hMCtotal.GetXaxis().SetTitleOffset(0.8) 
    hMCtotal.GetYaxis().SetLabelFont(42) 
    hMCtotal.GetYaxis().SetLabelOffset(0.007) 
    hMCtotal.GetYaxis().SetLabelSize(0.03 * yscale) 
    hMCtotal.GetYaxis().SetTitleFont(42) 
    hMCtotal.GetYaxis().SetTitleSize(0.035 * yscale) 
    hMCtotal.GetYaxis().SetTitleOffset(0.3)

    #add comparisons
    dataToObsH=hdata.Clone("myratio")  
    dataToObsH.Divide(htotal)  

#    dataToObs = rt.TGraphErrors(hdata.GetNbinsX()) 
#    RPoint = 0 
    for i in range(1,hdata.GetXaxis().GetNbins()+1): 
#       dataToObs.SetPoint(RPoint, hdata.GetBinCenter(i), hdata.GetBinContent(i)/htotal.GetBinContent(i)) 
       if (hdata.GetBinContent(i) != 0): 
          err=hdata.GetBinError(i)/hdata.GetBinContent(i);
          dataToObsH.SetBinError(i,err);
#          dataToObs.SetPointError(RPoint, hdata.GetBinWidth(i)/2, hdata.GetBinError(i)/hdata.GetBinContent(i)) 
#          RPoint += 1 
       else: 
          dataToObsH.SetBinError(i,0);  
          continue
          #          dataToObs.SetPointError(RPoint, hdata.GetBinWidth(i)/2, 0) 
#          RPoint += 1 ; continue   

#    dataToObs = rt.TGraphErrors(dataToObsH)  
    dataToObsH.Draw("e0psame0")
    dataToObsH.Draw("e1same")
#    dataToObs.SetMarkerColor(dataToObsH.GetMarkerColor())   
#    dataToObs.SetMarkerStyle(dataToObsH.GetMarkerStyle())
#    dataToObs.SetMarkerSize(dataToObsH.GetMarkerSize())
#    dataToObs.Draw("P 0 same");
    dataToObsH.SetLineWidth(2)  
    
    if(blind>0):
#        blinding_box2 = rt.TPave(hdata.GetBinLowEdge(hdata.FindBin(blind)), 0.4, hdata.GetXaxis().GetXmax(), 1.6, 0, "NB" )  
        blinding_box2 = rt.TPave(hMCtotal.GetBinLowEdge(blind), 0.2, hMCtotal.GetXaxis().GetXmax(), 1.8, 0, "NB" )  
        blinding_box2.SetFillColor(15)
        blinding_box2.SetFillStyle(3013)
        blinding_box2.Draw("same F");
    
    
    line = rt.TLine(hdata.GetXaxis().GetXmin(),1,hdata.GetXaxis().GetXmax(),1);
    line.SetLineColor(rt.kBlack)
    line.Draw("same")
    
    t2.Update()
    
    canvas.cd()
#    canvas.SetLogy(True)
    canvas.Update()
    canvas.RedrawAxis()
    #frame = canvas.GetFrame()
    #frame.Draw()

    #pdf_file = "{}_{}".format(dir1,dir2)
    pdf_file = "%s/%s_%s" % (outdir, dir1, dir2)
    if ( do_linzoom ):
       pdf_file = pdf_file + "-linzoom"
    elif ( do_liny ):
       pdf_file = pdf_file + "-liny"
    else:
       pdf_file = pdf_file + "-logy"

    pdf_file = pdf_file + ".pdf"
    canvas.SaveAs( pdf_file )

    canvas.Update()

########raw_input("Press Enter to end")
#input("Press Enter to end")



