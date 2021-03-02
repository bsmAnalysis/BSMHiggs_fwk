import ROOT as rt
#import tdrstyle
import array
import os
import sys
import math

from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double


year="2016" #to use only in 2016 DY MC

vh_tag="ZH"

## Input file without DY DR  corrections:
funcor = rt.TFile("plotter_"+vh_tag+"_"+year+"_2020_06_19_forLimits_v1.root","READ")

## Input file with DY DR corrections:      
fcor = rt.TFile("plotter_"+vh_tag+"_"+year+"_2020_06_19_forLimits_v2.root","READ")       

## Target ROOT file to store up and down templates in BDT (must be similar to v1, ie. w/o corrections)
flimit = rt.TFile("plotter_"+vh_tag+"_"+year+"_2020_06_19_forLimits.root","UPDATE")

def drawHist(cname, hcor, h_up, h_down):

    canvas = rt.TCanvas(cname,cname,700,700)
    hcor.Draw("E1")
    hcor.SetMarkerSize(0.5)
    
    h_up.Draw("HISTSAME")
    h_up.SetFillColor(0)
    h_up.SetLineColor(1)
    h_up.SetLineStyle(2)
    h_up.SetLineWidth(2)
    
    h_down.Draw("HISTSAME")
    h_down.SetFillColor(0)
    h_down.SetLineColor(1)
    h_down.SetLineStyle(3)
    h_down.SetLineWidth(2)
    
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()
    canvas.SaveAs(cname+".pdf")

    canvas.Update()
    
    return canvas


def makeUncHisto(hcor, huncor, h_up, h_down):
    
    ## Makee the up and down variations of hcor
    for i in range(-1,hcor.GetNbinsX()+1):
        deltah=abs(hcor.GetBinContent(i)-huncor.GetBinContent(i))
        err=math.sqrt(hcor.GetBinError(i)*hcor.GetBinError(i) +
                      huncor.GetBinError(i)*huncor.GetBinError(i))
        h_up.SetBinContent(i,hcor.GetBinContent(i)+deltah)
        h_up.SetBinError(i,err)
        h_down.SetBinContent(i,hcor.GetBinContent(i)-deltah)
        h_down.SetBinError(i,err)

    return h_up,h_down


## START MAIN:

histos = ['emu_A_CR_3b','emu_A_SR_3b','ee_A_CR_3b','mumu_A_CR_3b','ee_A_SR_3b','mumu_A_SR_3b',
          'emu_A_CR_4b','emu_A_SR_4b','ee_A_CR_4b','mumu_A_CR_4b','ee_A_SR_4b','mumu_A_SR_4b']

dirs = ['t#bar{t} + light_filt1','t#bar{t} + b#bar{b}_filt5','t#bar{t} + c#bar{c}_filt4','data',
        'Other Bkgds','Z#rightarrow ll','W#rightarrow l#nu','QCD',
        'Wh (12)','Wh (15)','Wh (20)','Wh (25)','Wh (30)','Wh (40)','Wh (50)','Wh (60)']

for dir in dirs:
    for histo in histos:
        hname=dir+"/"+histo+"_bdt"
        hname_shapes=dir+"/"+histo+"_bdt_shapes"

        if dir.__contains__("Z#rightarrow ll") and histo.__contains__("CR_3b"):
            hcor = fcor.Get(hname)
            huncor = funcor.Get(hname)
            
            if hcor==None: continue
            if huncor==None: continue    
            
            h_up=hcor.Clone(histo+"_bdt_dydRUp")
            h_up.Reset()
            h_down=hcor.Clone(histo+"_bdt_dydRDown")
            h_down.Reset()
            
            ## Create the up and down variations due to DRave modeling unc.    
            makeUncHisto(hcor, huncor, h_up, h_down)
            
            ## make the 2d versions BDT-vs-index and store as well:
            hcor_shapes = fcor.Get(hname_shapes)
            h_up_shapes=hcor_shapes.Clone(histo+"_bdt_shapes_dydRUp")
            h_down_shapes=hcor_shapes.Clone(histo+"_bdt_shapes_dydRDown")
            
            h_up_shapes.Reset()
            h_down_shapes.Reset()
            
            n=35 ## That is 35 bins in BDT
            for i in range(n):
                for j in range(i,n):
                    h_up_shapes.SetBinContent(i,j,h_up.GetBinContent(j))
                    h_down_shapes.SetBinContent(i,j,h_down.GetBinContent(j))
                    
        else: 
            hcor = flimit.Get(hname)
            hcor_shapes = flimit.Get(hname_shapes)
            if hcor==None:
                print("Histo is Null for that process ", hname) 
                continue
            else:   
                h_up = hcor.Clone(histo+"_bdt_dydRUp")
                h_down = hcor.Clone(histo+"_bdt_dydRDown")

            ## make the 2d versions BDT-vs-index and store as well:
            h_up_shapes=hcor_shapes.Clone(histo+"_bdt_shapes_dydRup")
            h_down_shapes=hcor_shapes.Clone(histo+"_bdt_shapes_dydRdown")
    
                
        ## Draw nominal BDT and up and down variations superimposed
        drawHist(dir+"_"+histo+"_bdt",hcor,h_up,h_down)
    
        ## STore up and down variations in original input file for limits:
        flimit.cd(dir)
        
        h_up.Write()
        h_up_shapes.Write()
        h_down.Write()
        h_down_shapes.Write()

flimit.Close()
        
raw_input("Press Enter to end")
