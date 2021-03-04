import ROOT as rt
import tdrstyle
import array
import os
import sys
import math

from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double

#set the tdr style
tdrstyle.setTDRStyle()

year="2016"

vh_tag="WH"
## Input file with top pt corrections:
fcor = rt.TFile("plotter_"+vh_tag+"_"+year+"_ttbar_2020_06_19_forLimits.root","READ")

## Input file without top pt corrections:
funcor = rt.TFile("plotter_"+vh_tag+"_"+year+"_ttbar_2020_06_19_v0_forLimits.root","READ")


## Target ROOT file to store up and down templates in BDT
flimit = rt.TFile("plotter_"+year+"_"+vh_tag+"_Sys_noSoftb_forLimits-replaced-signal-mc.root","UPDATE")

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

if vh_tag == "WH":
    histos = ['e_A_CR_3b','mu_A_CR_3b','e_A_SR_3b','mu_A_SR_3b',
              'e_A_CR_4b','mu_A_CR_4b','e_A_SR_4b','mu_A_SR_4b']
elif vh_tag == "ZH":
    histos = ['emu_A_CR_3b','emu_A_SR_3b','ee_A_CR_3b','mumu_A_CR_3b','ee_A_SR_3b','mumu_A_SR_3b',
              'emu_A_CR_4b','emu_A_SR_4b','ee_A_CR_4b','mumu_A_CR_4b','ee_A_SR_4b','mumu_A_SR_4b']
else:
    print("Error, cannot recognize the tag: {}".format(vh_tag))
    print(vh_tag)
    exit(-1)


dirs = ['t#bar{t} + light_filt1','t#bar{t} + b#bar{b}_filt5','t#bar{t} + c#bar{c}_filt4',
        'Other Bkgds','Z#rightarrow ll','W#rightarrow l#nu','QCD',
        'Wh (12)','Wh (15)','Wh (20)','Wh (25)','Wh (30)','Wh (40)','Wh (50)','Wh (60)']

for dir in dirs:
    for histo in histos:
        hname=dir+"/"+histo+"_bdt"
        hname_shapes=dir+"/"+histo+"_bdt_shapes"

        if dir.__contains__("filt"): ## modify only for ttbar MC
            hcor = fcor.Get(hname)
            huncor = funcor.Get(hname)
            
            h_up=hcor.Clone(histo+"_bdt_topptup")
            h_up.Reset()
            h_down=hcor.Clone(histo+"_bdt_topptdown")
            h_down.Reset()

            ## Create the up and down variations due to top pt reweighting unc.    
            makeUncHisto(hcor, huncor, h_up, h_down)

            ## make the 2d versions BDT-vs-index and store as well:
            hcor_shapes = fcor.Get(hname_shapes)
            h_up_shapes=hcor_shapes.Clone(histo+"_bdt_shapes_topptup")
            h_down_shapes=hcor_shapes.Clone(histo+"_bdt_shapes_topptdown")
    
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
                h_up = hcor.Clone(histo+"_bdt_topptup")
                h_down = hcor.Clone(histo+"_bdt_topptdown")

            ## make the 2d versions BDT-vs-index and store as well:
            h_up_shapes=hcor_shapes.Clone(histo+"_bdt_shapes_topptup")
            h_down_shapes=hcor_shapes.Clone(histo+"_bdt_shapes_topptdown")
    
                
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
