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

parser = argparse.ArgumentParser()
parser.add_argument( "prodmode", type=str, help="Production mode (wh or zh)")
parser.add_argument( "year", type=str, help="Year (2016, 2017 or 2018)")

args = parser.parse_args()

wz = args.prodmode

year = args.year

if wz=="wh":
    channels = ["mu_A_CR_3b", "mu_A_CR_4b", "mu_A_SR_3b", "mu_A_SR_4b", "e_A_CR_3b", "e_A_CR_4b", "e_A_SR_3b", "e_A_SR_4b"]
    e_mu = [""]
elif wz=="zh":
    channels = ["mumu_A_CR_3b", "mumu_A_CR_4b", "mumu_A_SR_3b", "mumu_A_SR_4b", "ee_A_CR_3b", "ee_A_CR_4b", "ee_A_SR_3b", "ee_A_SR_4b"]
    e_mu = ["e", "mu"]

verbose = True

inf = rt.TFile("haa4b_60_13TeV_{}.root".format(wz),"READ")     

dir = "shapes_fit_s/"  

for ch in e_mu:

  if ( verbose ) : print("\n\n verbose:  ch = ", ch , "\n")

  if wz=="zh":
      file = rt.TFile("fit{}_m60_zh_{}.root".format(year,ch),"update")
  elif wz=="wh":
      file = rt.TFile("fit{}_m60_wh.root".format(year),"update")

  for chan in channels:

      if (ch not in chan): continue

      data=file.Get(dir+chan+"/"+"data")    
      
      inh = inf.Get(chan+"/data_obs")
      edges = []
      np = inh.GetNbinsX()
      for i in range(1, np+1):
          edges.append(inh.GetXaxis().GetBinLowEdge(i))
      edges.append(inh.GetXaxis().GetBinUpEdge(i))
      print(edges)
          
      hdata = rt.TH1F("hdata", "data bdt", np, array('d',edges))
      
      px, py = Double(), Double() 
      pyerr = Double() 
      
      alpha = 1 - 0.6827;
      
      nPoints=data.GetN()
      for i in range(0,nPoints):
          data.GetPoint(i,px,py)

          N=px  
          if N==0:
              L=0
              U= rt.Math.gamma_quantile_c(alpha/2,N+1,1) 
              data.SetMarkerSize(0.5);
              data.SetMarkerStyle (20);
              data.SetPointEYlow(i,0)
              data.SetPointEYhigh(i,U-N)
            
          pyerr=data.GetErrorY(i)            
         
          hdata.SetBinContent(i+1,py)
          hdata.SetBinError(i+1,pyerr)

      mydir=file.GetDirectory(dir+chan)
      mydir.cd()
#    hdata.SetDirectory(mydir)    
      hdata.Write()

  file.Close()
