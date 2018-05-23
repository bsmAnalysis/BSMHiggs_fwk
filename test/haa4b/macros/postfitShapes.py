import ROOT as rt
import CMS_lumi, tdrstyle
import array
import os

from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double

#set the tdr style
tdrstyle.setTDRStyle()

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12

H_ref = 600; 
W_ref = 800; 
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
L = 0.12*W_ref
R = 0.04*W_ref

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

h =  rt.TH1F("h","h; m_{e^{+}e^{-}} (GeV); Events / 0.5 GeV",30, -0.3, 0.3)

#h.SetMaximum(260)

xAxis = h.GetXaxis()
xAxis.SetNdivisions(6,5,0)

yAxis = h.GetYaxis()
yAxis.SetNdivisions(6,5,0)
yAxis.SetTitleOffset(1)

#h.Draw()

file = rt.TFile("LIMIT_final/cards_SB13TeV_SM_Wh/0060/fitDiagnostics.root","READ")

dir = "shapes_prefit/E_SR_qcdB_Mt_3b/"
data = file.Get(dir+"data")
data.SetMarkerStyle(20)
data.SetMarkerColor(1)

total = file.Get(dir+"total_background")

sig = file.Get(dir+"wh")
sig.SetLineColor(rt.kMagenta)
sig.SetLineWidth(2)
sig.Scale(50)

qcd = file.Get(dir+"ddqcd")
qcd.SetFillColor(634)
qcd.SetLineColor(1)
singleto = file.Get(dir+"singleto")
singleto.SetFillColor(831)
singleto.SetLineColor(1)

#ttbargam = file.Get(dir+"ttbargam")
#ttbargam.SetFillColor(424)
#ttbargam.SetLineColor(1)

ttbarjet = file.Get(dir+"ttbarjet")
ttbarjet.SetFillColor(408)
ttbarjet.SetLineColor(1)

## vv = file.Get(dir+"vv")
## vv.SetFillColor(595)
## vv.SetLineColor(1)
## zvv = file.Get(dir+"zvv")
## zvv.SetFillColor(17)
## zvv.SetLineColor(1)

zll = file.Get(dir+"zll")
zll.SetFillColor(624)
zll.SetLineColor(1)
wlnu = file.Get(dir+"wlnu")
wlnu.SetFillColor(622)
wlnu.SetLineColor(1)
#qcd = file.Get(dir+"qcd")


MC = rt.THStack()
MC.Add(qcd)
MC.Add(singleto)
#MC.Add(ttbargam)
MC.Add(ttbarjet)
#MC.Add(vv)
#MC.Add(zvv)
MC.Add(zll)
MC.Add(wlnu)
#MC.Add(qcd)

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

MC.Draw("hist")
data.Draw("epsamex0")
sig.Draw("histsame")

t1.Update()

MC.SetMaximum(1.2*MC.GetMaximum())
MC.GetXaxis().SetTitle("BDT")

## canvas.cd()
## canvas.Update()
## canvas.RedrawAxis()
## frame = canvas.GetFrame()
## frame.Draw()


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

legend =  rt.TPad("legend_0","legend_0",x0_l,y0_l,x1_l, y1_l )
#legend.SetFillColor( rt.kGray )
legend.Draw()
legend.cd()

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

#array must be converted 
x_l = array.array("f",x_l)
ex_l = array.array("f",ex_l)
y_l = array.array("f",y_l)
ey_l = array.array("f",ey_l)

gr_l =  rt.TGraphErrors(1, x_l, y_l, ex_l, ey_l)

rt.gStyle.SetEndErrorSize(0)
gr_l.SetMarkerSize(0.9)
gr_l.Draw("0P")

latex.SetTextFont(42)
latex.SetTextAngle(0)
latex.SetTextColor(rt.kBlack)    
latex.SetTextSize(0.15)    
latex.SetTextAlign(12) 

box_ = rt.TBox()
xx_ = x_l[0]
yy_ = y_l[0]
latex.DrawLatex(xx_+1.*bwx_,yy_,"Data")

yy_ -= gap_
box_.SetLineStyle( rt.kSolid )
box_.SetLineWidth( 1 )
# box_.SetLineColor( kBlack )
box_.SetLineColor( wlnu.GetLineColor() )
box_.SetFillColor( wlnu.GetFillColor() )
box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 )
box_.SetFillStyle(0)
box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 )
#Draw Z->ee text
latex.DrawLatex(xx_+1.*bwx_,yy_,"W #rightarrow l #nu")

yy_ -= gap_
box_.SetLineColor( zll.GetLineColor() )
box_.SetFillColor( zll.GetFillColor() )
box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 )
box_.SetFillStyle(0)
box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 )
#Draw Z->ee text
latex.DrawLatex(xx_+1.*bwx_,yy_,"Z #rightarrow l l")


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


hdata = rt.TH1F("hdata","data bdt",30, 0, 30) 
htotal=total.Clone()

px, py = Double(), Double()
#x, y = Double(), Double()
pyerr = Double()

nPoints=data.GetN()
for i in range(0,nPoints):
    data.GetPoint(i,px,py)
    pyerr=data.GetErrorY(i)
    hdata.Fill(px,py)
    hdata.SetBinError(i,pyerr)

hratio=hdata.Clone("myratio")
hratio.Divide(htotal)

hratio.Draw("e2")
hratio.Draw("same 2 0")


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

line = rt.TLine(0.,1,30.,1);
line.SetLineColor(rt.kBlack)
line.Draw("same")

t2.Update()

canvas.cd()
canvas.Update()
canvas.RedrawAxis()
frame = canvas.GetFrame()
frame.Draw()

canvas.Update()

raw_input("Press Enter to end")
