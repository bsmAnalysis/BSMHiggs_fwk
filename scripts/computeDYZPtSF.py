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


"""
Gets the value of a given item
(if not available a default value is returned)
"""
def getByLabel(desc,key,defaultVal=None) :
    try :
        return desc[key]
    except KeyError:
        return defaultVal

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

def weightedAverage(ratio, hNLO, threshold, end=300):
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
    # Define the Canvas
    c = r.TCanvas(name,name,800,800)
    
    # Upper plot will be in pad1
    hLO.SetLineColor(r.kBlue)
    hLO.SetLineWidth(2)
    hNLO.SetLineColor(r.kRed)
    hNLO.SetLineWidth(2)
    
    ymin = min(hLO.GetMinimum(), hNLO.GetMinimum())
    ymax = 1.1*max(hLO.GetMaximum(), hNLO.GetMaximum())
    if ymin<0: ymin = 1.1*ymin
    else: ymin = 0.9*ymin
    hLO.GetYaxis().SetRangeUser(ymin,ymax)
    hLO.GetXaxis().SetRangeUser(0,300)
    hLO.GetYaxis().SetTitleOffset(1.55)
    hLO.GetYaxis().SetTitleSize(25)
    hLO.GetYaxis().SetTitleFont(43)
    hLO.GetYaxis().SetLabelFont(43)
    hLO.GetYaxis().SetLabelSize(30)
    hLO.GetXaxis().SetTitleSize(30)
    hLO.GetXaxis().SetTitleFont(43)
    hLO.GetXaxis().SetLabelFont(43)
    hLO.GetXaxis().SetLabelSize(30)
    pad1 = r.TPad("pad1","pad1", 0, 0.3, 1, 1)
    pad1.SetBottomMargin(0) # Upper and lower plot are joined
    pad1.SetGridx()         # Vertical grid
    pad1.SetGridy()         # Vertical grid
    pad1.Draw()
    pad1.cd()		    # pad1 becomes the current pad
    hLO.SetStats(0)	    # No statistics on upper plot
    hLO.Draw("ehist")
    hNLO.Draw("ehistsame")
    leg = r.TLegend(0.65,0.75,0.9,0.85)
    leg.AddEntry(hLO,"LO DY scaled","lp")
    leg.AddEntry(hNLO,"NLO DY scaled","lp")
    leg.Draw("same")
    SetOwnership( leg, 0 ) # 0 = release (not keep), 1 = keep

#    hLO.GetYaxis().SetLabelSize(0.)
#    axis = r.TGaxis(-5, 20, -5, 220, 20,220,510,"")
#    axis.SetLabelFont(43)
#    axis.SetLabelSize(15)
#    axis.Draw()

    c.cd()		    # Go back to the main canvas before defining pad2
    pad2 = r.TPad("pad2", "pad2", 0, 0, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.SetGridy()
    pad2.Draw()
    pad2.cd()

    # Define the ratio plot
    ratio = ratio_h.DrawCopy("ehist")
    ratio.GetXaxis().SetRangeUser(0,300)
    ratio.GetYaxis().SetRangeUser(0,2)
#    ratio = hNLO.Clone(name+"_clone")
    ratio.SetLineColor(r.kBlack)
    ratio.SetLineWidth(2)
#    ratio.SetMinimum(0.8)
#    ratio.SetMaximum(1.35)
    ratio.Sumw2()
    ratio.SetStats(0);      # No statistics on lower plot
#    ratio.Divide(hLO)
    ratio.SetMarkerStyle(21)

    ratio.SetTitle("")      # Remove the ratio title
    ratio.GetYaxis().SetTitle("NLO/LO ")
    ratio.GetYaxis().SetTitleSize(25)
    ratio.GetYaxis().SetTitleFont(43)
    ratio.GetYaxis().SetTitleOffset(1.55)
    ratio.GetYaxis().SetLabelFont(43)
    ratio.GetYaxis().SetLabelSize(25)
    ratio.GetXaxis().SetTitleSize(35)
    ratio.GetXaxis().SetTitleFont(43)
    ratio.GetXaxis().SetTitleOffset(4.)
    ratio.GetXaxis().SetLabelFont(43)
    ratio.GetXaxis().SetLabelSize(35)
    ratio.Draw("ehist");       # Draw the ratio plot
    r.gPad.Update()
    line = r.TLine(r.gPad.GetUxmin(), 1, r.gPad.GetUxmax(), 1)
    line.Draw("same")
    SetOwnership( line, 0 )
    
    return c


def produceZptSFs(inputLO, inputNLO, output_name):
    print(output_name+":")
    inFileLO = r.TFile(inputLO, 'READ')
    inFileNLO = r.TFile(inputNLO, 'READ')
    output = r.TFile(output_name, 'RECREATE')
    
    for hist in ['alljets']:
	histLO = inFileLO.Get(hist+"_jetsMulti")
	histNLO = inFileNLO.Get(hist+"_jetsMulti")
	if histLO:
	    canvas = drawHist(histLO,hist+'Multiplicity_LO')
	    canvas.Write()
	    canvas.SaveAs(os.path.split(output_name)[0]+'/'+hist+'Multiplicity_LO.pdf')
	if histNLO:
	    canvas = drawHist(histNLO,hist+'Multiplicity_NLO')
	    canvas.Write()
	    canvas.SaveAs(os.path.split(output_name)[0]+'/'+hist+'Multiplicity_NLO.pdf')
    
    for hist in ['alljets','3jets','4jets','5+jets','2b_3j_jets','2b_4j_jets','2b_geq5j_jets','3b_3j_jets','3b_4j_jets','3b_geq5j_jets','4b_4j_jets','4b_geq5j_jets','5b_geq5j_jets']:
	for ztype in ['']:
    	    histLO = inFileLO.Get(hist+ztype+"_ptw")
    	    histNLO = inFileNLO.Get(hist+ztype+"_ptw")
    	    if not (histLO and histNLO): continue
    	    if not (abs(histLO.Integral())>0 and abs(histNLO.Integral())>0): continue
    	
    	    histLO.SetDirectory(0)
    	    histNLO.SetDirectory(0)
    	    histLO.Scale(1./abs(histLO.Integral()))
    	    histNLO.Scale(1./abs(histNLO.Integral()))
    	    ratios_out = r.TH1F(hist+ztype+'_sf',hist+ztype+'_sf',histLO.GetXaxis().GetNbins(), histLO.GetXaxis().GetXmin(), histLO.GetXaxis().GetXmax())
    	    for i in range(1,histLO.GetXaxis().GetNbins()+1):
    	        denominator = histLO.GetBinContent(i)
    	        numerator = histNLO.GetBinContent(i)
    #	        if(abs(denominator)>0.):print("numerator: {}, denominator: {},numerator/denominator:{}".format(numerator,denominator,numerator/denominator))
    	        if(abs(denominator)>0.): 
    		    ratios_out.SetBinContent(i,numerator/denominator)
    		    denominator_error = histLO.GetBinError(i)
    		    numerator_error = histNLO.GetBinError(i)
    		    if(abs(numerator)>0.):
    		        ratios_out.SetBinError(i, math.sqrt((denominator_error/denominator)**2+(numerator_error/numerator)**2)*(numerator/denominator) )
    		    else:
    		        ratios_out.SetBinError(i,0)
    	        else: ratios_out.SetBinContent(i,0)
    	    print("  "+hist+ztype+'_sf')
	    ratio_out = weightedAverage(ratios_out,histNLO,150)
    	    if "lowPt" in inputLO: c = ratioPlot(histLO,histNLO, ratios_out, hist+ztype+" Low Mass DY")
    	    elif "highPt" in inputLO: c = ratioPlot(histLO,histNLO, ratios_out, hist+ztype+ " High Mass DY")
    	    c.Write()
    	    c.SaveAs(os.path.split(output_name)[0]+'/'+hist+ztype+ "_HighMassDY.pdf")
    	    ratios_out.Write()
    inFileLO.Close()
    inFileNLO.Close()
    output.Close()

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

samplesDB=''
inputdir=''
outdir=''

DtagsList = []
who = commands.getstatusoutput('whoami')[1]

for o,a in opts:
    if o in('-j'): samplesDB = a
    elif o in('-d'): inputdir = a
    elif o in('-o'): outdir = a

jsonFile = open(samplesDB,'r')
procList=json.load(jsonFile,encoding='utf-8').items()

inputdir += '/MC'

DYs = []
#run over sample to merge LO and NLO DY files
for proc in procList :
    #run over processes
    for desc in proc[1] :

        mctruthmode=getByLabel(desc,'mctruthmode',0)
        tag = getByLabel(desc,'tag','')
	# extract Z pt reweights only from LO and NLO DY samples
	if not tag == "Z#rightarrow ll": continue 
	
	order = "LO_DY"
        data = desc['data']
        for d in data :
            origdtag = getByLabel(d,'dtag','')
	    if "10to50" in origdtag: continue
            dtag = origdtag
            if(mctruthmode!=0) : dtag+='_filt'+str(mctruthmode)

            outfile = outdir +'/'+ dtag + '_' + 'ZPt.root'
	    if os.path.isfile(outfile): 
	        DYs.append(dtag)
		continue
            status, output = commands.getstatusoutput('ls '+inputdir+'/'+dtag+'_*.root')
            if status > 0 :
                print "!!!!! Warning: No root files for the dtag: " + origdtag
                continue
            segment = 0
            for file in glob.glob(inputdir+'/'+dtag+'_*.root'):
                out_temp = outdir +'/'+ dtag + '_' + str(segment) + '_zpt.root'
                commands.getstatusoutput('rootcp --recreate '+file+':*jets* '+out_temp)
#                commands.getstatusoutput('rootcp '+file+':*npartons* '+out_temp)
                segment += 1
            commands.getstatusoutput('hadd -f '+outfile+' '+outdir +'/'+ dtag + '_*' + '_zpt.root')
            commands.getstatusoutput('rm -rf '+outdir +'/'+ dtag + '_*' + '_zpt.root')
	    DYs.append(dtag)
            
LO_lowpt = ''
LO_highpt = ''
NLO_lowpt = ''
NLO_highpt = ''
for dtag in DYs:
    root_file = outdir +'/'+ dtag + '_' + 'ZPt.root'
    if 'amcNLO' in dtag:
        if '10to50' in dtag: NLO_lowpt += (' ' + root_file)
        else: NLO_highpt += (' ' + root_file)
    else:
        if '10to50' in dtag: LO_lowpt += (' ' + root_file)
        else: LO_highpt += (' ' + root_file)
if len(LO_lowpt) > 0: commands.getstatusoutput('hadd -f '+outdir +'/LODY_lowPt.root'+' '+LO_lowpt)
if len(LO_highpt) > 0: commands.getstatusoutput('hadd -f '+outdir +'/LODY_highPt.root'+' '+LO_highpt)
if len(NLO_lowpt) > 0: commands.getstatusoutput('hadd -f '+outdir +'/NLODY_lowPt.root'+' '+NLO_lowpt)
if len(NLO_highpt) > 0: commands.getstatusoutput('hadd -f '+outdir +'/NLODY_highPt.root'+' '+NLO_highpt)

if os.path.isfile(outdir +'/LODY_lowPt.root') and os.path.isfile(outdir +'/NLODY_lowPt.root'): produceZptSFs(outdir+'/LODY_lowPt.root',outdir+'/NLODY_lowPt.root',outdir+'/DYSF_lowPt.root')
if os.path.isfile(outdir +'/LODY_highPt.root') and os.path.isfile(outdir +'/NLODY_highPt.root'): produceZptSFs(outdir+'/LODY_highPt.root',outdir+'/NLODY_highPt.root',outdir+'/DYSF_highPt.root')
raw_input("Press [ENTER] to exit...")
