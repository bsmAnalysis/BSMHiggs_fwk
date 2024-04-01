import re
from sys import argv, stdout, stderr, exit
from optparse import OptionParser
from collections import defaultdict
import math
from decimal import Decimal

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
import ROOT
ROOT.gROOT.SetBatch(True)

SRonly = True
channelNameMapping = {"otherbkg": "Other Bkgds", "ttbarbba": "$t\\bar{t} + b\\bar{b}$", "ttbarcba": "$t\\bar{t} + c\\bar{c}$", "ttbarlig": "$t\\bar{t} + light$", "qcd": "qcd", "ddqcd": "ddqcd", "wlnu": "$W\\rightarrow l\\nu$", "zll": "$Z\\rightarrow ll$", "wh": "Wh", "total": "total"}

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("-u", "--uncertainties", default=False, action="store_true", help="Report the uncertainties from the fit(s) too")

(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

errors = False
if options.uncertainties: 
    errors = True

def sci_notation(number, sig_fig=1):
    if(number==""): return number
    value, error = number.split(' \pm ')
    value = float(value)
    error = float(error)
    if(value<1000): return number
    #order = math.log10(value/error)
    if abs(error) > 0: order = int(math.log10(value)) - int(math.log10(error))
    else: order = int(math.log10(value))
    if(order>1): 
        sig_fig = int(order)
	if(value > error*pow(10,sig_fig)): sig_fig += 1
    ret_string = "{0:.{1:d}e}".format(value, sig_fig)
    a,b = ret_string.split("e")
    b = int(b) #removed leading "+" and strips leading zeros too.
    e = error/pow(10,b)
    return "(" + a + " \pm " + str(eval("%.0e"%e)) + ") \\times 10^" + str(b)

file = ROOT.TFile.Open(args[0]);
prefit = file.Get("norm_prefit")
fit_s = file.Get("norm_fit_s")
fit_b = file.Get("norm_fit_b")
if prefit == None: stderr.write("Missing prefit in %s. Did you run FitDiagnostics in a recent-enough version of combine and with --saveNorm?\n" % file);
#if fit_s  == None: raise RuntimeError, "Missing fit_s in %s. Did you run FitDiagnostics with --saveNorm?" % file;
if fit_b  == None: raise RuntimeError, "Missing fit_b in %s. Did you run FitDiagnostics with --saveNorm?" % file;

outf = open("evtYields_postfit.tex", "w+")
outf.write("\\documentclass{article}\n\\usepackage{graphicx}\n\\usepackage{geometry}\n\\geometry{\n\tleft=10mm,\n\tright=10mm,\n\ttop=10mm,\n\tbottom=10mm\n}\n\\usepackage[utf8]{inputenc}\n\\usepackage{rotating}\n\\begin{document}\n\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Post b-fit event yields expected for backgroun processes}\n\\label{tab:table}\n\\resizebox{\\textwidth}{!}{\n ")
bFitTables = defaultdict(defaultdict)

#iter = fit_s.createIterator()
iter = fit_b.createIterator()
Headline = "%-30s %-30s         | pre-fit |         | signal+background Fit |         | bkg-only Fit |"%("| Channel |","| Process |") if (prefit and errors) else "%-30s %-30s         signal+background Fit         bkg-only Fit"%("Channel","Process")
#print Headline

while True:
    norm_b = iter.Next()
    if norm_b == None: break;
    norm_p = prefit.find(norm_b.GetName()) if prefit else None
    m = re.match(r"(\w+)/(\w+)", norm_b.GetName());
    if m == None: m = re.match(r"n_exp_(?:final_)?(?:bin)+(\w+)_proc_(\w+)", norm_b.GetName());
    if m == None: raise RuntimeError, "Non-conforming object name %s" % norm_b.GetName()
    #if norm_b == None: raise RuntimeError, "Missing normalization %s for background fit" % norm_s.GetName()
    if prefit and norm_p and errors:
#        print "%-30s %-30s %7.1f +/- %7.1f         %7.1f +/- %7.1f         %7.1f +/- %7.1f" % (m.group(1), m.group(2), norm_p.getVal(), norm_p.getError(), -1.0, -1.0, norm_b.getVal(), norm_b.getError())
	if m.group(1) not in bFitTables: bFitTables[m.group(1)] = defaultdict(str) 
	if m.group(2) == "wh": bFitTables[m.group(1)][m.group(2)] = " %.1f \pm %.1f" % (norm_p.getVal(), norm_p.getError())
	else:
	    bFitTables[m.group(1)][m.group(2)] = " %.1f \pm %.1f" % (norm_b.getVal(), norm_b.getError()) 
    else:
        if norm_p and prefit:
#            print "%-30s %-30s %7.1f %7.1f %7.1f" % (m.group(1), m.group(2), norm_p.getVal(),  -1.0,  norm_b.getVal())
	    if m.group(2) == "wh": bFitTables[m.group(1)][m.group(2)] = " %.1f " % (norm_p.getVal()) 
	    else: 
	        bFitTables[m.group(1)][m.group(2)] = " %.1f " % (norm_b.getVal()) 
        else:
#            print "%-30s %-30s %7.1f %7.1f" % (m.group(1), m.group(2), -1.0, norm_b.getVal())
	    bFitTables[m.group(1)][m.group(2)] = " %.1f " % (norm_b.getVal()) 

# here to filter channels 
if SRonly:
    todelete = []
    for channel, procs in bFitTables.iteritems():
        if "SR" not in channel or ("emu" in channel): todelete.append(channel)
    for c in todelete: del bFitTables[c]

outf.write("\\begin{tabular}{|c|" + "c|"*len(bFitTables) + "}\\hline\n")
_, max_d = max(bFitTables.items(), key = lambda x: len(set(x[1]))) # look for the channel with max number of procs
process = [p for p,v in max_d.iteritems()] # store these procs into a list
channels = [c for c,v in bFitTables.iteritems()] # store channels into a list
header = "channel &" + " & ".join(["$"+c.replace("_","\\_").replace("mu","\\mu")+"$" for c in channels]) + "\\\\\\hline\n" # header for the table
outf.write(header)
process.remove("wh")
# deal with total process
for c in channels:
    value = 0.
    err = 0.
    for p in process:
	if bFitTables[c][p] == "": continue
        v, e = bFitTables[c][p].split(' \pm ')
	value += float(v)
	err += pow(float(e),2)
    bFitTables[c]["total"] = " %.1f \pm %.1f" % (value, math.sqrt(err))

process.append("total")
process.append("wh")
for p in process:
    if p == "total": 
	outf.write("\\hline\n")
	line = channelNameMapping[p]
    elif p == "wh": 
        outf.write("\\hline\n") # save wh process at last
        line = channelNameMapping[p] + "60 (prefit)"
    else: line = channelNameMapping[p]
#    for c in channels: line += " & $ %s" % (bFitTables[c][p]) + "$"
    for c in channels: line += " & $ %s" % (sci_notation(bFitTables[c][p])) + "$"
    line += "\\\\\n"
    outf.write(line)
outf.write("\\hline\n")
outf.write("\\end{tabular}\n}\n\\end{center}\n\\end{sidewaystable}\n\\end{document}\n")

#for channel, procs in bFitTables.iteritems():
#    print(channel)
#    for proc, yields in procs.iteritems():
#        print("%s : %s" % (proc, yields))


outf.close()
