#! /usr/bin/env python

from ROOT import TFile, TH1D, TCanvas, TChain, THStack, TLegend, TLatex
import sys, os, math, argparse

import rootlogonTDR

from histograms import colorsMap, histoDefinition

finalStates = {1:'me'}

def setWeight(sample, njets, mass, luminosity):

    # DY
    if sample.find('Drell-Yan')>-1 or sample.find('Z+jets')>-1 or sample.find('Z+X')>-1:
        if njets == 0: weight=3
        if njets == 1: weight=3

    # top
    elif sample.find('top')>-1: 
        if njets==0: weight=1.2
        elif njets==1: weight=1.06

    elif sample.find('tt')>-1:
        if njets==0: weight=1.2*0.00212776
        elif njets==1: weight=1.06*0.00212776

    elif sample.find('wt')>-1  or sample.find('tw')>-1 or sample.find('Wt')>-1 or sample.find('tW')>-1:
        if njets==0: weight=1.2*0.0126537
        elif njets==1: weight=1.06*0.0126537

    # Wgamma
    elif sample.find('W#gamma')>-1:
        weight=0.026 #  after applying MC->Data scale factor of 1.6
        #weight=0.016 #original weight from madgraph

    # WZ
    elif sample.find('WZ')>-1:
        weight = 0.0007155

    # ZZ
    elif sample.find('ZZ')>-1:
        weight = 0.0002763*4


    # other samples
    else: weight=1.0

    # multiply by luminosity
    return weight*luminosity


def roundFloat(number,digits):
    stringedNumber = str(number)[:str(number).find('.')+1+digits]
    if number - float(stringedNumber) >= float('0.'+'0'*digits+'5'):
        stringedNumber = str(float(stringedNumber)+float('0.'+'0'*(digits-1)+'1'))
    return float(stringedNumber)


parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),description="analysis")
parser.add_argument('-s',dest='selection',action='store',required=False,help='selection')
parser.add_argument('-e',dest='era',action='store',required=False,help='data taking era')
parser.add_argument('-p',dest='plot',action='store',required=False,help='plot to display')
parser.add_argument('-m',dest='mass',action='store',required=False,help='higgs mass')
parser.add_argument('-n',dest='njets',action='store',required=False,help='jet bin')
parser.add_argument('-u',dest='uniqueBackground',action='store',required=False,help='collapse backgrounds')
args=parser.parse_args()

# lmm
selection =  'type!=3 && njets<2 && jetLowBtag < 2.1 && jet1Btag < 2.1 && lep1.Pt()>20 && lep2.Pt()>10 && lid3!=0 && gammaCharge==0 && gammaMass<12 && min(pmet,pTrackMet)>20 && gammaMt>20'
# mmmm
#selection =  'type==0 && njets<2 && jetLowBtag < 2.1 && jet1Btag < 2.1 && lep1.Pt()>20 && lep2.Pt()>10 && lid3!=0 && gammaCharge==0 && gammaMass<12 && min(pmet,pTrackMet)>20 && gammaMt>20'
# emm
#selection =  '(type==1 || type==2) && njets<2 && jetLowBtag < 2.1 && jet1Btag < 2.1 && lep1.Pt()>20 && lep2.Pt()>10 && lid3!=0 && gammaCharge==0 && gammaMass<12 && min(pmet,pTrackMet)>20 && gammaMt>20'


# plot
plot = 'gammaMass'
if args.plot in histoDefinition: plot = args.plot
plotVariable = histoDefinition[plot]['variable']

# signal mass
mass = 140
if args.mass: mass = int(args.mass)

# jet bin
njets = 0
if args.njets: njets = int(args.njets)

# backgroundColor
uniqueBackground = False
if args.uniqueBackground: uniqueBackground = bool(int(args.uniqueBackground))


# luminosity
luminosity = {'Run2011A':1.2,'Run2011B':1.9,'Run2011':4.0}

# data taking era
era = 'Run2011'
if args.era in ['Run2011','Run2011A','Run2011B']: era = args.era

# datasets
dataDir = '/data/smurf/mzanetti/data/'+era+'/Wgamma/'
dataTypes = {
    '1data': [dataDir+'data_2l.root'],
    '2tt': [dataDir+'ttbar_weighted.root'],
    '3tW': [dataDir+'ttbar_weighted.root'],
    '4W#gamma': [dataDir+'wglmm_weighted.root'],
    '5WZ': [dataDir+'/wz_weighted.root'],
    '6ZZ': [dataDir+'/zz_weighted.root'],
    }


trees = {}

# data
dataInputFile = TFile(dataTypes['1data'][0])
trees['1data'] = dataInputFile.Get('tree')

# load the MC samples
for dataType in dataTypes:
    if dataType.find('data')>-1 or dataType.find('signal')>-1: continue
    trees[dataType] = TChain()
    for sample in dataTypes[dataType]:
        trees[dataType].Add(sample+'/tree')


# set the dictionaries
dumpHistos = {}
tDumpHistos = {}
stacks = {}
lLegends = {}
yields = {}
yields['backgroundPerFlavor'] = {}
filledAlready = {}
filledAlreadyTotal = False
for fState in finalStates:
    stacks[fState] = THStack(finalStates[fState],'')
    yields[fState]= {}
    yields['backgroundPerFlavor'][fState] = 0
    lLegends[fState] = TLegend(0.5,0.4,0.75,0.7)
    lLegends[fState].SetFillColor(10)
    # background color
    filledAlready[fState] = False

# global objects
tStack = THStack('totalStack','')
yields['total'] = {}
lLegends['total'] = TLegend(0.5,0.4,0.75,0.7)
lLegends['total'].SetFillColor(10)


###################
# A N A L Y S I S #
###################
for tree in sorted(trees):

    print 'Processing ', tree[1:]

    # total histos
    tHistoName = 'h'+tree
    tDumpHistos[tHistoName] = TH1D(tHistoName,'',histoDefinition[plot]['binning'][0],histoDefinition[plot]['binning'][1],histoDefinition[plot]['binning'][2])


    # loop over final states
    for fState in sorted(finalStates):

        # no weight for data
        if tree.find('data')>-1: weight = '1'
        # weights (with lumi) for other samples
        else: weight = str(setWeight(tree, njets, mass, luminosity[era]))+'*weight'
        

        # define the histos
        histoName = 'h'+tree+str(fState)
        dumpHistos[histoName] = TH1D(histoName,'',histoDefinition[plot]['binning'][0],histoDefinition[plot]['binning'][1],histoDefinition[plot]['binning'][2])
        if tree.find('signal')==-1 and tree.find('data')==-1:
            if not uniqueBackground:
                dumpHistos[histoName].SetFillColor(colorsMap[tree[1:]])
                dumpHistos[histoName].SetLineColor(colorsMap[tree[1:]])
            else: 
                dumpHistos[histoName].SetFillColor(920)
                dumpHistos[histoName].SetLineColor(920)
        if tree.find('data')>-1:
            dumpHistos[histoName].SetMarkerStyle(20)
            dumpHistos[histoName].SetMarkerSize(1.3)

        
        # filling the histogram
        trees[tree].Draw(plotVariable+'>>'+histoName,'('+ selection +')*'+weight);
        # updating the yields
        yields[fState][tree] = dumpHistos[histoName].Integral(0,10000)

        # adding up flavor blind histos
        if fState==1:
            tDumpHistos[tHistoName] = dumpHistos[histoName].Clone(tHistoName)
        else:
            tDumpHistos[tHistoName].Add(dumpHistos[histoName])
            
        # filling the per-final state stack (and legends too)
        if tree.find('signal')==-1 and tree.find('data')==-1:
            stacks[fState].Add(dumpHistos[histoName])
            if not uniqueBackground: lLegends[fState].AddEntry(dumpHistos[histoName], tree[1:], 'f')
            elif not filledAlready[fState]:
                lLegends[fState].AddEntry(dumpHistos[histoName], 'background', 'f')
                filledAlready[fState] = True
        if tree.find('data')>-1: lLegends[fState].AddEntry(dumpHistos[histoName], 'data', 'p')

        
    # filling the final-state blind stack
    if tree.find('signal')==-1 and tree.find('data')==-1:
        tStack.Add(tDumpHistos[tHistoName])
        if not uniqueBackground: lLegends['total'].AddEntry(tDumpHistos[tHistoName], tree[1:], 'f')
        elif not filledAlreadyTotal:
            lLegends['total'].AddEntry(tDumpHistos[tHistoName], 'background', 'f')
            filledAlreadyTotal = True
    if tree.find('data')>-1: lLegends['total'].AddEntry(tDumpHistos[tHistoName],'data' , 'p')
        
#################
# D R A W I N G #
#################
tex = {}
cCanvases = {} 
cCanvases['total'] = TCanvas('total','total',800,600)
tDumpHistos['h1data'].SetMaximum(2*tDumpHistos['h1data'].GetMaximum())
tDumpHistos['h1data'].GetXaxis().SetTitle(histoDefinition[plot]['xLabel'])
tDumpHistos['h1data'].GetYaxis().SetTitle("Events / 0.5 GeV/c^{2}")
tDumpHistos['h1data'].GetYaxis().SetTitleOffset(1.0)
tDumpHistos['h1data'].Draw('ep')
tStack.Draw('fhist,same')
tDumpHistos['h1data'].Draw('ep,same')
lLegends['total'].Draw('same')
tex['total'] = TLatex(tDumpHistos['h1data'].GetBinLowEdge(tDumpHistos['h1data'].GetNbinsX()/2),
                     0.9*tDumpHistos['h1data'].GetMaximum(),
                     'CMS Preliminary');
tex['total'].Draw();
tex['totalL'] = TLatex(tDumpHistos['h1data'].GetBinLowEdge(tDumpHistos['h1data'].GetNbinsX()/2),
                       0.8*tDumpHistos['h1data'].GetMaximum(),
                       'L='+str(luminosity[era])+' fb^{-1}');
tex['totalL'].Draw();


outputFile = TFile('wGammaStar.root','RECREATE')
for h in tDumpHistos:
    tDumpHistos[h].Write()
for c in cCanvases:
    cCanvases[c].Write()
    cCanvases[c].SaveAs(plot + '.png')


###################
# P R I N T I N G #
###################
yields['background']=0
for tree in trees:
    yields['total'].update({tree:0})
    for fState in finalStates:
        yields['total'][tree]+=yields[fState][tree]
        # summing up all background for all flavors
        if tree.find('signal')==-1 and tree.find('data')==-1: yields['background']+=yields[fState][tree]

for fState in finalStates:
    print ''
    print finalStates[fState]+' final state:'
    for tree in sorted(trees):
        if tree.find('signal')==-1 and tree.find('data')==-1:
            yields['backgroundPerFlavor'][fState] += yields[fState][tree]
            print tree[1:], roundFloat(yields[fState][tree],1)
    print '--------------'
    print 'total background: ', roundFloat(yields['backgroundPerFlavor'][fState],1)
    print 'data: ', roundFloat(yields[fState]['1data'],1)

print ''
print 'Overall counts:'
for tree in sorted(trees):
    if tree.find('signal')==-1 and tree.find('data')==-1:
        print tree[1:], roundFloat(yields['total'][tree],1)
print '--------------'
print 'total background: ', roundFloat(yields['background'],1)
print 'data: ', roundFloat(yields['total']['1data'],1)

raw_input('type whatever to quit')
