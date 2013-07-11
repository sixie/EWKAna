#! /usr/bin/env python

from ROOT import TH1F, TH2F

histoDefinition = {

    # WGamma* plots
    'gammaDR' : {'variable':'gammaDR','xLabel':'#Delta R','yLabel':'Events / 0.1 rad','binning':[25,0,2.5]},
    'gammaMass' : {'variable':'gammaMass','xLabel':'m_{#mu#mu} [GeV/c^{2}]','yLabel':'Events / 0.5 GeV/c^{2}','binning':[24,0,12]},
    'gammaMet' : {'variable':'met','xLabel':'MET [GeV/c]','yLabel':'Events / 4 GeV/c','binning':[25,0,100]},
    'gammaMinMet' : {'variable':'min(pmet,pTrackMet)','xLabel':'MinMET [GeV/c]','yLabel':'Events / 4 GeV/c','binning':[25,0,100]},
    'gammaMt' : {'variable':'gammaMt','xLabel':'MT [GeV/c]','yLabel':'Events / 4 GeV/c','binning':[25,0,100]},
    'gammaLep3Pt' : {'variable':'min(39.9,lep3.Pt())','xLabel':'p_{T}^{2} [GeV/c]','yLabel':'Events / 2 GeV/c','binning':[20,0,40]},
    'gammaCharge' : {'variable':'abs(lq1+lq2)','xLabel':'|q_{l1}+q_{l2}|','yLabel':'Events','binning':[3,-0.5,2.5]},
    'gammaWWCharge' : {'variable':'abs(gammaWWCharge)','xLabel':'|q_{l1}+q_{l2}|','yLabel':'Events','binning':[3,-0.5,2.5]},


    # WW plots
    'dPhi' : {'variable':'dPhi','xLabel':'#Delta#Phi [rad]','yLabel':'Events / 0.1 rad','binning':[15,0,3.3]},  
    'metPhi' : {'variable':'metPhi','xLabel':'MET #Phi [rad]','yLabel':'Events / 0.1 rad','binning':[33,0,3.3]},
    'dPhiMetDiLep' : {'variable':'acos(cos(metPhi-lep2.Phi()))','xLabel':' #Delta#Phi ','yLabel':'Events / 0.1 rad','binning':[33,0,3.3]},  
    'trackMetPhi' : {'variable':'trackMetPhi','xLabel':'track MET #Phi [rad]','yLabel':'Events / 0.1 rad','binning':[33,0,3.3]},  
    'mt' : {'variable':'mt','xLabel':'m_{T} [GeV/c^{2}]','yLabel':'Events / 8 GeV/c^{2}','binning':[50,0,400]},  
    'mll' : {'variable':'dilep.M()','xLabel':'m_{ll} [GeV/c^{2}]','yLabel':'Events / 8 GeV/c^{2}','binning':[40,0,320]},
    'lep1Pt' : {'variable':'lep1.Pt()','xLabel':'p_{T}^{1} [GeV/c]','yLabel':'Events / 4 GeV/c','binning':[30,0,120]},
    'lep2Pt' : {'variable':'lep2.Pt()','xLabel':'p_{T}^{2} [GeV/c]','yLabel':'Events / 4 GeV/c','binning':[30,0,120]},
    'met' : {'variable':'min(pmet,pTrackMet)','xLabel':'MET [GeV/c]','yLabel':'Events / 4 GeV/c','binning':[40,0,160]},
    'mex' : {'variable':'met*cos(metPhi)','xLabel':'MEX [GeV/c]','yLabel':'Events / 4 GeV/c','binning':[80,-160,160]},
    'mey' : {'variable':'met*sin(metPhi)','xLabel':'MEY [GeV/c]','yLabel':'Events / 4 GeV/c','binning':[80,-160,160]},
    'dilepPt' : {'variable':'dilep.Pt()','xLabel':'p_T^{ll} [GeV/c^{2}]','yLabel':'Events / 8 GeV/c^{2}','binning':[30,0,120]},
    'dilepPtAndMet' : {'variable':'dilep.Pt()+met','xLabel':'p_T^{ll} [GeV/c^{2}]','yLabel':'Events / 8 GeV/c^{2}','binning':[40,0,240]},
    'ptDiff' : {'variable':'min(lep1.Pt()-lep2.Pt(),49)','xLabel':'p_T^1-p_T^2 [GeV/c^{2}]','yLabel':'Events / 8 GeV/c^{2}','binning':[50,0,50]},
    'njets' : {'variable':'njets','xLabel':'jet multiplicity','yLabel':'Events','binning':[4,-0.5,3.5]},  

    
    # top tagging
    'lowBtagDiscriminator' : {'variable':'max(-4.99,min(jetLowBtag,9.99))','xLabel':'bTag discriminator','yLabel':'Events / 0.5 units','binning':[30,-5,10]},
    'jet1BtagDiscriminator' : {'variable':'max(-4.99,min(jet1Btag,9.99))','xLabel':'bTag discriminator','yLabel':'Events / 0.5 units','binning':[30,-5,10]},
    }



colorsMap = {
    'signal':1,
    'WW': 422,

    'Z+jets': 418,
    'Z+X': 418,
    'Drell-Yan': 418,
    'Z#rightarrow#tau#tau': 910,

    'W+jets': 920,
    'W#gamma': 872,
    'DiBoson': 861,
    'di-boson': 861,
    'WZ': 861,
    'ZZ': 863,
    
    'top': 419,
    'tt': 419,
    'wt': 797,
    'tw': 797,    
    'tW': 797,    
    }

