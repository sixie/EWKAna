##################################
#2011 Analysis
#################################

#Make DY Bkg Control Region and Scale Factor Table
root -l ComputeDYBkgScaleFactor.C+'(2)'

#Make Rout/in Plot
root -l ComputeDataDrivenRoutin.C+'(2)'

#Make Shape Systematics Plots
root -l PlotShapeSystematics.C+


##################################
#New MVAIdIsoCombined Analysis
#################################

#Make DY Bkg Control Region and Scale Factor Table
root -l ComputeDYBkgScaleFactor.C+'(13)'

#Make Rout/in Plot
root -l ComputeDataDrivenRoutin.C+'(13)'
root -l PlotROutIn.C+

#Make Shape Systematics Plots
root -l PlotShapeSystematics.C+
