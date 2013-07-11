#Make WW Bkg Control Region and Scale Factor Table
root -l ComputeWWBkgScaleFactor.C+'(13,"/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/mitf-alljets/backgroundC_skim2.root","/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP//mitf-alljets/data_2l_skim2.root")'

#Make Mll plot illustrating WW control region
root -l MakeMllExtrapolationPlot.C+

#Make Shape Systematics Plots
root -l PlotShapeSystematics.C+
