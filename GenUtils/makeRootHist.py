#! /usr/bin/python

import sys,os
from math import sqrt

infile=sys.argv[1]
outfileName=sys.argv[2]

numbins=1000
minval=0.
maxval=500.
binsize=(maxval-minval)/float(numbins)

hist=[]
for i in range(numbins):
    hist.append(0)

outfile=open(outfileName,"w")
nStr="TH1D* "+outfileName.split('.')[0]+"() {\n"
nStr+= "  TH1D* hist=new TH1D(\"h"+outfileName.split('.')[0]+"\",\"\","+str(numbins)+","+str(minval)+","+str(maxval)+");\n"
outfile.write(nStr)

rfcatC="rfcat /data/blue/sixie/Madgraph/WW/"+infile
(fin,fout)=os.popen4(rfcatC)
for line in fout.readlines():
    line=line.rstrip()
    if (len(line.split()) > 0):
        pdgId=line.split()[0]
        if (pdgId=="24"):      # it's the higgs
            pt=sqrt(float(line.split()[6])**2+float(line.split()[7])**2)
#            pt=float(line.split()[10])
            whichbin=int((pt-minval)/binsize)
            if(whichbin<numbins and whichbin>=0):
                hist[whichbin]+=1


for i in range(numbins):
    prStr="  hist->SetBinContent("+str(i+1)+","+str(hist[i])+");\n"
    outfile.write(prStr);

endStr="  hist->Draw();\n"
endStr+="  return hist;\n\n}"
outfile.write(endStr);

outfile.close()
