#!/usr/bin/env python3

import ROOT
from matplotlib import pyplot as plt
#import matplotlib.pyplot as plt
import sys
from array import array

x1 = array ('d')
y1 = array ('d')
x2 = array ('d')
y2 = array ('d')
x3 = array ('f')
y3 = array ('f')

SumCharge=0
SumSamp=0
SumDCCharge=0
SumDCSamp=0
DCsCount=0

error=0
nline=0

h1 = ROOT.TH1F('interval','interval',30,-0.5,29.5)

with open(sys.argv[1],'r') as f:
    for line in f:
        parts = line.strip().split()

        if(parts[0] =="Raw:"):
            SumCharge += int(parts[3])
            SumSamp += int(parts[4])
            nline+=1
        else:
            DCsCount+=1
            SumDCCharge += int(parts[3])
            SumDCSamp += int(parts[4])

            if(DCsCount==2):
                print(parts[3],SumCharge,parts[4],SumSamp)
                x1.append(int(SumDCCharge))
                y1.append(int(SumCharge))
                x2.append(int(SumDCSamp))
                y2.append(int(SumSamp))
                SumCharge=0
                SumSamp=0
                SumDCCharge=0
                SumDCSamp=0
                DCsCount=0

            h1.Fill(int(nline))
            nline=0


npt = len(x1)
print("points are:",npt)

c1 = ROOT.TCanvas('c1','Charge correlation',600,700)
c1.Divide(1,2)
c1.cd(1)
c1.GetPad(1).SetGridy()
c1.GetPad(1).SetGridx()
gr1 = ROOT.TGraph(npt, x1, y1)
gr1.SetMarkerStyle(24)
gr1.SetMarkerSize(0.8)
gr1.SetMarkerColor(ROOT.kBlue)
gr1.SetTitle("DC: FEE firmware vs 100% streaming, run 72443, sector 2, Fee 11, ch 201")
gr1.GetYaxis().SetTitleOffset(1.3)
gr1.GetXaxis().SetTitle("DC from FEE firmware [adu]")
gr1.GetYaxis().SetTitle("DC from 100% streaming [adu]")
gr1.Draw('ap')

c1.cd(2)
c1.GetPad(2).SetGridy()
c1.GetPad(2).SetGridx()
gr2 = ROOT.TGraph(npt, x2, y2)
gr2.SetMarkerStyle(24)
gr2.SetMarkerSize(0.8)
gr2.SetMarkerColor(ROOT.kBlue)
gr2.SetTitle("Nsamp: FEE firmware vs 100% streaming, run 72443, sector 23, Fee 11, ch 201")
gr2.GetYaxis().SetTitleOffset(1.3)
gr2.GetXaxis().SetTitle("Nsamp from FEE firmware [adu]")
gr2.GetYaxis().SetTitle("Nsamp from 100% streaming [adu]")
gr2.Draw('ap')

c1.Update()
c2 = ROOT.TCanvas('c2','Interval',600,400)
h1.Draw()
c2.Update()

print("error:",error)
d= input("OK?")
