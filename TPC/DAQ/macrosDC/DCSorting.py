import ROOT
import sys
import time
import numpy

#type=sys.argv[1]
#run=sys.argv[2]

#filename=sys.argv[1]

RawList=[]

file0 = ROOT.TFile('TPC_ebdc02_0_physics-00072443-0050.evt_TPCRawDataTree_skip10.root');
#file0 = ROOT.TFile(filename);

tree1 = file0.Get('SampleTree')
sz = tree1.Draw('SumDig:NSampDig:BCO','fee==22&&Channel==201','goff');
#sz = tree1.Draw('SumDig:NSampDig:BCO:Channel','fee==22','goff');
print(sz)

SumDig=tree1.GetV1()
NSampDig=tree1.GetV2()
BCO=tree1.GetV3()

rollover=0
first_flag=1
prevBCO=0
linecount = 0

for i in range(sz):

    linecount +=1
    if(first_flag==1):
       prevBCO=BCO[i]
       first_flag=0

    if(BCO[i]<prevBCO and linecount>10):
        rollover +=1
        linecount=0

    actualBCO = BCO[i] + rollover * 1048576

    RawList.append(['Raw', actualBCO, SumDig[i], NSampDig[i]])

    prevBCO=BCO[i]


tree2 = file0.Get('SampleTreeDC')
#sz2 = tree2.Draw("DCValues[1]:DCSamples[1]:BCO:sampaMAXChannel:sampaAddress","fee==22&&sampaAddress==6&&sampaMAXChannel==15","goff");
sz2 = tree2.Draw("DCValues[1]:DCSamples[1]:BCO","fee==22&&sampaAddress==6&&sampaMAXChannel==15","goff");
print(sz2)

SumDig=tree2.GetV1()
NSampDig=tree2.GetV2()
BCO=tree2.GetV3()

rollover=0
first_flag=1
prevBCO=0
linecount = 0

for i in range(sz2):

    linecount +=1
    if(first_flag==1):
       prevBCO=BCO[i]
       first_flag=0

    if(BCO[i]<prevBCO and linecount>5):
        rollover +=1
        linecount=0

    actualBCO = BCO[i] + rollover * 1048576

    RawList.append(['DCs', actualBCO, SumDig[i], NSampDig[i]])

    prevBCO=BCO[i]

#Sorting
#RawList.sort(key=lambda x: x[1])
RawList2=sorted(RawList,key=lambda x: x[1])


for i in range(sz+sz2):
#    print(RawList[i])
    print(RawList2[i][0],RawList2[i][1],RawList2[i][2])
