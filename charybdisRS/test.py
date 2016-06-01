#!/usr/bin/env python

from xml.dom.minidom import parse
import sys, os

from ROOT import *
from array import array

lheFileName1 = sys.argv[1]
lheFileName2 = sys.argv[2]

lheFile1 = parse(lheFileName1)

gROOT.ProcessLine(".x rootlogon.C")

## Analysis histograms
hStAdd  = TH1F("hStAdd" , "Transverse momentum;S_{T} (GeV/c);Candidates per 150GeV/c", 100, 0, 15000)
hStRs  = TH1F("hStRs" , "Transverse momentum;S_{T} (GeV/c);Candidates per 150GeV/c", 1000, 0, 15000)
hMEtAdd = TH1F("hMEtAdd", "Missing E_{T};Missing transverse momentum;Events per 10GeV/c", 50, 0, 500)
hMEtRs = TH1F("hMEtRs", "Missing E_{T};Missing transverse momentum;Events per 10GeV/c", 50, 0, 500)


def fill(hist, value):
    valueUp = hist.GetXaxis().GetXmax()
    binWidth = hist.GetXaxis().GetBinWidth(hist.GetNbinsX()-1)
    hist.Fill(min(value, valueUp-binWidth*1e-9))

## Load init info
event = 0
for eventNode in lheFile1.getElementsByTagName("event"):

    eventTexts = eventNode.childNodes[0].nodeValue.strip().split('\n')

    info = eventTexts[0].strip().split()
    n, proc = int(info[0]), int(info[1])        
    weight, qscale = float(info[2]), float(info[3])
    aqed, aqcd = float(info[4]), float(info[5])

    event += 1

    st = 0.
    met = 0.
    nJets = 0
    nTquark, nBquark, nLquark = 0, 0, 0
    nGluons, nLepton, nPhoton, nOthers = 0, 0, 0, 0

    for i in range(n):
        pline = eventTexts[i+1].strip()
        if pline == "" or pline[0] == '#': continue
        pline = pline.split()

        pdgId   = int(pline[0])
        status  = int(pline[1])
        mother1 = int(pline[2])
        mother2 = int(pline[3])
        color1  = int(pline[4])
        color2  = int(pline[5])
        px = float(pline[ 6])
        py = float(pline[ 7])
        pz = float(pline[ 8])
        e  = float(pline[ 9])
        m  = float(pline[10])
        ctau = float(pline[11])
        spin = float(pline[12])

        if status != 1: continue

        absPdgId = abs(pdgId)
        candLVec = TLorentzVector(px, py, pz, e)
        pt = candLVec.Pt()


        if absPdgId == 6:
            nTquark += 1
        elif absPdgId == 5:
            nBquark += 1
            nJets += 1
        elif absPdgId == 21:
            nGluons += 1
            nJets += 1
        elif absPdgId == 22:
            nPhoton += 1
        elif absPdgId in (1,2,3,4):
            nLquark += 1
            nJets += 1
        elif absPdgId in (11, 13, 15):
            nLepton += 1
        elif absPdgId in (12, 14, 16):
            met += pt
        else:
            nOthers += 1

        if pt > 50 and absPdgId not in (11, 13, 15):
            st += pt
    if met > 50: st += met

    fill(hStAdd, st)
    fill(hMEtAdd, met)


lheFile2 = parse(lheFileName2)
def fill(hist, value):
    valueUp = hist.GetXaxis().GetXmax()
    binWidth = hist.GetXaxis().GetBinWidth(hist.GetNbinsX()-1)
    hist.Fill(min(value, valueUp-binWidth*1e-9))

## Load init info
event = 0
for eventNode in lheFile2.getElementsByTagName("event"):

    eventTexts = eventNode.childNodes[0].nodeValue.strip().split('\n')

    info = eventTexts[0].strip().split()
    n, proc = int(info[0]), int(info[1])
    weight, qscale = float(info[2]), float(info[3])
    aqed, aqcd = float(info[4]), float(info[5])

    event += 1

    st = 0.
    met = 0.
    nJets = 0
    nTquark, nBquark, nLquark = 0, 0, 0
    nGluons, nLepton, nPhoton, nOthers = 0, 0, 0, 0

    for i in range(n):
        pline = eventTexts[i+1].strip()
        if pline == "" or pline[0] == '#': continue
        pline = pline.split()

        pdgId   = int(pline[0])
        status  = int(pline[1])
        mother1 = int(pline[2])
        mother2 = int(pline[3])
        color1  = int(pline[4])
        color2  = int(pline[5])
        px = float(pline[ 6])
        py = float(pline[ 7])
        pz = float(pline[ 8])
        e  = float(pline[ 9])
        m  = float(pline[10])
        ctau = float(pline[11])
        spin = float(pline[12])

        if status != 1: continue

        absPdgId = abs(pdgId)
        candLVec = TLorentzVector(px, py, pz, e)
        pt = candLVec.Pt()


        if absPdgId == 6:
            nTquark += 1
        elif absPdgId == 5:
            nBquark += 1
            nJets += 1
        elif absPdgId == 21:
            nGluons += 1
            nJets += 1
        elif absPdgId == 22:
            nPhoton += 1
        elif absPdgId in (1,2,3,4):
            nLquark += 1
            nJets += 1
        elif absPdgId in (11, 13, 15):
            nLepton += 1
        elif absPdgId in (12, 14, 16):
            met += pt
        else:
            nOthers += 1

        if pt > 50 and absPdgId not in (11, 13, 15):
            st += pt
    if met > 50: st += met

    fill(hStRs, st)
    fill(hMEtRs, met)



cN = TCanvas("cN", "cN", 500, 500)
#hNFrame.SetMaximum(max(x.GetMaximum() for x in hNs)*1.2)
hStAdd.SetLineColor(kRed)
hStAdd.GetXaxis().SetNdivisions(505)
hStAdd.SetMinimum(0.1)

hStAdd.Draw()
legN = TLegend(0.65, 0.65, 0.92, 0.92)
legN.SetFillStyle(0)
legN.SetBorderSize(0)
legN.AddEntry(hStAdd, hStAdd.GetTitle(), "l")
hStRs.SetLineColor(kBlue)
hStRs.Draw("same")
legN.AddEntry(hStRs, hStRs.GetTitle(), "l")
legN.Draw()

