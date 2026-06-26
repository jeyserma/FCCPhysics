
import os, sys
import numpy as np
import math
from scipy.constants import c, micro, nano, pi, milli, micro

from PIL import Image
import glob


import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(0)




class GuineaPigReader:
    
    def __init__(self, fInName):
        #self.rdfevts = ROOT.RDataFrame("events", fInName)
        #self.fIn = ROOT.TFile(fInName)

        #self.metadata = self.fIn.Get('metadata')
        self.metadata = ROOT.TChain("metadata")
        for f in fInName:
            self.metadata.Add(f)
        #self.metadata.Add(fInName)

        self.metadata.SetBranchStatus("*", 0)
        self.metadata.SetBranchStatus("GPIntValues", 1)
        self.metadata.SetBranchStatus("GPFloatValues", 1)
        self.metadata.SetBranchStatus("GPDoubleValues", 1)
        self.metadata.SetBranchStatus("GPStringValues", 1)
        self.metadata.SetBranchStatus("GPIntKeys", 1)
        self.metadata.SetBranchStatus("GPFloatKeys", 1)
        self.metadata.SetBranchStatus("GPDoubleKeys", 1)
        self.metadata.SetBranchStatus("GPStringKeys", 1)

        self.metadata.GetEntry(0)
        self.metaKeysI = [str(k) for k in self.metadata.GPIntKeys]
        self.metaKeysF = [str(k) for k in self.metadata.GPFloatKeys]
        self.metaKeysD = [str(k) for k in self.metadata.GPDoubleKeys]
        self.metaKeysS = [str(k) for k in self.metadata.GPStringKeys]
        #self.metaValsI = [k for k in self.metadata.GPIntValues]
        #self.metaValsF = [k for k in self.metadata.GPFloatValues]
        #self.metaValsD = [k for k in self.metadata.GPDoubleValues]
        self.metaValsS = [k for k in self.metadata.GPStringValues]
        self.metadata_vals_average()
        #self.ff(fInName)

    def ff(self, fInName):
        print("TEST")
        sumValsI = None
        sumValsF = None
        sumValsD = None

        nfiles = 0
        import time
        for i, fname in enumerate(fInName):
            print(fname)
            t0 = time.time()
            f = ROOT.TFile.Open(fname, "READ")
            t1 = time.time()

            tree = f.Get("metadata")
            t2 = time.time()

            tree.GetEntry(0)
            t3 = time.time()

            _ = [v[0] for v in tree.GPDoubleValues]
            t4 = time.time()

            f.Close()
            del tree, f
            t5 = time.time()

            print(
                i,
                f"open {t1-t0:.3f}s",
                f"get {t2-t1:.3f}s",
                f"entry {t3-t2:.3f}s",
                f"copy {t4-t3:.3f}s",
                f"close {t5-t4:.3f}s",
            )

        quit()
        for ifn, fname in enumerate(fInName):
            print(fname, ifn)
            f = ROOT.TFile.Open(fname)
            t = f.Get("metadata")

            if not t or t.GetEntries() == 0:
                f.Close()
                continue

            t.GetEntry(0)

            valsI = t.GPIntValues
            valsF = t.GPFloatValues
            valsD = t.GPDoubleValues

            if sumValsI is None:
                sumValsI = [0.0] * len(valsI)
                sumValsF = [0.0] * len(valsF)
                sumValsD = [0.0] * len(valsD)

            for i in range(len(valsI)):
                sumValsI[i] += valsI[i][0]

            for i in range(len(valsF)):
                sumValsF[i] += valsF[i][0]

            for i in range(len(valsD)):
                sumValsD[i] += valsD[i][0]

            nfiles += 1
            f.Close()

            del t
            del f
            #ROOT.gROOT.GetListOfFiles().Print()

        nentries = len(fInName)
        self.metaValsI = [v / nentries for v in sumValsI]
        self.metaValsF = [v / nentries for v in sumValsF]
        self.metaValsD = [v / nentries for v in sumValsD]

    def metadata_vals_average(self):
        nentries = self.metadata.GetEntries()
        print(f"Average metadata over {nentries} events")

        sumValsI = None
        sumValsF = None
        sumValsD = None

        for ie in range(nentries):
            self.metadata.GetEntry(ie)

            valsI = self.metadata.GPIntValues
            valsF = self.metadata.GPFloatValues
            valsD = self.metadata.GPDoubleValues

            if ie == 0:
                sumValsI = [0.0] * len(valsI)
                sumValsF = [0.0] * len(valsF)
                sumValsD = [0.0] * len(valsD)

            for i in range(len(valsI)):
                sumValsI[i] += valsI[i][0]
            for i in range(len(valsF)):
                sumValsF[i] += valsF[i][0]
            for i in range(len(valsD)):
                sumValsD[i] += valsD[i][0]

            '''
            valsI = list(self.metadata.GPIntValues)
            valsF = list(self.metadata.GPFloatValues)
            valsD = list(self.metadata.GPDoubleValues)

            if ie == 0:
                sumValsI = [0.0] * len(valsI)
                sumValsF = [0.0] * len(valsF)
                sumValsD = [0.0] * len(valsD)

            for i, v in enumerate(valsI):
                sumValsI[i] += v[0]
            for i, v in enumerate(valsF):
                sumValsF[i] += v[0]
            for i, v in enumerate(valsD):
                sumValsD[i] += v[0]
            '''
        self.metaValsI = [v / nentries for v in sumValsI]
        self.metaValsF = [v / nentries for v in sumValsF]
        self.metaValsD = [v / nentries for v in sumValsD]

    def get_beam_particles_slice(self, slice):

        beam1p, beam2p = [], []
        tree = self.fIn.events
        for t in tree:
            for p in t.Beam1Slice:
                if p.time != slice:
                    continue
                x,y,z,px,py,pz,m = p.vertex.x, p.vertex.y, p.vertex.z, p.momentum.x, p.momentum.y, p.momentum.z, p.mass
                E = (m*m + px*px + py*py + pz*pz)**0.5
                p = {E,x,y,z}
                beam1p.append(p)

            for p in t.Beam2Slice:
                if p.time != slice:
                    continue
                x,y,z,px,py,pz,m = p.vertex.x, p.vertex.y, p.vertex.z, p.momentum.x, p.momentum.y, p.momentum.z, p.mass
                E = (m*m + px*px + py*py + pz*pz)**0.5
                p = {E,x,y,z}
                beam2p.append(p)

        return beam1p, beam2p

    def get_metdata(self, tag):
        if tag in self.metaKeysI:
            return self.metaValsI[self.metaKeysI.index(tag)]
        elif tag in self.metaKeysF:
            return self.metaValsF[self.metaKeysF.index(tag)]
        elif tag in self.metaKeysD:
            return self.metaValsD[self.metaKeysD.index(tag)]
        elif tag in self.metaKeysS:
            return self.metaValsS[self.metaKeysS.index(tag)]
        else:
            print("NOT FOUND")
        return
        if tag in self.metaKeysI:
            return self.metaValsI[self.metaKeysI.index(tag)][0]
        elif tag in self.metaKeysF:
            return self.metaValsF[self.metaKeysF.index(tag)][0]
        elif tag in self.metaKeysD:
            return self.metaValsD[self.metaKeysD.index(tag)][0]
        elif tag in self.metaKeysS:
            return self.metaValsS[self.metaKeysS.index(tag)][0]
        else:
            print("NOT FOUND")
        
def get_quantile(h, qval=0.95):
    import array
    q = array.array('d', [qval])
    xq = array.array('d', [0.0])
    h.GetQuantiles(1, xq, q)

    quant = xq[0]
    return quant