#!/usr/bin/python2

from ROOT import TGraph

import sys

graph = TGraph()

#scale to taggEff:
scaleFactor = .58 / .132

def saveFill(line):
    try:
        graph.SetPoint(graph.GetN(),float(line[0]),scaleFactor /  float(line[1].strip('\n')))
    except:
        print line


if len(sys.argv) !=2:
    print "need file to convert"
    exit(1)

filename = sys.argv[1]

with open(filename) as f:
    map(saveFill,map(lambda line: line.split('\t'),f.readlines()))

graph.SaveAs(filename+'.root')
