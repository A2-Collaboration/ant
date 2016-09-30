#!/usr/bin/python2

from ROOT import TGraph

import sys

graph = TGraph()

#scale to taggEff:
scaleFactor = .58 / .132

yFunction = {
        ''                    : lambda x: x,
        'p2Lscaled'           : lambda x: scaleFactor / x,
        'invert'              : lambda x : 1.0 / x
}

def saveFill(line):
    try:
        graph.SetPoint(graph.GetN(),float(line[0]),yFunction['p2Lscaled'](float(line[1].strip('\n'))))
    except:
        print '   Skipping line: ', line


if not (len(sys.argv)==2 or len(sys.argv)==3):
    print  'usage: ', sys.argv[0] , 'filename [functionName]'
    exit(1)

filename = sys.argv[1]
function = ''

if len(sys.argv) == 3: 
    function = sys.argv[2]
    if not function in yFunction.keys():
        print 'Function ', function , ' not in list of functions:  ', yFunction.keys()
        exit(2)


with open(filename) as f:
    map(saveFill,map(lambda line: line.split('\t'),f.readlines()))

graph.SaveAs(filename+'.root')
