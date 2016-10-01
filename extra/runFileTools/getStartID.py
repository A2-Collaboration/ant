#!/usr/bin/python2

from __future__ import print_function

import sys
import os
import dateutil.parser




def startTime(filename):
    out = os.popen("Ant-info " + filename).readlines()
    TIDline = filter(lambda x: x.count("TID") == 1,out)[0]
    return TIDline.split("'")[1]

def parseTime(isostring):
    return dateutil.parser.parse(isostring)

def startTimes(filelist):
    return map(startTime,filelist)

def main():
    if not len(sys.argv) == 2:
        print("Provide Filelist")

    try:
        with open(sys.argv[1]) as f:
            dates = map(parseTime,startTimes(f.readlines()))
    except:
        dates = [""]
        dates[0] = parseTime(startTime(sys.argv[1]))

    map(print,[d.strftime("%s") for d in dates])


if __name__ == "__main__":
    main()
