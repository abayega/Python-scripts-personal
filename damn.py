#!/usr/bin/env python

import re, os, sys, shelve, patternCount4

def openFiles():
    refdict = shelve.open('/home/abayega/Desktop/hmm/refqrydict.txt')
    #refdict = refdictfile['refdict']
  
    qrydict = shelve.open('/home/abayega/Desktop/hmm/refqrydict.txt')
    #qrydict = qrydictfile['qrydict']

    patternCount4.patternCount(qrydict['qrydict'], refdict['refdict'])
