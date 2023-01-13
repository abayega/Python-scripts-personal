#!/bin/bash

#This script takes a list of directories which contain many fasta files and loops over all files then returns sizes of DNA in bases

cat *.fast5 | awk '/!^>/{print length}' >sizes.txt


