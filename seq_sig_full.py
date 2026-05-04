#!/usr/bin/python3
#-*- coding : utf-8 -*-

import os, sys, re
import requests
from requests.adapters import HTTPAdapter, Retry

url = 'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28organism_id%3A2697049%29%20AND%20%28reviewed%3Atrue%29'
curl -H "Accept: text/plain; format=flatfile" "https://rest.uniprot.org/uniprotkb/P12345"

with open(sys.argv[1], 'r') as f # file name as argument
	accession_number = []
	full_sequence = []

	for line if f:
		row = line.rstrip().split('\t') # for column view

		accession_number = row[0]
		
