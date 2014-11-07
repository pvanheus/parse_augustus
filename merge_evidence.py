#!/usr/bin/env python

import argparse
import json
import glob

parser = argparse.ArgumentParser('merged evidence files into single file')
parser.add_argument('input_pattern', help='Pattern (glob) for files to use as input')
parser.add_argument('output_file', type=argparse.FileType('w'), help='Output file to write')
args = parser.parse_args()

evidence = dict()
for filename in glob.glob(args.input_pattern):
    input_file = open(filename)
    try:
        evidence_dict = json.load(input_file)
    except ValueError as e:
        message is "Failed to load JSON from {}: {}".format(input_file.name, str(e))
        raise ValueError(message)
    input_file.close()
    for gene, evidence_list in evidence_dict.items():
        assert not gene in evidence, "Gene {} from {} already present in evidence: duplicate key!".format(gene, input_file.name)
        evidence[gene] = evidence_list
json.dump(evidence, args.output_file)