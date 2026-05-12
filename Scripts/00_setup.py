# Install packages
import os
import shutil
import glob
import pandas as pd
import plotly.express as px

print("Packages installed.")

# Required to run Qiime in python notebook
import sys
_ = (sys.path.append("/usr/local/lib/python3.9/site-packages"))

# Functions 
## Container object holding metadata and file paths for sample
class Sample_metadata:
    def __init__(self, sampleid, fastq):
        self.id = sampleid
        self.f = fastq
    def add_rev(self, rev):
        self.r = rev

## Retrieve all gzipped FASTQ files
def get_files(results_loc):
    dir_abspath = os.path.abspath(results_loc)
    file_paths = glob.glob(str(dir_abspath) + "/*.fastq.gz")
    return file_paths

## Parse a FASTQ filename to extract the sample ID and read direction.
def parse_file_name(file_name):
    parts = file_name.split("_")
    #miseq.2
    if len(parts) == 5:
        sampleid, b, c, direction, e = parts
        return sampleid, direction
    #nextseq: 2423211-01-A5_S33_R1_001.fastq.gz
    if len(parts) == 4:
        sampleid, b, direction, d = parts
        bits = sampleid.split("-")
        try:
            sampleid = f"{bits[0]}-{bits[1]}"
        except IndexError:
            print(sampleid)
        return sampleid, direction
    #nextseq2: 015-03-B10_R1.fastq.gz or Gblock-F3_R1.fastq.gz or Gull-d-17-07-24-D11_R1.fastq.gz
    else:
        bits = file_name.split("-")
        if len(bits) == 2:
            sampleid, rem = bits
            a = rem.split(".")[0]
            direction = a.split("_")[1]
            return sampleid, direction
        else:
            rem = bits[-1]
            sampleid = bits[0]
            for bit in bits[1:-1]:
                sampleid = f"{sampleid}-{bit}"
            a = rem.split(".")[0]
            direction = a.split("_")[1]
            return sampleid, direction
        
## match forward and reverse FASTQ files into paired samples       
def sort_paths(file_paths, results_loc):
    samples = {}
    for file_path in file_paths:
        path,file_name = os.path.split(file_path)
        sampleid, direction = parse_file_name(file_name)
        if direction == "R1":
            samples[sampleid] = Sample_metadata(sampleid,file_path)
    for file_path in file_paths:
        path,file_name = os.path.split(file_path)
        sampleid, direction = parse_file_name(file_name)
        if direction == "R2":
            try:
                samples[sampleid].add_rev(file_path)
            except KeyError:
                print("Sample: "+str(sampleid)+" does not have a forward read file please check input")
    return samples
        
## Export a paired-end FASTQ manifest file for downstream analysis
def export(samples, results_loc):
    f = open(results_loc+"/pe-33-manifest", "w")
    f.writelines("sampleid\tforward-absolute-filepath\treverse-absolute-filepath\n")
    for key in samples.keys():
        try:
            line = str(samples[key].id)+"\t"+str(samples[key].f)+"\t"+str(samples[key].r)+"\n"
        except AttributeError:
            print("Sample: "+str(samples[key].id)+" does not have a reverse read file please check input")
        f.writelines(line)
    f.close()

print("Functions loaded.")