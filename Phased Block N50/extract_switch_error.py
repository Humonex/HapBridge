#!/bin/env python3
import sys
from collections import defaultdict
import gzip
import pysam
import re
import itertools
import misc
import time
import logging

#根据switch_error文件中的信息，对VCF文件进行修改并输出

logger = misc.create_logger("Modified_switch_error")
class PhasedBlock(object):
    def __init__(self, chrom, ps):
        self.chrom = chrom
        self.ps = ps
        self.variants = []
    def add(self, pos):
        self.variants.append(pos)
    def start(self):
        return self.variants[0]
    def end(self):
        return self.variants[-1]
    
class PhasedBlockSet(object):

    @staticmethod

    def from_vcf(fname, overlapping=False):
        blk_dict = {}
        CHROM, FORMAT, SAMPLE = 0, -2, -1
        logger.info("PhasedBlockSet.form_vcf")
        for rec in pysam.VariantFile(fname):
             vid = (rec.chrom, rec.pos)
             assert(len(rec.samples) == 1)
             if "PS" in rec.format and type(rec.samples[0]["PS"]) == int:
                 bid = (rec.chrom, rec.samples[0]["PS"])
                 if bid not in blk_dict:
                     blk_dict[bid] = PhasedBlock(bid[0], bid[1])
                 blk_dict[bid].add(rec.pos)
        logger.info("Load %d phased blocks from %s" % (len(blk_dict), fname))
    
        if not overlapping:
            prev = None
            trivial = set()
            for k, v in sorted(blk_dict.items()):
                if prev != None and prev[0][0] == k[0] and prev[1].end() > k[1]:
                    trivial.add(k)
                    logger.info("detected overlapping blocks: %s:%d" % k)
                else:
                    prev = (k, v)
            [blk_dict.pop(k) for k in trivial]
        blk_list = sorted(blk_dict.values(), key=lambda x: (x.chrom, x.ps))
        
        return blk_list 


def Disconnect_Block_with_flip_error(swith_error_file, ifname, truth_name, ofname):
    with open(swith_error_file, 'r') as file:
        lines = file.readlines()

    switch = {}
    for line in lines:
        columns = line.split()  # 分割每一行
        key = columns[0]  # 第一列作为键
        value = int(columns[1])  # 第二列作为值
        if key in switch:  # 如果键已经存在
            switch[key].append(value)
        else:
            switch[key] = [value]
    print(switch)
    blocks = PhasedBlockSet.from_vcf(ifname)
    blocks_vcf = pysam.VariantFile(ifname)


    ovcf = pysam.VariantFile(ofname, "w", header=blocks_vcf.header)

    for blk in blocks :
        if blk.chrom in switch:
            blk_switch = {}
            for s in switch[blk.chrom]:
                if blk.start() <= s <= blk.end():
                    if blk.ps in blk_switch:
                        blk_switch[blk.ps].append(s)
                    else:
                        blk_switch[blk.ps] = [s]
            for rec in blocks_vcf.fetch(blk.chrom, blk.start(), blk.end()):
                if "PS" in rec.format and type(rec.samples[0]["PS"]) == int:
                    ps_value = rec.samples[0]["PS"]
                    if ps_value in blk_switch:
                        for s in blk_switch[ps_value]:
                            if rec.pos >= s:
                                rec.samples[0]["PS"] = s
                ovcf.write(rec)


def Disconnect_Block_no_flip_error(swith_error_file, ifname, truth_name, ofname):
    with open(swith_error_file, 'r') as file:
        lines = file.readlines()

    switch = {}
    chromosome_list = [f'chr{i}' for i in range(1, 23)]
    print (chromosome_list)
    for line in lines:
        columns = line.split()  # 分割每一行
        key = columns[0]  # 第一列作为键
        value = int(columns[1])  # 第二列作为值
        if key in switch:  # 如果键已经存在
            switch[key].append(value)
        else:
            switch[key] = [value]
    print(switch)
    blocks = PhasedBlockSet.from_vcf(ifname)
    blocks_vcf = pysam.VariantFile(ifname)
    truth_vcf = pysam.VariantFile(truth_name)
    true_switch = {}
    for blk in blocks:
        truth_gt = {}
        if blk.chrom in chromosome_list:
            for truth_rec in truth_vcf.fetch(blk.chrom, blk.start(), blk.end()):
                if "PS" in truth_rec.format and type(truth_rec.samples[0]["PS"]) == int:
                    gt1 = truth_rec.samples[0]["GT"]
                    if truth_rec.pos in truth_gt:
                        truth_gt[truth_rec.pos].append(gt1[1])
                    else:
                        truth_gt[truth_rec.pos] = gt1[1]
        result = {}
        bv = blocks_vcf.fetch(blk.chrom, blk.start(), blk.end())
        for rec_true in bv:
            if "PS" in rec_true.format and type(rec_true.samples[0]["PS"]) == int:
                gt = rec_true.samples[0]["GT"]
                if rec_true.pos in truth_gt:
                    if gt[1] == truth_gt[rec_true.pos]:
                        result[rec_true.pos] = 2    
                    else:
                        result[rec_true.pos] = 1
        keys = list(result.keys())
        if blk.chrom in switch:
            for s in switch[blk.chrom]:
                if s in result:
                    current_index = keys.index(s)
                    if current_index < len(keys) -3 :
                            next_key = keys[current_index + 3]
                            pre_key =  keys[current_index - 3]
                            if blk.start() < next_key < blk.end():
                                next_value = result[next_key]
                                pre_value = result[pre_key]
                                if result[s] != next_value and result[s] == pre_value:
                                    if blk.chrom in true_switch:
                                        true_switch[blk.chrom].append(s)
                                    else:
                                        true_switch[blk.chrom] = [s]
    print("true_switch", true_switch)

    ovcf = pysam.VariantFile(ofname, "w", header=blocks_vcf.header)

    for blk in blocks :
        if blk.chrom in true_switch:
            blk_switch = {}
            for s in true_switch[blk.chrom]:
                if blk.start() <= s <= blk.end():
                    if blk.ps in blk_switch:
                        blk_switch[blk.ps].append(s)
                    else:
                        blk_switch[blk.ps] = [s]
            for rec in blocks_vcf.fetch(blk.chrom, blk.start(), blk.end()):
                if "PS" in rec.format and type(rec.samples[0]["PS"]) == int:
                    ps_value = rec.samples[0]["PS"]
                    if ps_value in blk_switch:
                        for s in blk_switch[ps_value]:
                            if rec.pos >= s:
                                rec.samples[0]["PS"] = s
                ovcf.write(rec)
        else:
            for rec in blocks_vcf.fetch(blk.chrom, blk.start(), blk.end()):
                ovcf.write(rec)

def main():
    import argparse
    parse = argparse.ArgumentParser(description='Modified the VCF file according to the switch error file')
    parse.add_argument('iswitch_error',help='The input Switch_error file.')
    parse.add_argument('ivcf',help='The input VCF file.')
    parse.add_argument('ivcf_truth',help='The input truth VCF file.')
    parse.add_argument('ovcf',help='The output VCF file.')
    args = parse.parse_args()
    Disconnect_Block_with_flip_error(args.iswitch_error, args.ivcf,args.ivcf_truth, args.ovcf)
    #Disconnect_Block_no_flip_error(args.iswitch_error, args.ivcf,args.ivcf_truth, args.ovcf)

if __name__ == '__main__':
    main()