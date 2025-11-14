#!/bin/env python3
import sys
from collections import defaultdict
import gzip
import pysam
import re
import itertools
import misc
import time

logger = misc.create_logger("bridge")


class SVPosition(object):
    def __init__(self, chrom, pos, reads):
        self.chrom = chrom
        self.position = pos
        self.reads = reads
        #self.gtype = gtype
        assert len(self.reads) == 2
    @staticmethod
    def from_cutesv(fname):
        svs = []
        CHROM, POS, REF, ALT, QUAL, FILTER, INFO, FORMAT, NULL = 0, 1, 3, 4, 5, 6, 7, 8, 9
        for line in misc.open_file(fname):
            if line.startswith("##"): continue
    
            if line.startswith("#"):
                its = line[1:].split() 
                assert CHROM == its.index("CHROM")
                assert POS == its.index("POS")
                assert REF == its.index("REF")
                assert ALT == its.index("ALT")
                assert QUAL == its.index("QUAL")
                assert FILTER == its.index("FILTER")
                assert FORMAT == its.index("FORMAT")
                assert NULL == its.index("NULL")
                continue
            
            its = line.split()
            if not its[NULL].startswith("0/1") or its[FILTER] != "PASS": continue
            for i in its[INFO].split(";"):
                ii = i.split("=")
                if len(ii) >= 2 and ii[0] == "RNAMES":
                     ns = ii[1].split(",")
                     svp = SVPosition(its[CHROM], int(its[POS])-1, [set(ns), set()])
                     
                     svs.append(svp)
        logger.info("Load %d SV from %s" % (len(svs), fname))
        return svs

    @staticmethod
    def from_straglr(fname):
        strs = defaultdict(lambda:defaultdict(list))
    
        CHROM, POS  = 0, 1
        
        for line in open(fname):
            if line.startswith("#"): continue
    
            its = line.split()
            if its[0] != "chr1": continue
    
            chrom = its[CHROM]
            pos = int(its[POS])
            gt = its[6].split(";")
     
            if len(gt) == 2:
                gt0 = float(gt[0][:gt[0].find('(')])
                gt1 = float(gt[1][:gt[1].find('(')])
                if its[10] == 'NA': continue 
                if abs(float(its[10]) - gt0) < abs(float(its[10]) - gt1):
                    strs[(chrom,pos)][0].append(its[7])
                else:
                    strs[(chrom,pos)][1].append(its[7])
        svs = []
        for (chrom, pos), reads in strs.items():
            r = len(reads[0]) / (len(reads[0]) + len(reads[1]))
            if r >= 0.3 and r <= 0.7:
                svs.append(SVPosition(chrom, pos, reads))
        svs.sort(key=lambda x: x.position)
        return svs


    

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

def test_sv_in_block():

    blocks = load_block("../whatshap/out.vcf")
    svs = SVPosition.from_cutesv("../cuteSV/SV.vcf")

    for b0, b1 in zip(blocks[:-1], blocks[1:]):
        print(b0[0], b1[0])
        # for sv in svs:
        #      if sv.position > b0[1][-1] and sv.position < b1[1][0]:
        #          print(sv.position, b0[1][-1], b1[1][0])


def test_sv_consist_bam():
    rtypes = load_rtype("../whatshap/out.bam")
    #svs = SVPosition.from_cutesv("../cuteSV/SV.vcf")
    svs = SVPosition.from_straglr("../str/chr1_straglr.tsv")
    blocks = load_block("../whatshap/out.vcf")
    for sv in svs:
        in_block = False
        for b in blocks:
            if sv.position > b[1][0] and sv.position < b[1][-1]:
                 in_block = True
                 break

        if not in_block: continue

        bs0 = defaultdict(int)
        for r in sv.reads[0]:
             if r in rtypes:
                 bs0[rtypes[r]] += 1
             else:
                 bs0["NA"] += 1
        bs1 = defaultdict(int)
        for r in sv.reads[1]:
             if r in rtypes:
                 bs1[rtypes[r]] += 1
             else:
                 bs1["NA"] += 1

        gt = list((set(bs0.keys()) | set(bs1.keys())) - set(["NA"]))
        counts = [0, 0]
        if len(gt) == 1:
            counts[0] += bs0[gt[0]]
            counts[1] += bs1[gt[0]]
        elif len(gt) == 2:
            counts[0] += bs0[gt[0]] + bs1[gt[1]]
            counts[1] += bs0[gt[1]] + bs1[gt[0]]
 
        
        s = min(counts) / sum(counts)
        #print("%.02f" % s, sv.position, bs0, bs1)


def sv_test():
    svs = load_sv("../cuteSV/SV.vcf")
    rtypes = load_rtype("../whatshap-indel/out.bam")
    for sv in svs:
        bs = defaultdict(int)
        for r in sv[1]:
            if r in rtypes:
                bs[rtypes[r]] += 1
            else:
                bs["NA"] += 1

        keys = set(bs.keys())  - set(["NA"])
        ppp = set([k[0] for k in keys])
        # print(len(ppp), ppp,  keys)        
        # print(sv[0], bs)
#        if len(bs) >=2 or True:
#            keys = list(bs.keys())
#            if keys[0][0] != keys[1][0]:
#                print(sv[0], bs)


mean = lambda x: sum(x) / len(x)

def kimean(scores):

    center = lambda x: sum(x) / len(x)
    split_set = lambda s, m:[[i for i in s if i < m],[i for i in s if i >= m]]

    split_set2 = lambda s, m: [[i for i in s if abs(i - m[0]) < abs(i-m[1])], [i for i in s if abs(i-m[0]) >= abs(i-m[1])]]

    sz = [0, len(scores)]
    sets = split_set(scores, center(scores))
    while sz[0] != len(sets[0]):
        sz = len(sets[0]), len(sets[1])
        sets = split_set2(scores, [center(sets[0]), center(sets[1])])

    return sets




class Bridger(object):
    def __init__(self):
        self.options = {}
        self.options["margin"] = 100000

    def load_read_haplotype(self, fname):
        bf = pysam.AlignmentFile(fname)

        rtypes = {}
        for al in bf:
            if al.has_tag("HP") and al.has_tag("PS"):
                rtypes[al.qname] = (al.get_tag("PS"), al.get_tag("HP"))

        return rtypes

 
    def get_gap_between_blocks(self, blk0, blk1):
        if(blk0.chrom == blk1.chrom):
            return max(blk0.start(), blk0.end() - self.options["margin"]), min(blk1.start() + self.options["margin"], blk1.end())        


    def extend_block_in_vcf(self, ifname, ofname, links, switchs = {}):
        '''根据links标记得两个相邻得blocks扩展vcf'''
        # rewrite VCF
        #print(links) 
        haplotype = {} # block: block, hp
    
        for blk1, (blk0, flip) in sorted(links.items(), key = lambda x: int(x[0][1])):
            #print(blk0, blk1)
            assert int(blk1[1]) > int(blk0[1]) and blk1[0] == blk0[0]
            swt = 0 if blk0 not in switchs else len(switchs[blk0])
            
            if blk0 in haplotype:
                 blk2, flip1 = haplotype[blk0]
                 haplotype[blk1] = (blk2, flip1 ^ flip ^ (swt % 2 != 0))
            else:
                 haplotype[blk1] = (blk0, flip ^ (swt % 2 != 0))
   
        ivcf = pysam.VariantFile(ifname)
        ovcf = pysam.VariantFile(ofname, "w", header=ivcf.header)
        for rec in ivcf.fetch():
            #assert len(rec.samples) == 1
            if "PS" in rec.format and type(rec.samples[0]["PS"]) == int:
                 blk = (rec.chrom, rec.samples[0]["PS"])
                 swt = [] if blk not in switchs else switchs[blk]
                 gt = 0
                 #if swt != None:
                 for s in swt:
                    if rec.pos > s:
                        gt += 1

                 if blk in haplotype:
                    nblk, flip = haplotype[blk]
                    assert blk[0] == nblk[0]  # chrom

                    rec.samples[0]["PS"] = nblk[1]
            
                    flip = flip ^ (gt % 2 == 1)
                    if flip:
                        gt = rec.samples[0]["GT"]
                        rec.samples[0]["GT"] = (gt[1], gt[0])
                        rec.samples[0].phased = True
                 else:
                    flip = (gt % 2 == 1)
                    if flip:
                        gt = rec.samples[0]["GT"]
                        rec.samples[0]["GT"] = (gt[1], gt[0])
                        rec.samples[0].phased = True
                                              
            ovcf.write(rec)
                 
    

class SVBridger(Bridger):
    def __init__(self):
        super().__init__()
    def bridge_by_str(self, vcf_fname, bam_fname, sv_fname, out_fname):
        #input_vcf_fname = "../whatshap/out.vcf"
        #input_vcf_fname = "./str-out-sv.vcf"
        #blocks = PhasedBlockSet.from_vcf(input_vcf_fname)
        
        #svs = SVPosition.from_cutesv("../cuteSV/SV.vcf")
        #svs = SVPosition.from_straglr("../str/chr1_straglr.tsv")
        rtypes = self.load_read_haplotype(bam_fname)
        blocks = PhasedBlockSet.from_vcf(vcf_fname)
        svs = SVPosition.from_straglr(sv_fname)
   
        links = {} 
        for blk0, blk1 in zip(blocks[:-1], blocks[1:]):
            if blk0.chrom != blk1.chrom: continue
            if blk0.start() != 1279441: continue
            start, end = self.get_gap_between_blocks(blk0, blk1)
           
            cands = []
            for sv in svs:
                if sv.chrom == blk0.chrom and sv.position >= start and sv.position <= end:
                    cands.append(sv)
            bs0 = defaultdict(int)
            bs1 = defaultdict(int)
            for sv in cands:
                print("--", sv.chrom, sv.position)
                for r in sv.reads[0]:
                    if r in rtypes:
                        bs0[rtypes[r]] += 1
                        print(sv.position, 0, r, rtypes[r])
                    else:
                        print(sv.position, 0, r, "NA")
                        bs0["NA"] += 1
                for r in sv.reads[1]:
                    if r in rtypes:
                        bs1[rtypes[r]] += 1
                        print(sv.position, 1, r, rtypes[r])
                    else:
                        print(sv.position, 0, r, "NA")
                        bs1["NA"] += 1
                 
            keys = (set(bs0.keys()) | set(bs1.keys()))  - set(["NA"]) 
            print(keys)
            ppp = set([k[0] for k in keys])
            print("ppp", ppp) 
            if len(ppp) >= 2:
                assert len(ppp) == 2
                pair = []
                for p in ppp:
                    if bs0[(p,1)] + bs1[(p,2)] > bs0[(p,2)] + bs1[(p,1)]:
                        pair.append((p,1))
                    elif bs0[(p,1)] + bs1[(p,2)] < bs0[(p,2)] + bs1[(p,1)]:
                        pair.append((p,2))
    
                if len(pair) == 2:
                    print(pair)
                    links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps) , pair[0][1] != pair[1][1] )
    
        # rewrite VCF
        self.extend_block_in_vcf(vcf_fname, out_fname, links)
    



class MethBridger(Bridger):
    def __init__(self, threads):
        super().__init__()
        self.options = {}
        self.options["threads"] = threads
        self.options["margin"] = 100000
        
        # read 有多少meth就可以判断单倍型
        self.options["read_count"] = 5
        self.options["read_rate"] = 0.6
        self.find_switch_time = 0
        self.bridge_two_blocks_time = 0
   

    def bridge2(self, vcf_fname, bam_fname, ofname, detect_switch=True):
        blocks = PhasedBlockSet.from_vcf(vcf_fname)
        print("start bridge2")  

        bf = pysam.AlignmentFile(bam_fname)
        jobs = [] # (job_type, job_data)
        if detect_switch:
            print("detect_switch")
            for blk in blocks:
                block_alignment = bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True)
                jobs.append([0, blk, len(list(block_alignment))])

        for blk0, blk1 in zip(blocks[:-1], blocks[1:]):

            if(blk0.chrom == blk1.chrom):
                start, end = max(blk0.start(), blk0.end() - self.options["margin"]), min(blk1.start() + self.options["margin"], blk1.end()) 
                block_alignment = bf.fetch(blk.chrom, start, end, multiple_iterators=True)

                jobs.append([1, blk0, blk1, len(list(block_alignment))])

        jobs.sort(key = lambda x: -x[-1] )

        # import multiprocessing as mp
        # with mp.Pool(self.options["threads"]) as pool:
        #     result = pool.starmap(self.run_job, zip(jobs, [bam_fname]*len(jobs)), chunksize=1)
        
        import time 
        start_time_pool = time.time()  
        import multiprocessing as mp
        with mp.Pool(self.options["threads"]) as pool:
            result = pool.starmap(self.run_job, zip(jobs, [bam_fname]*len(jobs)), chunksize=1)
        end_time_pool = time.time()
        run_time_pool = end_time_pool - start_time_pool
        print("----------------------------------------------------------")
        print("Thread pool run time: {:.2f} seconds(the bridge2 )".format(run_time_pool))
        print("----------------------------------------------------------")

        links, switchs = {}, {}
        for job, rz in zip(jobs, result):
            if job[0] == 0:
                swt, blk = rz, job[1]
                switchs[(blk.chrom, blk.ps)] = swt
            elif job[0] == 1:
                lnk, blk0, blk1 = rz, job[1], job[2]
                if lnk == 1:
                    links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), False)
                elif lnk == -1:
                    links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), True)
        # for job, rz in zip(jobs, result):
        #     if job[0] == 1:
        #         lnk, blk0, blk1 = rz, job[1], job[2]
        #         if lnk == 1:
        #             links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), False)
        #         elif lnk == -1:
        #             links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), True)


        self.extend_block_in_vcf(vcf_fname, ofname, links, switchs)

    # def run_job(self, job, bam_fname):
    #     logger.info("XXX start %d %d" % (job[0], job[1].ps))
    #     if job[0] == 0:
    #         result = self.find_switch_in_block(job[1], bam_fname)
    #     elif job[0] == 1:
    #         result =  self.bridge_two_blocks(job[1], job[2], bam_fname)

    #     logger.info("XXX end %d %d" % (job[0], job[1].ps))
    #     return result
    
    #统计job的运行时间
    def run_job(self, job, bam_fname):
        start_time = time.time()
        # if job[0]==0:
        #     logger.info("find_switch_in_block start %d %d %s" % (job[0], job[1].ps,job[1].chrom))
        # elif job[0]==1:
        #     logger.info("bridge_two_blocks start %d %d %s" % (job[0], job[1].ps,job[1].chrom))
        if job[0] == 0:
            # print("不执行find_switch_in_block")
            # result = None
            result = self.find_switch_in_block(job[1], bam_fname)
        elif job[0] == 1:
            #result = None
            result =  self.bridge_two_blocks(job[1], job[2], bam_fname)
        end_time = time.time()
        # if job[0]==0:
        #     logger.info("find_switch_in_block end %d %d %s" % (job[0], job[1].ps,job[1].chrom))
        # elif job[0]==1:
        #     logger.info("bridge_two_blocks end %d %d %s" % (job[0], job[1].ps,job[1].chrom))
        elapsed_time = end_time - start_time
        logger.info("Job %d %d %s took %d seconds to complete" % (job[0], job[1].ps,job[1].chrom, elapsed_time))
        return result

    def bridge(self, vcf_fname, bam_fname, ofname):
        blocks = PhasedBlockSet.from_vcf(vcf_fname)
  
        links = {}
        if self.options["threads"] > 1:
            import multiprocessing as mp
            from functools import partial
            logger.info("start bridge")
            with mp.Pool(self.options["threads"]) as pool:
                result = pool.starmap(self.bridge_two_blocks, zip(blocks[:-1],blocks[1:], [bam_fname]*(len(blocks)-1)))
        
            for lnk, blk0, blk1 in zip(result, blocks[:-1], blocks[1:]):
                if lnk == 1:
                    links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), False)
                elif lnk == -1:
                    links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), True)

        else:
            for blk0, blk1 in zip(blocks[:-1], blocks[1:]):

                if blk0.ps != 107169460: continue
                lnk = self.bridge_two_blocks(blk0, blk1, bam_fname)
                if lnk == 1:
                    links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), False)
                elif lnk == -1:
                    links[(blk1.chrom, blk1.ps)] = ((blk0.chrom, blk0.ps), True)

        self.extend_block_in_vcf(vcf_fname, ofname, links, self.find_switch(vcf_fname, bam_fname))

    def find_switch(self, vcf_fname, bam_fname):
        blocks = PhasedBlockSet.from_vcf(vcf_fname)
  
        switchs = {}
        if self.options["threads"] > 1:
            import multiprocessing as mp
            from functools import partial
            
            with mp.Pool(self.options["threads"]) as pool:
                result = pool.starmap(self.find_switch_in_block, zip(blocks, [bam_fname]*len(blocks)))
            for swt, blk in zip(result, blocks):
                if len(swt) > 0:
                    switchs[(blk.chrom, blk.ps)] = swt

        else:
            for blk in blocks:
               
                #if blk.ps != 2746360: continue
                #if blk.ps != 107169460: continue
                #if blk.ps != 153690139: continue
                #if blk.ps != 195119852: continue
                #if blk.ps != 246823403: continue
 
                #if blk.ps != 1357359: continue
                #if blk.ps != 610545: continue
                swt = self.find_switch_in_block(blk, bam_fname)
                if len(swt) > 0:
                    switchs[(blk.chrom, blk.ps)] = swt

        print("switchs", switchs)
        return switchs

    @staticmethod
    def are_two_blocks_adjacent(blk0, blk1):
        if blk0.chrom != blk1.chrom:
            return False
        # if not (blk0.start() <= blk0.end() and blk0.end() < blk1.start() and blk1.start() <= blk1.end()):
        #     print(blk0.start(), blk0.end(), blk0.end(), blk1.start(), blk1.start(), blk1.end())
        assert blk0.start() <= blk0.end() and blk0.end() < blk1.start() and blk1.start() <= blk1.end()
        return True

    def bridge_two_blocks(self, blk0, blk1, bam_fname):
        bf = pysam.AlignmentFile(bam_fname)
        gap = self.Gap(bf, blk0, blk1, self.options)
        return gap.bridge() 

    def find_switch_in_block(self, blk, bam_fname):
        bf = pysam.AlignmentFile(bam_fname)
        block = self.Block(bf, blk, self.options)
        return block.find_switch()

    class Gap(object):
        """Gap between two adjacent blocks"""
        def __init__(self, bf, blk0, blk1, opts):
            assert MethBridger.are_two_blocks_adjacent(blk0, blk1)
            self.bf = bf
            self.chrom = blk0.chrom
            self.blk0 = blk0
            self.blk1 = blk1
            self.opts = opts
            self.start, self.end = max(blk0.start(), blk0.end() - self.opts["margin"]), min(blk1.start() + self.opts["margin"], blk1.end())        
            #self.info("position=%s:%d-%d" % (self.blk0.chrom, self.start, self.end))

        # def info(self, msg):
        #     logger.info("Gap(%d-%d): %s" %(self.blk0.ps, self.blk1.ps, msg))
    
        def get_meths_linkscores(self, meths):
            #logger.info("get_meths_linkscores(Gap)")
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1
            links = defaultdict(lambda: [0, 0, 0, 0])
            block_alignment = self.bf.fetch(self.chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
               valid_meth = []
               for rloc, m in MethBridger.collect_meths_in_read(read):
                   if rloc in meths and meths[rloc].valid():
                       ps, hp, d = meths[rloc].test(m/256)
                       valid_meth.append((rloc, hp))
               for (p0, h0), (p1, h1) in itertools.combinations(valid_meth, 2):
                    assert p0 < p1 and (h0 == 1 or h0 == 2) and (h1 == 1 or h1 == 2)
                    links[(p0, p1)][h0-1 + (h1-1)*2] += 1
            return links
        
        def extract_consistent_meths(self, meths, links):
            print("extract_consistent_meths")
            def extract_top(links, cands, rate) :
                count = defaultdict(lambda: [0,0])
                for k, v in links.items():
                    if k[0] not in cands or k[1] not in cands: continue
                    c = v[0] + v[3], v[1] + v[2]
                    r = c[0] / sum(c)
                    if abs(r - (1-r)) >= 0.6 and sum(c) >= 10:
                        count[k[0]][0] += 1
                        count[k[1]][0] += 1
                    count[k[0]][1] += 1
                    count[k[1]][1] += 1
                scount = sorted(count.items(), key=lambda x: -x[1][0]/(x[1][1]+1))
            
                ext = set()
                for x in scount:
                    ext.add(x[0])
                    if len(ext) >= len(scount)*rate:
                            break
                return ext

            valid_locs = set([k[0] for k in links]) | set([k[1] for k in links])
            total = len(valid_locs)
            #valid_locs = extract_top(links, valid_locs, 0.5)
            i = 1 
            # while i < 3:
            #     valid_locs = extract_top(links, valid_locs, 0.5)
            #     i = i + 1

            while len(valid_locs) > total*0.15 and len(valid_locs) >= 10:
                valid_locs = extract_top(links, valid_locs, 0.5)
            return valid_locs

        def meth_options(self):
            return 8, 0.3, 0.3
        
        global zero_count
        global no_zero_count
        def bridge(self):
            meths0, meths1, unphased = self.collect_meths()

            meths3 = self.collect_middle_block_meths()
            meths_ASM = {}
            meths0_ASM = {}
            meths1_ASM = {}

            print("the length of the meths0, meths1, meths3", len(meths0), len(meths1), len(meths3))
            if len(meths0) or len(meths1) == 0:
                zero_count += 1
            else:
                no_zero_count += 1
       
            [v.verify_phased(*self.meth_options()) for v in itertools.chain(meths0.values(), meths1.values())]
            [v2.verify(*self.meth_options()) for v2 in itertools.chain(meths3.values())]


            # phasing unphased reads
            # 思路：首先将block中间的甲基化值搜集起来，然后对这些甲基化值进行一致性分析，挑选出一致性分数较高的甲基化值，然后在判断这些甲基化值的时候，加入判断条件
            self.extend_phased_reads(meths0, unphased)
            #在这个地方我需要调用一致性信息，在phase过程中就需要加入一致性信息，如果一致性信息差，就不把unphase的read加入到phased中
            self.extend_phased_reads(meths1, unphased)
            print("extend_phased_reads 后的the length of the meths0, meths1长度", len(meths0), len(meths1))
            #[v3.verify_phased(*self.meth_options()) for v3 in itertools.chain(meths0.values(), meths1.values())]
            #self.add_ASM(meths_ASM)
            links = self.get_meths_linkscores(meths3)
            valid_locs = self.extract_consistent_meths(meths3, links)
            for key in valid_locs:
                if key in meths3:
                    meths_ASM[key] = meths3[key]

            # links0 = self.get_meths_linkscores(meths0)
            # valid_locs0 = self.extract_consistent_meths(meths0, links0)
            # for key0 in valid_locs0:
            #     if key0 in meths0:
            #         meths0_ASM[key0] = meths0[key0]
            # links1 = self.get_meths_linkscores(meths1)
            # valid_locs1 = self.extract_consistent_meths(meths1, links1)
            # for key1 in valid_locs1:
            #     if key1 in meths1:
            #         meths1_ASM[key1] = meths1[key1]
            #meths0_ASM, meths1_ASM = self.statistics_meths_ASM(meths0, meths1, meths0_ASM, meths1_ASM)

            # print("length of the valid_locs", len(valid_locs))
            # for key in valid_locs:
            #     if key in meths0:
            #         meths0_ASM[key] = meths0[key]
            #     if key in meths1:
            #         meths1_ASM[key] = meths1[key]
            for key0 in meths0:
                if key0 in meths_ASM:
                    meths0_ASM[key0] = meths0[key0]
            #print("length of the meths0_ASM:", meths0_ASM)
            for key1 in meths1:
                if key1 in meths_ASM:
                    meths1_ASM[key1] = meths1[key1]
            #print("length of the meths1_ASM:", meths1_ASM)               
            logger.info("------------------------------------------------------------------")
            logger.info("get_link")
            #return self.get_link(meths0, meths1)
            #return self.get_link_ASM(meths0_ASM, meths1_ASM)
            #[v4.verify_phased(*self.meth_options()) for v4 in itertools.chain(meths0_ASM.values(), meths1_ASM.values())]
            return self.get_link_ASM(meths0_ASM, meths1_ASM)


        def extend_phased_reads(self, meths, unphased):
            phased = set()
            while True:
                phasing = self.phase_reads(meths, unphased - phased)
                if len(phasing) > 0:
                    self.add_phasing_reads(meths, phasing)
                    [v.verify_phased(*self.meth_options()) for v in meths.values()]
                    [phased.add(r[0]) for r in phasing]
                else:
                    break
           

        def collect_meths(self):
            # unzip self params
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1

            meths0, meths1 = {}, {}
            unphased = set()
    
            block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
                ps, hp = 0, 0
                if read.has_tag("PS") and read.has_tag("HP"):
                    ps = read.get_tag("PS")
                    hp = read.get_tag("HP")
                    
                    for rloc, m in MethBridger.collect_meths_in_read(read):
                        for meths, blk in ((meths0, blk0), (meths1, blk1)):
                            if blk.ps != ps: continue
    
                            if rloc not in meths:
                                 meths[rloc] = MethBridger.Meth(blk.chrom, rloc)
                            meths[rloc].add(read.qname, read.is_forward, m/256, ps, hp)
        
                else:
                    if read.modified_bases != None: 
                        unphased.add(read)
            #self.info("collect meths, left=%d, right=%d, unphased=%d" % (len(meths0), len(meths1), len(unphased))) 
            return meths0, meths1, unphased


            
        def collect_middle_block_meths(self):
            # unzip self params
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1
            meths3 = {}           
            blk1_start = blk1.start()
            blk0_end = blk0.end()
            block_alignment = self.bf.fetch(chrom, blk0_end, blk1_start, multiple_iterators=True)
            for read in block_alignment:
                ps, hp = 0, 0
                if read.has_tag("PS") and read.has_tag("HP"):
                    ps = read.get_tag("PS")
                    hp = read.get_tag("HP")
                    for rloc, m in MethBridger.collect_meths_in_read(read):
                        if rloc not in meths3:
                            meths3[rloc] = MethBridger.Meth(blk0.chrom, rloc)
                        meths3[rloc].add(read.qname, read.is_forward, m/256, ps, hp)
                else:         
                    for rloc, m in MethBridger.collect_meths_in_read(read):
                        if rloc not in meths3:
                            meths3[rloc] = MethBridger.Meth(blk0.chrom, rloc)
                        meths3[rloc].add(read.qname, read.is_forward, m/256, ps, hp)                    
            return meths3                    


        def phase_reads(self, meths, unphased):
            phasing = []
    
            remove_abnormal = lambda x: x if type(x) ==int else 0
            # params
            sup_count, sup_rate = self.opts["read_count"], self.opts["read_rate"] 
            for read in sorted(unphased, key=lambda x: remove_abnormal(x.reference_start)):
                checks = defaultdict(int)
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc in meths and meths[rloc].valid_phased():
                         (ps, hp, w) = meths[rloc].test_phased(m/256) 
                         if w >= 0.3: # TODO 检查  
                             checks[(ps, hp)] += 1
        
                sel = MethBridger.get_best_support(checks, sup_count, sup_rate)
                if sel != None:
                    phasing.append((read, sel)) #read和他read上每个甲基化位点的12值
            return phasing
        
        def add_phasing_reads(self, meths, phasing):
            for read, (ps, hp) in phasing:
           
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc not in meths:
                        meths[rloc] = MethBridger.Meth(self.chrom, rloc) 
                    meths[rloc].add(read.qname, read.is_forward, m/256, ps, hp)

        def add_ASM(self, meths):
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1
            meths3 = {}           
            blk1_start = blk1.start()
            blk0_end = blk0.end()
            block_alignment = self.bf.fetch(chrom, blk0_end, blk1_start, multiple_iterators=True)
            for read in block_alignment:
                valid_meth = []
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc in meths and meths[rloc].valid():
                        ps, hp, d = meths[rloc].test(m/256) 
                        valid_meth.append((rloc, hp, ps, d))
                if len(valid_meth) != 0:
                    print("valid_meth", valid_meth)           

        def statistics_meths_ASM(self, meths0, meths1, meths0_ASM, meths1_ASM):
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1    
            block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            num_meth0_ASM = len(meths0_ASM)
            num_meth1_ASM = len(meths1_ASM)
            for key, value  in meths0.items():
                if key in meths0_ASM:
                    if meths0[key].valid_phased():
                        num_meth0_ASM -= 1
                        meths0_ASM.pop(key)

            for key , value in meths1.items():
                if key in meths1_ASM:
                    if meths1[key].valid_phased():
                        num_meth1_ASM -= 1
                        meths1_ASM.pop(key)
            print("num_meth0_ASM, num_meth1_ASM", num_meth0_ASM, num_meth1_ASM)
            return meths0_ASM, meths1_ASM

        def get_link_ASM(self, meths0_ASM, meths1_ASM):
            # unzip self params
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1
    
            sup_count, sup_rate = self.opts["read_count"], self.opts["read_rate"]
    
            scores = defaultdict(int)
            block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
                sup0, sup1 = defaultdict(int), defaultdict(int)
                for rloc, m in MethBridger.collect_meths_in_read(read):

                    for meths, sup in ((meths0_ASM, sup0), (meths1_ASM, sup1)):
                        if rloc in meths and meths[rloc].valid_phased():
                            (ps, hp, d) = meths[rloc].test_phased(m/256)
                            if d >= 0.3:
                                sup[(ps, hp)] += 1
                sel0 = MethBridger.get_best_support(sup0, sup_count, sup_rate)
                sel1 = MethBridger.get_best_support(sup1, sup_count, sup_rate)
                if sel0 != None and sel1 != None: 
                    consist = sel0[1] == sel1[1]
                    scores[(sel0[0], sel1[0], consist)] += 1
            #print(scores)
            sel_scores = MethBridger.get_best_support(scores, 2, 0.5)
            if sel_scores != None:
                if sel_scores[2]: # consist
                    return 1
                else:
                    return -1
            else:
                return 0     




            
        
    
        def get_link(self, meths0, meths1):
            # unzip self params
            chrom, start, end, blk0, blk1 = self.chrom, self.start, self.end, self.blk0, self.blk1
    
            sup_count, sup_rate = self.opts["read_count"], self.opts["read_rate"]
    
            scores = defaultdict(int)
            block_alignment = self.bf.fetch(chrom, start, end, multiple_iterators=True)
            for read in block_alignment:
                sup0, sup1 = defaultdict(int), defaultdict(int)
                for rloc, m in MethBridger.collect_meths_in_read(read):
    
                    for meths, sup in ((meths0, sup0), (meths1, sup1)):
                        if rloc in meths and meths[rloc].valid_phased():
                            (ps, hp, d) = meths[rloc].test_phased(m/256)
                            if d >= 0.3:
                                sup[(ps, hp)] += 1
                sel0 = MethBridger.get_best_support(sup0, sup_count, sup_rate)
                sel1 = MethBridger.get_best_support(sup1, sup_count, sup_rate)
                if sel0 != None and sel1 != None:
                    consist = sel0[1] == sel1[1]
                    scores[(sel0[0], sel1[0], consist)] += 1
            #print(scores)
            sel_scores = MethBridger.get_best_support(scores, 2, 0.5)
            if sel_scores != None:
                if sel_scores[2]: # consist
                    return 1
                else:
                    return -1
            else:
                return 0          
    


    @staticmethod
    def collect_meths_in_read(read):
        meths = []
        if read.modified_bases != None:
            ref_loc = read.get_reference_positions(full_length=True)
            for m, locs in read.modified_bases.items():
                if m[0] == 'C' and m[2] == 'm':
                    for lc in locs:
                        if ref_loc[lc[0]] != None:  
                            rloc = ref_loc[lc[0]] + 1 if read.is_forward else ref_loc[lc[0]]
                            meths.append((rloc, lc[1]))
        return meths


    @staticmethod
    def get_best_support(sup, count, rate):
        '''support must be >= max(count, total*rate)
        '''
        
        sorted_sup = sorted(sup.items(), key = lambda x: -x[1])
        sum_sup = sum([i[1] for i in sorted_sup])
 
        if sum_sup > 0 and sorted_sup[0][1] >= max(sum_sup*rate, count):
            return sorted_sup[0][0]
        else:
           return None     



    class Meth(object):
        def __init__(self, chrom, position):
            self.chrom = chrom
            self.position = position
            self.infos = []
            self.centroids = None
            self.centroids_phased = None
            pass
         
        @staticmethod
        def mean(x):

            ava = sum(x) / len(x)
            #return ava
            xlen = len(x)
            off = xlen // 5
            x_sorted = sorted(x)
            assert x_sorted[-1] >= x_sorted[0] 
            if ava < 0.5:
                x_core = x_sorted[0:xlen-off]
                #print(x, ava)
                #print(x_core, sum(x_core)/len(x_core))
                assert ava >= sum(x_core) / len(x_core)
                return sum(x_core) / len(x_core)
            else:
                x_core = x_sorted[off:]
                #print(x, ava)
                #print(x_core, sum(x_core)/len(x_core))
                assert ava <= sum(x_core) / len(x_core)
     
                return sum(x_core) / len(x_core)
 
        def add(self, qname, is_forward, mvalue, ps=0, hp=0):
            assert len(qname) > 10 and mvalue <= 1
            self.infos.append((qname, is_forward, mvalue, ps, hp))

            # TODO
            self.ps = ps
            
        def verify(self, sup_count, sup_rate, distance):
            signals = [i[2] for i in self.infos]
            splited_signals = kimean(signals)
            #print("verify", self.position, splited_signals[0], splited_signals[1])
            
            self.centroids = self._verify_two_signal_sets(splited_signals[0], splited_signals[1], sup_count, sup_rate, distance)

        def verify_phased(self, sup_count, sup_rate, distance):
            signals0 = [i[2] for i in self.infos if i[4] == 1]
            signals1 = [i[2] for i in self.infos if i[4] == 2]
            #print("verify_phase", self.position, signals0, signals1)
            self.centroids_phased = self._verify_two_signal_sets(signals0, signals1, sup_count, sup_rate, distance)


        def _verify_two_signal_sets(self, signals0, signals1, sup_count, sup_rate, distance):
            '''验证两个集合是否是有效的两类，如果是则返回质心，否则返回None
               sup_count: 每一类的最小支持数目
               sup_rate:  每一类的最小支持比例
               distance: 中心的最小距离 [0,1]
            '''

            if len(signals0) > 0 and len(signals1) > 0:
                mean0, mean1 = self.mean(signals0), self.mean(signals1)
                size_total = len(signals0) + len(signals1)
                size_threshold = max(sup_count, sup_rate*size_total)
                #print("vvv", size_threshold, mean0, mean1, distance)
                if len(signals0) > size_threshold and len(signals1) >= size_threshold and abs(mean0-mean1) >= distance:
                    # 返回集合的质心
                    return [mean0, mean1]
            
            return None
              

            

     

        def valid(self):
            return self.centroids != None

        def valid_phased(self):
            return self.centroids_phased != None

        def test(self, m):
            assert self.valid() and m <= 1
            return self._test(m, self.centroids)
      
        def test_phased(self, m):
            assert self.valid_phased() and m <= 1
            return self._test(m, self.centroids_phased)
        
      
        def _test(self, m, centroids):
            hp1_ave, hp2_ave = centroids
            t = m - (hp1_ave + hp2_ave) / 2

            if abs(hp1_ave - m) < abs(m - hp2_ave):
                return (self.ps, 1, abs(t))
            else:
                return (self.ps, 2, abs(t))
         

    class Block(object):
        def __init__(self, bf, blk, opts):
            self.bf = bf
            self.blk = blk
            self.opts = opts

        def meth_options(self):
            return 15, 0.3, 0.4

        def collect_meths(self):
            #logger.info("collect_meths(Block)")
            bf, blk = self.bf, self.blk
            meths = {} # { Meth } 
            unphased = set()
            #print(blk.start(), blk.end())
            block_alignment = bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True)

            for read in block_alignment:
               hp = read.get_tag("HP") if read.has_tag("HP") else 0
               ps = read.get_tag("PS") if read.has_tag("PS") else 0
               
               for rloc, m in MethBridger.collect_meths_in_read(read):
                   if rloc not in meths:
                       meths[rloc] = MethBridger.Meth(blk.chrom, rloc)
                   meths[rloc].add(read.qname, read.is_forward, m/256, ps, hp)
               if ps == 0:
                   unphased.add(read)
            return meths, unphased

        def extract_consistent_meths(self, meths, links):
            #logger.info("extract_consistent_meths(Block)")
            def extract_top(links, cands, rate) :
                count = defaultdict(lambda: [0,0])
                for k, v in links.items():
                    if k[0] not in cands or k[1] not in cands: continue
                    c = v[0] + v[3], v[1] + v[2]
                    r = c[0] / sum(c)
                    if abs(r - (1-r)) >= 0.6 and sum(c) >= 10:
                        count[k[0]][0] += 1
                        count[k[1]][0] += 1
                    count[k[0]][1] += 1
                    count[k[1]][1] += 1
                scount = sorted(count.items(), key=lambda x: -x[1][0]/(x[1][1]+1))
            
                ext = set()
                for x in scount:
                    ext.add(x[0])
                    if len(ext) >= len(scount)*rate:
                          break
                return ext
    
            valid_locs = set([k[0] for k in links]) | set([k[1] for k in links])
            total = len(valid_locs)
    
            while len(valid_locs) > total*0.15 and len(valid_locs) >= 10:
                valid_locs = extract_top(links, valid_locs, 0.5)
            
            return valid_locs

        def get_meths_linkscores(self, bf, blk, meths):
            logger.info("get_meths_linkscores(Block)")
            links = defaultdict(lambda: [0, 0, 0, 0])
            block_alignment = bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True)
            # counter_1 = 0
            # counter_2 = 0
            # counter_3 = 0
            start_time = time.time()
            for read in block_alignment:
               valid_meth = []
               #counter_3 += 1
               for rloc, m in MethBridger.collect_meths_in_read(read):
                   if rloc in meths and meths[rloc].valid():
                       ps, hp, d = meths[rloc].test(m/256)
                       valid_meth.append((rloc, hp))
                       #counter_1 += 1
               for (p0, h0), (p1, h1) in itertools.combinations(valid_meth, 2):
                    assert p0 < p1 and (h0 == 1 or h0 == 2) and (h1 == 1 or h1 == 2)
                    links[(p0, p1)][h0-1 + (h1-1)*2] += 1
                    #counter_2 += 1
            # print("get_meths_linkscores最外面的循环执行了%d次" % counter_3)
            # print("get_meths_linkscores第二层的第一个循环执行了%d次" % counter_1)
            # print("get_meths_linkscores第二层的第二个循环执行了%d次" % counter_2)
            end_time = time.time()
            run_time = end_time - start_time
            print("get_meths_linkscores run time: {:.2f} seconds".format(run_time))
            return links

        def compare_phased_meths(self, meths, valid_locs):
            #logger.info("compare_phased_meths(Block)")
            count = [0, 0, 0, 0]
            consist = {} 
            for pos, m in meths.items():
                if m.valid_phased():
                    count[0] += 1
                
                if m.valid(): 
                    count[1] += 1
    
                if m.valid_phased() and m.valid():
                    count[2] += 1
    
                if m.valid_phased() and pos in valid_locs:
                    count[3] += 1
            
                if pos in valid_locs or m.valid_phased():
                    consist[pos] = (pos in valid_locs, m.valid_phased())
            
            for p, c in sorted(consist.items()):
                print("CC", p, c)
            print(count, count[3] / len(valid_locs))

     
        def find_switch(self):
            logger.info("find_switch(Block)")
            bf, blk = self.bf, self.blk
            #block_alignment = bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True)
            meths = {} # Meth
            meths, unphased = self.collect_meths()
    
            for p, v in sorted(meths.items()):
                v.verify(*self.meth_options())
                v.verify_phased(*self.meth_options())

            #self.extend_phased_reads(meths, unphased)
    
            improvement = self.test_switch(meths)     
            for p, c in sorted(improvement.items()):
                imp = c[0]/c[1] if c[1] > 0 else 0
                # if True or imp < 0:
                #     print(p, c, imp)    
                    
            links = self.get_meths_linkscores(bf, blk, meths)
            valid_locs = self.extract_consistent_meths(meths, links)
            sorted(valid_locs)
            # for v in sorted(valid_locs):
            #     print("valid", v)
            #self.compare_phased_meths(meths, valid_locs)
    
            #improvement = self.test_switch(meths, valid_locs)     
            return self.detect_switch(improvement)
        
        def detect_switch(self, improvement):
            ranges = []
            rr = []
            for p, c in sorted(improvement.items()):
                imp = c[0]/c[1] if c[1] > 0 else 0
                if imp < 0:
                    rr.append((p, c)) 
                else:
                    if len(rr) > 0:
                        ranges.append(rr)
                        rr = []

            switchs = []            
            for rr in ranges:
                if len(rr) >= 3:
                    sorted_rr = sorted(rr, key = lambda x: -x[1][0])
                    for ir in sorted_rr:
                        print(ir)
                        print("ir")
                        if ir[1][0] <= -10 and ir[1][0] / ir[1][1] <= -0.03:
                            switchs.append(ir[0])
                            break
            print("detect_switch", switchs)
            return switchs

                

        def test_switch(self, meths, valid_locs=None):
            #logger.info("test_switch(Block)")
            bf, blk = self.bf, self.blk

            improvement = defaultdict(lambda: [0, 0, 0, 0])
            block_alignment = bf.fetch(blk.chrom, blk.start(), blk.end(), multiple_iterators=True)
            for read in block_alignment:
                hap = {}
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc in meths and meths[rloc].valid_phased() and (valid_locs==None or rloc in valid_locs):
                        (ps, hp, d) = meths[rloc].test_phased(m/256)
                        if d >= 0.3: # TODO
                            hap[rloc] = hp
                if len(hap) < 4: continue
    
                c = sum([1 for i in hap.values() if i == 1]), sum([1 for i in hap.values() if i == 2])
                sorted_hap = sorted(hap.items()) # by position
                left_c = [0, 0]
                for p, v in sorted_hap[0:-1]:
                    if v == 1:
                        left_c[0] += 1
                    elif v == 2:
                        left_c[1] += 1
    
                    right_c = [c[0] - left_c[0], c[1] - left_c[1]]
                    if sum(left_c) < 2 or sum(right_c) < 2: continue 
                    new_c = [left_c[0] + right_c[1], left_c[1] + right_c[0]] # switch
                    improvement[p][0] += min(new_c) - min(c)
                    improvement[p][1] += sum(new_c)
                    improvement[p][2] += sum(left_c)
                    improvement[p][3] += sum(right_c)
            
            return improvement

        def extend_phased_reads(self, meths, unphased):
            #logger.info("extend_phased_reads(Block)")  
            phased = set()
            while True:
                phasing = self.phase_reads(meths, unphased - phased)
                if len(phasing) > 0:
                    #logger.info(str(len(phasing)))
                    self.add_phasing_reads(meths, phasing)
                    [v.verify_phased(*self.meth_options()) for v in meths.values()]
                    [phased.add(r[0]) for r in phasing]
                else:
                    break
           

        def phase_reads(self, meths, unphased):
            #logger.info("phase_reads(Block)")
            phasing = []
            remove_abnormal = lambda x: x if type(x) ==int else 0
            # params
            sup_count, sup_rate = self.opts["read_count"], self.opts["read_rate"]
            for read in sorted(unphased, key=lambda x: remove_abnormal(x.reference_start)):
                checks = defaultdict(int)
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc in meths and meths[rloc].valid_phased():
                         (ps, hp, w) = meths[rloc].test_phased(m/256) 
                         if w >= 0.3: # TODO 检查  
                             checks[(ps, hp)] += 1
        
                sel = MethBridger.get_best_support(checks, sup_count, sup_rate)
                if sel != None:
                    phasing.append((read, sel))
            return phasing
        
        def add_phasing_reads(self, meths, phasing):
            for read, (ps, hp) in phasing:
           
                for rloc, m in MethBridger.collect_meths_in_read(read):
                    if rloc not in meths:
                        meths[rloc] = MethBridger.Meth(self.chrom, rloc) 
                    meths[rloc].add(read.qname, read.is_forward, m/256, ps, hp)
        

def main():
    import argparse
    parse = argparse.ArgumentParser(description='Bridge phased blocks.')
    parse.add_argument('vcf', help='The VCF file.')
    parse.add_argument('bam', help='The BAM file.')
    parse.add_argument("out_vcf", help="The output VCF file")
    parse.add_argument("-t", "--threads", help="number of threads", default=1, type=int)
    parse.add_argument("--str", help="")
    args = parse.parse_args()

    logger.info("Input VCF: %s" % args.vcf)
    logger.info("Input BAM: %s" % args.bam)
    logger.info("Output VCF: %s" % args.out_vcf)

    bridger = MethBridger(args.threads)
    #bridger = SVBridger()
    #bridger.bridge_by_str(args.vcf, args.bam, args.str, args.out_vcf) 
    #bridger.check_blocks(args.vcf, args.bam, args.out_vcf) 
    bridger.bridge2(args.vcf, args.bam, args.out_vcf)
    print("---------------------------------------------------")
    #bridger.find_switch(args.vcf, args.bam)


if __name__ == '__main__':
    main()
