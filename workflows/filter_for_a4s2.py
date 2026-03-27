#!/usr/bin/env python3

import os,sys,re,argparse
import pysam

parser = argparse.ArgumentParser(description="Filter BAM files for read consensus")
parser.add_argument('--bam',   action='store', default='input.bam', help='path of the BAM file (input.bam)')
parser.add_argument('--out',   action='store', default='output.bam', help='path of the output BAM file (output.bam)')
parser.add_argument('--status',   action='store', default='status.txt', help='path of the status file (status.txt)')
parser.add_argument('--min-cov-per-strand',   action='store', default=2, type=int, help='minimum number of reads per strand in a read bundle to be considered sufficient for consensus (default: 2, giving us a4s2)')
args = parser.parse_args()

print_read_out = False
current_read_bundle = None
current_read_bundle_orientation = {'pos_strand' : 0,
                                   'neg_strand' : 0}
track_readbundle = []
# track_nonReplicate_readbundle = []
total_bundles = 0
counter = 0

bam = pysam.AlignmentFile(args.bam, "rb")
outbam = pysam.AlignmentFile(args.out, "wb", template=bam)

for read in bam.fetch(until_eof=True):
    if current_read_bundle is None:
        # the very first read bundle is being processed, so we need to set the current read bundle to the read bundle of this read
        current_read_bundle = read.get_tag("RB")
        total_bundles += 1
    elif ( read.get_tag("RB") != current_read_bundle ):
        # at this stage, we have now reached a new read bundle. we need to check if the previous read bundle was a4s2 and if so, write it out. then we need to reset the read bundle and orientation information for the new read bundle
        if ( print_read_out == True ):
            # write the reads in this bundle
            for i in track_readbundle:
                    outbam.write(i)
            counter += 1 # we have a new number of a4s2 read bundles
            
            # print a status check every 1000 a4s2 read bundles processed
            if ( counter % 1000 == 0 ):
                print("%d out of %d read bundles processed are a4s2; %f percent"%(counter, total_bundles, counter/total_bundles*100), end='\n',file=open(args.status,'a'))
            
            print_read_out = False # reset the print read out flag for the new read bundle

        track_readbundle = [] # empty out the read bundles; we will start tracking the new read bundle
        current_read_bundle = read.get_tag("RB") # reset the read bundle to the new read bundle
        # print("starting with read bundle %s"%(current_read_bundle))
        current_read_bundle_orientation = {'pos_strand' : 0,
                                        'neg_strand' : 0}
        total_bundles += 1 # we have a new read bundle, so we add to the total number of read bundles that we have processed
    
    # disqualify read
    if ( read.is_unmapped ) : continue # we don't want reads in a bundle that are unmapped
    if ( read.is_secondary ) : continue # we don't want reads in a bundle that belong to a secondary alignment
    if ( read.is_supplementary ) : continue # we don't want reads in a bundle that belong to a supplementary alignment
    if ( read.is_qcfail ) : continue # we don't want reads that fail quality control
    if ( read.mapping_quality < 20 ) : continue # we don't want reads with low mapping quality
    # if ( read.is_paired and read.is_proper_pair == False ) : continue # we don't want reads not mapped in a proper pair
    # isolate the CIGAR flag for this read -- we want to look for whether the read in this bundle comes from positive and negative. we will monitor if we have reached a4s2
    breakdown_flag = [int(j)*int(2**i) for i,j in enumerate(bin(read.flag)[2:][::-1]) if int(j) > 0]
    # monitor the first-in-strand read in the bundle and add to the count of reads in the bundle that come from the positive and negative strand. we will monitor if we have reached a4s2
    if ( 64 in breakdown_flag and 16 in breakdown_flag):
        current_read_bundle_orientation['neg_strand'] += 1
    elif ( 64 in breakdown_flag and 32 in breakdown_flag):
        current_read_bundle_orientation['pos_strand'] += 1
    
    # add the read
    track_readbundle.append(read)
    
    # now check if we have reached a4s2
    if ( current_read_bundle_orientation['pos_strand'] >= args.min_cov_per_strand and current_read_bundle_orientation['neg_strand'] >= args.min_cov_per_strand ) :
        print_read_out = True


bam.close()
outbam.close()
print("now sorting and indexing the output BAM file\n")
sorted_out = '.'.join(['sorted', args.out])
pysam.sort("-o", sorted_out, args.out)
pysam.index(sorted_out)
print("completed sorting and indexing the output BAM file\n")
print("done")

# print("now filtering the BAM to just one non-duplicate read per bundle")
# sorted_out_nondup = '.'.join(['sorted_nonDup', args.out])
# pysam.view("-h", "-F", "1024", "-o", sorted_out_nondup, sorted_out)
# pysam.index(sorted_out_nondup)
# print("completed filtering the BAM to just one non-duplicate read per bundle\n")