import argparse, os
from Bio import SeqIO
from Bio.SeqUtils import GC
import multiprocessing
from multiprocessing import Pool, Manager
import pandas as pd
import pyranges as pr
import numpy as np

def argument_parser():
    parser = argparse.ArgumentParser()    
    parser.add_argument('-i', '--inputFragPath', required=True, help="path of input fragment file")
    parser.add_argument('-o', '--output', required=True, help='outputDir')
    parser.add_argument('-c', '--cores', required=False, default = 500, help='Number of core will be used')
    
    return parser

def get_bioGC(q_chr, q_st, q_end):
    seq = bio_ref_chr[q_chr].seq[q_st:q_end]
    seq_str = str(seq)
    st_motif = seq_str[:4]
    ed_motif = seq_str[-4:]
    return GC(seq), st_motif.upper(), ed_motif.upper()

def run_part_v2(indices, frag_df_dict, OUTPUT_PATH):
    
    for idx in indices:
#       print(idx)
        frag_df = frag_df_dict[idx]
        gc_col={}
        st_motif_col={}
        end_motif_col={}
        for x in frag_df.index:
            c, s, e, _, _, _, _ = frag_df.loc[x]

            o_gc, o_sm, o_em = get_bioGC(c,s,e)

            gc_col[x]=o_gc
            st_motif_col[x]=o_sm
            end_motif_col[x]=o_em
            
        frag_df['GC']=pd.Series(gc_col)
        frag_df['start_motif']=pd.Series(st_motif_col)
        frag_df['end_motif']=pd.Series(end_motif_col)

        #PAULO    
        #frag_df.to_csv(OUTPUT_PATH+'_batch_%s'%idx)
        filename = '_batch_' + str(idx) + '.csv'
        frag_df.to_csv(filename)

def main():

    parser = argument_parser()
    options = parser.parse_args()

    frag_path = options.inputFragPath
    output_dir = options.output
    n_cpu = int(options.cores)
    sample_id = frag_path.split('/')[-1].split('.be')[0]
    
    bio_ref = SeqIO.parse('./hg19.fa', "fasta")
    global bio_ref_chr
    bio_ref_chr = {}
    for c in bio_ref:
        if not c.id in bio_ref_chr:
            bio_ref_chr[c.id]=c
        else:
            print(c.id)
#     print('ref load done')
    
    output_dir_sample = output_dir+'/'+sample_id
    #if not os.path.isdir(output_dir_sample):
        #PAULO
        #os.mkdir(output_dir_sample)
    OUTPUT_PATH = output_dir_sample+'/%s@'%sample_id
    
    test_bed = pr.read_bed(frag_path)
    test_bed.frag_width = test_bed.End-test_bed.Start
    test_bed = test_bed[test_bed.frag_width<1000]
    test_bed_df = test_bed.df
    
    fs = len(test_bed_df)
    one_size = int(fs/n_cpu)
    SPLITTED_FRAGS = {}
    for i in range(n_cpu):
        sub = test_bed_df.iloc[(i)*one_size:(i+1)*one_size]
        SPLITTED_FRAGS[i]=sub
        
    sub = test_bed_df.iloc[(i+1)*one_size:]
    SPLITTED_FRAGS[i+1]=sub
#     print('run_multi')
    pool = multiprocessing.Pool(n_cpu)
    m = Manager()
    all_paths = np.array_split(list(SPLITTED_FRAGS), n_cpu)

    pool.starmap(run_part_v2, [(XXX, {zz:SPLITTED_FRAGS[zz] for zz in XXX}, OUTPUT_PATH) for XXX in all_paths])
    pool.close()
    pool.join()

    
main()
