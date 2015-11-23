#!/usr/bin/python

"""Run bowtie2 alignment against a reference genome on FASTQ files from 
recently converted bcltofastq diretories.
"""

import os
import sys
import subprocess
import logging
import ConfigParser
import pandas as pd


def main():
    # Open config file.
    config = ConfigParser.ConfigParser()
    if len(sys.argv) == 2:
        config.readfp(open(sys.argv[1]))
    else:
        config.readfp(open('aligner.cfg'))
        
    # Create logger.
    logger = logging.getLogger(sys.argv[0])
    log_f_name = logging.FileHandler(config.get('Globals', 'LogFile'))
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    log_f_name.setFormatter(formatter)
    logger.addHandler(log_f_name)
    logger.setLevel(logging.INFO)
    logger.info('Starting...')

    # Read in values from .cfg file (align.cfg unless otherwise specified).
    base_dir = config.get('Globals', 'InputSeqDirectory')
    to_be_aligned_dirs = config.get('Globals', 'ToBeAlignedList')
    
    try:
        with open(to_be_aligned_dirs) as f:
            to_be_aligned = [z.strip() for z in f.readlines()]
    
    except IOError:
        logger.error(
            'Cannot find list of directories to be aligned.'
            'Check aligner.cfg')
        
    
    #try:
    for project in to_be_aligned:
        project_dir = os.path.join(base_dir, project)
        sample_names = [sample for sample in os.listdir(project_dir)
                        if os.path.isdir(os.path.join(project_dir, sample))]
        for sample in sample_names:
            paired_status, fastq = single_or_paired(project_dir, sample, logger)
            ref_genome = get_reference_genome(project_dir, sample, logger)
            index_basename = get_genome_index(ref_genome)
            if not paired_status:
                run_bowtie2_single(project_dir, sample, index_basename, fastq, logger)    
            else: 
                run_bowtie2_paired(project_dir, sample, index_basename, fastq[0], fastq[1], logger)   

    #except: pass
    
def single_or_paired(project_dir, sample, logger):
    paired = True
    sample_path = os.path.join(project_dir, sample)
    if os.path.exists(sample_path):
        files = os.listdir(sample_path)
        fastq = [f for f in files if f.endswith('fastq.gz') or f.endswith('fastq')]
        print fastq
        if len(fastq) == 2:
            logger.info('Found paired end fastq files, %s.' % fastq)
            return paired, fastq
        if len(fastq) == 1:
            logger.info('Found single end fastq file, %s.' % fastq)
            paired = False
            return paired, fastq
        else: 
            logger.info('%d Fastq Files found; Sample Dir %s; No alignment done.' % (len(fastq), sample_path))
            return None
            
            
def get_reference_genome(project_dir, sample, logger):
    sample_path = os.path.join(project_dir, sample)
    if os.path.exists(sample_path):
        tab = pd.read_csv(sample_path + '/SampleSheet.csv', header = 0)
        ref_genome = tab.loc[0,'SampleRef']
        logger.info('%s reference genome code is %s.' % (sample, ref_genome))
        return ref_genome
            
    else: 
        logger.info('Could not find sample or sample sheet for %s.' % sample_path)
        return None
        
def get_genome_index(ref_genome):
    gen_dict = {'mm9': '/mnt/hiseq2/genomes/mm9',
                'mm10': '/mnt/hiseq2/genomes/mm10',
                'hg19': '/mnt/hiseq2/genomes/hg19',
                'hg38': '/mnt/hiseq2/genomes/hg38',
                'LMBDA': '/Applications/bowtie2-2.2.6/example/index/lambda_virus'}  #for test, remove for production
    index_basename = gen_dict[ref_genome]
    return index_basename

    
    
def run_bowtie2_paired(project_dir, sample, index_base_name, fastq_r1, fastq_r2, logger):
    '''Run bowtie2 in paired end mode for fastq files. 
    Return path to bowtie2 output.
    '''
    logger.info('Running bowtie2 on paired-end reads; aligned to ref %s.' % index_base_name)
    filepath = os.path.join(project_dir, sample)
    args = [
        '/Applications/bowtie2-2.2.6/bowtie2',
        #'/home/seqproc/bowtie2-2.2.6/bowtie2',
        '-x', index_base_name,
        '-1', os.path.join(filepath, fastq_r1),
        '-2', os.path.join(filepath, fastq_r2),
        '--phred33',
        '--sensitive',
        '--threads 8',
        '-S', os.path.join(filepath, 'test_paired.out')
    ]
    
    bwt2_process = subprocess.call(args)

    if not bwt2_process:  #return code of 0 is success
        logger.info(
            'bowtie2 alignment completed successfully for %s' % filepath
        )
    else:
        logger.info(
            'Error in bowtie2 alignment? Return code: %d.' % bwt2_process
        )
    
    return bwt2_process


def run_bowtie2_single(project_dir, sample, index_base_name, fastq, logger):
    '''Run bowtie2 in single-end mode for fastq files. 
    Return path to bowtie2 output.
    '''
    logger.info('Running bowtie2 on single-end reads; aligned to ref %s.' % index_base_name)
    filepath = os.path.join(project_dir, sample)
    args = [
        '/Applications/bowtie2-2.2.6/bowtie2',  #temp path for testing
        #'/home/seqproc/bowtie2-2.2.6/bowtie2',
        '-x', index_base_name,
        '-U', os.path.join(filepath, fastq[0]),
        '--phred33',
        '--sensitive',
        '--threads 8',
        '-S', os.path.join(filepath, 'test_single.out')
    ]
    
    bwt2_process = subprocess.call(args)
    if not bwt2_process:  #return code of 0 is success
        logger.info(
            'bowtie2 alignment completed successfully for %s' % filepath
        )
    else:
        logger.info(
            'Error in bowtie2 alignment. Return code: %d.' % bwt2_process
        )
    
    return bwt2_process


if __name__ == '__main__':
    main()
