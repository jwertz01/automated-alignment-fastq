#!/usr/bin/python

"""Run alignment against a reference genome on FASTQ files from 
recently converted bcltofastq diretories.

The program reads the "SampleRef" and "Description" fields of the SampleSheet
CSV file.  

From these, it determines the reference genome to align against, and the type
of experiment (DNAseq or RNAseq).  

If DNAseq, then either BWA-MEM (default) or Bowtie2 is invoked on the reads. 
If RNASeq, then Tophat2 is invoked on the reads. 

The functions for calling Bowtie2 on the reads are present below but not used in 
the script as written.  The script will have to be edited to use Bowtie2 instead 
of BWA-MEM.  
"""

import os
import sys
import subprocess
import logging
import ConfigParser
import pandas as pd
import shlex


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
            'Cannot find list of directories to be aligned'
            'Check aligner.cfg'
        )
        raise
        

    for project in to_be_aligned:
        project_dir = os.path.join(base_dir, project)
        sample_names = [sample for sample in os.listdir(project_dir)
                        if os.path.isdir(os.path.join(project_dir, sample))]
        for sample in sample_names:
            logger.info(
                '***Working on "%s" in project directory: %s' % (sample, project_dir)
            )
            exp_type = dnaseq_or_rnaseq(project_dir, sample, logger)
            
            if exp_type == 'dnaseq':
                paired_status, fastq = single_or_paired(project_dir, sample, logger)
                ref_genome = get_reference_genome(project_dir, sample, logger)
                index_fa_name = get_genome_index(ref_genome, logger)
                if not paired_status:
                    run_bwamem_single(project_dir, sample, index_fa_name, fastq[0], logger)    
                else: 
                    run_bwamem_paired(project_dir, sample, index_fa_name, fastq[0], fastq[1], logger) 
                      
            if exp_type == 'rnaseq':
                paired_status, fastq = single_or_paired(project_dir, sample, logger)
                ref_genome = get_reference_genome(project_dir, sample, logger)
                index_basename = get_genome_index(ref_genome, logger)
                if not paired_status:
                    run_tophat2_single(project_dir, sample, index_basename, fastq[0], logger)    
                else: 
                    run_tophat2_paired(project_dir, sample, index_basename, fastq[0], fastq[1], logger) 
                      
            else:
                continue


def single_or_paired(project_dir, sample, logger):
    paired = True
    sample_path = os.path.join(project_dir, sample)
    if os.path.exists(sample_path):
        files = os.listdir(sample_path)
        fastq_list = [f for f in files if f.endswith(('fastq.gz', 'fq.gz', 'fastq', '.fq'))] #catch all types of fastq files
        if len(fastq_list) == 2:
            logger.info('***Found paired end fastq files in directory %s' % sample_path)
            logger.info('***Found fastq files, %s,%s' % (fastq_list[0], fastq_list[1]))
            return paired, fastq_list
        if len(fastq_list) == 1:
            logger.info('***Found fastq file in directory: %s' % sample_path)
            logger.info('***Single end fastq file: %s' % fastq_list[0])
            
            paired = False
            return paired, fastq_list
        else: 
            logger.info('***%d Fastq Files found; Sample Dir %s; No alignment done' % (len(fastq_list), sample_path))
            return None
            
            
def get_reference_genome(project_dir, sample, logger):
    sample_path = os.path.join(project_dir, sample)
    try:
        tab = pd.read_csv(sample_path + '/SampleSheet.csv', header = 0)
        ref_genome = tab.loc[0,'SampleRef']
        logger.info('***"%s" reference genome code is: %s' % (sample, ref_genome))
        return ref_genome.strip().lower()
            
    except IOError: 
        logger.info('***Could not find sample sheet or ref genome code for %s' % sample_path)
        return None
        
        
def get_genome_index(ref_genome, logger):
    try:
        gen_dict = {'mm9': '/mnt/hiseq2/genomes/mm9',
                'mm10': '/mnt/hiseq2/genomes/mm10',
                'hg19': '/mnt/hiseq2/genomes/hg19',
                'hg38': '/mnt/hiseq2/genomes/hg38',
                'lmbda': '/Applications/bwa-0.7.12/example/lambda_virus.fa', #for testing 
                'unk' : '/Applications/tophat-2.1.0.OSX_x86_64/test_data/test_ref'  #for testing
                } 
                
        index_basename = gen_dict[ref_genome]
        return index_basename
    except KeyError:
        logger.info("Invalid reference genome.")
        return None
    
    
def dnaseq_or_rnaseq(project_dir, sample, logger):
    sample_path = os.path.join(project_dir, sample)
    try:
        tab = pd.read_csv(sample_path + '/SampleSheet.csv', header = 0)
        exp_type = tab.loc[0,'Description']
        logger.info('***"%s" is a %s type of experiment' % (sample, exp_type))
        return exp_type.strip().lower()
            
    except IOError: 
        logger.info('***Could not find sample sheet in %s' % sample_path)
        return None
    
    except AttributeError:
        logger.info('***No experiment type field in the Sample Sheet in %s' % sample_path)
        return None


def run_bwamem_paired(project_dir, sample, index_fa_name, fastq_r1, fastq_r2, logger):
    '''Run BWA MEM in paired end mode for fastq files.  BWA-MEM doesn't offer an option 
    for writing output to file, so have to use context manager to do it in python
    '''
    logger.info('***Running BWA-MEM on paired-end reads; aligned to ref %s' % index_fa_name)
    filepath = os.path.join(project_dir, sample)
    
    args = [
            '/Applications/bwa-0.7.12/bwa',
            'mem',
            '-a',
            '-M',
            '-t', '4',
            index_fa_name,
            os.path.join(filepath, fastq_r1),
            os.path.join(filepath, fastq_r2)
        ]
    
    ## context manager to redirect stdout to file handle
    with open(os.path.join(filepath, 'bwamem_paired_align.sam'), 'w') as outfile:  
        bwamem_process = subprocess.call(args, stdout = outfile)
   
    if not bwamem_process:  #return code of 0 is success
        logger.info(
            '***BWA MEM alignment completed successfully for %s' % filepath
        )
    else:
        logger.info(
            '***Error in BWA-MEM paired alignment. Return code: %d' % bwamem_process
        )
    
    
def run_bwamem_single(project_dir, sample, index_fa_name, fastq_r1, logger):
    '''Run BWA MEM in single end mode for fastq files.
    '''
    logger.info('***Running BWA-MEM on single-end reads; aligned to ref %s' % index_fa_name)
    filepath = os.path.join(project_dir, sample)
    args = [
            '/Applications/bwa-0.7.12/bwa',
            'mem',
            '-a',
            '-M',
            '-t', '4',
            index_fa_name, 
            os.path.join(filepath, fastq_r1)
        ]
        
    ## context manager to redirect stdout to file handle
    with open(os.path.join(filepath, 'bwamem_align.sam'), 'w') as outfile:  
        bwamem_process = subprocess.call(args, stdout=outfile)
    
    if not bwamem_process:  #return code of 0 is success
        logger.info(
            '***BWA MEM alignment completed successfully for %s' % filepath
        )
    else:
        logger.info(
            '***Error in BWA-MEM single alignment. Return code: %d' % bwamem_process
        )
   
     
def run_tophat2_paired(project_dir, sample, index_basename, fastq_r1, fastq_r2, logger):
    '''Run tophat2 in paired-end mode for fastq files. 
    '''
    logger.info('***Running tophat2 on paired-end reads; aligned to ref %s' % index_basename)
    filepath = os.path.join(project_dir, sample)
    args = [
        '/Applications/tophat-2.1.0.OSX_x86_64/tophat2',  #temp path for testing
        #'/home/seqproc/tophat2/'
        '--num-threads','10',
        '--mate-inner-dist','200',
        '--max-multihits' ,'1',
        '--splice-mismatches', '1',
        index_basename,
        os.path.join(filepath, fastq_r1),
        os.path.join(filepath, fastq_r2)
    ]
    
    print subprocess.list2cmdline(args)
    top2_process = subprocess.call(args)
    
    if not top2_process:  #return code of 0 is success
        logger.info(
            '***Bowtie2 alignment completed successfully for %s' % filepath
        )
        
    else:
        logger.info(
            '***Error in bowtie2 alignment. Return code: %d' % top2_process
        )
    
    
def run_tophat2_single(project_dir, sample, index_basename, fastq_r1, logger):
    '''Run tophat2 in single-end mode for fastq files. 
    '''
    logger.info('***Running tophat2 on single-end reads; aligned to ref %s' % index_basename)
    filepath = os.path.join(project_dir, sample)
    args = [
        '/Applications/tophat2/tophat2',  #temp path for testing
        #'/home/seqproc/tophat2/
        '--num-threads','10',
        '--mate-inner-dist','200',
        '--max-multihits' ,'1',
        '--splice-mismatches', '1',
        index_basename,
        os.path.join(filepath, fastq_r1)
    ]
    
    top2_process = subprocess.call(args)
    
    if not top2_process:  #return code of 0 is success
        logger.info(
            '***Bowtie2 alignment completed successfully for %s' % filepath
        )
        
    else:
        logger.info(
            '***Error in bowtie2 alignment. Return code: %d' % top2_process
        )
    
    
def run_bowtie2_paired(project_dir, sample, index_base_name, fastq_r1, fastq_r2, logger):
    '''Run bowtie2 in paired end mode for fastq files. 
    '''
    logger.info('***Running bowtie2 on paired-end reads; aligned to ref %s' % index_base_name)
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
        '-S', os.path.join(filepath, 'test_paired.sam')
    ]
    
    bwt2_process = subprocess.call(args)

    if not bwt2_process:  #return code of 0 is success
        logger.info(
            '***Bowtie2 alignment completed successfully for %s' % filepath
        )
    else:
        logger.info(
            '***Error in bowtie2 alignment? Return code: %d' % bwt2_process
        )


def run_bowtie2_single(project_dir, sample, index_base_name, fastq_r1, logger):
    '''Run bowtie2 in single-end mode for fastq files. 
    '''
    logger.info('***Running bowtie2 on single-end reads; aligned to ref %s' % index_base_name)
    filepath = os.path.join(project_dir, sample)
    args = [
        '/Applications/bowtie2-2.2.6/bowtie2',  #temp path for testing
        #'/home/seqproc/bowtie2-2.2.6/bowtie2',
        '-x', index_base_name,
        '-U', os.path.join(filepath, fastq_r1),
        '--phred33',
        '--sensitive',
        '--threads 8',
        '-S', os.path.join(filepath, 'test_single.sam')
    ]
    
    bwt2_process = subprocess.call(args)
    
    if not bwt2_process:  #return code of 0 is success
        logger.info(
            '***Bowtie2 alignment completed successfully for %s' % filepath
        )
        
    else:
        logger.info(
            '***Error in bowtie2 alignment. Return code: %d' % bwt2_process
        )
 
    
if __name__ == '__main__':
    main()
