#!/home/seqproc/illumina_processing/env/bin/python

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
import smtplib
import ConfigParser
import pandas as pd
import glob


def main():
    genome_index_dict = {
        'bwamem': {
            'mm9': 'tbd',
            'mm10': 'tbd',
            'hg19': '/mnt/hiseq2/genomes/hg19/UCSC/hg19/Sequence/BWAIndex/version0.6.0/genome.fa',
            'hg38': 'tbd',
        },
        'tophat2': {
            'mm9': 'tbd',
            'mm10': 'tbd',
            'hg19': '/mnt/hiseq2/genomes/hg19/UCSC/hg19/Sequence/Bowtie2Index/genome',
            'hg38': 'tbd',
        }
    }
    config = get_config('pathway.cfg')
    logger = initialize_logger(config.get('Globals', 'LogFile'))

    lock_path = '/home/seqproc/illumina_processing/alignment_lock'
    if os.path.exists(lock_path):
        logger.info('Alignment is already running. Skipping.')
        sys.exit(0)
    open(lock_path, 'w').close()  # Set alignment-already-running lock.

    to_be_aligned_fname = config.get('Globals', 'ToBeAlignedList')
    proj_dirs_to_run = find_dirs_to_run(
        to_be_aligned_fname, lock_path, logger
    )
    proj_dirs_not_done = []
    for project_dir in proj_dirs_to_run:
        logger.info('Processing %s' % project_dir)
        proj_done = True
        for sample_path in glob.glob(os.path.join(project_dir, 'Sample*')):
            logger.info(
                'Processing "%s" in project directory "%s"' %
                (os.path.basename(sample_path), project_dir)
            )
            exp_type, ref_genome = read_in_sample_sheet(
                sample_path, ['dnaseq', 'rnaseq'],
                genome_index_dict['bwamem'].keys(),
                lock_path, config, logger
            )
            aligner = 'bwamem' if exp_type == 'dnaseq' else 'tophat2'
            index_path = get_genome_index(
                genome_index_dict, aligner, ref_genome, lock_path, logger
            )
            fastq_pairs = get_fastq_pairs(
                sample_path, lock_path, logger
            ) # list of (R1, R2) or (R1,) tuples
            success = run_alignments(
                fastq_pairs, aligner, sample_path, ref_genome,
                index_path, lock_path, logger
            )
            if not success:
                proj_done = False
        if not proj_done:
            proj_dirs_not_done.append(project_dir)
        logger.info('Done processing %s' % project_dir)

    with open(to_be_aligned_fname, 'w') as f:
        f.write('\n'.join(proj_dirs_not_done))
    if os.path.exists(lock_path):
        os.remove(lock_path)
    logger.info('Done processing.')


def raise_error(error_type, message, lock_path, logger):
    logger.error(message)
    if os.path.exists(lock_path):
        os.remove(lock_path)
    raise error_type(message)


def get_config(default_fname):
    config = ConfigParser.ConfigParser()
    if len(sys.argv) == 2:
        config.readfp(open(sys.argv[1]))
    else:
        config.readfp(open(default_fname))
    return config


def initialize_logger(log_filename):
    logger = logging.getLogger(sys.argv[0])
    log_fh = logging.FileHandler(log_filename)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    log_fh.setFormatter(formatter)
    logger.addHandler(log_fh)
    logger.setLevel(logging.INFO)
    logger.info('Starting...')
    return logger


def find_dirs_to_run(dir_list, lock_path, logger):
    try:
        with open(dir_list) as f:
            return [z.strip() for z in f.readlines()]
    except IOError:
        logger.error(
            'Cannot find list of directories to be aligned. '
            'Check ToBeAlignedList in config file.'
        )
        if os.path.exists(lock_path):
            os.remove(lock_path)
        raise


def get_fastq_pairs(sample_path, lock_path, logger):
    if not os.path.exists(sample_path):
        message = 'Path to fastq files (%s) doesn\'t exist.' % sample_path
        raise_error(IOError, message, lock_path, logger)

    fastq_list = [
        f for f in os.listdir(sample_path) if (
            f.endswith(('fastq.gz', 'fq.gz', 'fastq', '.fq'))
        )
    ]
    fastq_pairs = []
    for filename in fastq_list:
        if '_R1' in filename:
            r2_filename = filename.replace('_R1', '_R2')
            if r2_filename in fastq_list:
                fastq_pairs.append((filename, r2_filename))
            else:
                fastq_pairs.append((filename,))
        elif '_R2' in filename:
            r1_filename = filename.replace('_R2', '_R1')
            if r1_filename not in fastq_list:
                fastq_pairs.append((filename,))
        else:
            logger.warning(
                'File %s does not have R1 or R2 in name. '
                'Will do single-end alignment.'
            )
            fastq_pairs.append((filename,))
    logger.info('FASTQ pairs found: %s' % str(fastq_pairs))
    return fastq_pairs


def read_in_sample_sheet(
    sample_path, possible_exp_types, possible_ref_genomes, lock_path,
    config, logger
):
    try:
        table = pd.read_csv(
            os.path.join(sample_path, 'SampleSheet.csv'), header=0
        )
    except IOError: 
        message = 'Could not find sample sheet in %s' % sample_path
        logger.error(message)
        if os.path.exists(lock_path):
            os.remove(lock_path)
            send_email(
                config.get('Globals', 'DebugEmailRecipient'), message, config
            )
        raise

    exp_type = get_sample_sheet_value(
        table, 'experiment type', 'Description', possible_exp_types,
        lock_path, config, logger
    )
    ref_genome = get_sample_sheet_value(
        table, 'reference genome', 'SampleRef', possible_ref_genomes,
        lock_path, config, logger
    )
    return (exp_type, ref_genome)


def get_sample_sheet_value(
    table, category, sample_sheet_field, possible_values, lock_path,
    config, logger
):
    value = str(table.loc[0, sample_sheet_field]).strip().lower()
    logger.info('%s: %s' % (category, value))
    if value not in possible_values:
        message = (
            'Unrecognized %s: "%s". Check that %s field in sample sheet '
            'has %s.' % (
                category, value, sample_sheet_field,
                ' or '.join(['"%s"' % z for z in possible_values])
            )
        )
        send_email(
            config.get('Globals', 'DebugEmailRecipient'), message, config
        )
        raise_error(ValueError, message, lock_path, logger)
    return value


def get_genome_index(
    genome_index_dict, aligner, ref_genome, lock_path, logger
):
    # Return path to genome index.
    if aligner not in genome_index_dict:
        message = (
            'Unrecognized aligner. Check experiment type in sample sheet.'
        )
        raise_error(ValueError, message, lock_path, logger)
    if ref_genome not in genome_index_dict[aligner]:
        message = (
            'Unrecognized reference genome. '
            'Check reference genome in sample sheet.'
        )
        raise_error(ValueError, message, lock_path, logger)
    return genome_index_dict[aligner][ref_genome]


def run_alignments(
    fastq_pairs, aligner, sample_path, ref_genome, index_path,
    lock_path, logger
):
    if aligner not in ['bwamem', 'tophat2']:
        message = (
            'Unrecognized aligner. Check experiment type in sample sheet.'
        )
        raise_error(ValueError, message, lock_path, logger)

    all_successful = True
    for fastq_pair in fastq_pairs:
        if len(fastq_pair) in [1, 2]:
            is_paired = len(fastq_pair) == 2
            ret = 1
            logger.info(
                'Running %s on %s-end reads %s; aligning to ref %s' % (
                    aligner,
                    'paired' if is_paired else 'single',
                    ', '.join(fastq_pair),
                    index_path
                )
            )
            # Run alignment.
            outfile_path = get_outfile_path(
                fastq_pair, aligner, ref_genome, sample_path, logger
            )
            if aligner == 'bwamem':
                logger.info('Will output alignment to %s' % outfile_path)
                if os.path.exists(outfile_path):
                    logger.warning(
                        '%s already exists. Skipping. ' % outfile_path
                    )
                    continue
                ret = run_bwamem(
                    sample_path, index_path, fastq_pair, outfile_path,
                    is_paired, logger
                )
            elif aligner == 'tophat2':
                outdir_path = outfile_path.strip('.sam')
                logger.info('Will output alignment to %s' % outdir_path)
                if os.path.exists(outdir_path):
                    logger.warning(
                        '%s already exists. Skipping. ' % outdir_path
                    )
                    continue
                ret = run_tophat2(
                    sample_path, index_path, fastq_pair, outdir_path,
                    is_paired, logger
                )
            # Check return code.
            if not ret:  
                logger.info(
                    '%s alignment completed successfully.' % aligner
                )
            else:
                logger.warning(
                    'Error in %s alignment. Return code: %d' % (aligner, ret)
                )
                all_successful = False
        else:
            logger.warning(
                'Invalid fastq pair length. Ignoring pair. Pair: %s.' %
                fastq_pair
            )
    return all_successful


def get_outfile_path(fastq_pair, aligner, ref_genome, sample_path, logger):
    first_fname = fastq_pair[0].rstrip('.gz')
    outfile_name = None
    for substr in ['_R1', '_R2', '.']:
        if substr in first_fname:
            outfile_name = first_fname[:first_fname.rindex(substr)]
            break
    if not outfile_name:
        outfile_name = first_fname
    outfile_path = os.path.join(
        sample_path, '%s.%s_%s.sam' % (outfile_name, aligner, ref_genome)
    )
    return outfile_path


def send_email(recipient, content, config):
    """Send email with given content from EmailSender (in pathway.cfg) to
    recipient.
    """
    server = smtplib.SMTP('smtp.gmail.com:587')
    server.ehlo()
    server.starttls()
    sender = config.get('Globals', 'EmailSender')
    msg = '\r\n'.join([
        'From: %s' % sender, 'To: %s' % recipient,
        'Subject: Sample sheet error', '', content
    ])
    server.login(sender, config.get('Globals', 'EmailPassword'))
    server.sendmail(sender, recipient, msg)
    server.quit()


def run_bwamem(
    sample_path, index_path, fastqs, outfile_path, is_paired, logger
):
    """Run BWA MEM on fastq files. BWA-MEM doesn't offer an option for
    writing output to file, so have to use context manager to do it
    in python.
    """
    # bwa must be on path
    args = [
        'bwa', 'mem', '-a', '-M', '-t', '4', index_path,
        os.path.join(sample_path, fastqs[0])
    ]
    if is_paired:
        args.append(os.path.join(sample_path, fastqs[1]))

    # context manager to redirect stdout to file handle
    with open(outfile_path, 'w') as outfile:  
        return subprocess.call(args, stdout=outfile)


def run_tophat2(
    sample_path, index_path, fastqs, outdir_path, is_paired, logger
):
    # tophat2 must be on path
    args = [
        'tophat2', '--num-threads', '10', '--mate-inner-dist', '200',
        '--max-multihits', '1', '--splice-mismatches', '1',
        '--output-dir', outdir_path,
        index_path, os.path.join(sample_path, fastqs[0])
    ]
    if is_paired:
        args.append(os.path.join(sample_path, fastqs[1]))
    return subprocess.call(args)


def run_bowtie2(
    sample_path, index_path, fastqs, outfile_path, is_paired, logger
):
    """Run bowtie2 on fastq files."""
    # bowtie2 must be on path
    args = ['bowtie2', '-x', index_path]
    if is_paired:
        args += [
            '-1', os.path.join(sample_path, fastq_r1),
            '-2', os.path.join(sample_path, fastq_r2)
        ]
    else:
        args += ['-U', os.path.join(sample_path, fastq_r1)]
    args += [
        '--phred33', '--sensitive', '--threads', '8', '-S', outfile_path
    ]
    return subprocess.call(args)


if __name__ == '__main__':
    main()
