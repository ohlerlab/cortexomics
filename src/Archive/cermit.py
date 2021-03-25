from collections import Counter
from itertools import product
from Bio import SeqIO
import glob
import shutil
import os
import subprocess
from pybedtools import BedTool
import pandas as pd
import argparse
import sys

# parse options from command line
parser = argparse.ArgumentParser()

parser.add_argument('-bed', dest='bed',
  type=str,
  help="Input regions",
  default='/fast/AG_Akalin/wkopp/alison_sciatac_2018/scipipe_v4/state_calls_90_15_32_v3/DE_peaks_jamm/cluster_16.narrowPeak')

parser.add_argument('-fasta', dest='fasta',
  type=str,
  help="Input sequence")

parser.add_argument('-genome', dest='genome',
 type=str,
 help='Genome sequence',
 default='/fast/AG_Ohler/Scott/danRer11/noAlt_danRer11.fa')

parser.add_argument('-output', dest='output', 
  type=str,
  help='Output folder',
default='/fast/AG_Akalin/wkopp/alison_sciatac_2018/scipipe_v4/state_calls_90_15_32_v3/cermit_DE_peaks_jamm/cluster_16_results')

parser.add_argument('-scorefield', dest='scorefield', type=int, help='Score field index. Default=6', default=6)
parser.add_argument('-ntop', dest='ntop', type=int, help='Number of regions. Default=sys.maxsize', default=sys.maxsize)

parser.add_argument('-motifdb', dest='motifdb', type=str, help='Motif database for tomtom', default='/fast/AG_Akalin/wkopp/alison_sciatac_2018/extra/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme')
parser.add_argument('-no_tomtom', dest='no_tomtom', action="store_true", help='Compare motifs to database with tomtom')

# program args
parser.add_argument('-prog', dest='prog', type=str, help='cERMIT path',
                   default='/gnu/store/i6s2cavvb21mxkm46i13l2q6d0z30zag-cermit-1.11.1/bin/cERMIT')

parser.add_argument('-prefix', dest='prefix', type=str, help='Prefix for motif name. Default: "input". ', default='input')
parser.add_argument('-seed_length', dest='seed_length', type=int, help='seed length ', default = 6)
parser.add_argument('-min_motif_length', dest='min_motif_length', type=int, help='min motif length. Default=6', default = 6)
parser.add_argument('-sequence_length_threshold', dest='sequence_length_threshold', type=int, help='Max. sequence length (skip longer sequences). Default=10000', default = 10000)
parser.add_argument('-keep_tmp', dest='keep_tmp', action='store_true', help='Keep temporary files')
parser.add_argument('-max_motif_length', dest='max_motif_length', type=int, help='max motif length. Default=12', default = 12)
parser.add_argument('-required_core_length', dest='required_core_length', type=int, help='required_core_length. Default=5', default = 5)
parser.add_argument('-pssm_crop_threshold', dest='pssm_crop_threshold', type=float, help='pssm_crop_threshold. Default=.25', default = .25)
parser.add_argument('-cluster_pssm_threshold_fraction', dest='cluster_pssm_threshold_fraction', type=float, help='cluster_pssm_threshold_fraction. Default=.5', default = .5)
parser.add_argument('-cluster_sim_threshold', dest='cluster_sim_threshold', type=float, help='similarity threshold used ONLY when counting number of successful '
 'predictions (e.g. for benchmark evaluation on ChIP-chip data) '
 '. Default=.75', default = .75)
parser.add_argument('-fraction_of_top_sequences_to_consider', dest='fraction_of_top_sequences_to_consider', type=float, help='fraction_of_top_sequences_to_consider. Default=.5', default = 0.5)
parser.add_argument('-hypegeom_p_val_cutoff', dest='hypegeom_p_val_cutoff', type=float, help='hypegeom_p_val_cutoff. Default=1e-30', default = 1e-30)
parser.add_argument('-num_motif_clusters_to_output', dest='num_motif_clusters_to_output', type=int, help='num_motif_clusters_to_output. Default=10', default = 10)
parser.add_argument('-degen_threshold', dest='degen_threshold', type=float, help='degen_threshold. Default=75', default = .75)
parser.add_argument('-fraction_degen_pos_threshold', dest='fraction_degen_pos_threshold', type=float, help='fraction_degen_pos_threshold. Default=.75', default = .75)
parser.add_argument('-max_gene_set_size_percentage_threshold', dest='max_gene_set_size_percentage_threshold', type=float, help='max_gene_set_size_percentage_threshold. Default=.1', default = .1)
parser.add_argument('-min_gene_set_size_absolute_threshold', dest='min_gene_set_size_absolute_threshold', type=int, help='min_gene_set_size_absolute_threshold. Default=20', default = 20)
parser.add_argument('-min_gene_set_size_percentage_threshold', dest='min_gene_set_size_percentage_threshold', type=float, help='min_gene_set_size_percentage_threshold. Default=.01', default = 0.01)
parser.add_argument('-num_clusters_to_output_motif_occurrences_for', dest='num_clusters_to_output_motif_occurrences_for', type=int, help='num_clusters_to_output_motif_occurrences_for. Default=3', default = 3)
parser.add_argument('-species', dest='species', type=str, help='species. Default="human"', default = 'human')
parser.add_argument('-random_runs_filename', dest='random_runs_filename', type=str, help='random_runs_filename. Default="random_scores"', default = 'random_scores')
parser.add_argument('-num_random_runs', dest='num_random_runs', type=int, help='num_random_runs. Default=0', default = 0)
parser.add_argument('-significance_estimates_file', dest='significance_estimates_file', type=str, help='significance_estimates_file. Default="significance_estimates"', default = 'significance_estimates')
parser.add_argument('-strand_specific_analysis', dest='strand_specific_analysis', type=str, help='Whether to use the reverse complement. For DNA -> no, for RNA -> yes. Default="no"', default = 'no')
parser.add_argument('-bootstrap_Std_threshold', dest='bootstrap_Std_threshold', type=float, help='bootstrap_Std_threshold. Default=2.0', default = 2.0)
parser.add_argument('-bootstrap_motif_score', dest='bootstrap_motif_score', type=str, help='bootstrap_motif_score. Default="no"', default = 'no')
parser.add_argument('-num_bootstrap_draws', dest='num_bootstrap_draws', type=int, help='num_bootstrap_draws. Default=20', default = 20)
parser.add_argument('-use_regression_scoring', dest='use_regression_scoring', type=str, help='use_regression_scoring. Default="no"', default = 'no')
parser.add_argument('-fast_mode', dest='fast_mode', type=str, help='fast_mode. Default="no"', default = 'no')


args = parser.parse_args()

if __name__ == '__main__':
    shutil.rmtree(args.output, ignore_errors=True)

    tmppath = os.path.join(args.output, 'tmp')
    fastapath = os.path.join(tmppath, 'fasta')
    
    os.makedirs(os.path.join(tmppath, 'data', 'evidence'), exist_ok=True)
    os.makedirs(os.path.join(tmppath, 'data', 'seq'), exist_ok=True)
    os.makedirs(os.path.join(tmppath, 'covariates', 'input'), exist_ok=True)
    os.makedirs(os.path.join(tmppath, 'evidence_file_list'), exist_ok=True)
    os.makedirs(os.path.join(tmppath, 'seq_file_list'), exist_ok=True)
    os.makedirs(os.path.join(tmppath, 'cov_file_list'), exist_ok=True)
    os.makedirs(fastapath, exist_ok=True)
    
    inputfasta = os.path.join(fastapath, 'input.fasta')
    
    if args.fasta is not None:    
        shutil.copy(args.fasta, inputfasta)
    else:
        print('extract fasta from genome ...')
        # extract the input sequences
        df = pd.read_csv(args.bed, sep='\t', header=None)
        df = df.sort_values(by=args.scorefield, ascending=False)
        if args.ntop is not None:
            df = df.iloc[:min(args.ntop,df.shape[0]),]
        
        df.to_csv(inputfasta+'.bed', sep='\t', index=False, header=False)
        genome = args.genome
        cmd = 'bedtools getfasta -fi {genome} -fo {output} -bed {bed}'.format(
           genome=genome, output=inputfasta + '.tmp',
           bed=inputfasta+'.bed')
        
        subprocess.run([cmd], shell=True)
        
        bed = iter(BedTool(inputfasta+'.bed'))
        
        scorefield = args.scorefield
        with open(inputfasta + '.tmp', 'r') as fi:
            with open(inputfasta, 'w') as fo:
               for line in fi:
                   if line[0] == '>':
                       _ = fo.write(line[:-1] + '|{}\n'.format(next(bed).fields[scorefield]))
                   else:
                       _ = fo.write(line)
        
        os.remove(inputfasta + '.tmp')
    
    ALPHABET = ['A', 'C', 'G', 'T']
    with open(os.path.join(tmppath, 'oligos'), 'w') as f:
        f.write('\n'.join([''.join(x) for x in product(*([ALPHABET]*args.seed_length))]))
    
    print('convert fasta sequence to cERMIT format')
    # remove too long sequences, do zero-padding to get equal length and remove sequences with zeros.
    # for strand_specific ='no' add both sequence strands.
    #subprocess.run(['/fast/AG_Akalin/wkopp/alison_sciatac_2018/extra/cermitAll/fasta_format_reader.rb io_params_file={}'.format(paramsfile)], shell=True)
    #
    records = SeqIO.parse(inputfasta, 'fasta')
    
    # remove too long sequences and sequences with Ns
    max_sequence_length = args.sequence_length_threshold
    records_ = [r.upper() for r in records if ('N' not in r and 'n' not in r) and len(r) <= max_sequence_length]
    
    maxlen = max([len(x) for x in records_])
    print('processing {} sequences of length {} bp each.'.format(len(records_), maxlen))
    
    # padding to get equal length sequences
    for i, _ in enumerate(records_):
        records_[i] += 'X' * (maxlen - len(records_[i]))
    
    evfile = os.path.join(tmppath, 'data', 'evidence', 'evidence_input')
    seqfile = os.path.join(tmppath, 'data', 'seq', 'seq_input')
    with open(evfile, 'w') as fevi, open(seqfile, 'w') as fseq:
        for rec in records_:
            name, score = rec.id.split('|')
            fevi.write('{}\t{}\n'.format(name, score))
            fseq.write('{}$'.format(str(rec.seq)))
            if args.strand_specific_analysis == 'no':
                fseq.write('{}$'.format(str(rec.seq.reverse_complement())))
            
    
    print('extract covariates')
    #cwd = os.getcwd()
    #os.chdir('/fast/AG_Akalin/wkopp/alison_sciatac_2018/extra/cermitAll/')
    #subprocess.run(['/fast/AG_Akalin/wkopp/alison_sciatac_2018/extra/cermitAll/extract_seq_covariates.rb io_params_file={}'.format(paramsfile)], shell=True)
    #os.chdir(cwd)
    #
    with open(os.path.join(tmppath, 'covariates', 'input', 'length_MM1.var_names'), 'w') as f:
        f.write(' '.join(['length', 'AC_or_GT', 'CA_or_TG', 'TC_or_GA', 'CT_or_AG', 'TT_or_AA', 'CC_or_GG', 'CG', 'GC', 'AT', 'TA']))
        f.write('\n')
    
    with open(os.path.join(tmppath, 'covariates', 'input', 'length_MM1'), 'w') as f:
        f.write('{}\t{}\n'.format(len(records_), 11))
        for rec in records_:   
            cnts = Counter()
            for i in range(len(rec)-1):
                cnts[str(rec.seq[i:(i+2)])] +=1
       
            n = (sum([ cnts[''.join(x)] for x in product(*([ALPHABET] *2))]) + 1) * 2
            #n = len(rec) * 2
            #print(cnts)
            #print([''.join(x) for x in product(*([ALPHABET] *2))])
            #print([ cnts[''.join(x)] for x in product(*([ALPHABET] *2))])
            f.write('{} {} {} {} {} {} {} {} {} {} {}\n'.format(n, cnts['AC'] + cnts['GT'], cnts['CA'] + cnts['TG'],
                                                                cnts['TC'] + cnts['GA'], cnts['CT'] + cnts['AG'],
                                                                cnts['TT'] + cnts['AA'], cnts['CC'] + cnts['GG'],
                                                                cnts['CG'] *2, cnts['GC'] *2, cnts['AT'] *2, cnts['TA'] *2))
    
    covfile = os.path.join(tmppath, 'covariates', 'input', 'MM1')
    
    #print('generate auxiliary files')
    #subprocess.run(['/fast/AG_Akalin/wkopp/alison_sciatac_2018/extra/cermitAll/generate_filelist_files.rb io_params_file={}'.format(paramsfile)], shell=True)
    ##
    with open(os.path.join(tmppath, 'evidence_file_list', 'evidence_file_list_input'), 'w') as f:
        f.write(evfile)
    with open(os.path.join(tmppath, 'seq_file_list', 'sequence_file_list_input'), 'w') as f:
        f.write(seqfile)
    with open(os.path.join(tmppath, 'cov_file_list', 'covariates_file_list_input'), 'w') as f:
        f.write(covfile)
    # /fast/AG_Akalin/wkopp/cermitAll/cERMIT \
    # input/evidence_file_list/evidence_file_list_state9.bed \
    # input/seq_file_list/sequence_file_list_state9.bed \
    # out_seed_length=6:min_motif_length=6:use_regression_scoring=no:bootstrap_motif_score=no:min_gene_set_size_percentage_threshold=0.01_tmp \
    # ./oligos_size_6 \
    # 6 12 5 0.25 0.5 0.75 0.5 1.0e-30 10 0.75 0.75 0.1 20 0.01 3 \
    # human random_scores 0 significance_estimates no 2.0 no 20 no input/cov_file_list/covariates_file_list_state9.bed no
    
    cmd = "{prog} {evidence_file_list} {sequence_file_list} {out_dir} "
    cmd += "{in_seeds_file} {min_motif_len} {max_motif_len} "
    cmd += "{required_core_length} {pssm_crop_threshold} {cluster_pssm_threshold_fraction} "
    cmd += "{cluster_sim_threshold} {fraction_of_top_sequences_to_consider} {hypegeom_p_val_cutoff} "
    cmd += "{num_motif_clusters_to_output} {degen_threshold} {fraction_degen_pos_threshold} "
    cmd += "{max_gene_set_size_percentage_threshold} {min_gene_set_size_absolute_threshold} {min_gene_set_size_percentage_threshold} "
    cmd += "{num_clusters_to_output_motif_occurrences_for} {species} {random_runs_filename} "
    cmd += "{num_random_runs} {significance_estimates_file} {strand_specific_analysis} "
    cmd += "{bootstrap_Std_threshold} {bootstrap_motif_score} {num_bootstrap_draws} "
    cmd += "{use_regression_scoring} {covariates_file_list} {fast_mode} "
    
    tmpresultspath = os.path.join(tmppath, 'results')
    
    options = {'prog': args.prog, 
               'evidence_file_list': os.path.join(tmppath, 'evidence_file_list', 'evidence_file_list_input'),
               'sequence_file_list': os.path.join(tmppath, 'seq_file_list', 'sequence_file_list_input'),
               'covariates_file_list': os.path.join(tmppath, 'cov_file_list', 'covariates_file_list_input'),
               'out_dir': tmpresultspath,
               'in_seeds_file': os.path.join(tmppath, 'oligos'),
               'min_motif_len': args.min_motif_length,
               'max_motif_len': args.max_motif_length,
               'required_core_length': args.required_core_length,
               'pssm_crop_threshold': args.pssm_crop_threshold,
               'cluster_pssm_threshold_fraction': args.cluster_pssm_threshold_fraction,
               'cluster_sim_threshold': args.cluster_sim_threshold,
               'fraction_of_top_sequences_to_consider': args.fraction_of_top_sequences_to_consider,
               'hypegeom_p_val_cutoff': args.hypegeom_p_val_cutoff,
               'num_motif_clusters_to_output': args.num_motif_clusters_to_output,
               'degen_threshold': args.degen_threshold,
               'fraction_degen_pos_threshold': args.fraction_degen_pos_threshold,
               'max_gene_set_size_percentage_threshold': args.max_gene_set_size_percentage_threshold,
               'min_gene_set_size_absolute_threshold': args.min_gene_set_size_absolute_threshold,
               'min_gene_set_size_percentage_threshold': args.min_gene_set_size_percentage_threshold,
               'num_clusters_to_output_motif_occurrences_for': args.num_clusters_to_output_motif_occurrences_for,
               'species': args.species,
               'random_runs_filename': args.random_runs_filename,
               'num_random_runs': args.num_random_runs,
               'significance_estimates_file': args.significance_estimates_file,
               'strand_specific_analysis': args.strand_specific_analysis,
               'bootstrap_Std_threshold': args.bootstrap_Std_threshold,
               'bootstrap_motif_score': args.bootstrap_motif_score,
               'num_bootstrap_draws': args.num_bootstrap_draws,
               'use_regression_scoring': args.use_regression_scoring,
               'fast_mode': args.fast_mode }
    
    
    cmd = cmd.format(**options)
    print(cmd)
    print('running cERMIT ...')
    subprocess.run([cmd], shell=True)
    
    oldprefix = 'input'
    newprefix = args.prefix

    motifs = glob.glob(os.path.join(tmpresultspath, '*_TRANSFAC_like.txt'))
    for motif in motifs:
        shutil.move(motif, motif.replace(oldprefix, newprefix))
        oldmotifpdf = motif.replace('.txt', '.pdf')
        shutil.move(oldmotifpdf, oldmotifpdf.replace('_TRANSFAC_like', '').replace(oldprefix, newprefix))
        

    motifs = glob.glob(os.path.join(tmpresultspath, '*_TRANSFAC_like.txt'))
    print('discovered {} motifs'.format(len(motifs)))

    if len(motifs) > 0:
        print('converting motifs to meme format')
        
        with open(os.path.join(tmpresultspath, 'cermitmotifs.meme'), 'w') as fo:
            fo.write("MEME version 5\n\nALPHABET= ACGT\n\nstrand: + -\n\n")
            for motif in motifs:
               fo.write('MOTIF {}\n'.format(os.path.basename(motif).split('_TRANSFAC')[0]))
               with open(motif) as fi:
                   content = fi.readlines()
                   newlines = []
                   for line in content[2:]:
                       # cermit output is weired, which makes the parsing difficult
                       newlines.append('\t'.join(str(float(e[:-1])) for e in line.split('\t')[1:-1]))
               fo.write('letter-probability matrix: alength= 4 w= {}\n'.format(len(newlines)))
               fo.write('\n'.join(newlines))
               fo.write('\n\n')
           
    
    if os.path.exists(os.path.join(tmpresultspath, 'cermitmotifs.meme')) and not args.no_tomtom:
        print('comparing motifs to database')
        cmd = 'tomtom -verbosity 1 {} {} -oc {}'.format(os.path.join(tmpresultspath, 'cermitmotifs.meme'),
                                                        args.motifdb,
                                                        os.path.join(tmpresultspath, 'tomtom_out'))
        print(cmd)
        subprocess.run([cmd], shell=True)
    
        if os.path.exists(os.path.join(tmpresultspath, 'tomtom_out')):
            shutil.move(os.path.join(tmpresultspath, 'tomtom_out'), os.path.join(args.output))
    
    for motif in glob.glob(os.path.join(tmpresultspath, '*.pdf')):
        shutil.move(motif, args.output)

    if os.path.exists(os.path.join(tmpresultspath, 'input_summary.txt')):
        shutil.move(os.path.join(tmpresultspath, 'input_summary.txt'), os.path.join(args.output, 'summary.txt'))

    if os.path.exists(os.path.join(tmpresultspath, 'cermitmotifs.meme')):
        shutil.move(os.path.join(tmpresultspath, 'cermitmotifs.meme'), os.path.join(args.output))
    
    if not args.keep_tmp:
        shutil.rmtree(tmppath, ignore_errors=True)

    print('finished {}'.format(args.output))
