#!/usr/bin/env python3
"""Author:  Cesare de Filippo
Contact: casare_filippo(at)eva.mpg.de

Report the estimates of spurious alignments as well as the the count of spurious and true alignments as described in:
de Filippo et al. Quantifying and reducing spurious alignments for the analysis of ultra-short ancient DNA sequences. BMC Biology 2018 16:121
https://doi.org/10.1186/s12915-018-0581-9
https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-018-0581-9
"""

##################################################
# Libraries required.
import argparse
import textwrap as _textwrap
import time
import math
from collections import OrderedDict
# Libraries to be installed before
import pysam
from scipy.stats import beta


class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return _textwrap.wrap(text, 90)
##################################################


parser = argparse.ArgumentParser(
    description="This is a script to report the count of spurious and true alignments as described in de Filippo, Meyer and Pruefer (2018)\n'Quantifying and reducing spurious alignments for the analysis of ultra-short ancient DNA sequences.'\n",
    formatter_class=LineWrapRawTextHelpFormatter)

parser.add_argument('-i', dest='input_file',
                    help='Input file in bam format.', required=True)
parser.add_argument('-s', dest='sites',
                    help='File with informative sites in bed format. It must be sorted, bgziped and tabixed.', required=True)
parser.add_argument('-l', dest='minLength',
                    help='Minimum sequence length to be considered. Default is 20 bp.', default=20, type=int)
parser.add_argument('-L', dest='maxLength',
                    help='Maximum sequence length to be reported in the output. Longer sequences are treated as of this length. Default is 40.', type=int, default=40)
parser.add_argument('-q', dest='BaseQual',
                    help='Minimum base quality (BQ) value. Default is 30.',   type=int, default=30)
parser.add_argument('-m', dest='MapQual',
                    help='Minimum mapping quality (MQ) value. Default is 1.', type=int, default=1)
parser.add_argument('-d', dest='deam',
                    help="Consider putatively deaminated sequences at the last i,j (5',3') terminal positions.", type=str, default="0,0", required=False)
parser.add_argument('-I', dest='Indels',
                    help="Remove insertions and deletions (indels).", default=False, required=False, action='store_true')
parser.add_argument('-D', dest='DoubleStrand',
                    help="The library is double stranded. Single stranded is the default.", default=False, required=False, action='store_true')
parser.add_argument('-T', dest='Transversions',
                    help="Use only transversions from the informative sites (-s option).", default=False, required=False, action='store_true')

args = parser.parse_args()


# The arguments to be used later in the script.
input_file = args.input_file
infosites = args.sites
BQ_cutoff = args.BaseQual
MQ_cutoff = args.MapQual
minLength = args.minLength
maxLength = args.maxLength
rm_Indels = args.Indels
DoubleStrand = args.DoubleStrand
Transversions = args.Transversions
terminal_deam = args.deam.split(',')
# Arguments for the deamination filter and to check the consistency of the -d option.
deam_filter_skip = True
if (terminal_deam[0] != '0' or terminal_deam[1] != '0'):
    deam_filter_skip = False
try:
    int(terminal_deam[0])
    int(terminal_deam[1])
except ValueError:
    print("-d {} incosistent. It must be two integers separeted by comma as follow:\n-d 1,1".format(args.deam))
    exit()

###########################
## Define some functions ##
###########################

#############################################
# Function to align the reference sequence (refseq) and the sequence of interest (myseq), and return the positions and basequelities.


def alnseq(refseq, myseq, CIGAR, basequalities, start):
    """Function to align the referece sequence (refseq) and the sequence of interest (myseq), and return the positions and basequelities."""
    # a1 = refseq; a2 = myseq.
    a1 = a2 = ''
    # i1 and i2 will be the ith element modified for indels for a1 and a2 respectively.
    i1 = i2 = 0
    # The positions and the basequalities (bq) to be 'adjusted' in case of indels.
    positions = []
    bq = []
    pos = int(start)-1
    for (cigar_op, cigar_len) in CIGAR:
        if cigar_op == 0:
            a1 += refseq[i1:(cigar_len+i1)]
            a2 += myseq[i2:(cigar_len+i2)]
            bq += basequalities[i2:(cigar_len+i2)]
            positions += list(range(pos+1, pos+cigar_len+1))
            i1 += cigar_len
            i2 += cigar_len
        if cigar_op == 1:  # It's an insertion relative to the reference
            a2 += myseq[i2:(cigar_len+i2)]
            bq += basequalities[i2:(cigar_len+i2)]
            positions += list(range(pos+1, pos+cigar_len+1))
            i2 += cigar_len
            a1 += '-' * cigar_len
        if cigar_op == 2:  # It's a deletion relative to the reference
            a2 += '-' * cigar_len
            # 90 is an arbitrary base quality for the reference sequence.
            bq += [90] * cigar_len
            positions += [pos+j for j in range(1, cigar_len+1)]
            a1 += refseq[i1:(cigar_len+i1)]
            i1 += cigar_len
        pos = positions[-1]
    # must return refseq, myseq, pos, bq
    return (a1, a2, positions, bq)
#############################################

#############################################

# Function to deternime whether a sequence is deaminated condition on the last 'i,j' terminal positions


def is_deaminated(refseq, myseq, terminal_positions=terminal_deam, isReverse=False, isDoubleStranded=False):
    """Function to deternime whether a sequence is deaminated condition on the last 'i,j' terminal positions"""
    refseq = refseq.upper()
    output = False
    nt = 'CT'
    if (isDoubleStranded == False):
        if (isReverse):
            nt = 'GA'
        for p in range(0, int(terminal_positions[0])):
            if (refseq[p]+myseq[p] in nt):
                output = True
        if (output == False):
            for p in range(1, int(terminal_positions[1])+1):
                if(refseq[-p]+myseq[-p] in nt):
                    output = True
    else:
        for p in range(0, int(terminal_positions[0])):
            if (refseq[p]+myseq[p] in 'CT'):
                output = True
        if (output == False):
            for p in range(1, int(terminal_positions[1])+1):
                if(refseq[-p]+myseq[-p] in 'GT'):
                    output = True

    return(output)


#############################################


#############################################
# Function to determine the state of the read and whether it passed the Strand Filter.
def count(Allele, Reference, Modified, isReverse=False):
    ret_type = 'Oth'
    pass_SF = True
    if Modified == Allele:
        ret_type = 'Mod'
    elif Reference == Allele:
        ret_type = 'Ref'
    # if Reference+Modified == 'CG' or Reference+Modified == 'GC':  # The SF removes CG or GC sites
    #    output['pass_SF'] = False
    # with the 2 if below, CG/GC are covered whether isReverse is true or false
    if 'C' in Reference+Modified and isReverse == False:
        pass_SF = False
    elif 'G' in Reference+Modified and isReverse == True:
        pass_SF = False
    return (ret_type, pass_SF)
#############################################


def binom_interval(success, total, confint=0.95):
    quantile = (1 - confint) / 2.
    lower = beta.ppf(quantile, success, total - success + 1)
    upper = beta.ppf(1 - quantile, success + 1, total - success)
    if (math.isnan(lower)):
        lower = 0
    if (math.isnan(upper)):
        upper = 1
    return (lower, upper)

#############################################
# End of the functions!!!!
#############################################


def main():
    # The results in a dictionary to be printed at the end of the script.
    output_table = OrderedDict({i: {'Ref': 0, 'Mod': 0, 'Oth': 0,
                                    'Ref_SF': 0, 'Mod_SF': 0, 'Oth_SF': 0} for i in range(minLength, maxLength+1)})

    # Max divergence allowed in bwa using the ancient paramenters and used to correct the estimates of spurious alignments.
    MaxDivBWA = {'20': 2, '21': 2,
                 '22': 3, '23': 3,  '24': 3,  '25': 3,  '26': 3,  '27': 3,  '28': 3,  '29': 3,  '30': 3,  '31': 3,  '32': 3,  '33': 3,  '34': 3,  '35': 3,  '36': 3,  '37': 3,  '38': 3,  '39': 3,  '40': 3,  '41': 3,
                 '42': 4, '43': 4, '44': 4, '45': 4, '46': 4, '47': 4, '48': 4, '49': 4, '50': 4, '51': 4, '52': 4, '53': 4, '54': 4, '55': 4, '56': 4, '57': 4, '58': 4, '59': 4, '60': 4}

    start = time.time()
    r = list(range(minLength, maxLength+1))
    with pysam.AlignmentFile(input_file, "rb", check_sq=False) as samfile, pysam.TabixFile(infosites) as tabixfile:
        for chrom in [str(k) for k in range(1, 23)] + ['X']:
            for read in samfile.fetch(chrom, until_eof=True): # until_eof=True prevent pysam to complain if there is no index file.
                Cigar = read.cigarstring
                if (rm_Indels):
                    if 'I' in Cigar or 'D' in Cigar:
                        continue
                # Filter out softclip, hardclip and for MapQuality cutoff
                if 'S' not in Cigar and 'H' not in Cigar and read.mapping_quality >= MQ_cutoff:
                    pos = read.get_reference_positions(full_length=False)
                    site_position = 0
                    passTvFilter = True
                    try:
                        for s in tabixfile.fetch(chrom, pos[0], pos[-1], parser=pysam.asBed()):
                            site_position = int(s[1])
                            reference = s[3]
                            modified = s[4]
                            if(Transversions == True):
                                passTvFilter = not reference + \
                                    modified in ['CT', 'TC', 'GA', 'AG']
                    except ValueError:
                        break
                    if site_position != 0 and passTvFilter:
                        refseq = read.get_reference_sequence()
                        myseq = read.query_sequence
                        bq = read.query_qualities
                        L = min(len(myseq), maxLength)
                        if site_position in pos and L >= minLength:
                            p = site_position
                            # Sequences have different length, we need to align them
                            # update the variables accordingly
                            if len(myseq) != len(refseq):
                                (refseq, myseq, pos, bq) = alnseq(
                                    refseq, myseq, CIGAR=read.cigartuples, basequalities=bq, start=pos[0])
                            # increment counters
                            if deam_filter_skip or is_deaminated(refseq, myseq, terminal_deam, read.is_reverse, isDoubleStranded=DoubleStrand):
                                Allele = myseq[pos.index(p)]
                                BQ = bq[pos.index(p)]
                                if Allele in ['A', 'C', 'G', 'T']:
                                    if BQ >= BQ_cutoff:
                                        ret_type, pass_SF = count(
                                            Allele, Reference=reference, Modified=modified, isReverse=read.is_reverse)
                                        output_table[L][ret_type] += 1
                                        if pass_SF:
                                            output_table[L][ret_type + '_SF'] += 1

    print('bp\tRef\tMod\tOth\tRef_SF\tMod_S\tOth_SF\tSpuriousAln_(95%CI)\tSpuriousAln_SF(95%CI)')
    SpAl = []
    TrAl = []
    for i, elem in sorted(output_table.items()):
        print(i, end='')
        d = MaxDivBWA[str(i)]/i
        for j in ['Ref', 'Mod', 'Oth', 'Ref_SF', 'Mod_SF', 'Oth_SF']:
            print('\t'+str(elem[j]), end='')
        for j in ['', '_SF']:
            TrueAln = float(elem['Ref'+j])
            SpuriousAln = float(elem['Mod'+j]+elem['Oth'+j])
            TrueAln1 = max(TrueAln - SpuriousAln*d/(3-d), 0)
            SpuriousAln1 = SpuriousAln / (1-d/3)
            SpAl.append(SpuriousAln1)
            TrAl.append(TrueAln1)
            spu = round(SpuriousAln1 / (SpuriousAln1 + TrueAln1), 4)
            ci = binom_interval(SpuriousAln1, SpuriousAln1+TrueAln1)
            print(
                '\t'+str(spu)+' ('+str(round(ci[0], 4))+','+str(round(ci[1], 4))+')', end='')
        print()

    # Split the SpAl and TrAl and then print the cutoffs using the cumulative estimates.
    spal = SpAl[0::2]
    spal_sf = SpAl[1::2]
    tral = TrAl[0::2]
    tral_sf = TrAl[1::2]
    j001 = j01 = j1 = j001sf = j01sf = j1sf = True
    for i in range(0, len(spal)):
        cum_spal = sum(spal[i:])/(sum(spal[i:])+sum(tral[i:]))
        cum_spal_sf = sum(spal_sf[i:])/(sum(spal_sf[i:])+sum(tral_sf[i:]))
        if(cum_spal < 0.001 and j001):
            print('# 0.1% cutoff is', r[i], 'bp')
            j001 = False
        if(cum_spal < 0.01 and j01):
            print('# 1% cutoff is', r[i], 'bp')
            j01 = False
        if(cum_spal < 0.1 and j1):
            print('# 10% cutoff is', r[i], 'bp')
            j1 = False
        if(cum_spal_sf < 0.001 and j001sf):
            print('# 0.1% cutoff with SF is', r[i], 'bp')
            j001sf = False
        if(cum_spal_sf < 0.01 and j01sf):
            print('# 1% cutoff with SF is', r[i], 'bp')
            j01sf = False
        if(cum_spal_sf < 0.1 and j1sf):
            print('# 10% cutoff with SF is', r[i], 'bp')
            j1sf = False

    end = time.time()
    print("#...done in", round((end - start)/60, 3), "minute(s)!")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        exit()
    exit()
