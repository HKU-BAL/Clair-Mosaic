import os
import shlex
import gc
import subprocess
import concurrent.futures

from collections import Counter
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

import shared.to.param as param
from shared.to.vcf import VcfReader, VcfWriter
from shared.to.utils import str2bool, str_none, reference_sequence_from, subprocess_popen
from shared.to.interval_tree import bed_tree_from, is_region_in

LOW_AF_SNV = 0.2
LOW_AF_INDEL = 0.2
flanking = 100


def delete_lines_after(target_str, delimiter):
    lines = target_str.split('\n')
    index = 0
    for i, line in enumerate(lines):
        if delimiter in line:
            index = i
            break
    processed_lines = lines[:index+1]
    processed_str = '\n'.join(processed_lines) + '\n'
    return processed_str


def get_base_list(columns):
    pileup_bases = columns[4]

    base_idx = 0
    base_list = []
    read_end_set = set()
    read_start_set = set()
    while base_idx < len(pileup_bases):
        base = pileup_bases[base_idx]
        if base == '+' or base == '-':
            base_idx += 1
            advance = 0
            while True:
                num = pileup_bases[base_idx]
                if num.isdigit():
                    advance = advance * 10 + int(num)
                    base_idx += 1
                else:
                    break
            base_list[-1][1] = base + pileup_bases[base_idx: base_idx + advance]  # add indel seq
            base_idx += advance - 1

        elif base in "ACGTNacgtn#*":
            base_list.append([base, ""])
        elif base == '^':  # start of read, next base is mq, update mq info
            base_idx += 1
            read_start_set.add(len(base_list) - 1)
        # skip $, the end of read
        if base == "$":
            read_end_set.add(len(base_list) - 1)
        base_idx += 1
    read_start_end_set = read_start_set if len(read_start_set) > len(read_end_set) else read_end_set
    upper_base_counter = Counter([''.join(item).upper() for item in base_list])
    return upper_base_counter, base_list, read_start_end_set


def postfilter_per_pos(args):
    pos = args.pos
    ctg_name = args.ctg_name
    ref_base = args.ref_base
    alt_base = args.alt_base
    tumor_bam_fn = args.tumor_bam_fn
    ref_fn = args.ref_fn
    samtools = args.samtools
    min_bq = args.min_bq
    min_mq = args.min_mq
    is_snp = len(ref_base) == 1 and len(alt_base) == 1
    is_ins = len(ref_base) == 1 and len(alt_base) > 1
    is_del = len(ref_base) > 1 and len(alt_base) == 1
    pass_hetero_both_side = True
    hetero_germline_set = set()
    homo_germline_set = set()

    args.af = args.af if args.af is not None else 1.0

    if not os.path.exists(tumor_bam_fn):
        tumor_bam_fn += ctg_name + '.bam'

    if args.hetero_info is not None and args.hetero_info != "":
        hetero_germline_set = set([tuple(item.split('-')) for item in args.hetero_info.split(',')])
    if args.homo_info is not None and args.homo_info != "":
        homo_germline_set = set([tuple(item.split('-')) for item in args.homo_info.split(',')])

    flanking = args.flanking

    ctg_range = "{}:{}-{}".format(ctg_name, max(pos - flanking, 1), pos + flanking + 1)
    samtools_command = "{} mpileup --min-MQ {} --min-BQ {} --excl-flags 2316 -r {} --output-QNAME --output-extra HP ".format(
        samtools, min_mq, min_bq, ctg_range)

    tumor_samtools_command = samtools_command + tumor_bam_fn

    reference_sequence = reference_sequence_from(
        samtools_execute_command=samtools,
        fasta_file_path=ref_fn,
        regions=[ctg_range]
    )

    # tumor
    pos_dict = defaultdict(defaultdict)
    pos_counter_dict = defaultdict(defaultdict)
    hap_dict = defaultdict(int)
    ALL_HAP_LIST = [0, 0, 0]
    HAP_LIST = [0, 0, 0]
    alt_base_read_name_set = set()
    homo_germline_pos_set = set([int(item[0]) for item in homo_germline_set])
    hetero_germline_pos_set = set([int(item[0]) for item in hetero_germline_set])

    samtools_mpileup_tumor_process = subprocess_popen(shlex.split(tumor_samtools_command), stderr=subprocess.PIPE)
    for row in samtools_mpileup_tumor_process.stdout:
        columns = row.split('\t')

        read_name_list = columns[6].split(',')
        base_counter, base_list, read_start_end_set = get_base_list(columns)

        base_list = [[''.join(item[0]).upper()] + [item[1]] for item in base_list]

        p = int(columns[1])
        ctg = columns[0]
        ctg = p if args.ctg_name is not None else (ctg, p)
        if ctg in hetero_germline_pos_set or p == pos:
            phasing_info = columns[7].strip('\n').split(',')
            for hap_idx, hap in enumerate(phasing_info):
                if hap in '12' and read_name_list[hap_idx] not in hap_dict:
                    hap_dict[read_name_list[hap_idx]] = int(hap)

        pos_dict[p] = dict(zip(read_name_list, base_list))
        center_ref_base = reference_sequence[p - max(pos - flanking, 1)]

        if p == pos:
            for rn in read_name_list:
                ALL_HAP_LIST[hap_dict[rn]] += 1

            if is_snp:
                alt_base_read_name_set = set(
                    [key for key, value in zip(read_name_list, base_list) if ''.join(value) == alt_base])
            elif is_ins:
                alt_base_read_name_set = set(
                    [key for key, value in zip(read_name_list, base_list) if
                     ''.join(value).replace('+', '').upper() == alt_base and '+' in ''.join(value)])
            elif is_del:
                alt_base_read_name_set = set(
                    [key for key, value in zip(read_name_list, base_list) if
                     len(ref_base) == len(value[1]) and '-' in value[1]])

            for rn in alt_base_read_name_set:
                HAP_LIST[hap_dict[rn]] += 1

        if len(base_counter) == 1 and base_counter[center_ref_base] > 0:
            continue
        pos_counter_dict[p] = base_counter

    alt_hap_counter = Counter([hap_dict[key] for key in alt_base_read_name_set])

    hp0, hp1, hp2 = alt_hap_counter[0], alt_hap_counter[1], alt_hap_counter[2]
    MAX = max(hp1, hp2)
    MIN = min(hp1, hp2)
    af = float(args.af)
    if is_snp and af < LOW_AF_SNV:
        if hp1 * hp2 > 0 and (MIN > args.min_alt_coverage or MAX / MIN <= 10):
            pass_hetero_both_side = False
    elif not is_snp and af < LOW_AF_INDEL:
        if hp1 * hp2 > 0 and (MIN > args.min_alt_coverage or MAX / MIN <= 10):
            pass_hetero_both_side = False

    all_hp0, all_hp1, all_hp2 = ALL_HAP_LIST
    hp0, hp1, hp2 = HAP_LIST
    phaseable = all_hp1 * all_hp2 > 0 and hp1 * hp2 == 0 and (int(hp1) > args.min_alt_coverage or int(hp2) > args.min_alt_coverage)

    pass_hard_filter = pass_hetero_both_side

    print(' '.join([ctg_name, str(pos), str(pass_hard_filter), str(phaseable), str(pass_hetero_both_side)]))


def update_filter_info(args, key, row_str, fail_set, phasable_set, fail_pass_hetero_both_side_set):
    ctg_name = key[0] if args.ctg_name is None else args.ctg_name
    pos = key[1] if args.ctg_name is None else key
    k = (ctg_name, pos)
    columns = row_str.split('\t')

    is_candidate_filtered = 0
    phaseable = k in phasable_set

    if phaseable:
        if columns[7] == '.':
            columns[7] = 'H'
        else:
            columns[7] = 'H;' + columns[7]

    if k in fail_set:
        columns[5] = '0.0000'
        columns[6] = "LowQual"
        is_candidate_filtered = 1
    if k in fail_pass_hetero_both_side_set:
        columns[6] += ";"
        columns[6] += "MultiHap"

    row_str = '\t'.join(columns)

    return row_str, is_candidate_filtered


def postfilter(args):
    ctg_name = args.ctg_name
    threads = args.threads
    threads_low = max(1, int(threads * 4 / 5))
    enable_postfilter = args.enable_postfilter
    pileup_vcf_fn = args.pileup_vcf_fn
    germline_vcf_fn = args.germline_vcf_fn
    flanking = args.flanking
    output_dir = args.output_dir
    is_indel = args.is_indel
    cmrg_bed_fn = args.cmrg_bed_fn
    if not os.path.exists(output_dir):
        subprocess.run("mkdir -p {}".format(output_dir), shell=True)

    pileup_output_vcf_fn = args.output_vcf_fn
    if not enable_postfilter:
        subprocess.run("ln -sf {} {}".format(pileup_vcf_fn, pileup_output_vcf_fn), shell=True)
        return

    germine_input_vcf_reader = VcfReader(vcf_fn=germline_vcf_fn,
                                         ctg_name=ctg_name,
                                         show_ref=False,
                                         keep_row_str=False,
                                         filter_tag="PASS",
                                         save_header=False,
                                         skip_genotype=False)
    germine_input_vcf_reader.read_vcf()
    germline_input_variant_dict = germine_input_vcf_reader.variant_dict

    germline_gt_dict = defaultdict(list)
    # only keep homo germline
    for key in list(germline_input_variant_dict.keys()):
        ctg = args.ctg_name if args.ctg_name is not None else key[0]
        pos = key if args.ctg_name is not None else key[1]
        alt_base = germline_input_variant_dict[key].alternate_bases[0]
        if sum(germline_input_variant_dict[key].genotype) == 1:
            germline_gt_dict[ctg].append((pos, 1, alt_base))
        elif sum(germline_input_variant_dict[key].genotype) == 2:
            germline_gt_dict[ctg].append((pos, 2, alt_base))

    for k, v in germline_gt_dict.items():
        germline_gt_dict[k] = list(sorted(v, key=lambda x: x[0]))

    input_vcf_reader = VcfReader(vcf_fn=pileup_vcf_fn,
                                 ctg_name=ctg_name,
                                 show_ref=args.show_ref,
                                 keep_row_str=True,
                                 discard_indel=False if is_indel else True,
                                 filter_tag=args.input_filter_tag,
                                 save_header=True,
                                 keep_af=True)
    input_vcf_reader.read_vcf()
    pileup_variant_dict = input_vcf_reader.variant_dict

    bed_tree = bed_tree_from(bed_file_path=cmrg_bed_fn, contig_name=ctg_name)

    input_variant_dict = defaultdict()
    for k, v in pileup_variant_dict.items():
        if v.filter != "PASS":
            continue
        if args.test_pos and k != args.test_pos:
            continue
        pass_bed_region = len(bed_tree) == 0 or is_region_in(tree=bed_tree,
                                                             contig_name=ctg_name,
                                                             region_start=v.pos - 1,
                                                             region_end=v.pos)
        if not pass_bed_region:
            continue
        input_variant_dict[k] = v

    output_vcf_header = input_vcf_reader.header
    last_format_line = '##FORMAT=<ID=TU,Number=1,Type=Integer,Description="Count of T in the input BAM">'
    output_vcf_header = delete_lines_after(output_vcf_header, last_format_line)
    p_vcf_writer = VcfWriter(vcf_fn=pileup_output_vcf_fn,
                             ctg_name=ctg_name,
                             ref_fn=args.ref_fn,
                             header=output_vcf_header,
                             show_ref_calls=True)

    if not is_indel:
        pf_info_output_path = os.path.join(output_dir, "PF_INFO_SNV")
    else:
        pf_info_output_path = os.path.join(output_dir, "PF_INFO_INDEL")
    with open(pf_info_output_path, 'w') as f:
        for key, POS in input_variant_dict.items():
            ctg_name = args.ctg_name if args.ctg_name is not None else key[0]
            pos = key if args.ctg_name is not None else key[1]
            hetero_flanking_list = []
            homo_flanking_list = []
            germline_gt_list = germline_gt_dict[ctg_name] if ctg_name in germline_gt_dict else []
            for p, gt, alt_base in germline_gt_list:
                if p > pos + flanking:
                    break
                if p > pos - flanking and p != pos:
                    if gt == 1:
                        hetero_flanking_list.append('-'.join([str(p), str(alt_base)]))
                    else:
                        homo_flanking_list.append('-'.join([str(p), str(alt_base)]))
            POS.extra_infos = [set(hetero_flanking_list), set(homo_flanking_list)]
            info_list = [ctg_name, str(pos), POS.reference_bases, POS.alternate_bases[0], str(POS.af), str(POS.qual), ','.join(hetero_flanking_list), ','.join(homo_flanking_list)]
            f.write(' '.join(info_list) + '\n')

    file_directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    main_entry = os.path.join(file_directory, "clair_mosaic_to.py")

    parallel_command = "{} -C ' ' -j{} {} {} haplotype_filtering".format(args.parallel, threads_low, args.pypy3, main_entry)
    parallel_command += " --ctg_name {1}"
    parallel_command += " --pos {2}"
    parallel_command += " --ref_base {3}"
    parallel_command += " --alt_base {4}"
    parallel_command += " --af {5}"
    parallel_command += " --qual {6}"
    parallel_command += " --hetero_info {7}"
    parallel_command += " --homo_info {8}"
    parallel_command += " --samtools " + str(args.samtools)
    parallel_command += " --tumor_bam_fn " + str(args.tumor_bam_fn)
    parallel_command += " --ref_fn " + str(args.ref_fn)
    parallel_command += " --debug " if args.debug else ""
    parallel_command += " --flanking " + str(args.flanking) if args.flanking is not None else ""
    parallel_command += " :::: " + str(pf_info_output_path)

    postfilter_process = subprocess_popen(shlex.split(parallel_command))

    total_num = 0
    fail_set = set()
    phasable_set = set()
    fail_pass_hetero_both_side_set = set()

    for row in postfilter_process.stdout:
        columns = row.rstrip().split()
        if len(columns) < 4:
            continue
        total_num += 1
        ctg_name, pos, pass_hard_filter, phasable, pass_hetero_both_side = columns[:5]
        pos = int(pos)
        pass_hard_filter = str2bool(pass_hard_filter)
        phasable = str2bool(phasable)
        pass_hetero_both_side = str2bool(pass_hetero_both_side)

        if not pass_hard_filter:
            fail_set.add((ctg_name, pos))
        if phasable:
            phasable_set.add((ctg_name, pos))
        if not pass_hetero_both_side:
            fail_pass_hetero_both_side_set.add((ctg_name, pos))

        if total_num > 0 and total_num % 1000 == 0:
            print("[INFO] Processing in {}, total processed positions: {}".format(ctg_name, total_num))

    postfilter_process.stdout.close()
    postfilter_process.wait()

    for key in sorted(pileup_variant_dict.keys()):
        row_str = pileup_variant_dict[key].row_str.rstrip()
        row_str, is_candidate_filtered = update_filter_info(args, key, row_str, fail_set, phasable_set, fail_pass_hetero_both_side_set)
        p_vcf_writer.vcf_writer.write(row_str + '\n')

    p_vcf_writer.close()

    print("[INFO] Total input calls: {}, filtered by multi haplotypes: {}".format(len(input_variant_dict), len(fail_pass_hetero_both_side_set)))


def main():
    parser = ArgumentParser(description="Post-filtering for long-read data")

    parser.add_argument('--tumor_bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--pileup_vcf_fn', type=str, default=None,
                        help="Pileup VCF input")

    parser.add_argument('--output_vcf_fn', type=str, default=None,
                        help="Output vcf file")

    parser.add_argument('--germline_vcf_fn', type=str, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Output VCF directory")

    parser.add_argument('--cmrg_bed_fn', type=str, default=None,
                        help="CMRG regions. Default: 'CMRGv1.0.GRCh38.bed'")

    parser.add_argument('--python', type=str, default="python3",
                        help="Absolute path to the 'python3', default: %(default)s")

    parser.add_argument('--threads', type=int, default=4,
                        help="Max #threads to be used")

    parser.add_argument('--input_filter_tag', type=str_none, default=None,
                        help='Filter variants with tag from the input VCF')

    parser.add_argument('--pypy3', type=str, default="pypy3",
                        help="Absolute path of pypy3, pypy3 >= 3.6 is required")

    parser.add_argument('--parallel', type=str, default="parallel",
                        help="Absolute path of parallel, parallel >= 20191122 is required")

    parser.add_argument('--show_ref', action='store_true',
                        help="Show reference calls (0/0) in VCF file")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Absolute path to the 'samtools', samtools version >= 1.10 is required. Default: %(default)s")

    # options for advanced users
    parser.add_argument('--enable_postfilter', type=str2bool, default=True,
                        help="EXPERIMENTAL: Apply haplotype filtering to the variant calls")

    parser.add_argument('--min_mq', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered, default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered, default: %(default)d")

    parser.add_argument('--min_alt_coverage', type=int, default=2,
                        help="Minimum number of reads supporting an alternative allele required for a somatic variant to be called. Default: %(default)d")

    ## filtering for Indel candidates
    parser.add_argument('--is_indel', action='store_true',
                        help=SUPPRESS)

    ## test using one position
    parser.add_argument('--test_pos', type=int, default=None,
                        help=SUPPRESS)

    ## flakning window size to process
    parser.add_argument('--flanking', type=int, default=100,
                        help=SUPPRESS)

    parser.add_argument('--debug', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--hetero_info', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--homo_info', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--pos', type=int, default=None,
                        help=SUPPRESS)

    parser.add_argument('--ref_base', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--alt_base', type=str, default=None,
                        help=SUPPRESS)

    parser.add_argument('--af', type=float, default=None,
                        help=SUPPRESS)

    parser.add_argument('--qual', type=float, default=None,
                        help=SUPPRESS)

    global args
    args = parser.parse_args()

    if args.pos is None:
        postfilter(args)
    else:
        postfilter_per_pos(args)


if __name__ == "__main__":
    main()