# BSD 3-Clause License
#
# Copyright 2023 The University of Hong Kong, Department of Computer Science
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os
import random
import shlex
import subprocess

from argparse import ArgumentParser, SUPPRESS
from subprocess import run as subprocess_run

from src.tn.utils import str2bool, str_none

random.seed(0)
cov_suffix = ".cov.mosdepth.summary.txt"


def get_coverage_from_bam(args, bam_fn, is_sample=False, output_prefix=None):
    ctg_name = args.ctg_name
    dry_run = args.dry_run

    mosdepth = args.mosdepth
    contig_option = "" if ctg_name is None else "-c {}".format(ctg_name)

    prefix = 'sample' if is_sample else 'control'
    if output_prefix is None:
        output_prefix = os.path.join(args.output_dir, 'tmp', 'raw_{}'.format(prefix))

    mos_depth_command = "{} -t {} {} -n -x --quantize 0:15:150: {}.cov {}".format(mosdepth,
                                                                                  args.samtools_threads,
                                                                                  contig_option,
                                                                                  output_prefix,
                                                                                  bam_fn)

    print("[INFO] Calculating coverge for {} BAM using mosdepth...".format(prefix))

    print('[INFO] Will run the following commands:')
    print(mos_depth_command)
    if not dry_run:
        subprocess.run(mos_depth_command, shell=True)
    else:
        return 1
    coverage_log = os.path.join(output_prefix + cov_suffix)
    if ctg_name is None:
        last_row = open(coverage_log).readlines()[-1]
        coverage = float(float(last_row.split()[3]))
    else:
        all_rows = open(coverage_log).readlines()
        ctg_row = [row for row in all_rows if row.split()[0] == ctg_name]
        if len(ctg_row) == 0:
            print('[ERROR] no contig coverage found for contig {}'.format(ctg_name))
        coverage = float(ctg_row[0].split()[3])
    return coverage


def random_sample(population, k, seed=0):
    random.seed(seed)
    return random.sample(population, k)


def gen_contaminated_bam(args):
    sample_bam_fn = args.sample_bam_fn
    control_bam_fn = args.control_bam_fn
    output_dir = args.output_dir
    ctg_name = args.ctg_name
    samtools_execute_command = args.samtools
    samtools_threads = args.samtools_threads
    control_purity = args.control_purity
    sample_purity = args.sample_purity

    if not os.path.exists(os.path.join(output_dir, "tmp")):
        rc = subprocess.run('mkdir -p {}/tmp'.format(output_dir), shell=True)
    control_bam_coverage = args.control_bam_coverage if args.control_bam_coverage else get_coverage_from_bam(args,
                                                                                                          control_bam_fn,
                                                                                                          False)
    sample_bam_coverage = args.sample_bam_coverage if args.sample_bam_coverage else get_coverage_from_bam(args,
                                                                                                       sample_bam_fn,
                                                                                                       True)

    print("[INFO] Control/Sample BAM coverage: {}/{}".format(control_bam_coverage,
                                                           sample_bam_coverage))

    # add control to control
    if control_purity is not None:
        contam_coverage = control_bam_coverage * (1 - control_purity)
        rest_control_coverage = control_bam_coverage - contam_coverage

        need_downsample = True
        if contam_coverage >= float(sample_bam_coverage):
            print("[WARNING] Contam BAM coverage is higher than sample BAM, use all sample BAM!")
            control_subsample_pro = "1.00"
            need_downsample = False
            sample_subsample_bam = sample_bam_fn
        else:
            control_subsample_pro = "%.3f" % (rest_control_coverage / float(control_bam_coverage))
            sample_subsample_bam = os.path.join(args.output_dir, 'tmp', 'sample_subsample.bam')

        sample_subsample_pro = "%.3f" % (contam_coverage / float(sample_bam_coverage))

        print("[INFO] Control/Sample subsample proportion: {}/{}".format(control_subsample_pro,
                                                                       sample_subsample_pro))

        control_subsample_bam = os.path.join(args.output_dir, 'tmp', 'control_rest.bam')

        contig_option = "" if ctg_name is None else ctg_name

        # sample subsample cmd
        t_s_cmd = "{} view -@{} -bh -s {} -o {} {} {}".format(samtools_execute_command,
                                                              samtools_threads,
                                                              sample_subsample_pro,
                                                              sample_subsample_bam,
                                                              sample_bam_fn,
                                                              contig_option,
                                                              )
        # control subsample cmd
        n_s_cmd = "{} view -@{} -bh -s {} -o {} {} {}".format(samtools_execute_command,
                                                              samtools_threads,
                                                              control_subsample_pro,
                                                              control_subsample_bam,
                                                              control_bam_fn,
                                                              contig_option)

        n_index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, sample_subsample_bam)
        t_index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, control_subsample_bam)

        print('[INFO] Will run the following commands:')
        print(t_s_cmd)
        print(n_s_cmd)
        print(n_index_cmd)
        print(t_index_cmd)

        if not args.dry_run:
            subprocess.run(n_s_cmd, shell=True)
            subprocess.run(n_index_cmd, shell=True)
            if need_downsample:
                subprocess.run(t_s_cmd, shell=True)
                subprocess.run(t_index_cmd, shell=True)

        control_output_bam = os.path.join(output_dir, "control_purity_{}.bam".format(control_purity))

        print("[INFO] Merging sample BAM into control BAM as contamination...")
        merge_cmd = "{} merge -f -@{} {} {} {}".format(samtools_execute_command, samtools_threads, control_output_bam,
                                                       sample_subsample_bam, control_subsample_bam)
        index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, control_output_bam)

        print('[INFO] Will run the following commands:')
        print(merge_cmd)
        print(index_cmd)
        if not args.dry_run:
            subprocess_run(merge_cmd, shell=True)
            subprocess.run(index_cmd, shell=True)

        if args.cal_output_bam_coverage:
            get_coverage_from_bam(args, control_output_bam, False, os.path.join(args.output_dir, 'cov'))

        if args.remove_intermediate_dir:
            tmp_file_path = os.path.join(args.output_dir, 'tmp')
            if os.path.exist(tmp_file_path):
                rc = subprocess.run("rm -rf {}".format(tmp_file_path), shell=True)

        print("[INFO] Finishing merging, output file: {}".format(control_output_bam))

    # add control to sample
    if sample_purity is not None:
        contam_coverage = sample_bam_coverage * (1 - sample_purity)
        rest_sample_coverage = sample_bam_coverage - contam_coverage

        need_downsample = True
        if contam_coverage > float(control_bam_coverage):
            print("[WARNING] Contaim coverage is higher than control, use all control BAM!")
            control_subsample_pro = "1.00"
            need_downsample = False
            control_subsample_bam = control_bam_fn
        else:
            control_subsample_pro = "%.3f" % (contam_coverage / float(control_bam_coverage))
            control_subsample_bam = os.path.join(args.output_dir, 'tmp', 'control_subsample.bam')
        sample_subsample_pro = "%.3f" % (rest_sample_coverage / float(sample_bam_coverage))

        print("[INFO] Control/Sample subsample proportion: {}/{}".format(control_subsample_pro,
                                                                       sample_subsample_pro))

        sample_subsample_bam = os.path.join(args.output_dir, 'tmp', 'sample_rest.bam')

        contig_option = "" if ctg_name is None else ctg_name

        # sample subsample cmd
        t_s_cmd = "{} view -@{} -bh -s {} -o {} {} {}".format(samtools_execute_command,
                                                              samtools_threads,
                                                              sample_subsample_pro,
                                                              sample_subsample_bam,
                                                              sample_bam_fn,
                                                              contig_option,
                                                              )
        # control subsample cmd
        n_s_cmd = "{} view -@{} -bh -s {} -o {} {} {}".format(samtools_execute_command,
                                                              samtools_threads,
                                                              control_subsample_pro,
                                                              control_subsample_bam,
                                                              control_bam_fn,
                                                              contig_option)

        t_index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, sample_subsample_bam)
        n_index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, control_subsample_bam)

        print('[INFO] Will run the following commands:')
        print(t_s_cmd)
        print(n_s_cmd)
        print(t_index_cmd)
        print(n_index_cmd)

        if not args.dry_run:
            subprocess.run(t_s_cmd, shell=True)
            subprocess.run(t_index_cmd, shell=True)
            if need_downsample:
                subprocess.run(n_s_cmd, shell=True)
                subprocess.run(n_index_cmd, shell=True)

        sample_output_bam = os.path.join(output_dir, "sample_purity_{}.bam".format(sample_purity))

        print("[INFO] Merging control BAM into sample BAM as contamination...")
        merge_cmd = "{} merge -f -@{} {} {} {}".format(samtools_execute_command, samtools_threads, sample_output_bam,
                                                       sample_subsample_bam, control_subsample_bam)
        index_cmd = '{} index -@{} {}'.format(samtools_execute_command, samtools_threads, sample_output_bam)

        print('[INFO] Will run the following commands:')
        print(merge_cmd)
        print(index_cmd)
        if not args.dry_run:
            subprocess_run(merge_cmd, shell=True)
            subprocess.run(index_cmd, shell=True)

        if args.cal_output_bam_coverage:
            get_coverage_from_bam(args, sample_output_bam, True, os.path.join(args.output_dir, 'cov'))

        if args.remove_intermediate_dir:
            tmp_file_path = os.path.join(args.output_dir, 'tmp')
            if os.path.exist(tmp_file_path):
                rc = subprocess.run("rm -rf {}".format(tmp_file_path), shell=True)

        print("[INFO] Finishing merging, output file: {}".format(sample_output_bam))


def main():
    parser = ArgumentParser(description="Generate contaminated BAM for benchmarking")

    parser.add_argument('--control_bam_fn', type=str, default=None,
                        help="Sorted control BAM file input")

    parser.add_argument('--sample_bam_fn', type=str, default=None,
                        help="Sorted sample BAM file input")

    parser.add_argument('--control_bam_coverage', type=float, default=None,
                        help="Control BAM coverage calculated using mosdepth")

    parser.add_argument('--sample_bam_coverage', type=float, default=None,
                        help="Sample BAM coverage calculated using mosdepth")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Sorted chunked BAM file output path")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required")

    parser.add_argument('--samtools_threads', type=int, default=32,
                        help="Samtools threads to read input BAM")

    parser.add_argument('--control_purity', type=float, default=None,
                        help="Control purity")

    parser.add_argument('--sample_purity', type=float, default=None,
                        help="Sample purity")

    parser.add_argument('--ctg_name', type=str_none, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--mosdepth', type=str, default="mosdepth",
                        help="Path to the 'mosdepth'")

    # options for advanced users
    parser.add_argument('--dry_run', type=str2bool, default=0,
                        help="EXPERIMENTAL: Only print the synthetic log, debug only")

    parser.add_argument('--cal_output_bam_coverage', type=str2bool, default=0,
                        help="EXPERIMENTAL: calculate output contaminated BAM using mosdepth")

    parser.add_argument('--remove_intermediate_dir', action='store_true',
                        help="Remove intermediate directory. Default: False")

    args = parser.parse_args()

    gen_contaminated_bam(args)


if __name__ == "__main__":
    main()
