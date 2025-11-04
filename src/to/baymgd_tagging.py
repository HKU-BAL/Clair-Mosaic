from argparse import ArgumentParser, SUPPRESS
import numpy as np
from scipy.stats import beta, truncnorm
import shlex
import os
import gzip
from shared.to.utils import str_none, subprocess_popen


class BayesianMosaicClassifier:
    def __init__(self):
        self.base_prior = {
            'mosaic': 0.05,
            'germline': 0.95
        }

        self.likelihood_params = {
            'AF_obs': {
                'mosaic': {'dist': 'beta', 'alpha': 3, 'beta': 30},
                'germline': {'dist': 'truncnorm', 'mu': 0.5, 'sigma': 0.18},
            },
            'gnomAD_AF': {
                'mosaic': {'dist': 'beta', 'alpha': 1.2, 'beta': 200},
                'germline': {'dist': 'beta', 'alpha': 2, 'beta': 2},
            },
            'QUAL': {
                'mosaic': {'dist': 'lognorm', 'mu': 3.0, 'sigma': 1.5},
                'germline': {'dist': 'lognorm', 'mu': 3.0, 'sigma': 1.5},
            }
        }

    def _trunc_normal(self, mu, sigma):
        a = (0 - mu) / sigma
        b = (1 - mu) / sigma
        return lambda x: truncnorm.pdf(x, a, b, loc=mu, scale=sigma)

    def _safe_likelihood(self, value, hypothesis, evidence_type):
        params = self.likelihood_params[evidence_type][hypothesis]
        dist = params['dist']

        try:
            if dist == 'beta':
                value = np.clip(value, 1e-10, 1 - 1e-10)
                return beta.pdf(value, params['alpha'], params['beta'])
            elif params['dist'] == 'truncnorm':
                trunc_fn = self._trunc_normal(params['mu'], params['sigma'])
                return trunc_fn(np.clip(value, 0, 1))
            elif dist == 'lognorm':
                value = max(value, 1e-10)
                sigma = max(params['sigma'], 1e-10)
                return (1 / (value * sigma * np.sqrt(2 * np.pi)) *
                        np.exp(-(np.log(value) - params['mu']) ** 2 / (2 * sigma ** 2)))
            return 1.0
        except:
            return 1.0

    def _compute_base_posterior(self, AF_obs, gnomAD_AF, QUAL):
        prior = self.base_prior

        AF_obs = np.clip(AF_obs, 0, 1)
        gnomAD_AF = np.clip(gnomAD_AF, 0, 1)
        QUAL = max(QUAL, 0)

        unnorm = {}
        for hyp in prior:
            lk_AF = max(self._safe_likelihood(AF_obs, hyp, 'AF_obs'), 1e-300)
            lk_g = max(self._safe_likelihood(gnomAD_AF, hyp, 'gnomAD_AF'), 1e-300)
            lk_Q = max(self._safe_likelihood(QUAL, hyp, 'QUAL'), 1e-300)
            hyp_prior = max(prior[hyp], 1e-300)

            log_posterior = np.log(lk_AF) + np.log(lk_g) + np.log(lk_Q) + np.log(hyp_prior)
            unnorm[hyp] = log_posterior

        max_log = max(unnorm.values())
        log_total = max_log + np.log(sum(np.exp(v - max_log) for v in unnorm.values()))

        return {k: float(f"{np.exp(v - log_total):.4g}") for k, v in unnorm.items()}

    def compute_posterior(self, AF_obs, gnomAD_AF, QUAL):
        raw_pp = self._compute_base_posterior(AF_obs, gnomAD_AF, QUAL)

        def calc_phred(p_correct, p_wrong):
            return max(0, -10 * np.log10(max(1e-10, p_wrong) / max(1e-10, p_correct)) + 2.0)

        return {
            'pp_mosaic': float(f"{raw_pp['mosaic']:.4g}"),
            'pp_germline': float(f"{raw_pp['germline']:.4g}"),
            'phred_mosaic': int(round(calc_phred(raw_pp['mosaic'], raw_pp['germline']))),
            'phred_germline': int(round(calc_phred(raw_pp['germline'], raw_pp['mosaic']))),
        }


def parse_input_vcf_line(line, ctg_name):
    if line.startswith('#'): return None

    fields = line.strip().split('\t')
    chr = fields[0]
    if chr != ctg_name: return None

    info = fields[7]
    fmt = fields[8].split(':')
    vals = fields[9].split(':')
    fmt_dict = dict(zip(fmt, vals))

    return {
        'chrom': fields[0],
        'pos': int(fields[1]),
        'ref': fields[3],
        'alt': fields[4],
        'qual': float(fields[5]),
        'filter': fields[6],
        'af': float(fmt_dict.get('AF', 0)),
        'gq': float(fmt_dict.get('GQ', 0)),
        'raw_line': line,
        'info': info,
    }


def load_gnomad_chrom(gnomad_dir, chrom, is_indel=False):
    gnomad_file = f"gnomad.r2.1.af-ge-0.00001-full.chr{chrom}.sites.vcf.gz"
    gnomad_path = os.path.join(gnomad_dir, gnomad_file)

    if not os.path.exists(gnomad_path):
        print(f"Warning: gnomAD file not found: {gnomad_path}")
        return {}

    gnomad_dict = {}
    try:
        with gzip.open(gnomad_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                pos = int(fields[1])
                ref = fields[3]
                alts = fields[4].split(',')

                info = dict(item.split('=') for item in fields[7].split(';') if '=' in item)
                af_list = info.get('AF', '0').split(',')

                for i, alt in enumerate(alts):
                    af = float(af_list[i]) if i < len(af_list) else 1e-6
                    key = (pos, ref, alt)
                    if not is_indel:
                        if (len(ref) == 1 and len(alt) == 1 ):
                            gnomad_dict[key] = af
                        else:
                            continue
                    else:
                        if not (len(ref) == 1 and len(alt) == 1):
                            gnomad_dict[key] = af
                        else:
                            continue

        print(f"Loaded gnomAD data for chr{chrom} with {len(gnomad_dict)} variants")
    except Exception as e:
        print(f"Error loading gnomAD file: {str(e)}")

    return gnomad_dict


def baymgd_tag(args):
    input_path = args.pileup_vcf_fn
    gnomad_dir = args.gnomad_dir
    output_path = args.output_vcf_fn
    ctg_name = args.ctg_name

    print("[INFO] Processing in {}...".format(ctg_name))

    # Read input VCF
    input_vars = []
    headers = []
    vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (input_path)))
    vcf_fo = vcf_fp.stdout
    with vcf_fo as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line)
                continue
            if var := parse_input_vcf_line(line, ctg_name):
                input_vars.append(var)

    gnomad_cache = {}

    # Process variants
    classifier = BayesianMosaicClassifier()

    with open(output_path, 'w') as f_out:
        header_lines = []
        last_filter_index = -1
        for i, line in enumerate(headers):
            header_lines.append(line)
            if line.startswith('##FILTER='):
                last_filter_index = i

        # Write all lines up to and including the last FILTER line
        for line in header_lines[:last_filter_index + 1]:
            f_out.write(line)

        # Insert our new INFO headers after the last FILTER line
        f_out.write('##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="gnomAD population allele frequency">\n')
        f_out.write('##INFO=<ID=phred_pp_mosaic,Number=1,Type=Integer,Description="Phred-scaled posterior probability score for mosaic variants">\n')
        f_out.write('##INFO=<ID=phred_pp_germline,Number=1,Type=Integer,Description="Phred-scaled posterior probability score for germline variants">\n')

        # Write the remaining lines (after the last FILTER line)
        for line in header_lines[last_filter_index + 1:]:
            f_out.write(line)

        # Process each variant
        for var in input_vars:
            chrom = var['chrom']
            if chrom not in gnomad_cache:
                gnomad_cache[chrom] = load_gnomad_chrom(gnomad_dir, chrom.replace('chr', ''), is_indel=args.is_indel)

            filter = var['filter']

            if filter == 'PASS':
                key = (var['pos'], var['ref'], var['alt'])
                gnomad_af = gnomad_cache[chrom].get(key, 1e-6)

                pp = classifier.compute_posterior(var['af'], gnomad_af, var['gq'])

                # Build new INFO
                new_info = var['info'].split(';')
                if key in gnomad_cache[chrom]:
                    new_info.append(f'gnomAD_AF={gnomad_af:.4g}')
                new_info.extend([
                    f'phred_pp_mosaic={pp["phred_mosaic"]}',
                    f'phred_pp_germline={pp["phred_germline"]}',
                ])

                # Write output line
                fields = var['raw_line'].strip().split('\t')
                fields[7] = ';'.join(new_info)
                f_out.write('\t'.join(fields) + '\n')
            elif args.disable_print_nonmosaic_calls and filter == 'NonMosaic':
                continue
            elif not args.show_ref and filter == 'RefCall':
                continue
            else:
                # Write output line
                fields = var['raw_line'].strip().split('\t')
                f_out.write('\t'.join(fields) + '\n')

    f_out.close()


def main():
    parser = ArgumentParser(description="BayMGD tagging for pileup data")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--pileup_vcf_fn', type=str, default=None,
                        help="Pileup VCF input")

    parser.add_argument('--gnomad_dir', type=str, default=None,
                        help="The path of the gnomAD population database used for BayMGD tagging.")

    parser.add_argument('--output_vcf_fn', type=str, default=None,
                        help="Output VCF file")

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

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Absolute path to the 'samtools', samtools version >= 1.10 is required. Default: %(default)s")

    parser.add_argument('--show_ref', action='store_true',
                        help="Show reference calls (0/0) in VCF file")

    parser.add_argument("--disable_print_nonmosaic_calls", action='store_true',
                        help="Disable print non-mosaic calls. Default: enable non-mosaic calls printing")

    ## filtering for Indel candidates
    parser.add_argument('--is_indel', action='store_true',
                        help=SUPPRESS)

    global args
    args = parser.parse_args()

    baymgd_tag(args)


if __name__ == "__main__":
    main()
