import os
import shlex
import hashlib

from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

from shared.tn.vcf import VcfReader, VcfWriter, Position
from shared.tn.utils import str2bool, str_none, reference_sequence_from, subprocess_popen

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]

def calculate_file_md5(file_path):
    md5_hash = hashlib.md5()

    with open(file_path, 'rb') as file:
        for chunk in iter(lambda: file.read(4096), b''):
            md5_hash.update(chunk)

    return md5_hash.hexdigest()


def insert_after_line(original_str, target_line, insert_str):
    lines = original_str.split('\n')
    for i, line in enumerate(lines):
        if line == target_line:
            lines.insert(i+1, insert_str.rstrip('\n'))
            break
    return '\n'.join(lines)


class VcfReader_Database(object):
    def __init__(self, vcf_fn,
                 ctg_name=None,
                 direct_open=False,
                 keep_row_str=False,
                 save_header=False):
        self.vcf_fn = vcf_fn
        self.ctg_name = ctg_name
        self.variant_dict = defaultdict(Position)
        self.direct_open = direct_open
        self.keep_row_str = keep_row_str
        self.header = ""
        self.save_header = save_header

    def read_vcf(self):
        if self.vcf_fn is None or not os.path.exists(self.vcf_fn):
            return

        if self.direct_open:
            vcf_fp = open(self.vcf_fn)
            vcf_fo = vcf_fp
        else:
            vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (self.vcf_fn)))
            vcf_fo = vcf_fp.stdout
        for row in vcf_fo:
            columns = row.strip().split()
            if columns[0][0] == "#":
                if self.save_header:
                    self.header += row
                continue
            chromosome, position = columns[0], columns[1]
            if self.ctg_name is not None and chromosome != self.ctg_name:
                continue
            reference, alternate = columns[3], columns[4]
            position = int(position)
            row_str = row if self.keep_row_str else False
            key = (chromosome, position) if self.ctg_name is None else position

            self.variant_dict[key] = Position(ctg_name=chromosome,
                                              pos=position,
                                              ref_base=reference,
                                              alt_base=alternate,
                                              row_str=row_str
                                              )


def mosaicbase_tag(args):
    mosaicbase_tag_vcf_header_info = ''
    ctg_name = args.ctg_name
    pileup_vcf_fn = args.pileup_vcf_fn
    mosaicbase_fn = args.mosaicbase_fn
    pileup_output_vcf_fn = args.output_vcf_fn

    input_vcf_reader = VcfReader(vcf_fn=pileup_vcf_fn,
                                 ctg_name=ctg_name,
                                 show_ref=args.show_ref,
                                 keep_row_str=True,
                                 filter_tag=args.input_filter_tag,
                                 save_header=True)
    input_vcf_reader.read_vcf()
    pileup_variant_dict = input_vcf_reader.variant_dict

    input_variant_dict_set_contig = defaultdict(set)
    input_variant_dict_id_set_contig = defaultdict(set)

    for k, v in pileup_variant_dict.items():
        if ctg_name is None:
            if args.show_ref:
                input_variant_dict_set_contig[k[0]].add(str(k[1]))
                input_variant_dict_id_set_contig[k[0]].add(str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])
            elif args.show_germline:
                if v.filter == "PASS" or v.filter == "Germline":
                    input_variant_dict_set_contig[k[0]].add(str(k[1]))
                    input_variant_dict_id_set_contig[k[0]].add(str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])
            elif not args.disable_print_nonmosaic_calls:
                if v.filter == "PASS" or v.filter == "NonMosaic":
                    input_variant_dict_set_contig[k[0]].add(str(k[1]))
                    input_variant_dict_id_set_contig[k[0]].add(str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])
            else:
                if v.filter == "PASS":
                    input_variant_dict_set_contig[k[0]].add(str(k[1]))
                    input_variant_dict_id_set_contig[k[0]].add(str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])
        else:
            if args.show_ref:
                input_variant_dict_set_contig[ctg_name].add(str(k))
                input_variant_dict_id_set_contig[ctg_name].add(str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])
            elif args.show_germline:
                if v.filter == "PASS" or v.filter == "Germline":
                    input_variant_dict_set_contig[ctg_name].add(str(k))
                    input_variant_dict_id_set_contig[ctg_name].add(str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])
            elif not args.disable_print_nonmosaic_calls:
                if v.filter == "PASS" or v.filter == "NonMosaic":
                    input_variant_dict_set_contig[ctg_name].add(str(k))
                    input_variant_dict_id_set_contig[ctg_name].add(str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])
            else:
                if v.filter == "PASS":
                    input_variant_dict_set_contig[ctg_name].add(str(k))
                    input_variant_dict_id_set_contig[ctg_name].add(str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])

    total_input = 0
    for k, v in input_variant_dict_id_set_contig.items():
        total_input += len(input_variant_dict_id_set_contig[k])

    print("[INFO] Processing in {}...".format(ctg_name))

    input_keys_list = list(input_variant_dict_id_set_contig.keys())

    input_inter_mosaicbase_variant_dict_id_set_contig = defaultdict(set)

    mosaicbase_vcf_reader = VcfReader_Database(vcf_fn=str(mosaicbase_fn),
                                          ctg_name=ctg_name,
                                          keep_row_str=True,
                                          save_header=True)
    mosaicbase_vcf_reader.read_vcf()
    mosaicbase_input_variant_dict = mosaicbase_vcf_reader.variant_dict

    mosaicbase_variant_dict_id_set_contig = defaultdict(set)

    for k, v in mosaicbase_input_variant_dict.items():
        if ctg_name is None:
            if k[0] not in input_keys_list:
                continue
            mosaicbase_variant_dict_id_set_contig[k[0]].add(str(k[1]) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])
        else:
            mosaicbase_variant_dict_id_set_contig[ctg_name].add(str(k) + '\t' + v.reference_bases + '\t' + v.alternate_bases[0])

    for k, v in input_variant_dict_id_set_contig.items():
        input_inter_mosaicbase_variant_dict_id_set_contig[k] = input_variant_dict_id_set_contig[k].intersection(mosaicbase_variant_dict_id_set_contig[k])

    total_tagging_by_mosaicbase = 0

    for k, v in input_inter_mosaicbase_variant_dict_id_set_contig.items():
        total_tagging_by_mosaicbase += len(input_inter_mosaicbase_variant_dict_id_set_contig[k])

    print("[INFO] Processing in {}: tagged by {}: {}".format(ctg_name, str(mosaicbase_fn), total_tagging_by_mosaicbase))

    mosaicbase_vcf_md5 = calculate_file_md5(str(mosaicbase_fn))
    mosaicbase_vcf_header_info_base = '##INFO=<ID=MosaicBase,Number=0,Type=Flag,Description="file={},md5={},variant tagged by MosaicBase">'.format(str(mosaicbase_fn), mosaicbase_vcf_md5)
    mosaicbase_vcf_header_info_gene = '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene information from MosaicBase">'
    mosaicbase_vcf_header_info_disease = '##INFO=<ID=DISEASE,Number=1,Type=String,Description="Disease information from MosaicBase">'
    mosaicbase_vcf_header_info = mosaicbase_vcf_header_info_base + '\n' + mosaicbase_vcf_header_info_gene + '\n' + mosaicbase_vcf_header_info_disease + '\n'
    mosaicbase_tag_vcf_header_info += mosaicbase_vcf_header_info

    contigs_order = major_contigs_order + input_keys_list
    contigs_order_list = sorted(input_keys_list, key=lambda x: contigs_order.index(x))

    original_input_vcf_header = input_vcf_reader.header
    last_filter_line = '##FILTER=<ID=NonMosaic,Description="Non-mosaic variant tagged by panels of normals">'
    if mosaicbase_tag_vcf_header_info != '':
        new_mosaicbase_tag_vcf_header = insert_after_line(original_input_vcf_header, last_filter_line, mosaicbase_tag_vcf_header_info)
    else:
        new_mosaicbase_tag_vcf_header = original_input_vcf_header

    with open(pileup_output_vcf_fn, 'w') as output:
        output.write(new_mosaicbase_tag_vcf_header)
        for contig in contigs_order_list:
            all_pos_info = sorted(input_variant_dict_id_set_contig[contig], key=lambda x: int(x.split('\t')[0]))
            for pos_info in all_pos_info:
                pos = int(pos_info.split('\t')[0])
                key = (contig, pos) if ctg_name is None else pos
                ori_row_str = pileup_variant_dict[key].row_str
                columns = ori_row_str.split('\t')
                if pos_info in input_inter_mosaicbase_variant_dict_id_set_contig[contig]:
                    mosaicbase_row_str = mosaicbase_input_variant_dict[key].row_str
                    mosaicbase_columns = mosaicbase_row_str.split('\t')
                    mosaicbase_info_str = mosaicbase_columns[7].strip('\n')
                    info_str_ori = columns[7]
                    info_str_db = ""
                    info_str_db += mosaicbase_info_str
                    columns[7] = info_str_ori + ";" + "MosaicBase" + ";" + info_str_db
                row_str = '\t'.join(columns)
                output.write(row_str)

    output.close()


def main():
    parser = ArgumentParser(description="MosaicBase tagging for pileup data")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--pileup_vcf_fn', type=str, default=None,
                        help="Pileup VCF input")

    parser.add_argument('--mosaicbase_fn', type=str, default=None,
                        help="The path of the MosaicBase used for tagging mosaic variants. Default: if not specified, provided 'MosaicBase.GRCh38.simplified.vcf.gz' will be included")

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

    parser.add_argument('--show_germline', action='store_true',
                        help="Show germline calls in VCF file")

    parser.add_argument("--disable_print_nonmosaic_calls", action='store_true',
                        help="Disable print non-mosaic calls. Default: enable non-mosaic calls printing")

    global args
    args = parser.parse_args()

    mosaicbase_tag(args)


if __name__ == "__main__":
    main()
