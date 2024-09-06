import sys
import gzip
import argparse
import os
from dataclasses import dataclass


def parse_args(args):
    parser = argparse.ArgumentParser(
        prog="tidy_vcf",
        description="Given a vcf file, produces a tidy versions of sites and genotypes data, in efforts to make it easier to calculate summary statitics and visualizations.",
    )

    parser.add_argument(
        "-s",
        "--sites",
        nargs="?",
        type=str,
        required=False,
        default=None,
        help="Input file listing tab delimited sites to output. Each line must be a chromosome/scaffold name and the 1-indexed position of the site.",
    )

    parser.add_argument(
        "-t",
        "--thin",
        type=int,
        required=False,
        default=0,
        help="Alternative to --sites, where a sites will be selected no less than THIN bases apart.",
    )

    parser.add_argument(
        "-v",
        "--vcf",
        nargs="?",
        type=str,
        required=True,
        help="VCF input file. Should end in .vcf or .gz (plain text or compressed).",
    )

    parser.add_argument(
        "-o",
        "--sites_out",
        nargs="?",
        type=str,
        required=True,
        help="Name of output file for tidy sites data. Cannot be used simultaneously with --thin",
    )

    parser.add_argument(
        "-g",
        "--genotype_out",
        nargs="?",
        type=str,
        required=True,
        help="Name of output file for tidy genotype data.",
    )

    return parser.parse_args()


def openfile(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")


@dataclass
class Header:
    """
    Class to store information about the header of a VCF file.

    Attributes:
    header: vcf_file
    """

    vcf_file: str

    def __post_init__(self):
        self.build_vcf_dicts()
        self.content_header = f"CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"

    def build_vcf_dicts(self):
        info_dict = {}
        format_dict = {}
        filter_dict = {}
        with openfile(self.vcf_file) as fl:
            for line in fl:
                if not line.startswith("#"):
                    break
                if line.startswith("##FORMAT"):
                    stat = line.split(",")[0].split("=")[-1]
                    format_dict.update({stat: "NA"})
                if line.startswith("##INFO"):
                    stat = line.split(",")[0].split("=")[-1]
                    info_dict.update({stat: "NA"})
                if line.startswith("##FILTER"):
                    stat = line.split(",")[0].split("=")[-1]
                    filter_dict.update({stat: "NA"})
                if line.startswith("#CHROM"):
                    self.individual_ids = line.strip().split()[9:]
        self.info = info_dict
        self.format = format_dict
        self.filter = filter_dict


@dataclass
class VCFLine:
    """
    Class to store information about a an individual VCF line.

    Attributes:
    vcf_line: str
    """

    vcf_line: str

    def __post_init__(self):
        self._process_line()

    def _process_line(self):
        # save all the parts as attributes
        line_list = self.vcf_line.split()
        self.gt_formats = line_list[9:]

        (
            self.CHROM,
            self.POS,
            self.ID,
            self.REF,
            self.ALT,
            self.QUAL,
            self.FILTER,
            self.INFO,
            self.FORMAT,
        ) = line_list[0:9]

        self.content_string = f"{self.CHROM}\t{self.POS}\t{self.ID}\t{self.REF}\t{self.ALT}\t{self.QUAL}\t{self.FILTER}"

    def process_site(self, header: Header):
        assert len(self.gt_formats) == len(header.individual_ids)
        self.build_info(header)
        self.build_formats(header)

    def build_info(self, header: Header):
        info_dict = header.info.copy()
        if self.INFO != ".":
            info_pairs = [i.split("=") for i in self.INFO.split(";")]
            stats = [i[0] for i in info_pairs]
            values = [i[1] for i in info_pairs]
        else:
            stats = []
            values = []
        for key, value in info_dict.items():
            if key in stats:
                idx = stats.index(key)
                info_dict.update({key: str(values[idx])})
        infos = "\t".join(info_dict.values())
        self.info_string = f"{self.content_string}\t{infos}"

    def build_formats(self, header: Header):
        format_list = []
        # list of the format names
        format_elements = self.FORMAT.split(":")
        for gt in self.gt_formats:
            ind_formats = gt.split(":")
            format_dict = header.format.copy()
            for element in format_elements:
                idx = format_elements.index(element)
                format_dict.update({element: ind_formats[idx]})
            format_list.append(format_dict)
        format_strings = [
            ind + "\t" + "\t".join(f.values())
            for ind, f in zip(header.individual_ids, format_list)
        ]
        self.format_strings = [f"{self.content_string}\t{f}" for f in format_strings]


@dataclass
class TidyVCF:
    """
    Collect information about a VCF file and parse.

    Attributes:
    vcf_file: str
        Path to the VCF file.
    """

    vcf_file: str
    sites_out: str
    genotype_out: str
    sites_file: str = None
    thin: int = 0

    def __post_init__(self):
        self.previous_chrom: str = ""
        self.previous_pos: int = 0
        if self.sites_file:
            self.parse_sites()

    def parse_sites(self):
        sites_dict = set()
        with open(self.sites_file) as fl:
            for line in fl:
                chrom, site = line.strip().split()
                site_str = f"{chrom}\t{site}"
                sites_dict.add(site_str)
        self.sites_dict = sites_dict

    def sites_pass(self, vcfline: VCFLine):
        state = False
        if self.sites_file is not None:
            if f"{vcfline.CHROM}\t{vcfline.POS}" in self.sites_dict:
                state = True
        if self.thin:
            if (
                vcfline.CHROM == self.previous_chrom
                and int(vcfline.POS) - int(self.previous_pos) > self.thin
            ):
                state = True
            if vcfline.CHROM != self.previous_chrom:
                state = True
        if self.sites_file is None and self.thin < 1:
            state = True
        return state

    def parse_file(self, vcf, header: Header):
        sites_out = open(self.sites_out, "w")
        genotype_out = open(self.genotype_out, "w")
        self.previous_chrom = ""
        self.previous_pos = 0
        first = True
        with openfile(vcf) as fl:
            for line in fl:
                ln = line.strip()
                if not ln:
                    continue
                if line.startswith("#"):
                    print(ln, file=sites_out)
                    print(ln, file=genotype_out)
                if not line.startswith("#") and first:
                    first = False
                    info_keys_string = "\t".join(header.info.keys())
                    format_keys_string = "\t".join(header.format.keys())
                    print(f"{header.content_header}\t{info_keys_string}", file=sites_out)
                    print(
                        f"{header.content_header}\tIND\t{format_keys_string}",
                        file=genotype_out,
                    )
                if not line.startswith("#") and not first:
                    vcfline = VCFLine(ln)
                    passing = self.sites_pass(vcfline)
                    self.previous_chrom = vcfline.CHROM
                    self.previous_pos = vcfline.POS
                    if passing:
                        vcfline.process_site(header)
                        print(vcfline.info_string, file=sites_out)
                        for gt in vcfline.format_strings:
                            print(gt, file=genotype_out)

def main():
    args = parse_args(sys.argv[1:])
    if args.sites and args.thin < 1:
        pass
    elif args.thin >= 1 and not args.sites:
        pass
    elif args.thin >= 1 and args.sites:
        raise ValueError(f"--sites (-s) OR --thin (-t) should be used, not both.")
    else:
        print(
            "neither --sites (-s) or --thin (-t) were given. Using all VCF sites. This might take a while."
        )

    h = Header(args.vcf)
    v = TidyVCF(
        vcf_file=args.vcf,
        sites_out=args.sites_out,
        genotype_out=args.genotype_out,
        sites_file=args.sites,
        thin=args.thin,
    )
    v.parse_file(args.vcf, h)
