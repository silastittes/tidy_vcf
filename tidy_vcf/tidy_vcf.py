## TODO: get passing bussiness and specfic sites working


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
                    info_keys_string = "\t".join(h.info.keys())
                    format_keys_string = "\t".join(h.format.keys())
                    print(f"{h.content_header}\t{info_keys_string}", file=sites_out)
                    print(
                        f"{h.content_header}\tIND\t{format_keys_string}",
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


if __name__ == "__main__":
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


### testing ###
line = "1       89211   .       G       .       19.64   QUAL_filter     BaseQRankSum=-1.383;DP=65;ExcessHet=3.01;MLEAC=.;MLEAF=.;MQ=23.2;MQRankSum=-0.674;ReadPosRankSum=0.21   GT:DP:RGQ       0/0:2:6   0/0:5:15        0/0:1:3 0/0:2:6 0/0:6:18        0/0:6:8 0/0:5:12        0/0:12:36       0/0:5:15        0/0:3:9 0/0:8:21        0/0:9:21        0/0:3:9"
vcf_content = (
    """##FORMAT=<ID=AD,Number=R,Type=Integer,Description="">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="">
##FORMAT=<ID=GT,Number=1,Type=String,Description="">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="">
##FORMAT=<ID=PID,Number=1,Type=String,Description="">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="">
##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="">
##INFO=<ID=AC,Number=A,Type=Integer,Description="">
##INFO=<ID=AF,Number=A,Type=Float,Description="">
##INFO=<ID=AN,Number=1,Type=Integer,Description="">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="">
##INFO=<ID=DP,Number=1,Type=Integer,Description="">
##INFO=<ID=END,Number=1,Type=Integer,Description="">
##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="">
##INFO=<ID=FS,Number=1,Type=Float,Description="">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="">
##INFO=<ID=MQ,Number=1,Type=Float,Description="">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="">
##INFO=<ID=QD,Number=1,Type=Float,Description="">
##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="">
##INFO=<ID=SOR,Number=1,Type=Float,Description="">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Oblat_F_10      Oblat_F_11      Oblat_F_2       Oblat_F_5       Oblat_F_6       Oblat_F_7       Oblat_F_8       Oblat_F_9 Oblat_M_01      Oblat_M_12      Oblat_M_13      Oblat_M_3       Oblat_M_4
"""
    + line
)

VCF = "test.vcf"


class TestClass:

    def make_test_vcf(self):
        with open(VCF, "w") as fl:
            print(vcf_content, file=fl)

    def delete_test_vcf(self):
        os.remove(VCF)

    def test_header(self):
        self.make_test_vcf()
        h = Header(VCF)
        IDs = [
            "Oblat_F_10",
            "Oblat_F_11",
            "Oblat_F_2",
            "Oblat_F_5",
            "Oblat_F_6",
            "Oblat_F_7",
            "Oblat_F_8",
            "Oblat_F_9",
            "Oblat_M_01",
            "Oblat_M_12",
            "Oblat_M_13",
            "Oblat_M_3",
            "Oblat_M_4",
        ]

        all_keys = [
            "AC",
            "AF",
            "AN",
            "BaseQRankSum",
            "DP",
            "END",
            "ExcessHet",
            "FS",
            "InbreedingCoeff",
            "MLEAC",
            "MLEAF",
            "MQ",
            "MQRankSum",
            "QD",
            "RAW_MQandDP",
            "ReadPosRankSum",
            "SOR",
        ]

        assert h.individual_ids == IDs
        info_keys = list(h.info.keys())
        assert sum([int(i == j) for i, j in zip(info_keys, all_keys)]) == len(all_keys)
        self.delete_test_vcf()

    def test_vcf_line(self):
        self.make_test_vcf()
        h = Header(VCF)
        v = VCFLine(line)
        v.build_formats(h)
        v.process_site(h)
        assert (
            v.info_string
            == "1\t89211\t.\tG\t.\t19.64\tQUAL_filter\tNA\tNA\tNA\t-1.383\t65\tNA\t3.01\tNA\tNA\t.\t.\t23.2\t-0.674\tNA\tNA\t0.21\tNA"
        )
        assert (
            v.format_strings[0]
            == "1\t89211\t.\tG\t.\t19.64\tQUAL_filter\tOblat_F_10\tNA\t2\tNA\t0/0\tNA\tNA\tNA\tNA\tNA\t6\tNA"
        )
        self.delete_test_vcf()
