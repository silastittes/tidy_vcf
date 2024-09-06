from tidy_vcf import *

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
