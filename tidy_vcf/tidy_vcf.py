import gzip
import argparse
import random

def main():
    parser = argparse.ArgumentParser(
        
    prog = "tidy_vcf",    

    description="Given a vcf file, produces a tidy versions of sites and genotypes data, in efforts to make it easier to calculate summary statitics and visualizations.")

    parser.add_argument('-s', '--sites', nargs="?", type=str, required = False,
                help='Input file listing tab delimited sites to output. Each line must be a chromosome/scaffold name and the 1-indexed position of the site.')

    parser.add_argument('-t', '--thin', type=int, required = False,
        help = 'Alternative to --sites, where a sites will be selected no less than THIN bases apart.')

    parser.add_argument('-v', '--vcf', nargs="?", type=str, required = True,
                help='VCF input file. Should end in .vcf or .gz (plain text or compressed).')

    parser.add_argument('-o', '--sites_out', nargs="?", type=str, required = True,
                help='Name of output file for tidy sites data. Cannot be used simultaneously with --thin')

    parser.add_argument('-g', '--genotype_out', nargs="?", type=str, required = True,
                help='Name of output file for tidy genotype data.')
    args = parser.parse_args()

    if args.sites and not args.thin:
        pass
    elif args.thin and not args.sites:
        pass
    elif args.thin and args.sites:
        raise ValueError(f'--sites (-s) OR --thin (-t) should be used, not both.')
    else:
        print("neither --sites (-s) or --thin (-t) were given. Using all VCF sites. This might take a while.")

    sites_out = open(args.sites_out, "w") 
    genotype_file = open(args.genotype_out, "w") 

    def vcf_dict(vcf_file):
        info_dict = {}
        format_dict = {}
        with openfile(vcf_file) as fl:
            for line in fl:
                if line[0:2] != "##":
                    break
                else:
                    if line[0:8] == "##FORMAT":
                        stat = line.split(",")[0].split('=')[-1]
                        format_dict[stat] = ""
                    if line[0:6] == "##INFO":
                        stat = line.split(",")[0].split('=')[-1]
                        info_dict[stat] = ""
        return {'info':info_dict , 'format':format_dict}



    def build_info(info_str, info_dictionary, print_header = False):
        info_dict = info_dictionary.copy()
        if info_str != ".":
            info_pairs = [i.split('=') for i in info_str.split(';')] 
            stats = [i[0] for i in info_pairs]          
            values = [i[1] for i in info_pairs]
        else:
            stats = []
            values = [] 
        for key, value in info_dict.items():
            if key in stats:
                idx = stats.index(key)  
                info_dict[key] += f"{values[idx]}"
            else: 
                info_dict[key] += "NA"
        #return info_dict
        if print_header:
            return '\t'.join(info_dict.keys())
        else:
            return '\t'.join(info_dict.values())

    #def build_GTs(GTs, format_dict, print_header = False):
    #    gts = [g.split(":") for g in GTs]
    #    if print_header:
    #        print('\t'.join(format_dict.keys()))
        

    def get_inds(vcf_line):
        return vcf_line.split("\t")[9:]
      
    def df_content(CHROM, POS, ID, REF, ALT, QUAL, FILTER):
        return f"{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}"

    def parse_GTs(vcf_line):
        GTs = get_inds(vcf_line)
        return [g.split(":") for g in GTs]

    def parse_sites(sites_file):
        sites_dict = {}
        with open(sites_file) as fl:
            for line in fl:
                chrom, site = line.strip().split()
                site_str = f'{chrom}\t{site}'
                if site_str not in sites_dict:
                    sites_dict[site_str] = ""
        return sites_dict

    def sites_pass(CHROM, POS, sites_dict = None, previous_chrom = '', previous_pos = ''):
        state = False
        if args.sites:
            if f'{CHROM}\t{POS}' in sites_dict:
                state =True
        if args.thin:
            if CHROM == previous_chrom and int(POS) - int(previous_pos) > args.thin:
                state = True
            if CHROM != previous_chrom:
                state = True
                #print(CHROM == previous_chrom)
                #print(f'prev is {int(previous_pos)}, current is {int(POS)}')
                #print(f'prev chrom is {previous_chrom}, current is {CHROM}')
        return state


    def openfile(filename):
        if filename.endswith(".gz"):
            return gzip.open(filename, "rt")
        else:
            return open(filename, "r")

    def process_site(line, inds, vcf_keys):
        vcf_list = line.strip().split("\t")
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = vcf_list[0:9]
        GTs = vcf_list[9:]
        gts = [g.split(":") for g in GTs]
        site_row = df_content(CHROM, POS, ID, REF, ALT, QUAL, FILTER)
        #print("LENGTHS")
        #print(len(inds), len(gts))
        row_iter = zip(inds, gts)
        info_data = build_info(INFO, vcf_keys['info'])
        
        print(f"{site_row}\t{info_data}", file = sites_out)
        for ind, gt in row_iter:
            gt_final = '\t'.join(gt)
            print(f"{ind}\t{site_row}\t{gt_final}", file = genotype_file)
      
    def parse_file(vcf, sites_file = None):
        previous_chrom = ''
        previous_pos = 0
        vcf_keys = vcf_dict(args.vcf)
        first = True
        if args.sites:
            sites_dict = parse_sites(sites_file)
        with openfile(vcf) as fl:
            for line in fl:
                if line[0:2] == "##":
                    ln = line.strip()
                    print(ln, file = sites_out)
                    print(ln, file = genotype_file)
                elif line[0:6] == "#CHROM":
                    inds = get_inds(line.strip())
                elif line[0] != "#":
                    vcf_list = line.strip().split("\t")
                    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = vcf_list[0:9]
                    if first:
                        GTs = vcf_list[9:]
                        #info_headers = '\t'.join([sc.split("=")[0] for sc in INFO.split(";")])
                        info_headers = build_info(INFO, vcf_keys['info'], True)
                        format_l = '\t'.join(FORMAT.split(":"))
                        print(f"CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t{info_headers}", file = sites_out)
                        print(f"IND\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t{format_l}", file = genotype_file)
                        process_site(line, inds, vcf_keys)
                        first = False
                    else:
                        passing = False
                        if args.sites:
                            passing = sites_pass(CHROM, POS, sites_dict)
                        elif args.thin: 
                            passing = sites_pass(CHROM, POS, None, previous_chrom, previous_pos)
                            previous_chrom = CHROM
                            previous_pos = POS
                        else:
                            passing = True
                        if passing:
                            process_site(line, inds, vcf_keys)

    parse_file(args.vcf, args.sites)     

if __name__ == "__main__":
#    # execute only if run as a script
    main()

#print(vcf_keys)

#print(vcf_keys['info'])
#ss = "AC=7;AF=0.130;AN=54;BaseQRankSum=-4.342;ClippingRankSum=0.000;DP=102;ExcessHet=0.0834"
#build_info(ss, vcf_keys['info'], True)
#print(build_info(ss, vcf_keys['format']))
#ss = "AC=7;AF=0.130;AN=54;BaseQRankSum=-4.342;ClippingRankSum=0.000;DP=102;ExcessHet=0.0834"
#print(build_info(ss, vcf_keys['info'], True))
#print(build_info(ss, vcf_keys['format']))
