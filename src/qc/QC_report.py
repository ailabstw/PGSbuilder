#!/usr/bin/python3
import os, sys, argparse, subprocess, json
from fpdf import FPDF
import pandas as pd

def parse_args() -> argparse.Namespace:
    """
    Returns:
        arguments
    """
    parser = argparse.ArgumentParser()

    # log files
    files = parser.add_argument_group('File Arguments')
    files.add_argument('--bfile', required=True, help='the bfile prefix of the input file')
    files.add_argument('--dedup_record', required=True, help='the SNP list after dedup and extracting valid chromosomes')
    files.add_argument('--snp_record', required=True, help='the json file records the information of SNP QC')
    files.add_argument('--ind_record', required=True, help='the json file records the information of sample QC')
    files.add_argument('--kinship_record', required=True, help='the sample list not passing kinship filter')
    files.add_argument('--population_record', required=True, help='the sample list passing population stratification')
    
    # figures
    figs = parser.add_argument_group('Figure Arguments')
    figs.add_argument('--fig_maf', required=True, help='the figure of MAF')
    figs.add_argument('--fig_vmiss', required=True, help='the figure of varint missing rate')
    figs.add_argument('--fig_hwe', required=True, help='the prefix of figures of H-W, e.g. [prefix].hwe.hist.png')
    figs.add_argument('--fig_het_vs_smiss', required=True, help='the figure of heterozygosity to sample missing rate')
    figs.add_argument('--fig_population_prefix', required=True, help='the prefix of figures of population stratification')

    # parameters
    params = parser.add_argument_group('Parameters')
    params.add_argument("--maf", type=float, required=True)
    params.add_argument("--vmiss", type=float, required=True)
    params.add_argument("--hwe_all", type=float, required=True)
    params.add_argument("--hwe_strict", type=str, required=True, help='"true" or "false"')
    params.add_argument("--hwe_control", type=float, required=True)
    params.add_argument("--hwe_case", type=float, required=True)
    params.add_argument("--flipscan", type=str, required=True, help='"true" or "false"')
    params.add_argument("--heter", type=float, required=True)
    params.add_argument("--smiss", type=float, required=True)
    params.add_argument("--sex", type=str, required=True, help='"true" or "false"')
    params.add_argument("--kinship", type=float, required=True)
    params.add_argument("--prune_window", type=float, required=True)
    params.add_argument("--prune_step", type=float, required=True)
    params.add_argument("--prune_threshold", type=float, required=True)
    params.add_argument("--population", type=str, required=True)
    params.add_argument("--population_sd", type=float, required=True)

    # others
    parser.add_argument("--logo_file", type=str, required=False, default="")
    parser.add_argument("--work_dir", type=str, required=True)
    
    args = parser.parse_args()

    return args


class PDF(FPDF):
    def page_header(self, Title):
        # Arial bold 15
        self.set_font('Times', 'B', 18)
        # Calculate width of title and position
        w = self.get_string_width(Title) + 6
        self.set_xy((210 - w) / 2, 15)
        # Colors of frame, background and text
        self.set_draw_color(256, 256, 256)
        self.set_fill_color(256, 256, 256)
        self.set_text_color(22, 143, 153)
        # Thickness of frame (1 mm)
        self.set_line_width(1)
        # Title
        self.cell(w, 9, Title, 1, 1, 'C', 1)
        # Line break
        self.ln(10)
        
    def logo(self, fig):
        if fig != '':
            self.image(fig, x = 5, y = 5, w = 50, h = 12)
        
    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Arial italic 8
        self.set_font('Times', 'I', 8)
        # Text color in gray
        self.set_text_color(51, 51, 51)
        # Page number
        self.cell(0, 10, 'Page ' + str(self.page_no()), 0, 0, 'C')

    def chapter_title(self, num, label):
        # Arial 12
        self.set_font('Times', 'B', 13)
        # Background color
        self.set_fill_color(22, 143, 153)
        self.set_text_color(256, 256, 256)
        # Title
        self.cell(0, 6, ' Section %s : %s' % (num, label), 0, 1, 'L', 1)
        # Line break
        self.ln(4)

    def chapter_body(self, txt):
        # Times 12
        self.set_font('Times', '', 12)
        # Output justified text
        self.multi_cell(0, 5, txt)

    def content(self, text_format_list):
        for text, form in text_format_list:
            if form == 'bold':
                self.set_font('Times', 'B', 12)
                self.write(5, text)
            else:
                self.set_font('Times', '', 12)
                self.write(5, text)


def Subprocess(command):
    return subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True).communicate()[0].split()[0]


def main():
    # Arguments
    args = parse_args()
    basename = args.bfile.split('/')[-1]
    inter_space = 3
    os.chdir(args.work_dir)
    
    # record
    records = dict()
    record_groups = ['input', 'filter', 'snp_qc', 'ind_qc', 'fig']
    for i in record_groups:
        records[i] = dict()
    
    # New page
    pdf = PDF()
    Title = "GWAS QC Report"
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(args.logo_file)
    pdf.set_line_width(0.5)
    

    ############## Introduction ################################
    # title
    pdf.chapter_title("1", "Introduction")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # get the number of samples and SNPs
    records['input']['name'] = basename
    records['input']['snp_num'] = int(Subprocess('wc -l {}.bim'.format(args.bfile)))
    records['input']['ind_num'] = int(Subprocess('wc -l {}.fam'.format(args.bfile)))
    records['input']['male_num'] = int(Subprocess("awk '$5 == 1' {}.fam | wc -l".format(args.bfile)))
    records['input']['female_num'] = int(Subprocess("awk '$5 == 2' {}.fam | wc -l".format(args.bfile)))

    # text and format
    text_list = [
        ('This report encompasses the quality control (QC) summary for the ', 'normal'),
        (records['input']['name'], 'bold'),
        (' data. A total of ', 'normal'),
        ('{} samples ({} males, {} females)'.format(records['input']['ind_num'], records['input']['male_num'], records['input']['female_num']), 'bold'),
        (' were genotyped for ', 'normal'),
        ('{} SNPs'.format(records['input']['snp_num']), 'bold'),
        ('. Quality control was performed across variants and samples. The results and the filtering criteria are discussed below.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)


    ################## cutoffs setting ########################### 
    # title
    pdf.chapter_title("2", "User-defined/Default cutoffs")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)
    
    # text and format
    text_list = [
        ('-- Minor Allele Frequency Cutoff:', 'normal'),
        ('{0} (lower bound)'.format(args.maf), 'bold'),
        ('-- Genotype Missingness Cutoff:', 'normal'),
        ('{0} (upper bound)'.format(args.vmiss), 'bold'),
        ('-- Hardy Weinberg Equilibrium Cutoff (all):', 'normal'),
        ('{0} (lower bound)'.format(args.hwe_all), 'bold'),
        ('-- Hardy Weinberg Equilibrium Cutoff (control):', 'normal'),
        ('{0} (lower bound)'.format(args.hwe_control), 'bold'),
        ('-- Hardy Weinberg Equilibrium Cutoff (case):', 'normal'),
        ('{0} (lower bound)'.format(args.hwe_case), 'bold'),
        ('-- Flip Scan:', 'normal'),
        ('{0}'.format(args.flipscan), 'bold'),
        ('-- Heterozygosity Cutoff:', 'normal'),
        ('{0} SD (upper bound)'.format(args.heter), 'bold'),
        ('-- Individual Missingness Cutoff:', 'normal'),
        ('{0} (upper bound)'.format(args.smiss), 'bold'),
        ('-- Check Sex:', 'normal'),
        ('{0}'.format(args.sex), 'bold'),
        ('-- Kinship Cutoff:', 'normal'),
        ('{0} (upper bound)'.format(args.kinship), 'bold'),
        ('-- Pruning Window:', 'normal'),
        ('{0}'.format(args.prune_window), 'bold'),
        ('-- Pruning Step:', 'normal'),
        ('{0}'.format(args.prune_step), 'bold'),
        ('-- Pruning Correlation Cutoff:', 'normal'),
        ('{0} (upper bound)'.format(args.prune_threshold), 'bold'),
        ('-- Population:', 'normal'),
        ('{0}'.format(args.population), 'bold'),
        ('-- Standard Deviation of Target Population:', 'normal'),
        ('{0}'.format(args.population_sd), 'bold')
    ]

    # content
    pdf.chapter_body("QC was processed using the following cutoffs:")
    pdf.ln(2)
    for text, form in text_list:
        if form == 'bold':
            pdf.set_font('Times', 'B', 12)
            pdf.cell(0, 5, text)
            pdf.ln(5)
        else:
            pdf.cell(7)
            pdf.set_font('Times', '', 12)
            pdf.cell(90, 5, text)
    pdf.ln(10)
    
    # record
    records['filter']['maf'] = args.maf
    records['filter']['vmiss'] = args.vmiss
    records['filter']['hwe_all'] = args.hwe_all
    records['filter']['hwe_case'] = args.hwe_case
    records['filter']['hwe_control'] = args.hwe_control
    records['filter']['flipscan'] = args.flipscan
    records['filter']['heter'] = args.heter
    records['filter']['smiss'] = args.smiss
    records['filter']['sex'] = args.sex
    records['filter']['kinship'] = args.kinship
    records['filter']['prune_window'] = args.prune_window
    records['filter']['prune_step'] = args.prune_step
    records['filter']['prune_threshold'] = args.prune_threshold
    records['filter']['population'] = args.population
    records['filter']['population_sd'] = args.population_sd


    ################## SNPs summary ################################
    # title
    pdf.chapter_title("3.1", "Summary of variant level QC") 
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # load record json file
    snp_filters = ['MAF', 'Missing', 'HWE_all', 'HWE_control', 'HWE_case', 'Flip']
    snp_record = json.load(open(args.snp_record, 'r'))
    snp_df = pd.DataFrame(snp_record['Filters'])
    snp_df = snp_df.set_index('name', drop=True)
    snp_df = snp_df.rename(index={
        'MAF_Filter': 'MAF',
        'GENO_Filter': 'Missing',
        'HWE_Filter_ALL': 'HWE_all',
        'HWE_Filter_CONTROL': 'HWE_control',
        'HWE_Filter_CASE': 'HWE_case',
        'FLIP_Filter': 'Flip'
    })
    for f in snp_filters:
        if f not in snp_df.index:
            snp_df.loc[f] = 0

    records['snp_qc']['Final'] = snp_record['Final_num']
    records['snp_qc']['Dup'] = records['input']['snp_num'] - int(Subprocess('wc -l {}'.format(args.dedup_record)))

    # text and format
    text_list = [
        (str(records['snp_qc']['Final']), 'bold'),
        (' out of ', 'normal'),
        (str(records['input']['snp_num']), 'bold'),
        (' pass variant QC procedures. ', 'normal'),
        ('The table below shows the number of unqualified variants in each criterion.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)
    
    # table
    header = ['', 'Dup & Chr']
    content = ['# SNPs', str(records['snp_qc']['Dup'])]
    for f in snp_filters:
        header.append(f)
        content.append(str(snp_df.loc[f, 'remove']))
        records['snp_qc'][f] = int(snp_df.loc[f, 'remove'])

    ## column width
    col_width = [max(pdf.get_string_width(str(header[i])), pdf.get_string_width(str(content[i]))) for i in range(len(header))]
    
    ## header
    table_width = sum(col_width) + inter_space * len(header) * 2
    text_size = 12
    pdf.set_font('Times', 'B', text_size)
    pdf.cell((190 - table_width)/2)
    pdf.set_draw_color(22, 143, 153)
    pdf.set_text_color(22, 143, 153)
    for i in range(len(header)):
        if i == 0:
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=header[i], align='C', border='')
        else:
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=header[i], align='C', border='L')
    pdf.ln(6)

    ## content
    pdf.cell((190 - table_width)/2)
    pdf.set_text_color(51, 51, 51)
    for i in range(len(content)):
        if i == 0:
            pdf.set_font('Times', 'B', text_size)
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=content[i], align='C', border='T')
        else:
            pdf.set_font('Times', '', text_size)
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=content[i], align='C', border='LT')
    pdf.ln(15)


    ################## Sample summary #################################
    # title
    pdf.chapter_title("3.2", "Summary of individual level QC") 
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # load record json file
    ind_filters = ['Sex', 'Missing', 'Heterozygosity']
    ind_record = json.load(open(args.ind_record, 'r'))
    ind_df = pd.DataFrame(ind_record['Filters'])
    ind_df = ind_df.set_index('name', drop=True)
    ind_df = ind_df.rename(index={
        'SEXCHECK_Filter': 'Sex',
        'MIND_Filter': 'Missing',
        'HET_Filter': 'Heterozygosity'
    })
    for f in ind_filters:
        if f not in ind_df.index:
            ind_df.loc[f] = 0

    records['ind_qc']['Kinship'] = int(Subprocess('wc -l {}'.format(args.kinship_record))) - 1
    records['ind_qc']['Final'] = int(Subprocess('wc -l {}'.format(args.population_record)))
    records['ind_qc']['Population'] = ind_record['Final_num'] - records['ind_qc']['Kinship'] - records['ind_qc']['Final']

    # text and format
    text_list = [
        (str(records['ind_qc']['Final']), 'bold'),
        (' out of ', 'normal'),
        (str(records['input']['ind_num']), 'bold'),
        (' pass individual QC procedures. ', 'normal'),
        ('The table below shows the number of unqualified samples in each criterion.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)
    
    # table
    header = ['']
    content = ['# Samples']
    for f in ind_filters:
        header.append(f)
        content.append(str(ind_df.loc[f, 'remove']))
        records['ind_qc'][f] = int(ind_df.loc[f, 'remove'])
    header += ['Kinship', 'Population']
    content += [str(records['ind_qc']['Kinship']), str(records['ind_qc']['Population'])]

    ## column width
    col_width = [max(pdf.get_string_width(str(header[i])), pdf.get_string_width(str(content[i]))) for i in range(len(header))]
    
    ## header
    table_width = sum(col_width) + inter_space * len(header) * 2
    text_size = 12
    pdf.set_font('Times', 'B', text_size)
    pdf.cell((190 - table_width)/2)
    pdf.set_draw_color(22, 143, 153)
    pdf.set_text_color(22, 143, 153)
    for i in range(len(header)):
        if i == 0:
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=header[i], align='C', border='')
        else:
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=header[i], align='C', border='L')
    pdf.ln(6)

    ## content
    pdf.cell((190 - table_width)/2)
    pdf.set_text_color(51, 51, 51)
    for i in range(len(content)):
        if i == 0:
            pdf.set_font('Times', 'B', text_size)
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=content[i], align='C', border='T')
        else:
            pdf.set_font('Times', '', text_size)
            pdf.cell(col_width[i] + 2*inter_space, 5, txt=content[i], align='C', border='LT')


    ################## SNPs QC #################################
    Title = 'Variant level QC'
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(args.logo_file)

    ##### Dup & Chr
    # title
    pdf.chapter_title("4.1", "Variant level QC --- Duplication and Invalid Chromosome")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)
    
    # text and format
    text_list = [
        ('This is our first step of GWAS QC. ', 'normal'),
        ('For duplicated variants, the first one is left and the others are removed. ', 'normal'),
        ('For the chromosome, autosomal chromosomes from Chr1 to Chr22 and X chromosome are valid in our QC. ', 'normal'),
        ('Y and other extra chromosme variants are removed.', 'normal'),
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)


    ##### Flip scan
    # title
    pdf.chapter_title("4.2", "Variant level QC --- Flip Scan")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # text and format
    text_list = [
        ('The genotype of controls and cases are checked for strand inconsistencies. ', 'normal'),
        ('Variants which cannot be fixed by flipping are removed.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)


    ##### MAF
    # title
    pdf.chapter_title("4.3", "Variant level QC --- Minor Allele Frequency (MAF)")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # text and format
    text_list = [
        ('MAF is the frequency of the least often occurring allele at a specific location. ', 'normal'),
        ('Most studies are underpowered to detect associations with SNPs with a low MAF and therefore exclude these SNPs. ', 'normal'),
        ('In this study, the MAF threshold is set at ', 'normal'),
        (str(args.maf), 'bold'),
        (', and ', 'normal'),
        (str(snp_df.loc['MAF', 'remove']), 'bold'),
        (' SNPs with MAF lower than {0} are removed. '.format(args.maf), 'normal'),
        ('Plot shown below illustrates the histogram of MAF distribution, and the region colored in red indicates the number of removed SNPs.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)

    # figure
    pdf.image(args.fig_maf, x=15, w=180, h=120)
    records['fig']['maf'] = args.fig_maf


    ##### Variant Missingness
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(args.logo_file)

    # title
    pdf.chapter_title("4.4", "Variant level QC --- Missingness")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # text and format
    text_list = [
        ('Variant level missingness is the number of individuals in the sample for whom information on a specific SNP is missing. ', 'normal'),
        ('SNPs with a high level of missingness can potentially lead to bias. ', 'normal'),
        ('The threshold for SNP level missingness is set at ', 'normal'),
        (str(args.vmiss), 'bold'),
        (', and ', 'normal'),
        (str(snp_df.loc['Missing', 'remove']), 'bold'),
        (' SNPs are filtered out.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)

    # figure
    pdf.image(args.fig_vmiss, x=15, w=180, h=120)
    records['fig']['vmiss'] = args.fig_vmiss


    ##### Hardy Weinberg Equilibrium
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(args.logo_file)

    # title
    pdf.chapter_title("4.5", "Variant level QC --- Hardy Weinberg Equilibrium (HWE)")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # text and format
    text_list = [
        ('HWE concerns the relationship between allele and genotype frequencies. ', 'normal'),
        ('Violation of the HWE can be attributed to inbreeding, contamination, or other forms of non-random mating. ', 'normal'),
        ('The P-value threshold of Hardy-Weinberg equilibrium test is set at ', 'normal')
    ]
    if args.hwe_strict == 'true':
        text_list += [('{} and {} for controls and cases respectively'.format(args.hwe_control, args.hwe_case), 'bold')]
    else:
        text_list += [('{}'.format(args.hwe_all), 'bold')]
    text_list += [
        (', resulting in ', 'normal'),
        (str(snp_df.loc[['HWE_all', 'HWE_control', 'HWE_case'], 'remove'].sum()), 'bold'),
        (' SNPs that violate the HWE law.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)

    # figure
    if args.hwe_strict == 'true':
        pdf.image('{}.hwe_control.hist.png'.format(args.fig_hwe), x=25, w=150, h=100)
        pdf.image('{}.hwe_case.hist.png'.format(args.fig_hwe), x=25, y=180, w=150, h=100)
        records['fig']['hwe'] = ['{}.hwe_control.hist.png'.format(args.fig_hwe), '{}.hwe_case.hist.png'.format(args.fig_hwe)]
    else:
        pdf.image('{}.hwe.hist.png'.format(args.fig_hwe), x=15, w=180, h=120)
        records['fig']['hwe'] = ['{}.hwe.hist.png'.format(args.fig_hwe)]


    ################## Sample QC #################################
    Title = "Individual level QC"
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(args.logo_file)


    ##### Sex
    # title
    pdf.chapter_title("5.1", "Individual level QC --- Sex Discrepancy")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # text and format
    text_list = [
        ('Checks for sex discrepancies between the record and X chromosome heterozygosity/homozygosity rates. ', 'normal'),
        (str(ind_df.loc['Sex', 'remove']), 'bold'),
        (' sample(s) with discordant sex information are removed.', 'normal')
    ]
    
    # content
    pdf.content(text_list)
    pdf.ln(10)


    ##### Individual Missingness
    # title
    pdf.chapter_title("5.2", "Individual level QC --- Missingness")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # text and format
    text_list = [
        ('Individual level missingness is the number of SNPs that are missing for a specific individual. ', 'normal'),
        ('High levels of missingness indicates a technical problem or poor quality within the data. ', 'normal'),
        ('The threshold of missingness is set at ', 'normal'),
        (str(args.smiss), 'bold'),
        ('. In this dataset, ', 'normal'),
        (str(ind_df.loc['Missing', 'remove']), 'bold'),
        (' individual(s) with high rates of genotype missingness were excluded.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)


    ##### Heterozygosity
    # title
    pdf.chapter_title("5.3", "Individual level QC --- Heterozygosity")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # text and format
    text_list = [
        ('Heterozygosity is the carrying of two different alleles of a specific SNP. ', 'normal'),
        ('The heterozygosity rate of an individual is the proportion of heterozygous genotypes. ', 'normal'),
        ('High levels of heterozygosity within an individual might be an indication of low sample quality, ', 'normal'),
        ('whereas low levels of heterozygosity may be due to inbreeding. ', 'normal'),
        ('In this step, ', 'normal'),
        (str(ind_df.loc['Heterozygosity', 'remove']), 'bold'),
        (' individuals who exceed ± ', 'normal'),
        (str(args.heter), 'bold'),
        (" SD from the samples' heterozygosity rate mean are excluded. ", 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)

    # text and format
    text_list = [
        ('Heterozygosity and missingness are plotted on the next page based on all GWAS samples analysed in this study. ', 'normal'),
        ('The red horizontal line indicates {0} standard deviations from the mean of heterozygosity, '.format(args.heter), 'normal'),
        (' and the red vertical line indicates the missingness threshold of {0}'.format(args.smiss), 'normal')
    ]
    
    # content
    pdf.content(text_list)
    pdf.ln(10)


    ##### Relationship Inference
    # titile
    pdf.chapter_title("5.4", "Individual level QC --- Kinship Coefficient")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # text and format
    text_list = [
        ('The kinship coefficient represents the degree of relation between a pair of individuals. ', 'normal'),
        ('Using this criterion with the threshold of ', 'normal'),
        (str(args.kinship), 'bold'),
        (', ', 'normal'),
        (str(records['ind_qc']['Kinship']), 'bold'),
        (' individuals are removed.', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)

    # plot of het-vs-imiss
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(args.logo_file)
    pdf.image(args.fig_het_vs_smiss, x=15, w=180, h=120)
    records['fig']['het_vs_smiss'] = args.fig_het_vs_smiss


    ##### Population
    pdf.add_page()
    pdf.page_header(Title)
    pdf.logo(args.logo_file)

    # titile
    pdf.chapter_title("5.5", "Individual level QC --- Population Stratification")
    pdf.set_text_color(51, 51, 51)
    pdf.set_draw_color(51, 51, 51)

    # text and format
    text_list = [
        ('Population stratification is the presence of a systematic difference in allele frequencies between subpopulaitons. ', 'normal'),
        (str(records['ind_qc']['Population']), 'bold'),
        (' samples with different ethnic backgroud ', 'normal'),
        ('(SD > {})'.format(args.population_sd), 'bold'),
        (' are removed. ', 'normal'),
        ('The following plot depicts the distribution of the data on the first two dimension computed by principal component analysis (PCA).', 'normal')
    ]

    # content
    pdf.content(text_list)
    pdf.ln(10)

    # figure
    pdf.image('{}.C1_C2.png'.format(args.fig_population_prefix), x=15, w=180, h=120)
    pdf.image('{}.C1_C3.png'.format(args.fig_population_prefix), x=15, y=200, w=90, h=60)
    pdf.image('{}.C2_C3.png'.format(args.fig_population_prefix), x=105, y=200, w=90, h=60)
    records['fig']['population'] = '{}.C1_C2.png'.format(args.fig_population_prefix)


    ################## Report #################################
    # PDF
    pdf.output('{0}.qc_report.pdf'.format(basename), 'F')
    
    # record json
    json.dump(records, open('{0}.qc_record.json'.format(basename), 'w'))


if __name__ == "__main__":
    main()
