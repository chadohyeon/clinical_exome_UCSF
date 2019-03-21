#!/usr/bin/python

import os,sys,argparse
from os.path import expanduser

#infile='/opt/vep/.vep/inputs/100917_proband_het.vcf'
#outpath='/opt/vep/.vep/outputs/'
with open ("config_texts/dbNSFP_cols.txt") as rf:
    dbNSFP_cols=[i.strip() for i in rf.readlines()]

def main(infile, outpath, number_forks, buffer_size):
    outfile = outpath + infile.split('/')[-1].replace('.vcf', '.vep_annot.vcf')
    # Set the run
    if '/Users/dcha' in expanduser("~"):
        vep_path = 'docker run -it -v ~/program/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep ./vep'
        custom_path = '/opt/vep/.vep/custom/'
        plugin_path = '/opt/vep/.vep/Plugins/'
    elif '/home/ec2-user' in expanduser("~") or '/root' in expanduser("~"):
        vep_path = '/opt/program/ensembl-vep/vep'
        custom_path = '/root/.vep/custom/'
        plugin_path = '/root/.vep/Plugins/'
    else:
        print(expanduser("~"))
        sys.exit(0)   # exit without ahy errors

    cmd = [vep_path,
            #'--cache',
            '--assembly GRCh38 --offline',
            #'--tab',
            '--vcf',
            #'--field '+fields,
            '--fork', number_forks,
            '--force_overwrite',
            '--buffer_size', buffer_size,
            '-i', infile,
            '-o', outfile,
            '--no_stats',     #Don't generate a stats file. Provides marginal gains in run time.
            '--polyphen b',
            '--sift b',
            '--ccds',
            #'--hgvs', # it adds 50% run time as it checks the fasta file
            '--numbers',
            '--canonical',
            '--protein',
            '--biotype',
            '--uniprot',
            '--tsl',
            '--appris',
            '--max_af',
            '--af',
            '--af_1kg',
            '--af_gnomad',
            '--pubmed',
            '--pick --pick_order canonical,appris,tsl,biotype,ccds,rank,length',
            '--no_intergenic',
            '--symbol'
            ]
    # Add plugins
    cmd = cmd + ['--plugin ExACpLI',
                 '--plugin LoFtool',
                 '--plugin dbNSFP,{0}dbNSFP3.5a.gz,'.format(plugin_path)+','.join(dbNSFP_cols),
                ]
    #if use_loftee:
    #cmd = cmd + ['--plugin LoF,loftee_path:{0},conservation_file:{0}loftee.sql,human_ancestor_fa:{0}human_ancestor.fa.gz\,gerp_bigwig:{0}gerp_conservation_scores.homo_sapiens.GRCh38.bw'.format(plugin_path)]
   
    # Add cytoband annotation
    cmd = cmd + [','.join(['--custom ' + custom_path + 'cytoBand.bed.gz', 'cytoBand', 'bed', 'overlap', '0'])]
    # Add additional allele frequency annotations
    cmd = cmd + [','.join(['--custom ' + custom_path + 'gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz', 'gnomADg', 'vcf', 'exact', '0',
                           'AC', 'AF', 'AF_Male', 'AF_Female', 'AF_raw', 'AF_POPMAX', 'AF_AFR', 'AF_AMR', 'AF_ASJ', 'AF_EAS', 'AF_FIN', 'AF_NFE', 'AF_OTH', 'AF_SAS']),
                 ','.join(['--custom ' + custom_path + 'ExAC.r0.3.nonpsych.sites.vcf.gz',  'ExAC_nonpsych', 'vcf', 'exact', '0',
                           'AC', 'AF', 'AN', 'AC_Adj', 'AN_Adj', 'AC_AFR',  'AC_AMR',  'AC_EAS', 'AC_FIN', 'AC_NFE', 'AC_OTH', 'AC_SAS', 'AN_AFR', 'AN_AMR', 'AN_EAS', 'AN_FIN', 'AN_NFE', 'AN_OTH', 'AN_SAS']),
                 ]

    # Add public databases annotations
    cmd = cmd + [#','.join(['--custom ' + custom_path + 'HGMD-PUBLIC_Ensembl_95.bed.gz', 'HGMD', 'bed', 'exact', '0']),
                 ','.join(['--custom ' + custom_path + 'MIM_morbid_Ensembl_95.bed.gz', 'OMIM\(Gene\)', 'bed', 'overlap', '0']),
                 ','.join(['--custom ' + custom_path + 'ClinVar_Variation_Ensembl_95.bed.gz', 'ClinVar\(SNV\)', 'bed', 'exact', '0']),
                 #','.join(['--custom ' + custom_path + 'ClinVar_StructuralVariation_Ensembl_95.bed.gz', 'ClinVar\(SV\)', 'bed', 'overlap', '0']),
                 ','.join(['--custom ' + custom_path + 'Orphanet_Ensembl_95.bed.gz', 'Orphanet\(Gene\)', 'bed', 'overlap', '0']),
                 #','.join(['--custom ' + custom_path + 'NHGRI-EBI_GWAS_catalog_Ensembl_95.bed.gz', 'GWAS_catalog', 'bed', 'exact', '0']),
                 ','.join(['--custom ' + custom_path + 'DDG2P_Ensembl_95.bed.gz', 'DDG2P\(Gene\)', 'bed', 'overlap', '0']),
                 ','.join(['--custom ' + custom_path + 'COSMIC_Ensembl_95.bed.gz', 'COSMIC', 'bed', 'exact', '0']),
                 ','.join(['--custom ' + custom_path + 'Cancer_Gene_Census_Ensembl_95.bed.gz', 'Cancer_Gene_Census\(Gene\)', 'bed', 'overlap', '0']),
                 ','.join(['--custom ' + custom_path + 'dbGaP_Ensembl_95.bed.gz', 'dbGaP', 'bed', 'exact', '0']),
                 #','.join(['--custom ' + custom_path + 'MAGIC_Ensembl_95.bed.gz', 'MAGIC', 'bed', 'exact', '0']),
                 #','.join(['--custom ' + custom_path + 'Teslovich_Ensembl_95.bed.gz', 'Teslovich', 'bed', 'exact', '0']),
                 #','.join(['--custom ' + custom_path + 'GIANT_Ensembl_95.bed.gz', 'GIANT', 'bed', 'exact', '0']),
                 #','.join(['--custom ' + custom_path + 'DGVa_Ensembl_95.bed.gz', 'DGVa\(SV\)', 'bed', 'overlap', '0']),
                 ]
    # Add conservations scores annotations
    cmd = cmd + [','.join(['--custom ' + custom_path + 'hg38.phastCons100way.bw', 'phastCons100way', 'bigwig', 'overlap', '0']),
                 ','.join(['--custom ' + custom_path + 'hg38.phyloP100way.bw', 'phyloP100way', 'bigwig', 'overlap', '0']),
                 ','.join(['--custom ' + custom_path + 'hg38.phastCons20way.bw', 'phastCons20way', 'bigwig', 'overlap', '0']),
                 ','.join(['--custom ' + custom_path + 'hg38.phyloP20way.bw', 'phyloP20way', 'bigwig', 'overlap', '0']),
                 ','.join(['--custom ' + custom_path + 'hg38.phastCons7way.bw', 'phastCons7way', 'bigwig', 'overlap', ' 0']),
                 ','.join(['--custom ' + custom_path + 'hg38.phyloP7way.bw', 'phyloP7way', 'bigwig', 'overlap', '0']),
                # ','.join(['-custom ' + custom_path + 'phastCons46way.primates.hg19ToHg38.bw','phastCons46wayPr','bigwig','overlap','0']),
                # ','.join(['-custom ' + custom_path + 'phyloP46way.primates.hg19ToHg38.bw','phyloP46wayPr','bigwig','overlap','0']),
                # ','.join(['-custom ' + custom_path + 'phastCons46way.vertebrate.hg19ToHg38.bw','phastCons46wayVt','bigwig','overlap','0']),
                # ','.join(['-custom ' + custom_path + 'phyloP46way.vertebrate.hg19ToHg38.bw','phyloP46wayVt','bigwig','overlap','0']),
                # ','.join(['-custom ' + custom_path + 'LINSIGHT.hg19ToHg38.bw','LINSIGHT','bigwig','overlap','0']),
                ]
    cmd = ' '.join(cmd)
    os.system(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--infile', required=True, type=str, help='Input File')
    parser.add_argument('-o', '--outpath', required=True, type=str, help='Output Path')
    #parser.add_argument('--use_loftee', dest='use_loftee', action='store_true')
    parser.add_argument('-f', '--number_forks', required=False, type=str, help='Number of forks', default='1')
    parser.add_argument('-b', '--buffer_size', required=False, type=str, help='Buffer size', default='5000')
    args = parser.parse_args()
    main(args.infile, args.outpath, args.number_forks, args.buffer_size)

