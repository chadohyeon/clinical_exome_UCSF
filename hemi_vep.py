#!/usr/bin/python
import vcf_main_vep
import sys

vcf_input="/Users/dcha/program/vep_data/files/100917_proband_het_vep_annot.vcf"
info_annot_file="config_texts/INFO_annot_list_vep.txt"
gene_summary_annot_file="config_texts/gene_summary_annot_list.txt"
family_file="config_texts/ClinEx_pedigree.txt"
gene_additional_annot_list='config_texts/dbNSFP3.5_gene.complete'
gene_annot_columns='config_texts/gene_summary_annot_list_dbNSFP.txt'
par_configure="config_texts/par_region.txt"

def hemizygous(vcf, qc_text, par_text):
    chrx={i:vcf[i] for i in vcf.keys() if vcf[i]['#CHROM']=='chrX'\
                                       and vcf[i]['father']['GT']=='0/0'\
                                       and vcf[i]['mother']['GT']=='0/1'\
                                       and vcf[i]['proband']['GT']!='0/0'}

    # Step2-1: Only containing exonic variants without synonymous SNV or unknown, and splicing variants
    chrx_exome=vcf_main_vep.get_exonic_variants(chrx)

    # Step2-2: QC based on given qc criteria configuration file.
    chrx_exome_qc=vcf_main_vep.get(qc_main(chrx_exome))
    
    # Step2-3: Considering zygosity of variants within pseudo-autosomal regions
    chrx_exome_qc_par={}
    with open (par_text) as _hg38_par:
        hg38_par={i.split(':')[0]:i.split(':')[1].strip().split('-') for i in _hg38_par.readlines()}
    for i in chrx_exome_qc.keys():
        if chrx_exome_qc[i]['proband']['GT']=='0/1':
            if int(hg38_par['PAR1'][0])<int(chrx_exome_qc[i]['POS'])<int(hg38_par['PAR1'][1]):
                chrx_exome_qc_par[i]=chrx_exome_qc[i]
            elif int(hg38_par['PAR2'][0])<int(chrx_exome_qc[i]['POS'])<int(hg38_par['PAR2'][1]):
                chrx_exome_qc_par[i] = chrx_exome_qc[i]
        else: chrx_exome_qc_par[i] = chrx_exome_qc[i]
    del chrx, chrx_exome, chrx_exome_qc
    return chrx_exome_qc_par


def gene_summary(family_dict, gene_summary_annot_file, output_path):
    with open(output_path + family_dict['familyID']+".vep_annot.hemi.variant_summary.txt") as rf:
        with open(output_path + family_dict['familyID']+".vep_annot.hemi.gene_summary.txt", 'w') as wf:
            header=rf.readline().strip().split('\t')
            variants=[i.strip().split('\t') for i in rf.readlines()]
            genes = {}
            for i in variants:
                if i[header.index('SYMBOL')] not in genes.keys(): genes[i[header.index('SYMBOL')]]=1
                else: genes[i[header.index('SYMBOL')]]+=1
            with open(gene_summary_annot_file) as rf:
                header_to_write_basic=[lines.strip() for lines in rf.readlines()]
            header_to_write_variant_info=['Variants_number', 'Variants_detail','LOFTEE','Priority']
            header_to_write_family_info=['FatherID', 'MotherID', 'ProbandID', 'ProbandSex']
            wf.write('\t'.join(['temp_pos']+header_to_write_basic + header_to_write_variant_info + header_to_write_family_info)+'\n')
            count=0
            gene_count=0
            priority="G2"
            li_to_write = []
            variant_info = []
            loftee=[]
            for i in variants:
                if count==0:
                    li_to_write.append(i[header.index('POS')])
                    current_gene=i[header.index('SYMBOL')]
                    gene_count=genes[current_gene]
                    for j in header_to_write_basic: li_to_write.append(i[header.index(j)])
                    li_to_write.append(str(gene_count))
                    variant_info.append(variant_type_for_gene_summary(i, header))
                    loftee.append(i[header.index("LOFTEE")])
                    count += 1
                    if gene_count==1:
                        for k in variant_info:
                            if k.split("///")[1]=="HIGH": priority="G1"
                        li_to_write.append(variant_info[0].split('///')[0])
                        li_to_write.append(loftee[0])
                        li_to_write.append(priority)
                        li_to_write+=[i[header.index('FatherID')], i[header.index('MotherID')], i[header.index('ProbandID')], i[header.index('ProbandSex')]]
                        wf.write('\t'.join(li_to_write) + '\n')
                        variant_info = []
                        loftee=[]
                        priority="G2"
                        li_to_write = []
                        count = 0
                elif count==gene_count-1:
                    variant_info.append(variant_type_for_gene_summary(i, header))
                    loftee.append(i[header.index("LOFTEE")])
                    for k in variant_info:
                        if k.split('///')[1]=='HIGH': priority='G1'
                    li_to_write.append(' |'.join([var.split('///')[0] for var in variant_info]))
                    li_to_write.append(' |'.join(loftee))
                    li_to_write.append(priority)

                    li_to_write+=[i[header.index('FatherID')], i[header.index('MotherID')], i[header.index('ProbandID')], i[header.index('ProbandSex')]]
                    wf.write('\t'.join(li_to_write)+'\n')
                    variant_info=[]
                    loftee=[]
                    priority = "G2"
                    li_to_write=[]
                    count=0
                else:
                    variant_info.append(variant_type_for_gene_summary(i,header))
                    loftee.append(i[header.index("LOFTEE")])
                    count+=1

def variant_type_for_gene_summary(i,header):
    gt_pair = (i[header.index('Father_GT')], i[header.index('Mother_GT')], i[header.index('Proband_GT')])
    variant_type = i[header.index('Consequence')]
    impact=i[header.index('IMPACT')]
    if gt_pair == ('0/0', '0/1', '0/1'):
        return '0:1:1;'+variant_type+'///'+impact
    elif gt_pair == ('0/0', '0/1', '1/1'):
        return '0:1:2;'+variant_type+'///'+impact

if __name__=="__main__":
    vcf_input = sys.argv[1]
    outpath = sys.argv[2]
    vcf_dict,fam_dict=vcf_main_vep.vcf_to_dict(vcf_input, info_annot_file, family_file)
    if fam_dict['probandSex']=='male':
        hemi_dict=hemizygous(vcf_dict, qc_configure, par_configure)
        vcf_main_vep.variant_summary(hemi_dict, fam_dict, info_annot_file, "hemi", outpath)
        gene_summary(fam_dict, gene_summary_annot_file, outpath)
        vcf_main_vep.gene_additional_annotation(outpath + fam_dict['familyID'] + ".vep_annot.hemi.gene_summary.txt",
                                                gene_additional_annot_list, gene_annot_columns)
        print("   Hemizygous ... Done")
    else:
        print("   Hemizygous analysis was not done for the female proband")
