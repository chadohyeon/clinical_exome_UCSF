#!/usr/bin/python
import vcf_main_vep
import sys

vcf_input="/Users/dcha/program/vep_data/files/100917_proband_het_vep_annot.vcf"
info_annot_file="config_texts/INFO_annot_list_vep.txt"
gene_summary_annot_file="config_texts/gene_summary_annot_list.txt"
family_file="config_texts/ClinEx_pedigree.txt"
gene_additional_annot_list='config_texts/dbNSFP3.5_gene.complete'
gene_annot_columns='config_texts/gene_summary_annot_list_dbNSFP.txt'

def compound_heterozygous(vcf):

    # Step1: Listing all of the proband heterozygous variants
    het={i:vcf[i] for i in vcf.keys() if vcf[i]['proband']['GT']=='0/1'}

    # Step2-1: Only containing exonic variants without synonymous SNV or unknown, and splicing variants
    het_exome=vcf_main_vep.get_exonic_variants(het)

    # Step2-2: QC based on given qc criteria configuration file.
    het_exome_qc=vcf_main_vep.qc_main(het_exome)

    #het_exome_qc_AF=vcf_main_vep.global_AF_filter(het_exome_qc, float(qc_criteria['AF']))           # This code is for AF filter, but won't be used
    # Step2-3: Only including variants within genes that both include F:M:P=1:0:1 and F:M:P=0:1:1, excluding genes with parental homozygous variants.
    genes={}
    for i in het_exome_qc.keys():
        gene_name=het_exome_qc[i]['INFO']['SYMBOL']
        if gene_name not in genes.keys(): genes[gene_name]=1
        else: genes[gene_name]+=1
    dup_genes=[i for i in genes.keys() if genes[i]>1]
    het_exome_qc_dup={i:het_exome_qc[i] for i in het_exome_qc.keys() if het_exome_qc[i]['INFO']['SYMBOL'] in dup_genes}
    gene_zygosity={i:{"1:0:1":0, "0:1:1":0, "parental_hom":0} for i in dup_genes}
    for i in het_exome_qc_dup.keys():
        gene=het_exome_qc_dup[i]['INFO']['SYMBOL']
        gt_pair = (het_exome_qc_dup[i]['father']['GT'], het_exome_qc_dup[i]['mother']['GT'], het_exome_qc_dup[i]['proband']['GT'])
        if gt_pair == ("0/1", "0/0", "0/1"):
            gene_zygosity[gene]["1:0:1"] += 1
        elif gt_pair == ("0/0", "0/1", "0/1"):
            gene_zygosity[gene]["0:1:1"] += 1
        elif gt_pair == ("1/1", "0/0", "0/1") or gt_pair == ("0/0", "1/1", "0/1"):
            gene_zygosity[gene]["parental_hom"] +=1
    comp_genes=[i for i in gene_zygosity.keys() if gene_zygosity[i]["1:0:1"]>=1 and gene_zygosity[i]["0:1:1"]>=1 and gene_zygosity[i]["parental_hom"]==0]
    het_exome_qc_dup_comp={}
    for i in het_exome_qc_dup.keys():
        gt_pair=(het_exome_qc_dup[i]['father']['GT'], het_exome_qc_dup[i]['mother']['GT'], het_exome_qc_dup[i]['proband']['GT'])
        if het_exome_qc_dup[i]['INFO']['SYMBOL'] in comp_genes:
            if gt_pair==("0/1", "0/0", "0/1") or gt_pair==("0/0", "0/1", "0/1"):
                het_exome_qc_dup_comp[i]=het_exome_qc_dup[i]
   # del het, het_exome, het_exome_qc, het_exome_qc_dup, genes, dup_genes,gene_zygosity, comp_genes
    return het_exome_qc_dup_comp


def gene_summary(family_dict, gene_summary_annot_file, output_path):
    with open(output_path + family_dict['familyID']+".vep_annot.comp_het.variant_summary.txt") as rf:
        with open(output_path + family_dict['familyID']+".vep_annot.comp_het.gene_summary.txt", 'w') as wf:
            header=rf.readline().strip().split('\t')
            variants=[i.strip().split('\t') for i in rf.readlines()]
            genes = {}
            for i in variants:
                if i[header.index('SYMBOL')] not in genes.keys(): genes[i[header.index('SYMBOL')]]=1
                else: genes[i[header.index('SYMBOL')]]+=1
            with open(gene_summary_annot_file) as rf:
                header_to_write_basic=[lines.strip() for lines in rf.readlines()]
            header_to_write_variant_info=['Variants_number', 'Variants_detail','LOFTEE', 'Priority']
            header_to_write_family_info=['FatherID', 'MotherID', 'ProbandID', 'ProbandSex']
            wf.write('\t'.join(['temp_pos']+header_to_write_basic + header_to_write_variant_info + header_to_write_family_info)+'\n')
            count=0
            gene_count=0
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
                    count+=1
                elif count==gene_count-1:
                    variant_info.append(variant_type_for_gene_summary(i, header))
                    loftee.append(i[header.index("LOFTEE")])
                    zygosity_dict={"1:0:1":[], "0:1:1":[]}
                    for j in variant_info:
                        zygosity_dict[j.split(";")[0]].append(j.split("///")[1])
                    paternal_LoF=False
                    maternal_LoF=False
                    if 'HIGH' in zygosity_dict['1:0:1']: paternal_LoF=True
                    if 'HIGH' in zygosity_dict['0:1:1']: maternal_LoF=True
                    if paternal_LoF and maternal_LoF: priority="G1"
                    elif paternal_LoF: priority="G2"
                    elif maternal_LoF: priority="G2"
                    else: priority="G3"
                    li_to_write.append(' |'.join([var.split('///')[0] for var in variant_info]))
                    li_to_write.append(' |'.join(loftee))
                    li_to_write.append(priority)
                    li_to_write+=[i[header.index('FatherID')], i[header.index('MotherID')], i[header.index('ProbandID')], i[header.index('ProbandSex')]]
                    wf.write('\t'.join(li_to_write)+'\n')
                    variant_info=[]
                    loftee=[]
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
    if gt_pair == ('0/1', '0/0', '0/1'):
        return '1:0:1;'+variant_type+'///'+impact
    elif gt_pair == ('0/0', '0/1', '0/1'):
        return '0:1:1;'+variant_type+'///'+impact

if __name__=="__main__":
    vcf_input = sys.argv[1]
    outpath = sys.argv[2]
    vcf_dict,fam_dict=vcf_main_vep.vcf_to_dict(vcf_input, info_annot_file, family_file)
    print("Analyzing family ID "+fam_dict['familyID'])
    comp_het_dict=compound_heterozygous(vcf_dict)
    vcf_main_vep.variant_summary(comp_het_dict, fam_dict, info_annot_file, "comp_het", outpath)
    gene_summary(fam_dict, gene_summary_annot_file, outpath)
    vcf_main_vep.gene_additional_annotation(outpath + fam_dict['familyID']+".vep_annot.comp_het.gene_summary.txt", gene_additional_annot_list, gene_annot_columns)
    print("   Compound heterozygous ... Done")
