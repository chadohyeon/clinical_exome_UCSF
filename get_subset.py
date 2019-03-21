import vcf_main_vep
import sys, os

family_file="config_texts/ClinEx_pedigree.txt"

def get_subset(fn, family_file, output_path):
    # Converting a joint-vcf file into a dictionary
    # 'info_annot_file' contains the annotated columns of interest in INFO column
    with open (fn) as _raw_vcf:
        wf = open(output_path + fn.split('/')[-1].replace(".vcf", ".proband_het_hom_and_X.vcf"), 'w')
        raw_vcf=_raw_vcf.readlines()
        double_hashtag=True
        for i in raw_vcf:
            if double_hashtag:
                if i.startswith("##"): wf.write(i)
                else:
                    double_hashtag=False
                    wf.write(i)
                    header=i.strip().split('\t')
                    family_dict=vcf_main_vep.family_check(header[-3:], family_file)
            else:
                x=i.strip().split('\t')
                if x[0] not in ['chr'+str(num) for num in list(range(1,23))+['X','Y']]: continue


                d_format = x[header.index(family_dict['fatherID'])].split(":")
                m_format = x[header.index(family_dict['motherID'])].split(":")
                p_format = x[header.index(family_dict['probandID'])].split(":")
                member_keys = x[header.index('FORMAT')].split(':')

                if not (len(d_format)==len(m_format)==len(p_format)==len(member_keys)): continue
                if not {'GT', 'GQ', 'AD', 'DP'}.issubset(set(member_keys)): continue

                if x[0]=='chrX':
                    if p_format[0]=='0/1': wf.write(i)
                    elif p_format[0]=='1/1': wf.write(i)
                else:
                    if p_format[0]=='0/1': wf.write(i)
                    elif p_format[0]=='1/1' and d_format[0]!='1/1' and m_format[0]!='1/1': wf.write(i)
        wf.close()

if __name__=="__main__":
    vcf_input = sys.argv[1]
    outpath = sys.argv[2]
    get_subset(vcf_input, family_file, outpath)
