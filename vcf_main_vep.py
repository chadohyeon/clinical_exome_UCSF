#!/usr/bin/python
import os

info_annot_basic="config_texts/INFO_annot_list_basic.txt"
qc_text="config_texts/qc_configure.txt"

def family_check(li, family_file):
    # Reading the ID columns in the joint-vcf and matching it to the pedigree data
    family_dict={}
    family_id_not_add=True
    with open(family_file) as _raw_family:
        header=_raw_family.readline().strip().split('\t')
        raw_family=[i.strip().split('\t') for i in _raw_family.readlines()]
        for i in raw_family:
            for j in li:
                if j in i[header.index('GMI_ID')] or j in i[header.index('SampleID')] or j in i[header.index('Identifier')]:
                    if family_id_not_add:
                        family_dict['familyID']=i[header.index('FamilyID')]
                        family_id_not_add=False
                    if i[header.index('Role')]=='father': family_dict['fatherID']=j
                    elif i[header.index('Role')]=='mother': family_dict['motherID']=j
                    elif i[header.index('Role')]=='proband':
                        family_dict['probandID']=j
                        family_dict['probandSex']=i[header.index('Sex')]
    # Example
    return family_dict


def vcf_to_dict(fn, info_annot_file, family_file):
    # Converting a joint-vcf file into a dictionary
    # 'info_annot_file' contains the annotated columns of interest in INFO column
    with open (fn.replace('.vep_annot','.loftee')) as _loftee:
        loftee_vcf={}
        loftee=_loftee.readlines()
        for i in loftee:
            if i.startswith('##'):
                if "Consequence annotations from Ensembl VEP." in i:
                    loftee_vep_header=i.replace('">\n','').split("Format: ")[1].split("|")
                else: continue
            elif i.startswith('#CHROM'): header=i.replace('\n','').split('\t')
            else:
                x=i.strip().split('\t')
                key=x[header.index('#CHROM')]+':'+x[header.index('POS')]+':'+x[header.index('REF')]+'>'+x[header.index('ALT')]
                loftee_raw={k.split('=')[0]:k.split('=')[1] for k in i.split('\t')[header.index('INFO')].split(';') if '=' in k}['CSQ'].split('|')
                loftee_csq={loftee_vep_header[num]:loftee_raw[num] for num in range(len(loftee_vep_header))}
                if loftee_csq['LoF']=="": loftee_csq['LoF']="."
                if loftee_csq['LoF_filter']=="": loftee_csq['LoF_filter']="."
                loftee_vcf[key]=loftee_csq
        del loftee
 
    with open (fn) as _raw_vcf:
        raw_vcf=_raw_vcf.readlines()
        vcf_dict={}
        double_hashtag=True
        for i in raw_vcf:
            if double_hashtag:
                if i.startswith("##"):
                    if "Consequence annotations from Ensembl VEP." in i:
                        vep_header=i.replace('">\n','').split("Format: ")[1].split("|")
                    else: continue
                else:
                    double_hashtag=False
                    family_dict=family_check(header[-3:], family_file)
            else:
                x=i.strip().split('\t')
                key=x[header.index('#CHROM')]+':'+x[header.index('POS')]+':'+x[header.index('REF')]+'>'+x[header.index('ALT')]
                inner_dict={j:x[header.index(j)] for j in ['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER']}
                info_raw={k.split('=')[0]:k.split('=')[1] for k in x[header.index('INFO')].split(';') if '=' in k}

                with open(info_annot_basic) as _info_basic:
                    info_annot_list=[lines.strip() for lines in _info_basic.readlines()]
                info_dict={l:info_raw[l] for l in info_annot_list if l in info_raw.keys()}

                if 'CSQ' not in info_raw.keys(): continue
                csq_split=[]
                for ann in info_raw['CSQ'].split('|'):
                    if ann=='':
                        ann='.'
                    csq_split.append(ann)
                vep_raw=list(zip(vep_header,csq_split))
                
                vep_raw+=[("LOFTEE",loftee_vcf[key]['LoF']),('LOFTEE_filter',loftee_vcf[key]['LoF_filter'])]

                with open(info_annot_file) as _vep_annot:
                    vep_annot_list=[lines.strip() for lines in _vep_annot.readlines()]
                for j in vep_raw:
                    if j[0] in vep_annot_list: info_dict[j[0]]=j[1]
                inner_dict['INFO']=info_dict

                member_keys = x[header.index('FORMAT')].split(':')

                for m in ['father', 'mother', 'proband']:
                    indiv_dict={}
                    member_values=x[header.index(family_dict[m+'ID'])].split(':')
                    indiv_dict['GT'] = member_values[member_keys.index('GT')]
                    if indiv_dict == './.': continue
                    indiv_dict['GQ'] = member_values[member_keys.index('GQ')]
                    DP_calc, AB_calc=calculateAB(member_values[member_keys.index('AD')])
                    indiv_dict["DP"] = DP_calc
                    indiv_dict["AB"] = AB_calc
                    inner_dict[m]=indiv_dict
                inner_dict["familyID"]=family_dict["familyID"]
                vcf_dict[key]=inner_dict
        return vcf_dict, family_dict


def calculateAB(ad):
    ref_count, alt_count = ad.split(',')[0], ad.split(',')[1]
    if ref_count == '.':
        ref_count = 0
    else:
        pass
    if alt_count == '.':
        alt_count = 0
    else:
        pass
    dp = int(ref_count) + int(alt_count)
    if dp == 0:
        ab = 0
    else:
        ab = int(alt_count) / float(dp)
    return str(dp), str(ab)

def get_exonic_variants(vcf)                     
    vcf_output={}
    for i in vcf.keys():
        if vcf[i]['INFO']['IMPACT']=='HIGH' or vcf[i]['INFO']['IMPACT']=='MODERATE':
            vcf_output[i]=vcf[i]
    return vcf_output

def qc_for_each_member(variant,qc_criteria, GQ_check=False):
    #QC based on each member's GQ and AB based on given qc criteria
    try:
        if GQ_check:
            if float(variant['GQ'] < float(qc_criteria['GQ']) : return False
        if variant['GT']=='0/1':
            if float(variant['AB']) >= float(qc_criteria['het_AB_LL']):
                if float(variant['AB']) <= float(qc_criteria['het_AB_UL']):
                    return True
        elif variant['GT']=='0/0':
            if float(variant['AB']) <= float(qc_criteria['wt_AB']): return True
        else:
            if float(variant['AB']) >= float(qc_criteria['hom_AB']): return True
    except: return False
    return False

def qc_GQ_MEAN(fam_variant, qc_criteria):
    try:
        GQ_SUM = float(fam_variant['father']['GQ']) + float(fam_variant['mother']['GQ']) + float(fam_variant['proband']['GQ'])
        GQ_MEAN = GQ_SUM/3.0
        if GQ_MEAN >= qc_criteria['GQ_MEAN']: return True
        else: return False
    except: return False

def qc_main(vcf)
    vcf_output={}
    with open (qc_text) as _qc_criteria:
        qc_criteria={i.split(':')[0]:i.split(':')[1].strip() for i in _qc_criteria.readlines()}
    for i in vcf.keys():
        if float(vcf[i]['QUAL'])>=float(qc_criteria['QUAL']):
            if vcf[i]['FILTER']==qc_criteria["FILTER"]:
                if qc_for_each_member(vcf[i]['proband'], qc_criteria):
                    if qc_for_each_member(vcf[i]['father'], qc_criteria):
                        if qc_for_each_member(vcf[i]['mother'], qc_criteria):
                            if qc_GQ_MEAN(vcf[i]):
                                vcf_output[i]=vcf[i]
    return vcf_output                
                     
def variant_summary(vcf, family_dict, info_annot_file, mode, output_path):
       # 'mode' should be comp_het or hemi
    with open(info_annot_basic) as _info_basic:
        with open(info_annot_file) as _info:
            info_annot_list = [lines.strip() for lines in _info_basic.readlines()]
            info_annot_list += [lines.strip() for lines in _info.readlines()]
            info_annot_list.insert(info_annot_list.index('OMIM(Gene)'),'ACMG_recommendations')
    with open("unsorted_variants.txt", "w") as wf:
        for i in vcf.keys():
            one_line=""
            for j in ['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER']:
                one_line+=vcf[i][j]+'\t'
            for k in info_annot_list:
                if k=="Consequence":
                    if ',' in vcf[i]['INFO'][k]: vcf[i]['INFO'][k]=vcf[i]['INFO'][k].split(',')[0]
                    elif '&' in vcf[i]['INFO'][k]: vcf[i]['INFO'][k]=vcf[i]['INFO'][k].split('&')[0]
                if k=='ACMG_recommendations':
                    with open('config_texts/ACMG_genes.txt') as acmg:
                        acmg_col = '.'
                        for acmg_line in [lines.strip() for lines in acmg.readlines()]:
                            acmg_li = acmg_line.split('\t')
                            if vcf[i]['INFO']['SYMBOL'] == acmg_li[0]:
                                acmg_col = acmg_li[1] + '///' + acmg_li[2]
                        vcf[i]['INFO'][k]=acmg_col
                try:    one_line+=vcf[i]['INFO'][k]+'\t'
                except: one_line+='.\t'
            for l in ['GT', 'GQ', 'DP', 'AB']:
                for m in ['father', 'mother', 'proband']:
                    one_line+=vcf[i][m][l]+'\t'
            one_line+="\t".join([family_dict['fatherID'], family_dict['motherID'], family_dict['probandID'], family_dict['probandSex']])
            one_line+='\n'
            wf.write(one_line)
    output_file_name=output_path + family_dict['familyID'] + ".vep_annot." + mode + ".variant_summary.txt"
    with open(output_file_name, "w") as wf:
        wf.write("\t".join(['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER']+info_annot_list+\
                           ['Father_GT','Mother_GT','Proband_GT','Father_GQ','Mother_GQ','Proband_GQ',\
                            'Father_DP','Mother_DP','Proband_DP','Father_AB','Mother_AB','Proband_AB',\
                            'FatherID','MotherID','ProbandID', 'ProbandSex'])+'\n') 
    os.system("cat unsorted_variants.txt | sort -k1,1V -k2,2n >> {0}".format(output_file_name))
    os.system("rm unsorted_variants.txt")

def gene_additional_annotation(previous_gene_summary, gene_additional_annot_list, annot_columns):
    total_dict={}
    mode=previous_gene_summary.split('vep_annot.')[1].split('.gene_summary')[0]
    with open(annot_columns) as _annot_col:
        annot_col=['temp_pos']+[line.strip() for line in _annot_col.readlines()]
    with open(gene_additional_annot_list) as _annot_db:
        db_header=_annot_db.readline().strip().split('\t')
        db_table=[line.strip().split('\t') for line in _annot_db.readlines()]
    with open (previous_gene_summary) as _prev_gene:
        header=_prev_gene.readline().strip().split('\t')
        table=[line.strip().split('\t') for line in _prev_gene.readlines()]
    for x in table:
        gene_id=x[header.index('Gene')]
        total_dict[gene_id]={header[j]:x[j] for j in range(len(header))}
    for y in db_table:
        if y[db_header.index('Ensembl_gene')] in total_dict.keys():
            for z in [k for k in annot_col if k not in header]:
                total_dict[y[db_header.index('Ensembl_gene')]][z]=y[db_header.index(z)]
    with open ('unsorted_genes.txt', 'w') as wf:
        for gene in total_dict.keys():
            li_to_write=[]
            for col in annot_col:
                if col in total_dict[gene].keys():
                    li_to_write.append(total_dict[gene][col])
                else:
                    li_to_write.append('.')
            li_to_write.append(mode)
            wf.write('\t'.join(li_to_write)+'\n')
    output_file_name=previous_gene_summary.replace('.txt','_extended.txt')
    with open ('temp_sorted_genes.txt', 'w') as wf:
        annot_col.append("Mode")
        wf.write('\t'.join(annot_col)+'\n')
    os.system('cat unsorted_genes.txt | sort -k2,2V -k1,1n >> temp_sorted_genes.txt')
    with open(output_file_name, 'w') as wf:
        with open('temp_sorted_genes.txt') as rf:
            for i in rf.readlines():
                wf.write('\t'.join(i.split('\t')[1:]))
    os.system('rm unsorted_genes.txt')
    os.system('rm temp_sorted_genes.txt')
    os.system('rm '+previous_gene_summary)
