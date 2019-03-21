#!/bin/bash
source ~/.bash_profile

input_path=/root/.vep/inputs/
annot_path=/root/.vep/annot/
output_path=/root/.vep/outputs/

interm_path=/root/.vep/interm/

input_vcfs=$(ls $input_path) 
for vcf in $input_vcfs;do
	echo "$(date +"%T") 1) Getting a subset only including proband het and hom variants"
	mkdir $interm_path
	python get_subset.py ${input_path}${vcf} $interm_path
	echo "$(date +"%T") 2) Annotating .vcf with vep"
	bash loftee_only.sh ${interm_path}${vcf//.vcf/.proband_het_hom_and_X.vcf} ${annot_path}${vcf//.vcf/.proband_het_hom_and_X.loftee.vcf} 4 25000
	python vep_annot.py -i ${interm_path}${vcf//.vcf/.proband_het_hom_and_X.vcf} \
			-o ${annot_path} \
			-f 4 \
			-b 25000
	echo "$(date +"%T") 3) Generating comp het, hemi, hom and inherited het .txt files"
	rm -r ${interm_path}
	python comp_het_vep.py ${annot_path}${vcf//.vcf/.proband_het_hom_and_X.vep_annot.vcf} ${output_path}
	python hemi_vep.py ${annot_path}${vcf//.vcf/.proband_het_hom_and_X.vep_annot.vcf} ${output_path}
	python hom_vep.py ${annot_path}${vcf//.vcf/.proband_het_hom_and_X.vep_annot.vcf} ${output_path}
	python inherited_het_vep.py ${annot_path}${vcf//.vcf/.proband_het_hom_and_X.vep_annot.vcf} ${output_path}
done
echo "$(date +"%T") Jobs finished"
