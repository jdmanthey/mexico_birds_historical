
# interactive job

interactive -p nocona -c 4 -m 8G

source activate bcftools

workdir=/lustre/scratch/jmanthey/08_mexico/

cd ${workdir}/05_pca

# cat the vcf files for each species
grep "^#" Atl_structure_CM051598.1.recode.vcf > Atl_structure.vcf
for i in $( ls Atl*recode.vcf ); do
	grep -v "^#" $i >> Atl_structure.vcf
done


grep "^#" Car_structure_CM027537.1.recode.vcf > Car_structure.vcf
for i in $( ls Car*recode.vcf ); do
	grep -v "^#" $i >> Car_structure.vcf
done


grep "^#" Cat_structure_NC_046258.1.recode.vcf > Cat_structure.vcf
for i in $( ls Cat*recode.vcf ); do
	grep -v "^#" $i >> Cat_structure.vcf
done


grep "^#" Myi_structure_CM027537.1.recode.vcf > Myi_structure.vcf
for i in $( ls Myi*recode.vcf ); do
	grep -v "^#" $i >> Myi_structure.vcf
done

# remove the unneeded files
rm *recode*


##########################################
##########################################
### run plink for pca
### and run ADMIXTURE
##########################################
##########################################

# make chromosome map
grep -v "#" Atl_structure.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map_atl.txt
grep -v "#" Car_structure.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map_car.txt
grep -v "#" Cat_structure.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map_cat.txt
grep -v "#" Myi_structure.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map_myi.txt

#plink output format
vcftools --vcf Atl_structure.vcf  --plink --chrom-map chrom_map_atl.txt --out structure_atl 
vcftools --vcf Car_structure.vcf  --plink --chrom-map chrom_map_car.txt --out structure_car
vcftools --vcf Cat_structure.vcf  --plink --chrom-map chrom_map_cat.txt --out structure_cat 
vcftools --vcf Myi_structure.vcf  --plink --chrom-map chrom_map_myi.txt --out structure_myi 

# convert  with plink
plink --file structure_atl --recode12 --allow-extra-chr --out structure_atl_plink
plink --file structure_car --recode12 --allow-extra-chr --out structure_car_plink
plink --file structure_cat --recode12 --allow-extra-chr --out structure_cat_plink
plink --file structure_myi --recode12 --allow-extra-chr --out structure_myi_plink

# run pca on dataset
plink --file structure_atl_plink --pca --allow-extra-chr --out structure_atl_plink_pca
plink --file structure_car_plink --pca --allow-extra-chr --out structure_car_plink_pca
plink --file structure_cat_plink --pca --allow-extra-chr --out structure_cat_plink_pca
plink --file structure_myi_plink --pca --allow-extra-chr --out structure_myi_plink_pca

# convert with plink for admixture dataset
plink --file structure_atl --recode12 --allow-extra-chr --out admixture_atl_plink
plink --file structure_car --recode12 --allow-extra-chr --out admixture_car_plink
plink --file structure_cat --recode12 --allow-extra-chr --out admixture_cat_plink
plink --file structure_myi --recode12 --allow-extra-chr --out admixture_myi_plink

# run admixture
for K in 1 2 3 4 5 6 7; do admixture --cv admixture_atl_plink.ped $K  | tee log_atl_${K}.out; done
for K in 1 2 3 4 5 6 7; do admixture --cv admixture_car_plink.ped $K  | tee log_car_${K}.out; done
for K in 1 2 3 4 5 6 7; do admixture --cv admixture_cat_plink.ped $K  | tee log_cat_${K}.out; done
for K in 1 2 3 4 5 6 7; do admixture --cv admixture_myi_plink.ped $K  | tee log_myi_${K}.out; done








