# interactive job

source activate bcftools

# 06_stats_atl

grep "After filtering" slurm-19155079_1.out | grep "Sites" | cut -d ' ' -f4 > 06_stats_atl.snps

grep "After filtering" slurm-19155079_1.out | grep "Sites" | cut -d ' ' -f9 > 06_stats_atl.sites

# 06_stats_car

grep "After filtering" slurm-19155079_2.out | grep "Sites" | cut -d ' ' -f4 > 06_stats_car.snps

grep "After filtering" slurm-19155079_2.out | grep "Sites" | cut -d ' ' -f9 > 06_stats_car.sites

# 06_stats_cat

grep "After filtering" slurm-19155079_3.out | grep "Sites" | cut -d ' ' -f4 > 06_stats_cat.snps

grep "After filtering" slurm-19155079_3.out | grep "Sites" | cut -d ' ' -f9 > 06_stats_cat.sites

# 06_stats_myi

grep "After filtering" slurm-19155079_4.out | grep "Sites" | cut -d ' ' -f4 > 06_stats_myi.snps

grep "After filtering" slurm-19155079_4.out | grep "Sites" | cut -d ' ' -f9 > 06_stats_myi.sites

# 06b_stats_atl

grep "After filtering" slurm-19155079_5.out | grep "Sites" | cut -d ' ' -f4 > 06b_stats_atl.snps

grep "After filtering" slurm-19155079_5.out | grep "Sites" | cut -d ' ' -f9 > 06b_stats_atl.sites

# 06b_stats_car

grep "After filtering" slurm-19155079_6.out | grep "Sites" | cut -d ' ' -f4 > 06b_stats_car.snps

grep "After filtering" slurm-19155079_6.out | grep "Sites" | cut -d ' ' -f9 > 06b_stats_car.sites

# 06b_stats_cat

grep "After filtering" slurm-19155079_7.out | grep "Sites" | cut -d ' ' -f4 > 06b_stats_cat.snps

grep "After filtering" slurm-19155079_7.out | grep "Sites" | cut -d ' ' -f9 > 06b_stats_cat.sites

# 06b_stats_myi

grep "After filtering" slurm-19155079_8.out | grep "Sites" | cut -d ' ' -f4 > 06b_stats_myi.snps

grep "After filtering" slurm-19155079_8.out | grep "Sites" | cut -d ' ' -f9 > 06b_stats_myi.sites



# in R

R

x_files <- list.files(pattern="*snps")
for(a in 1:length(x_files)) {
	a_rep <- scan(x_files[a])
	writeLines(paste0(x_files[a], "  ", sum(a_rep)))
}

x_files <- list.files(pattern="*sites")
for(a in 1:length(x_files)) {
	a_rep <- scan(x_files[a])
	writeLines(paste0(x_files[a], "  ", sum(a_rep)))
}















