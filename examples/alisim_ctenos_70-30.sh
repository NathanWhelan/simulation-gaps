#!/usr/bin/env bash

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=2
#SBATCH --hint=multithread
#SBATCH --partition=whelan
#SBATCH --ntasks=10
#SBATCH --mem=32GB
#SBATCH --job-name=alisim
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Whelan.105@osu.edu

# cd to dir from which I submit job
cd $SLURM_SUBMIT_DIR

#commands

##WAG
iqtree3 --alisim D16_choano_empirical-gaps_conflict_cteno-sister_70-30_WAG+C10+G1 -s Dataset16_Choano_Uncertain_slowestHalf_cteno-out_gene1-out.phy -te ../Dataset16_Choano_Uncertain_slowestHalf.tre -m WAG+C10+G{1.0} -pre D16_choano_empirical-gaps_conflict_cteno-sister_70-30_WAG+C10+G1 --site-freq SAMPLING --site-rate SAMPLING -nt 10 --seed 1 -blfix

##VT
iqtree3 --alisim D16_choano_empirical-gaps_conflict_cteno-sister_70-30_VT+C10+G05 -s Dataset16_Choano_Uncertain_slowestHalf_cteno-out_gene2-out.phy -te ../Dataset16_Choano_Uncertain_slowestHalf.tre -m VT+C10+G{0.5} -pre D16_choano_empirical-gaps_conflict_cteno-sister_70-30_VT+C10+G05 --site-freq SAMPLING --site-rate SAMPLING -nt 10 --seed 2 -blfix

##JTTDCMut
iqtree3 --alisim D16_choano_empirical-gaps_conflict_cteno-sister_70-30_JTTDCMut+C10+G02 -s Dataset16_Choano_Uncertain_slowestHalf_cteno-out_gene3-out.phy -te ../Dataset16_Choano_Uncertain_slowestHalf.tre -m JTTDCMut+C10+G{0.2} -pre D16_choano_empirical-gaps_conflict_cteno-sister_70-30_JTTDCMut+C10+G02 --site-freq SAMPLING --site-rate SAMPLING -nt 10 --seed 2 -blfix

##mtZOA
iqtree3 --alisim D16_choano_empirical-gaps_conflict_cteno-sister_70-30_mtZOA+C10+G3 -s Dataset16_Choano_Uncertain_slowestHalf_cteno-out_gene4-out.phy -te ../Dataset16_Choano_Uncertain_slowestHalf.tre -m mtZOA+C10+G{3.0} -pre D16_choano_empirical-gaps_conflict_cteno-sister_70-30_mtZOA+C10+G3 --site-freq SAMPLING --site-rate SAMPLING -nt 10 --seed 4 -blfix

##NQYEAST
iqtree3 --alisim D16_choano_empirical-gaps_conflict_cteno-sister_70-30_NQYEAST+C10+G08 -s Dataset16_Choano_Uncertain_slowestHalf_cteno-out_gene5-out.phy -te ../Dataset16_Choano_Uncertain_slowestHalf.tre -m NQ.yeast+C10+G{0.8} -pre D16_choano_empirical-gaps_conflict_cteno-sister_70-30_NQ.yeast+C10+G08 --site-freq SAMPLING --site-rate SAMPLING -nt 10 --seed 6 -blfix