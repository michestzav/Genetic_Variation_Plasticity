#!/bin/bash

module load anaconda/2021.11
conda activate /default/software/gemma

#==============================================================================
# Record start time and RAM usage
#==============================================================================

start_time=$(date +"%Y-%m-%d %H:%M:%S")
start_ram=$(free -m | awk 'NR==2{printf "%.2fGB", $3/1024}')

gemma="gemma"


#==============================================================================
# 1. Trait BLUPs
#==============================================================================

cd /Missing_data09_Plasticity/Traits_blups

### RUN GWAS-GEMMA with LMM
# Pheno is organized in R into FAM files

PHENO_ORDER=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21")

PHENO_NAMES=(
"BLUP_log_Dry_mass"
"BLUP_log_Area"
"BLUP_LMA_g_m2"
"BLUP_log_gsw"
"BLUP_PS2"
"BLUP_ETR"
"BLUP_Fs"
"BLUP_log_Fm"
"BLUP_GC_U"
"BLUP_GC_L"
"BLUP_PL_U"
"BLUP_PL_L"
"BLUP_SD_upper"
"BLUP_SD_lower"
"BLUP_SDR"
"BLUP_Ocurrence_upper"
"BLUP_d13C"
"BLUP_C"
"BLUP_d15N"
"BLUP_N"
"BLUP_CN"
)

#==============================================================================
# Create kinship matrix
#==============================================================================

gemma -bfile ALL_Admixed_FAI -gk 1 -o ALL_Admixed_kinship

#==============================================================================
# Run GEMMA-LMM for each phenotype
#==============================================================================

for index in ${!PHENO_ORDER[*]}; do
    $gemma \
        -bfile ALL_Admixed_FAI \
        -n ${PHENO_ORDER[$index]} \
        -k ALL_Admixed_kinship.cXX.txt \
        -lmm 4 \
        -c AdmixtureProp_COV.txt \
        -o gemma_lmm_blups_${PHENO_NAMES[$index]}
done


#==============================================================================
# 2. Plasticity slopes
#==============================================================================

cd /Missing_data09_Plasticity/Slopes

### RUN GWAS-GEMMA with LMM
# Pheno is organized in R into FAM files

PHENO_ORDER=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22")

PHENO_NAMES=(
"slope_Dry_mass"
"slope_Area"
"slope_LMA_g_m2"
"slope_gsw"
"slope_PS2"
"slope_ETR"
"slope_Fs"
"slope_Fm"
"slope_GC_U"
"slope_GC_L"
"slope_PL_U"
"slope_PL_L"
"slope_SD_upper"
"slope_SD_lower"
"slope_SDR"
"slope_Ocurrence_upper"
"slope_d13C"
"slope_C"
"slope_d15N"
"slope_N"
"slope_CN"
"slope_Ocurrence_upper_2"
)

#==============================================================================
# Run GEMMA-LMM for each phenotype
#==============================================================================

for index in ${!PHENO_ORDER[*]}; do
    $gemma \
        -bfile ALL_Admixed_FAI \
        -n ${PHENO_ORDER[$index]} \
        -k ALL_Admixed_kinship.cXX.txt \
        -lmm 4 \
        -c AdmixtureProp_COV.txt \
        -o gemma_lmm_slope_${PHENO_NAMES[$index]}
done


#==============================================================================
# 3. Scaled plasticity slopes
#==============================================================================

cd /Missing_data09_Plasticity/Scale_slopes

### RUN GWAS-GEMMA with LMM
# Pheno is organized in R into FAM files

PHENO_ORDER=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19")

PHENO_NAMES=(
"scaled_Dry_mass"
"scaled_Area_m2"
"scaled_LMA_g_m2"
"scaled_gsw"
"scaled_d13C"
"scaled_PS2"
"scaled_ETR"
"scaled_Fm"
"scaled_Fs"
"scaled_PL_L"
"scaled_PL_U"
"scaled_SD_lower"
"scaled_SD_upper"
"scaled_SDR"
"scaled_Ocurrence_upper"
"scaled_C"
"scaled_N"
"scaled_d15N"
"scaled_CN"
)

#==============================================================================
# Run GEMMA-LMM for each phenotype
#==============================================================================

# --- run only phenotypes 1..19 ---
for n in $(seq 1 19); do
  i=$((n-1))
  echo "Running GEMMA -n $n  ->  ${PHENO_NAMES[$i]}"

  "$gemma" -bfile ALL_Admixed_FAI \
    -n "$n" \
    -k ALL_Admixed_kinship.cXX.txt \
    -lmm 4 \
    -c AdmixtureProp_COV.txt \
    -o "gemma_lmm_scaled_${PHENO_NAMES[$i]}"
done