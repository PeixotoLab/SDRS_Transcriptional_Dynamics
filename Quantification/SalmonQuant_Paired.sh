#!/bin/sh

#SBATCH --partition=peixoto
#SBATCH --job-name=SalmonAlign_HC5
#SBATCH --error=SalmonAlign_HC5.err
#SBATCH --output=SalmonAlign_HC5.out
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1

module load salmon

salmon quant \
--threads 6 \
--index SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 FASTQ/WTHC5_PFC_1_R1.fastq.gz --mates2 FASTQ/WTHC5_PFC_1_R2.fastq.gz \
--output Output/WTHC5_PFC_1_quant \
--numBootstraps 30 \
&

salmon quant \
--threads 6 \
--index SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 FASTQ/WTHC5_PFC_3_R1.fastq.gz --mates2 FASTQ/WTHC5_PFC_3_R2.fastq.gz \
--output Output/WTHC5_PFC_3_quant \
--numBootstraps 30 \
&

salmon quant \
--threads 6 \
--index SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 FASTQ/WTHC5_PFC_4_R1.fastq.gz --mates2 FASTQ/WTHC5_PFC_4_R2.fastq.gz \
--output Output/WTHC5_PFC_4_quant \
--numBootstraps 30 \
&

salmon quant \
--threads 6 \
--index SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 FASTQ/WTHC5_PFC_5_R1.fastq.gz --mates2 FASTQ/WTHC5_PFC_5_R2.fastq.gz \
--output Output/WTHC5_PFC_5_quant \
--numBootstraps 30 \
&

salmon quant \
--threads 6 \
--index SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
--mates1 FASTQ/WTHC5_PFC_6_R1.fastq.gz --mates2 FASTQ/WTHC5_PFC_6_R2.fastq.gz \
--output Output/WTHC5_PFC_6_quant \
--numBootstraps 30 \
&

wait