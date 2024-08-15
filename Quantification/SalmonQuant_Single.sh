#!/bin/sh

#SBATCH --partition=peixoto
#SBATCH --job-name=SalmonAlign_HC3
#SBATCH --error=SalmonAlign_HC3.err
#SBATCH --output=SalmonAlign_HC3.out
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1

module load salmon

salmon quant \
--threads 6 \
--index SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
-r FASTQ/WTHC3_N2.fastq.gz \
--output Output/WTHC3_1_quant \
--numBootstraps 30 \
&

salmon quant \
--threads 6 \
--index SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
-r FASTQ/WTHC3_N3.fastq.gz \
--output Output/WTHC3_2_quant \
--numBootstraps 30 \
&

salmon quant \
--threads 6 \
--index SalmonQuant/Gencode_vM28/SalmonIndexFiles/gencode.vM28-salmon-index-v1.0.0-mouse-withdecoys \
--libType A \
-r FASTQ/WTHC3_N4.fastq.gz \
--output Output/WTHC3_3_quant \
--numBootstraps 30 \
&

wait