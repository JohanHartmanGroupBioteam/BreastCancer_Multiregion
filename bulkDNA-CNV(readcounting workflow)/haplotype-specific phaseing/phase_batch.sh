path=`pwd`
bampath=/home/qiaoy/proj_nobackup/1.MultiReg/1.WGS/1.preprocess/01.DuplicatesMarked

cd $bampath

for i in *Germline*.md.bam ; do

sample=${i%%.*}

patient=${i%%-*}

echo "$patient ; $sample"

mkdir -p $path/$patient
touch ${path}/${patient}/${chr}_phase.sh

for chr in chr{1..22} chrX chrY; do

echo "#!/bin/bash -l

#SBATCH -A sens2017114
#SBATCH -p node
#SBATCH -n 2
#SBATCH -t 30:00
#SBATCH -J phase_${patient}_${chr}
#SBATCH --output=slurm_${patient}_${chr}.out
#SBATCH --error=slurm_${patient}_${chr}.err


## how to set environment
## module load bioinfo-tools GATK/4.3.0.0 bcftools/1.17 samtools/1.20
## module load conda
## source activate gatk

# ml bioinfo-tools SHAPEIT/v4.2.2 WhatsHap/2.3-20240529-be88057 gcc/12.3.0 bcftools/1.17 


pathbam=/home/qiaoy/proj_nobackup/1.MultiReg/1.WGS/1.preprocess/01.DuplicatesMarked
pathbamWGBS=/proj/sens2017114/data/MultiRegion/WGBS/results/deduplicated
pathbamRNA=/home/qiaoy/proj_nobackup/1.MultiReg/2.RNA/4.results/0.mapped
ref=/sw/data/uppnex/ToolBox/hg38bundle/Homo_sapiens_assembly38.fasta
dbsnp=/sw/data/uppnex/ToolBox/hg38bundle/dbsnp_146.hg38.vcf.gz
hapmap=/sw/data/ToolBox/hg38bundle/hapmap_3.3.hg38.vcf.gz
Mills=/sw/data/ToolBox/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
funcotator=/castor/project/proj_nobackup/wharf/qiaoy/qiaoy-sens2017114/1.MultiReg/1.WGS/9.genome/1.Human/funcotator/funcotator_dataSourc
phasepath=/home/qiaoy/proj_nobackup/1.MultiReg/1.WGS/2.variantcalling/10.phase
toolpath="/home/qiaoy/proj_nobackup/0.tools/4.shell"

## haplotype

mkdir -p \${phasepath}/$patient
cd \${phasepath}/$patient

<<COMMENT

gatk HaplotypeCaller \\
-R \$ref \\
-I \$pathbam/${patient}*Germline*bqsr.bam \\
--native-pair-hmm-threads 16 \\
--intervals $chr \\
-O ${chr}.htp.vcf.gz

## This is 2D with pre-trained architecture
gatk CNNScoreVariants \\
-I \$pathbam/${patient}*Germline*bqsr.bam \\
--variant ${chr}.htp.vcf.gz \\
--reference \$ref \\
--intervals $chr \\
--output ${chr}.htp.cnn.vcf.gz \\
--tensor-type read_tensor

## filter
gatk FilterVariantTranches \\
--variant ${chr}.htp.cnn.vcf.gz \\
--resource \${hapmap} \\
--resource \${Mills} \\
--info-key CNN_2D \\
--intervals $chr \\
--snp-tranche 99.95 --indel-tranche 99.4 \\
--output ${chr}.htp.filtered.vcf.gz


bcftools view -m2 -M2 -f PASS -i 'GT=\"het\" & FORMAT/DP > 15' -Oz -o ${chr}.HQ.vcf.gz ${chr}.htp.filtered.vcf.gz
tabix -f -p vcf ${chr}.HQ.vcf.gz

COMMENT

## whatshap
#module unload
#ml bioinfo-tools WhatsHap/2.3-20240529-be88057 samtools
#ml gcc/12.3.0

ref=/sw/data/uppnex/ToolBox/hg38bundle/Homo_sapiens_assembly38.fasta

<<COMMENT

whatshap phase -o ${chr}.wtsp.vcf.gz --reference=\${ref} --tag=PS \\
--ignore-read-groups \\
${chr}.HQ.vcf.gz \$pathbam/${patient}*.md.bam \$pathbamWGBS/${patient}*.sorted.bam \$pathbamRNA/${patient}*/*.sorted.bam
tabix -f -p vcf ${chr}.wtsp.vcf.gz

COMMENT

## manual phase
zcat ${chr}.wtsp.vcf.gz | \\
awk -v OFS=\"\t\" '/^#/{print} \\
!/^#/{if (index(\$9, \"PS\") == 0 ){ \$9=\$9\":PS\"; \$10=\$10\":\"\$2; gsub(/\//, \"|\", \$10); print;} else print }' | \\
bgzip -c > ${chr}.wtsp.MP.vcf.gz
tabix -f -p vcf ${chr}.wtsp.MP.vcf.gz


## shapeit
#module unload
#ml bioinfo-tools SHAPEIT/v4.2.2

map="/home/qiaoy/proj_nobackup/1.MultiReg/1.WGS/9.genome/1.Human/phaseGRCh38/map/genetic_maps.b38"
reference="/home/qiaoy/proj_nobackup/1.MultiReg/1.WGS/9.genome/1.Human/phaseGRCh38/reference"


## remoce chr prefix
zcat ${chr}.wtsp.MP.vcf.gz | sed 's/##contig=<ID=chr\([0-9XYM]\+\),/##contig=<ID=\1,/' | \\
awk '/^#/{print} /^chr/{gsub(\"chr\",\"\"); print }' | bgzip > ${chr}.wtsp.temp.vcf.gz
tabix -f -p vcf ${chr}.wtsp.temp.vcf.gz

<<COMMENT

shapeit4 \\
--thread 8 \\
--input ${chr}.wtsp.temp.vcf.gz \\
--map \$map/${chr}.b38.gmap.gz \\
--region ${chr/chr/} \\
--use-PS 0 \\
--reference \$reference/ALL.${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \\
--output ${chr}.shapeit.vcf.gz \\
--bingraph ${chr}.shapeit.bingraph

# <<COMMENT
# COMMENT

## adding chr
zcat ${chr}.shapeit.vcf.gz | sed 's/##contig=<ID=\([0-9XYM]\+\)>/##contig=<ID=chr\1>/' | \\
awk -v OFS=\"\t\" '/^#/{print} !/^#/{\$1=\"chr\"\$1; print}' | bgzip > ${chr}.shapeit.chr.vcf.gz

mv ${chr}.shapeit.chr.vcf.gz ${chr}.shapeit.vcf.gz
tabix -f -p vcf ${chr}.shapeit.vcf.gz

rm ${chr}.shapeit.chr.vcf*

COMMENT

#### manual change PS tag ####
module load python/3.12.7
python3 \$toolpath/find_PS_to_change.py ${chr}.shapeit.vcf.gz ${chr}.wtsp.temp.vcf.gz ${chr}.test.vcf


( zcat ${chr}.wtsp.vcf.gz | grep \"#\" ; awk -v OFS=\"\t\" 'NR > 1 {\$1=\"chr\"\$1; print }' ${chr}.test.vcf ) | bgzip > ${chr}.wtsp.Manul.vcf.gz
tabix -f -p vcf ${chr}.wtsp.Manul.vcf.gz

rm ${chr}.test.vcf ${chr}.wtsp.temp.vcf.gz* ${chr}.wtsp.MP.vcf.gz*

" > ${path}/${patient}/${chr}_phase.sh

done

done


