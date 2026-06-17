path=/home/qiaoy/proj_nobackup/1.MultiReg/1.WGS/2.variantcalling/14.haplotag/haplotag_sh
pathbam="/home/qiaoy/proj_nobackup/1.MultiReg/1.WGS/2.variantcalling/14.haplotag"

for i in P{1..7}; do 

cd $pathbam/$i

patient=${i##*/}

for sample in P*; do 

echo "$patient $sample"

touch $path/$patient/${sample}.Count.sh

echo "#!/bin/bash -l

#SBATCH -A sens2017114
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J count${sample}
#SBATCH --output=slurm_count${sample}.out
#SBATCH --error=slurm_count${sample}.err


## ml bioinfo-tools samtools/1.20

htagpath=/home/qiaoy/proj_nobackup/1.MultiReg/1.WGS/2.variantcalling/14.haplotag
toolpath=/home/qiaoy/proj_nobackup/0.tools/4.shell/

cd \$htagpath/$patient/$sample


for bam in *WGS.bam; do

sh \$toolpath/count_halpotag_DNA.sh \$bam 50000 \${bam/.bam/.count}

done


for count in *WGS.count; do 

sh \$toolpath/covert_reads_bulk.sh \$count \${count/.count/.polished.count}

done 

" > $path/$patient/${sample}.Count.sh




done

done

