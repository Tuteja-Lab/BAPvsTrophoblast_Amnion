## Pre-processing 

### SRA datasets

```bash
while read a b; do
python3 enaBrowserTools-1.6/python3/enaDataGet.py -f fastq -as  ~/.aspera_setting.ini $a;
done<assets/needed-sra.txt 
```
### Genome/annotation

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v37.primary_assembly.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```


## Mapping

Script to index (`star-index.sh`):

```bash
#!/bin/bash
ml star
GTF=gencode.v37.primary_assembly.annotation.gtf
FASTA=GRCh38.primary_assembly.genome.fa
STAR --runThreadN 36 \
     --runMode genomeGenerate \
     --genomeDir GRCh38 \
     --sjdbOverhang 100 \
     --sjdbGTFfile ${GTF} \
     --genomeFastaFiles ${FASTA}
```

 Script to Map/Count (`star-map.sh`):

 ```bash
 #!/bin/bash
ml star
index=/work/LAS/geetu-lab/arnstrm/amnion-rnaseq/lo-etal/GRCh38
read1=$1
cpus=36
out=$(basename ${read1} | sed 's/.fastq.gz//g')
STAR \
--runThreadN 36 \
--genomeDir $index \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.3 \
--outFileNamePrefix ${out}_ \
--readFilesCommand zcat \
--readFilesIn ${read1}
```

prepare jobs:

```bash
find $(pwd) -name ".fastq.gz" > input.fofn
```

Slurm job script (`star.slurm`):

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=8:00:00
#SBATCH --job-name=star
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=arnstrm@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
IFS=$'\n' read -d '' -r -a LINES < input.fofn
LINE=${LINES[$SLURM_ARRAY_TASK_ID]}
echo -e "$(date +"%D  %r")\tthis command is ${LINE}"
./star-map.sh ${LINE}
```

Submit jobs:

```bash
sbatch --array=0-72 star.slurm
```

## Counts

```bash
for f in *_ReadsPerGene.out.tab; do
g=$(echo $f |sed 's/_ReadsPerGene.out.tab//g');
echo -e "AA_geneids\t$g" > ${g}_counts.txt;
tail -n +5 $f | cut -f1,4 >> ${g}_counts.txt;
done
join_files.sh *_counts.txt | grep -v "_counts.txt" | sort -k1,1 > counts-full.txt
# filter to retain only coding-genes

# rename SRR ids with info
```


