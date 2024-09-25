![logo](https://github.com/jelber2/hapmers/blob/main/hapmers.svg)
# hapmers manuscript code

sbatch script for basecalling WGS_HG002_EZ1_25kb data

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/basecall.sbatch`

```bash
#!/bin/bash

#SBATCH -c 24
#SBATCH --mem 200G
#SBATCH -p gpu
#SBATCH -n 1
#SBATCH --gpus-per-task=a100:6
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
/msc/home/jelber43/git/dorado-0.6.0-linux-x64/bin/dorado basecaller sup ./ > WGS_HG002_EZ1_25kb.pod5.bam
```


```bash
sbatch basecall.sbatch
```



sbatch script for herro-pre WGS_HG002_EZ1_25kb data

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/herro-pre.sbatch`

```bash
#!/bin/bash

#SBATCH -c 248
#SBATCH --mem 1900G
#SBATCH -p gpu
#SBATCH -n 1
#SBATCH --gpus-per-task=a100:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
eval "$(/msc/home/jelber43/.local/bin/micromamba shell hook --shell bash)"
micromamba activate herro
samtools view -h -@4 WGS_HG002_EZ1_25kb.pod5.bam |samtools fastq -@ 4 | pigz -p 128 > WGS_HG002_EZ1_25kb.pod5.bam.fastq.gz
seqkit seq -n -i WGS_HG002_EZ1_25kb.pod5.bam.fastq.gz > WGS_HG002_EZ1_25kb.pod5.bam.fastq.read.ids
/msc/home/jelber43/git/herro/scripts/create_batched_alignments.sh WGS_HG002_EZ1_25kb.pod5.bam.fastq.gz WGS_HG002_EZ1_25kb.pod5.bam.fastq.read.ids 248 batches/
```


```bash
sbatch herro-pre.sbatch
```

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/herro.sbatch`

```bash
#!/bin/bash

#SBATCH -c 12
#SBATCH --mem 300G
#SBATCH -p gpu
#SBATCH -n 1
#SBATCH --gpus-per-task=a100:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

srun --container-image=$HOME/scripts/herroNEW.sqush --no-container-mount-home --container-mounts=/msc:/msc /msc/home/jelber43/WGS_HG002_EZ1_25kb/herroNew.sh
```



`/msc/home/jelber43/WGS_HG002_EZ1_25kb/herroNew.sh`

```bash
#! /bin/bash
cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
eval "$(/msc/home/jelber43/.local/bin/micromamba shell hook --shell bash)"
micromamba activate herro
samtools fastq -@ 4 WGS_HG002_EZ1_25kb.pod5.bam > WGS_HG002_EZ1_25kb.pod5.bam.fastq
time /home/herro/target/release/herro inference --read-alns batches -t 12 -d 0 -m /home/herro/model_v0.1.pt -b 64 WGS_HG002_EZ1_25kb.pod5.bam.fastq WGS_HG002_EZ1_25kb.herroNew.fasta
```



This is a Snakemake file (snakefile) to run Brutal Rewrite and Peregrine_2021

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/herro-pg_asm1x.smk`

```yaml
#
shell.executable("/bin/bash")

# set input path here
IDS, = glob_wildcards("{id}.herroNew.fasta")

# set number of threads here
THREADS=101
THREADS2=range(1,THREADS+1)
THREADS3=expand(["{threads}"], threads=THREADS2)
THREADS4=[str(item).zfill(3) for item in THREADS3]
THREADS5=[str(item).zfill(2) for item in THREADS3]
MEMORY_SIZE2="-Xmx300g"

scattergather:
    split=128

rule all:
        input: expand(["{id}/pg_asm2/2_{threads2}.fa.gz"], threads2=THREADS5, allow_missing=True, id=IDS)

localrules: brutal_rewrite,all,cat,compress_fasta1,compress_overlap1,partition,pg_build_idx1,pg_build_sdb1,pg_ovlp1,pg_ovlp_ec1,reads1,shred

# run brutal rewrite
rule brutal_rewrite:
        input: "{id}.herroNew.fasta"
        output: "br/{id}.br.fasta"
        params: THREADS
        shell: """
        /msc/home/jelber43/git/br/target/release/br -t {params} -k 19 -i {input} -m graph -o {output}
        """

# add parition step
rule partition:
    input: "br/{id}.br.fasta"
    output: scatter.split("partition/{{id}}.herro_{scatteritem}.fasta.gz")
    params:
        split = "128",
        memory = MEMORY_SIZE2
    shell: '''
        partition.sh {params.memory} in={input} out=partition/{wildcards.id}.herro_%.fasta.gz ways={params.split}
        for i in `seq 1 128`
        do
          j=$((i-1))
          mv partition/{wildcards.id}.herro_${{j}}.fasta.gz partition/{wildcards.id}.herro_${{i}}-of-128.fasta.gz
        done
        '''

# add shredding step
rule shred:
    input: "partition/{id}.herro_{scatteritem}.fasta.gz"
    output: "shred/{id}.herro.shred_{scatteritem}.fasta.gz"
    params: MEMORY_SIZE2
    shell: '''
        shred.sh {params} in={input} out={output} median=1000000 variance=2500
        '''

rule cat:
    input: gather.split("shred/{{id}}.herro.shred_{scatteritem}.fasta.gz")
    output: "shred/{id}.herro.shred.fasta.gz"
    shell: """cat {input} > {output}"""

# now run peregrine-2021

# make a file listing the location of the raw HiFi reads
rule reads1:
        input: "shred/{id}.herro.shred.fasta.gz"
        output: "{id}/pg_asm1/reads1"
        shell: """ls `pwd`/{input} > {output}"""

# make a sequence database for the raw HiFi reads
rule pg_build_sdb1:
        input: "{id}/pg_asm1/reads1"
        output:
            database="{id}/pg_asm1/1.seqdb",
            index="{id}/pg_asm1/1.idx"
        params:
            prefix= "{id}/pg_asm1/1"
        shell: """/msc/home/jelber43/git/peregrine-2021/target/release/pg_build_sdb --input {input} --out_prefix {params.prefix}"""


# make a shimmer index for the raw reads
rule pg_build_idx1:
        input:
            sqdb="{id}/pg_asm1/1.seqdb",
            idx="{id}/pg_asm1/1.idx"
        output: expand(["{id}/pg_asm1/1-{threads2}-of-{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        params:
            k="56",
            r="6",
            w="80",
            prefix= "{id}/pg_asm1/1",
            chunks= THREADS
        shell: """/msc/home/jelber43/git/peregrine-2021/target/release/pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}"""

# overlap the raw reads
rule pg_ovlp1:
        input:
            sqdb="{id}/pg_asm1/1.seqdb",
            idx="{id}/pg_asm1/1.idx",
            dat=expand(["{id}/pg_asm1/1-{threads2}-of-{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        output: expand(["{id}/pg_asm1/overlap1.{threads2}"], threads2=THREADS5, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm1/1",
            prefix2= "{id}/pg_asm1/overlap1",
            chunks= THREADS
        shell: """/msc/home/jelber43/git/peregrine-2021/target/release/pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}"""

# error-correct the raw reads
rule pg_ovlp_ec1:
        input:
            sqdb="{id}/pg_asm1/1.seqdb",
            idx="{id}/pg_asm1/1.idx",
            ovlp=expand(["{id}/pg_asm1/overlap1.{threads2}"], threads2=THREADS5, allow_missing=True)
        output: expand(["{id}/pg_asm2/2_{threads2}.fa"], threads2=THREADS5, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm1/overlap1",
            prefix2= "{id}/pg_asm2/2",
            chunks= THREADS
        shell: """/msc/home/jelber43/git/peregrine-2021/target/release/pg_ovlp_ec {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}"""

# compress corrected reads
rule compress_fasta1:
        input: expand(["{id}/pg_asm2/2_{threads2}.fa"], threads2=THREADS5, allow_missing=True)
        output: expand(["{id}/pg_asm2/2_{threads2}.fa.gz"], threads2=THREADS5, allow_missing=True)
        params: THREADS
        shell: """pigz -p {params} {input}"""

# compress overlap1 files
rule compress_overlap1:
        input:
            reads=expand(["{id}/pg_asm2/2_{threads2}.fa.gz"], threads2=THREADS5, allow_missing=True),
            overlap=expand(["{id}/pg_asm1/overlap1.{threads2}"], threads2=THREADS5, allow_missing=True)
        output: expand(["{id}/pg_asm1/overlap1.{threads2}.gz"], threads2=THREADS5, allow_missing=True)
        params: THREADS
        shell: """pigz -p {params} {input.overlap}"""
```



This is an sbatch script to actually run the Snakefile and then DeChat

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/br-pg_asm1x-dechat.sbatch`

```bash
#!/bin/bash

#SBATCH -c 101
#SBATCH --mem 750G
#SBATCH -p cpu
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

date
cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
eval "$(/msc/home/jelber43/.local/bin/micromamba shell hook --shell bash)"
micromamba activate herro-pg_asm
time snakemake --local-cores 1 -j 1 -c 1 --snakefile herro-pg_asm1x.smk --printshellcmds all > herro-pg_asm1x.smk.log 2>&1
date
micromamba deactivate
time cat WGS_HG002_EZ1_25kb/pg_asm2/2_*.fa.gz > WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.fasta.gz
micromamba activate dechat
time dechat -r 1 -i WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.fasta.gz -o WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.dechat -t 248 > dechat.log 2>&1
date
```

```bash
cd /msc/home/jelber43/WGS_HG002_EZ1_25kb

sbatch br-pg_asm1x-dechat.sbatch
```





## 2024-06-05



Get MATERNAL and PATERNAL genomes from HG002

```bash
cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
srun -n 1 -p cpu -c 4 --mem 100g --pty bash -l
micromamba activate minimap2

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.fasta.gz

gunzip hg002v1.0.1.fasta.gz

seqtk seq -l0 hg002v1.0.1.fasta|paste - - |fgrep "MATERNAL" |tr '\t' '\n'|seqtk seq -l60 > hg002v1.0.1.maternal.fasta &
seqtk seq -l0 hg002v1.0.1.fasta|paste - - |fgrep "PATERNAL" |tr '\t' '\n'|seqtk seq -l60 > hg002v1.0.1.paternal.fasta &
```



Map SUP (raw) reads to MATERNAL AND PATERNAL haplotypes

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/sup-to-haps.sbatch`

```bash
#!/bin/bash

#SBATCH -c 248
#SBATCH -p cpu
#SBATCH --mem 300G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
eval "$(/msc/home/jelber43/.local/bin/micromamba shell hook --shell bash)"
micromamba activate minimap2.28

minimap2 --eqx --secondary=no -t 248 -Y -c -ax lr:hq hg002v1.0.1.maternal.fasta <(samtools fasta -@4 WGS_HG002_EZ1_25kb.pod5.bam) | samtools view -h -q60 -F 0x100 -@4 - |paftools.js sam2paf - |> WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.paf

minimap2 --eqx --secondary=no -t 248 -Y -c -ax lr:hq hg002v1.0.1.paternal.fasta <(samtools fasta -@4 WGS_HG002_EZ1_25kb.pod5.bam) | samtools view -h -q60 -F 0x100 -@4 - |paftools.js sam2paf - |> WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.paf
```




Map Herro reads to MATERNAL AND PATERNAL haplotypes

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/herroNew-to-haps.sbatch`

```bash
#!/bin/bash

#SBATCH -c 128
#SBATCH -p cpu
#SBATCH --mem 300G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
eval "$(/msc/home/jelber43/.local/bin/micromamba shell hook --shell bash)"
micromamba activate minimap2.28

minimap2 --eqx --secondary=no -t 128 -Y -c -ax lr:hq hg002v1.0.1.maternal.fasta WGS_HG002_EZ1_25kb.herroNew.fasta | samtools view -h -q60 -F 0x100 -@4 - |paftools.js sam2paf - > WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.paf

minimap2 --eqx --secondary=no -t 128 -Y -c -ax lr:hq hg002v1.0.1.paternal.fasta WGS_HG002_EZ1_25kb.herroNew.fasta | samtools view -h -q60 -F 0x100 -@4 - |paftools.js sam2paf - > WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.paf
```



Map Brutal Rewrite reads to MATERNAL AND PATERNAL haplotypes

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/brNew-to-haps.sbatch`

```bash
#!/bin/bash

#SBATCH -c 128
#SBATCH -p cpu
#SBATCH --mem 300G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
eval "$(/msc/home/jelber43/.local/bin/micromamba shell hook --shell bash)"
micromamba activate minimap2.28

minimap2 --eqx --secondary=no -t 128 -Y -c -ax lr:hq hg002v1.0.1.maternal.fasta br/WGS_HG002_EZ1_25kb.br.fasta | samtools view -h -q60 -F 0x100 -@4 - |paftools.js sam2paf - > WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.paf

minimap2 --eqx --secondary=no -t 128 -Y -c -ax lr:hq hg002v1.0.1.paternal.fasta br/WGS_HG002_EZ1_25kb.br.fasta | samtools view -h -q60 -F 0x100 -@4 - |paftools.js sam2paf - > WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.paf
```



Map Peregrine_2021 reads to MATERNAL AND PATERNAL haplotypes

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/pgNew-to-haps.sbatch`
```bash
#!/bin/bash

#SBATCH -c 248
#SBATCH -p cpu
#SBATCH --mem 300G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
eval "$(/msc/home/jelber43/.local/bin/micromamba shell hook --shell bash)"
micromamba activate minimap2.28

minimap2 --eqx --secondary=no -t 248 -Y -c -ax lr:hq hg002v1.0.1.maternal.fasta WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.fasta.gz | samtools view -h -q60 -F 0x100 -@4 - |paftools.js sam2paf - > WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.paf

minimap2 --eqx --secondary=no -t 248 -Y -c -ax lr:hq hg002v1.0.1.paternal.fasta WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.fasta.gz | samtools view -h -q60 -F 0x100 -@4 - |paftools.js sam2paf - > WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.paf
```



Map DeChat reads to MATERNAL AND PATERNAL haplotypes

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/dechatNew-to-haps.sbatch`
```bash
#!/bin/bash

#SBATCH -c 248
#SBATCH -p cpu
#SBATCH --mem 300G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
eval "$(/msc/home/jelber43/.local/bin/micromamba shell hook --shell bash)"
micromamba activate minimap2.28

minimap2 --eqx --secondary=no -t 248 -Y -c -ax lr:hq hg002v1.0.1.maternal.fasta recorrected.fa | samtools view -h -q60 -F 0x100 -@4 - |paftools.js sam2paf - > WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.paf

minimap2 --eqx --secondary=no -t 248 -Y -c -ax lr:hq hg002v1.0.1.paternal.fasta recorrected.fa | samtools view -h -q60 -F 0x100 -@4 - |paftools.js sam2paf - > WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.paf
```



This is the Julia Script for getting readmers and hapmers

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/phase-switching-0mers.jl`

```julia
using BioSequences
using FASTX
using Kmers
using CSV
using DataFrames
using Base.Threads
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--maternal_paf_file", "-m"
            help = "path to maternal paf file"
            required = true
        "--paternal_paf_file", "-p"
            help = "path to paternal paf file"
            required = true
        "--nanopore_fasta", "-n"
            help = "path to nanopore fasta file (SUP.fasta)"
            required = true
        "--maternal_fasta", "-i"
            help = "path to maternal fasta file (hg002v1.0.1.maternal.fasta)"
            required = true
        "--paternal_fasta", "-j"
            help = "path to paternal fasta file (hg002v1.0.1.paternal.fasta)"
            required = true
        "--kmer_length", "-k"
            arg_type = Int64
            help = "for example 15"
            required = true
        "--output_file", "-o"
            help = "Output file name (WGS_HG002_EZ1_10kb_SUP.kmer-matches.txt)"
            required = true
    end

    return parse_args(s)
end

function phase_switches(maternal_paf_file::String, 
                        paternal_paf_file::String,
                        nanopore_fasta::String,
                        maternal_fasta::String,
                        paternal_fasta::String,
                        kmer_length::Int64)

    println()
    println("Reading maternal paf file...")
    maternal_paf = CSV.read(maternal_paf_file, 
                            DataFrame; 
                            delim='\t',
                            header=false,
                            types=[String, UInt32, UInt32, String, String, UInt32, UInt32])

    rename!(maternal_paf, [:query_id, :query_start, :query_stop, :strand, :ref_id, :ref_start, :ref_stop])
    println("Done reading maternal paf file.")

    println()
    println("Reading paternal paf file...")
    paternal_paf = CSV.read(paternal_paf_file, 
                            DataFrame; 
                            delim='\t',
                            header=false,
                            types=[String, UInt32, UInt32, String, String, UInt32, UInt32])

    rename!(paternal_paf, [:query_id, :query_start, :query_stop, :strand, :ref_id, :ref_start, :ref_stop])
    println("Done reading paternal paf file.")

    # Function to extract the chromosome part using regular expressions
    function extract_chr_id(id::String)
        return match.(r"chr\w+_", id)
    end

    println()
    println("Extracting maternal chr ids...")
    # Apply the function to create a new column for matching
    maternal_paf[!, :chr_id] = extract_chr_id.(maternal_paf.ref_id)
    maternal_paf[!, :chr_id] = [m.match for m in maternal_paf.chr_id]
    println("Done extracting maternal chr ids.")

    println()
    println("Extracting paternal chr ids...")
    paternal_paf[!, :chr_id] = extract_chr_id.(paternal_paf.ref_id)
    paternal_paf[!, :chr_id] = [m.match for m in paternal_paf.chr_id]
    println("Done extracting paternal chr ids.")

    println()
    println("Performing leftjoin on maternal and paternal paf files...")
    pafs = leftjoin(maternal_paf, paternal_paf, on = [:query_id, :chr_id, :query_start, :query_stop], makeunique=true)

    dropmissing!(pafs)
    println("Done performing leftjoin on maternal and paternal paf files.")


    paternal_paf=nothing
    maternal_paf=nothing

    println()
    println("Reading read FASTA Index for preallocating dictionary size.")
    nanopore_fasta_index = string(nanopore_fasta, ".fai")
    total_records = length(FASTX.FASTA.Index(nanopore_fasta_index).lengths)
    println("There are $total_records read FASTA records")
    nanopore_dict = Dict{String, String}()
    reader = open(FASTA.Reader, nanopore_fasta)
    record = FASTA.Record()
    counter = Atomic{Int}(0)
    reader_lock = ReentrantLock()
    nanopore_lock = ReentrantLock()

    println()
    println("Reading read FASTA records for adding to dictionary.")
    Threads.@threads for i in 1:total_records
        local_record = FASTA.Record()
        lock(reader_lock)
        read!(reader, local_record)
        unlock(reader_lock)

        lock(nanopore_lock)
        nanopore_dict[identifier(local_record)] = sequence(local_record)
        unlock(nanopore_lock)

        atomic_add!(counter, 1)
        print("\rRead $(counter[]) FASTA records.")
        flush(stdout)
    end
    println()
    println("Done reading read FASTA records for adding to dictionary.")
    close(reader)

    println()
    println("Reading maternal FASTA Index for preallocating dictionary size.")
    maternal_fasta_index = string(maternal_fasta, ".fai")
    total_records_maternal = length(FASTX.FASTA.Index(maternal_fasta_index).lengths)
    println("There are $total_records_maternal FASTA records")

    maternal_dict = Dict{String, String}()
    reader = open(FASTA.Reader, maternal_fasta)
    record = FASTA.Record()
    counter = Atomic{Int}(0)
    reader_lock = ReentrantLock()
    maternal_lock = ReentrantLock()

    println()
    println("Reading maternal FASTA records for adding to dictionary.")
    Threads.@threads for i in 1:total_records_maternal
        local_record = FASTA.Record()
        lock(reader_lock)
        read!(reader, local_record)
        unlock(reader_lock)

        lock(maternal_lock)
        maternal_dict[identifier(local_record)] = sequence(local_record)
        unlock(maternal_lock)

        atomic_add!(counter, 1)
        print("\rRead $(counter[]) FASTA records.")
        flush(stdout)
    end
    println()
    println("Done reading FASTA records for adding to dictionary.")
    close(reader)

    println()
    println("Reading paternal FASTA Index for preallocating dictionary size.")
    paternal_fasta_index = string(paternal_fasta, ".fai")
    total_records_paternal = length(FASTX.FASTA.Index(paternal_fasta_index).lengths)
    println("There are $total_records_paternal FASTA records")

    paternal_dict = Dict{String, String}()
    reader = open(FASTA.Reader, paternal_fasta)
    record = FASTA.Record()
    counter = Atomic{Int}(0)
    reader_lock = ReentrantLock()
    paternal_lock = ReentrantLock()

    println()
    println("Reading paternal FASTA records for adding to dictionary.")
    Threads.@threads for i in 1:total_records_paternal
        local_record = FASTA.Record()
        lock(reader_lock)
        read!(reader, local_record)
        unlock(reader_lock)

        lock(paternal_lock)
        paternal_dict[identifier(local_record)] = sequence(local_record)
        unlock(paternal_lock)


        atomic_add!(counter, 1)
        print("\rRead $(counter[]) FASTA records.")
        flush(stdout)
    end
    println()
    println("Done reading paternal FASTA records for adding to dictionary.")
    close(reader)

    switches = DataFrame(read_id=Vector{String}(undef, nrow(pafs)),
                         maternal_hapmers_matching=Vector{UInt64}(undef, nrow(pafs)),
                         paternal_hapmers_matching=Vector{UInt64}(undef, nrow(pafs)),
                         ZeroMers=Vector{UInt64}(undef, nrow(pafs)))

    counter = Atomic{Int}(0)

    Threads.@threads for i in 1:nrow(pafs)
        switches_lock = ReentrantLock()

        # nanopore SUP sequences
        nanopore_start = Int(pafs.query_start[i]) + 1 
        nanopore_stop = Int(pafs.query_stop[i])
        nanopore_seq = LongSequence{DNAAlphabet{2}}(nanopore_dict[pafs.query_id[i]][nanopore_start:nanopore_stop])
        if pafs.strand[i] == "-"
            reverse_complement!(nanopore_seq)
        end

        # maternal sequences
        maternal_start = Int(pafs.ref_start[i]) + 1
        maternal_stop = Int(pafs.ref_stop[i])
        maternal_seq = LongSequence{DNAAlphabet{2}}(maternal_dict[pafs.ref_id[i]][maternal_start:maternal_stop])
        
        # maternal sequences
        paternal_start = Int(pafs.ref_start_1[i]) + 1
        paternal_stop = Int(pafs.ref_stop_1[i])
        paternal_seq = LongSequence{DNAAlphabet{2}}(paternal_dict[pafs.ref_id_1[i]][paternal_start:paternal_stop])

        # make an empty vector of kmer_length with dna alphabet 2
        nanopore_kmers = Vector{Kmer{DNAAlphabet{2}, kmer_length, 1}}()
    
        # add the k-mers to the vector
        for i in EveryKmer(nanopore_seq, Val(kmer_length))
            push!(nanopore_kmers, i[2])
        end
    
        # make an empty vector of kmer_length with dna alphabet 2
        maternal_kmers = Vector{Kmer{DNAAlphabet{2}, kmer_length, 1}}()
    
        # add the k-mers to the vector
        for i in EveryKmer(maternal_seq, Val(kmer_length))
            push!(maternal_kmers, i[2])
        end
        
        # make an empty vector of kmer_length with dna alphabet 2
        paternal_kmers = Vector{Kmer{DNAAlphabet{2}, kmer_length, 1}}()
    
        # add the k-mers to the vector
        for i in EveryKmer(paternal_seq, Val(kmer_length))
            push!(paternal_kmers, i[2])
        end
   
        # Find unique elements in each vector
        maternal_hapmers = setdiff(maternal_kmers, paternal_kmers)
        paternal_hapmers = setdiff(paternal_kmers, maternal_kmers)
    
        nanopore_vs_paternal_hapmers = intersect(paternal_hapmers, nanopore_kmers)
        nanopore_vs_maternal_hapmers = intersect(maternal_hapmers, nanopore_kmers)

        lock(switches_lock)
        switches[i, :] = [pafs.query_id[i], length(nanopore_vs_maternal_hapmers), length(nanopore_vs_paternal_hapmers), length(maternal_hapmers)]
        unlock(switches_lock)
        atomic_add!(counter, 1)
        print("\rProcessed $(counter[]) records for hapmer analysis.")
        flush(stdout)
    end
    return switches
end


function main()
    parsed_args = parse_commandline()
    maternal_paf_file = parsed_args["maternal_paf_file"]
    paternal_paf_file = parsed_args["paternal_paf_file"]
    nanopore_fasta = parsed_args["nanopore_fasta"]
    maternal_fasta = parsed_args["maternal_fasta"]
    paternal_fasta = parsed_args["paternal_fasta"]
    kmer_length = parsed_args["kmer_length"]
    output_file = parsed_args["output_file"]

    switches = phase_switches(maternal_paf_file,
                              paternal_paf_file,
                              nanopore_fasta,
                              maternal_fasta,
                              paternal_fasta,
                              kmer_length)

    CSV.write(output_file, switches; delim='\t', newline='\n')
end


main()
```




Sort PAF files by alignment length and retain only one alignment per read, remove split reads, and also keep only the first 1,000,000 bases of a read

```bash
srun -c 96 --mem=1000G -p cpu --pty bash -l
sort --parallel=96 -S 1000G -uk1,1 -rnk11,11 WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.paf | sort --parallel=96 -S 1000G -uk1,1 | cut -f 1,3-6,8-9 | fgrep -v ":" > WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf

sort --parallel=96 -S 1000G -uk1,1 -rnk11,11 WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.paf | sort --parallel=96 -S 1000G -uk1,1 | cut -f 1,3-6,8-9 | fgrep -v ":" > WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf

sort --parallel=96 -S 1000G -uk1,1 -rnk11,11 WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.paf | sort --parallel=96 -S 1000G -uk1,1| cut -f 1,3-6,8-9 | fgrep -v ":" > WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf

sort --parallel=96 -S 1000G -uk1,1 -rnk11,11 WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.paf | sort --parallel=96 -S 1000G -uk1,1 | cut -f 1,3-6,8-9 | fgrep -v ":" > WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf

sort --parallel=96 -S 1000G -uk1,1 -rnk11,11 WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.paf | sort --parallel=96 -S 1000G -uk1,1 | cut -f 1,3-6,8-9 | fgrep -v ":" > WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf

sort --parallel=96 -S 1000G -uk1,1 -rnk11,11 WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.paf | sort --parallel=96 -S 1000G -uk1,1 | cut -f 1,3-6,8-9 | fgrep -v ":" > WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf

sort --parallel=96 -S 1000G -uk1,1 -rnk11,11 WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.paf | sort --parallel=96 -S 1000G -uk1,1 | cut -f 1,3-6,8-9 | fgrep -v ":" |fgrep "_0-" -> WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf

sort --parallel=96 -S 1000G -uk1,1 -rnk11,11 WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.paf | sort --parallel=96 -S 1000G -uk1,1 | cut -f 1,3-6,8-9 | fgrep -v ":" | fgrep "_0-" > WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf

sort --parallel=96 -S 1000G -uk1,1 -rnk11,11 WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.paf | sort --parallel=96 -S 1000G -uk1,1 | cut -f 1,3-6,8-9 | fgrep -v ":" | fgrep "_0-" > WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf

sort --parallel=96 -S 1000G -uk1,1 -rnk11,11 WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.paf | sort --parallel=96 -S 1000G -uk1,1 | cut -f 1,3-6,8-9 | fgrep -v ":" | fgrep "_0-" > WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf
```



Total number of reads in the different datasets

```bash
wc -l *New_against_hg002v1.0.1.?aternal.fasta.mapQ60.primary.cols1345689.unique.paf WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.?aternal.fasta.mapQ60.primary.cols1345689.unique.paf
   1733152 WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf
   1705385 WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf
   1724727 WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf
   1697341 WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf
   1733156 WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf
   1705413 WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf
   1724488 WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf
   1697130 WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf
   4029514 WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf
   3962042 WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf
```




Remove the trailing sequence length from Peregrine_2021 and DeChat read alignments

```bash
srun -c 96 -p cpu --mem=1900G --pty bash -l
cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
perl -pe "s/_0-\d+//g" WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf > WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf2

perl -pe "s/_0-\d+//g" WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf > WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf2

perl -pe "s/_0-\d+//g" WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf > WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf2

perl -pe "s/_0-\d+//g" WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf > WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf2

julia
```



This code is to get the shared read ids among all: raw (SUP), Herro, Brutal Rewrite, Peregrine_2021, and DeChat datasets

```julia
using BioSequences
using FASTX
using Kmers
using CSV
using DataFrames
using Base.Threads

function get_pafs(maternal_paf_file::String, 
                  paternal_paf_file::String)

    println()
    println("Reading maternal paf file...")
    maternal_paf = CSV.read(maternal_paf_file, 
                            DataFrame; 
                            delim='\t',
                            header=false,
                            types=[String, UInt32, UInt32, String, String, UInt32, UInt32])

    rename!(maternal_paf, [:query_id, :query_start, :query_stop, :strand, :ref_id, :ref_start, :ref_stop])
    println("Done reading maternal paf file.")

    println()
    println("Reading paternal paf file...")
    paternal_paf = CSV.read(paternal_paf_file, 
                            DataFrame; 
                            delim='\t',
                            header=false,
                            types=[String, UInt32, UInt32, String, String, UInt32, UInt32])

    rename!(paternal_paf, [:query_id, :query_start, :query_stop, :strand, :ref_id, :ref_start, :ref_stop])
    println("Done reading paternal paf file.")

    # Function to extract the chromosome part using regular expressions
    function extract_chr_id(id::String)
        return match.(r"chr\w+_", id)
    end

    println()
    println("Extracting maternal chr ids...")
    # Apply the function to create a new column for matching
    maternal_paf[!, :chr_id] = extract_chr_id.(maternal_paf.ref_id)
    maternal_paf[!, :chr_id] = [m.match for m in maternal_paf.chr_id]
    println("Done extracting maternal chr ids.")

    println()
    println("Extracting paternal chr ids...")
    paternal_paf[!, :chr_id] = extract_chr_id.(paternal_paf.ref_id)
    paternal_paf[!, :chr_id] = [m.match for m in paternal_paf.chr_id]
    println("Done extracting paternal chr ids.")

    println()
    println("Performing leftjoin on maternal and paternal paf files...")
    pafs = leftjoin(maternal_paf, paternal_paf, on = [:query_id, :chr_id, :query_start, :query_stop], makeunique=true)

    dropmissing!(pafs)
    println("Done performing leftjoin on maternal and paternal paf files.")


    paternal_paf=nothing
    maternal_paf=nothing
    return pafs
end

sup_paf = get_pafs("WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf",
                   "WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf")

herro_paf = get_pafs("WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf",
                     "WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf")

br_paf = get_pafs("WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf",
                  "WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf")

pg_paf = get_pafs("WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf2",
                  "WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf2")

dechat_paf = get_pafs("WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf2",
                      "WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf2")


sup_herro = intersect(sup_paf[:,1], herro_paf[:,1])
sup_herro_br = intersect(sup_herro, br_paf[:,1])
sup_herro_br_pg = intersect(sup_herro, pg_paf[:,1])
sup_herro_br_pg_dechat = intersect(sup_herro, dechat_paf[:,1])


sup_herro_br_pg_dechat_DataFrame = DataFrame(query_id=sup_herro_br_pg_dechat)

CSV.write("sup_herro_br_pg_dechat_DataFrameNew.txt", sup_herro_br_pg_dechat_DataFrame; delim='\t', newline='\n')
```




Get the shredded suffix to filter the reads with filterbyname.sh from BBTools/BBMap

```bash
micromamba activate minimap2
grep -f sup_herro_br_pg_dechat_DataFrameNew.txt \
 <(cut -f 1 WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf) > pg_dechat_reads.to.filterNew.txt

filterbyname.sh ow=t include=t in=recorrected.fa names=pg_dechat_reads.to.filterNew.txt out=WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.dechat.same_reads.fasta > WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.dechat.same_reads.fasta.log 2>&1 &
filterbyname.sh ow=t include=t in=WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.fasta.gz names=pg_dechat_reads.to.filterNew.txt out=WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.same_reads.fasta > WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.same_reads.fasta.log 2>&1 &

filterbyname.sh ow=t include=t in=WGS_HG002_EZ1_25kb_SUP.fasta names=sup_herro_br_pg_dechat_DataFrameNew.txt out=WGS_HG002_EZ1_25kb_SUPNew.same_reads.fasta > WGS_HG002_EZ1_25kb_SUPNew.same_reads.fasta.log 2>&1 &

filterbyname.sh ow=t include=t in=WGS_HG002_EZ1_25kb.herroNew.fasta names=sup_herro_br_pg_dechat_DataFrameNew.txt out=WGS_HG002_EZ1_25kb.herroNew.same_reads.fasta > WGS_HG002_EZ1_25kb.herroNew.same_reads.fasta.log 2>&1 &

filterbyname.sh ow=t include=t in=br/WGS_HG002_EZ1_25kb.br.fasta names=sup_herro_br_pg_dechat_DataFrameNew.txt out=WGS_HG002_EZ1_25kb.brNew.same_reads.fasta > WGS_HG002_EZ1_25kb.brNew.same_reads.fasta.log 2>&1 &


# remove trailing read lengths
perl -pi -e "s/_0-\d+.+//g" WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.same_reads.fasta

perl -pi -e "s/_0-\d+.+//g" WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.dechat.same_reads.fasta

# make FASTA index
samtools faidx WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.same_reads.fasta &
samtools faidx WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.dechat.same_reads.fasta &
samtools faidx WGS_HG002_EZ1_25kb_SUPNew.same_reads.fasta &
samtools faidx WGS_HG002_EZ1_25kb.brNew.same_reads.fasta &
samtools faidx WGS_HG002_EZ1_25kb.herroNew.same_reads.fasta &

# Filter PAF files for same_read_ids 
fgrep -f sup_herro_br_pg_dechat_DataFrameNew.txt WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf > WGS_HG002_EZ1_25kb_SUPNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf4 &

fgrep -f sup_herro_br_pg_dechat_DataFrameNew.txt WGS_HG002_EZ1_25kb_SUP_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf > WGS_HG002_EZ1_25kb_SUPNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf4 &

fgrep -f sup_herro_br_pg_dechat_DataFrameNew.txt WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf > WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf4 &

fgrep -f sup_herro_br_pg_dechat_DataFrameNew.txt WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf > WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf4 &


fgrep -f sup_herro_br_pg_dechat_DataFrameNew.txt WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf > WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf4 &

fgrep -f sup_herro_br_pg_dechat_DataFrameNew.txt WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf > WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf4 &

fgrep -f sup_herro_br_pg_dechat_DataFrameNew.txt WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf2 > WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf4 &

fgrep -f sup_herro_br_pg_dechat_DataFrameNew.txt WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf2 > WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf4 &

fgrep -f sup_herro_br_pg_dechat_DataFrameNew.txt WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf2 > WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf4 &

fgrep -f sup_herro_br_pg_dechat_DataFrameNew.txt WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf2 > WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf4 &
```



Actually calculate phase-switching for SUP read alignments

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/SUP-phase-switching3-0mers.sbatch`

```bash
#!/bin/bash

#SBATCH -c 48
#SBATCH -p cpu
#SBATCH --mem 300G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb

time julia --threads 48 \
phase-switching-0mers.jl \
-m WGS_HG002_EZ1_25kb_SUPNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf4 \
-p WGS_HG002_EZ1_25kb_SUPNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf4 \
-n WGS_HG002_EZ1_25kb_SUPNew.same_reads.fasta \
-i hg002v1.0.1.maternal.fasta \
-j hg002v1.0.1.paternal.fasta \
-k 15 \
-o WGS_HG002_EZ1_25kb_SUPNew.kmer_matches.same_reads_0mers.txt
```



Actually calculate phase-switching for Herro read alignments

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/herro-phase-switching3-0mers.sbatch`
```bash
#!/bin/bash

#SBATCH -c 48
#SBATCH -p cpu
#SBATCH --mem 300G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb

time julia --threads 48 \
phase-switching-0mers.jl \
-m WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf4 \
-p WGS_HG002_EZ1_25kb_herroNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf4 \
-n WGS_HG002_EZ1_25kb.herroNew.same_reads.fasta \
-i hg002v1.0.1.maternal.fasta \
-j hg002v1.0.1.paternal.fasta \
-k 15 \
-o WGS_HG002_EZ1_25kb_herroNew.kmer_matches.same_reads_0mers.txt
```



Actually calculate phase-switching for DeChat read alignments

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/dechat-phase-switching3-0mers.sbatch`
```bash
#!/bin/bash

#SBATCH -c 48
#SBATCH -p cpu
#SBATCH --mem 300G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb

time julia --threads 48 \
phase-switching-0mers.jl \
-m WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf4 \
-p WGS_HG002_EZ1_25kb_dechatNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf4 \
-n WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.dechat.same_reads.fasta \
-i hg002v1.0.1.maternal.fasta \
-j hg002v1.0.1.paternal.fasta \
-k 15 \
-o WGS_HG002_EZ1_25kb_dechatNew.kmer_matches.same_reads_0mers.txt
```



Actually calculate phase-switching for Brutal Rewrite read alignments

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/br-phase-switching3-0mers.sbatch`
```bash
#!/bin/bash

#SBATCH -c 48
#SBATCH -p cpu
#SBATCH --mem 300G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb

time julia --threads 48 \
phase-switching-0mers.jl \
-m WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf4 \
-p WGS_HG002_EZ1_25kb_brNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf4 \
-n WGS_HG002_EZ1_25kb.brNew.same_reads.fasta \
-i hg002v1.0.1.maternal.fasta \
-j hg002v1.0.1.paternal.fasta \
-k 15 \
-o WGS_HG002_EZ1_25kb_brNew.kmer_matches.same_reads_0mers.txt
```



Actually calculate phase-switching for Peregrine_2021 read alignments

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/pg-phase-switching3-0mers.sbatch`
```bash
#!/bin/bash

#SBATCH -c 48
#SBATCH -p cpu
#SBATCH --mem 300G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb

time julia --threads 48 \
phase-switching-0mers.jl \
-m WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.maternal.fasta.mapQ60.primary.cols1345689.unique.paf4 \
-p WGS_HG002_EZ1_25kb_pgNew_against_hg002v1.0.1.paternal.fasta.mapQ60.primary.cols1345689.unique.paf4 \
-n WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.same_reads.fasta \
-i hg002v1.0.1.maternal.fasta \
-j hg002v1.0.1.paternal.fasta \
-k 15 \
-o WGS_HG002_EZ1_25kb_pgNew.kmer_matches.same_reads_0mers.txt
```



get matching same_reads_ids again

```bash
comm -12 <(cut -f 1 WGS_HG002_EZ1_25kb_dechatNew.kmer_matches.same_reads_0mers.txt|sort)  \
<(comm -12 <(comm -12 <(cut -f 1 WGS_HG002_EZ1_25kb_SUPNew.kmer_matches.same_reads_0mers.txt|sort) <(cut -f 1 WGS_HG002_EZ1_25kb_herroNew.kmer_matches.same_reads_0mers.txt|sort)) \
         <(comm -12 <(cut -f 1 WGS_HG002_EZ1_25kb_brNew.kmer_matches.same_reads_0mers.txt|sort) <(cut -f 1 WGS_HG002_EZ1_25kb_pgNew.kmer_matches.same_reads_0mers.txt|sort))) > Newsame_reads_ids2_0mers.txt

fgrep -f Newsame_reads_ids2_0mers.txt WGS_HG002_EZ1_25kb_SUPNew.kmer_matches.same_reads_0mers.txt > WGS_HG002_EZ1_25kb_SUPNew.kmer_matches.same_reads2_0mers.txt
fgrep -f Newsame_reads_ids2_0mers.txt WGS_HG002_EZ1_25kb_herroNew.kmer_matches.same_reads_0mers.txt > WGS_HG002_EZ1_25kb_herroNew.kmer_matches.same_reads2_0mers.txt
fgrep -f Newsame_reads_ids2_0mers.txt WGS_HG002_EZ1_25kb_brNew.kmer_matches.same_reads_0mers.txt > WGS_HG002_EZ1_25kb_brNew.kmer_matches.same_reads2_0mers.txt
fgrep -f Newsame_reads_ids2_0mers.txt WGS_HG002_EZ1_25kb_pgNew.kmer_matches.same_reads_0mers.txt > WGS_HG002_EZ1_25kb_pgNew.kmer_matches.same_reads2_0mers.txt
fgrep -f Newsame_reads_ids2_0mers.txt WGS_HG002_EZ1_25kb_dechatNew.kmer_matches.same_reads_0mers.txt > WGS_HG002_EZ1_25kb_dechatNew.kmer_matches.same_reads2_0mers.txt

On the MUW Cluster

```bash
cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
julia
```



Use Julia to plot histograms

```julia
using CSV
using DataFrames
using Plots

test = CSV.read("WGS_HG002_EZ1_25kb_SUPNew.kmer_matches.same_reads2_0mers.txt", DataFrame; delim="\t",header=true, types=[String,Int64,Int64,Int64,Int64])
test2 = test
test2.maternal_hapmers_percent = test2.maternal_hapmers_matching ./ (test2.paternal_hapmers_matching .+ test2.maternal_hapmers_matching) .* 100
test2.paternal_hapmers_percent = test2.paternal_hapmers_matching ./ (test2.paternal_hapmers_matching .+ test2.maternal_hapmers_matching) .* 100
test2 = filter(row -> all(x -> !(x isa Number && isnan(x)), row), test2)
count = nrow(filter(row -> row.maternal_hapmers == 0 && row.paternal_hapmers == 0 && row.maternal_hapmers_matching == 0 && row.paternal_hapmers_matching == 0, test))

nrow(test)
nrow(test2)
nrow(test)-nrow(test2)
count
nrow(test)-nrow(test2)-count

# test
# 1145743 reads with mapping quality 60

# test2
# 1009418 reads (88.10160742854201  %) with hapmers matching for at least one haplotype (are these reads with heterozygous variants)
#  136325 reads (11.898392571457997 %) with no-matching hapmers (are these reads in homozygous regions?) for either haplotype
#  128842 reads (11.245279264197992 %) with no hapmers, reads in homozygous regions
#    7483 reads (0.6531133072600051 %) with presumably readmers as errors

test3 = CSV.read("WGS_HG002_EZ1_25kb_herroNew.kmer_matches.same_reads2_0mers.txt", DataFrame; delim="\t",header=true, types=[String,Int64,Int64,Int64,Int64])
test4 = test3
test4.maternal_hapmers_percent = test4.maternal_hapmers_matching ./ (test4.paternal_hapmers_matching .+ test4.maternal_hapmers_matching) .* 100
test4.paternal_hapmers_percent = test4.paternal_hapmers_matching ./ (test4.paternal_hapmers_matching .+ test4.maternal_hapmers_matching) .* 100
test4 = filter(row -> all(x -> !(x isa Number && isnan(x)), row), test4)
count2 = nrow(filter(row -> row.maternal_hapmers == 0 && row.paternal_hapmers == 0 && row.maternal_hapmers_matching == 0 && row.paternal_hapmers_matching == 0, test3))

nrow(test3)
nrow(test4)
nrow(test3)-nrow(test4)
count2
nrow(test3)-nrow(test4)-count2


# test3
# 1145743 reads with mapping quality 60

# test4
# 1009411 reads (88.10099647128544   %)  with hapmers matching for at least one haplotype (are these reads with heterozygous variants)
#  136332 reads (11.899003528714555  %) with no-matching hapmers (are these reads in homozygous regions?) for either haplotype
#  134100 reads (11.70419544348078   %) with no hapmers, reads in homozygous regions
#    2232 reads (0.19480808523377406 %) with presumably readmers as errors


test5 = CSV.read("WGS_HG002_EZ1_25kb_brNew.kmer_matches.same_reads2_0mers.txt", DataFrame; delim="\t",header=true, types=[String,Int64,Int64,Int64,Int64])
test6 = test5
test6.maternal_hapmers_percent = test6.maternal_hapmers_matching ./ (test6.paternal_hapmers_matching .+ test6.maternal_hapmers_matching) .* 100
test6.paternal_hapmers_percent = test6.paternal_hapmers_matching ./ (test6.paternal_hapmers_matching .+ test6.maternal_hapmers_matching) .* 100
test6 = filter(row -> all(x -> !(x isa Number && isnan(x)), row), test6)
count3 = nrow(filter(row -> row.maternal_hapmers == 0 && row.paternal_hapmers == 0 && row.maternal_hapmers_matching == 0 && row.paternal_hapmers_matching == 0, test5))

nrow(test5)
nrow(test6)
nrow(test5)-nrow(test6)
count3
nrow(test5)-nrow(test6)-count3


# test5
# 1145743 reads with mapping quality 60

# test6
# 1009411 reads (88.10099647128544   %)  with hapmers matching for at least one haplotype (are these reads with heterozygous variants)
#  136332 reads (11.899003528714555  %) with no-matching hapmers (are these reads in homozygous regions?) for either haplotype
#  134100 reads (11.70419544348078   %) with no hapmers, reads in homozygous regions
#    2232 reads (0.19480808523377406 %) with presumably readmers as errors


test7 = CSV.read("WGS_HG002_EZ1_25kb_pgNew.kmer_matches.same_reads2_0mers.txt", DataFrame; delim="\t",header=true, types=[String,Int64,Int64,Int64,Int64])
test8 = test7
test8.maternal_hapmers_percent = test8.maternal_hapmers_matching ./ (test8.paternal_hapmers_matching .+ test8.maternal_hapmers_matching) .* 100
test8.paternal_hapmers_percent = test8.paternal_hapmers_matching ./ (test8.paternal_hapmers_matching .+ test8.maternal_hapmers_matching) .* 100
test8 = filter(row -> all(x -> !(x isa Number && isnan(x)), row), test8)
count4 = nrow(filter(row -> row.maternal_hapmers == 0 && row.paternal_hapmers == 0 && row.maternal_hapmers_matching == 0 && row.paternal_hapmers_matching == 0, test7))



nrow(test7)
nrow(test8)
nrow(test7)-nrow(test8)
count4
nrow(test7)-nrow(test8)-count4


# test7
# 1145743 reads with mapping quality 60

# test8
# 1007110 reads (87.90016609309418   %)  with hapmers matching for at least one haplotype (are these reads with heterozygous variants)
#  138633 reads (12.099833906905824  %) with no-matching hapmers (are these reads in homozygous regions?) for either haplotype
#  136549 reads (11.917943203667838  %) with no hapmers, reads in homozygous regions
#    2084 reads (0.18189070323798618 %) with presumably readmers as errors


test9 = CSV.read("WGS_HG002_EZ1_25kb_dechatNew.kmer_matches.same_reads2_0mers.txt", DataFrame; delim="\t",header=true, types=[String,Int64,Int64,Int64,Int64])
test10 = test9
test10.maternal_hapmers_percent = test10.maternal_hapmers_matching ./ (test10.paternal_hapmers_matching .+ test10.maternal_hapmers_matching) .* 100
test10.paternal_hapmers_percent = test10.paternal_hapmers_matching ./ (test10.paternal_hapmers_matching .+ test10.maternal_hapmers_matching) .* 100
test10 = filter(row -> all(x -> !(x isa Number && isnan(x)), row), test10)
count5 = nrow(filter(row -> row.maternal_hapmers == 0 && row.paternal_hapmers == 0 && row.maternal_hapmers_matching == 0 && row.paternal_hapmers_matching == 0, test9))


nrow(test9)
nrow(test10)
nrow(test9)-nrow(test10)
count5
nrow(test9)-nrow(test10)-count5


# test7
# 1145743 reads with mapping quality 60

# test8
# 1007150 reads (87.90365727741736   %)  with hapmers matching for at least one haplotype (are these reads with heterozygous variants)
#  138593 reads (12.096342722582638  %) with no-matching hapmers (are these reads in homozygous regions?) for either haplotype
#  136549 reads (11.917943203667838  %) with no hapmers, reads in homozygous regions
#    2044 reads (0.17839951891480026 %) with presumably readmers as errors


hist1d = histogram(abs.(test2.maternal_hapmers_percent - test2.paternal_hapmers_percent), bins=100, title="RAW", xlabel="", ylabel="Num. Reads", legend=false, ylim=(0, 6000))
hist2d = histogram(abs.(test4.maternal_hapmers_percent - test4.paternal_hapmers_percent), bins=100, title="Herro", xlabel="", ylabel="Num. Reads", legend=false, ylim=(0, 6000))
hist3d = histogram(abs.(test6.maternal_hapmers_percent - test6.paternal_hapmers_percent), bins=100, title="Brutal Rewrite", xlabel="", ylabel="Num. Reads", legend=false, ylim=(0, 6000))
hist4d = histogram(abs.(test8.maternal_hapmers_percent - test8.paternal_hapmers_percent), bins=100, title="Peregrine 2021", xlabel="", ylabel="Num. Reads", legend=false, ylim=(0, 6000))
hist5d = histogram(abs.(test10.maternal_hapmers_percent - test10.paternal_hapmers_percent), bins=100, title="DeChat", xlabel="                                     Absolute value of difference in percentage matching hapmers", ylabel="Num. Reads", legend=false, ylim=(0, 6000))

combined_plot3 = plot(hist1d, hist2d, hist3d, hist4d, hist5d, layout=(3, 2))

savefig(combined_plot3, "hapmers-histogram-Newsame_reads2.svg")

hist1a = histogram(abs.(test2.maternal_hapmers_percent - test2.paternal_hapmers_percent), bins=100, title="RAW", xlabel="", ylabel="Num. Reads", legend=false)
hist2a = histogram(abs.(test4.maternal_hapmers_percent - test4.paternal_hapmers_percent), bins=100, title="Herro", xlabel="", ylabel="Num. Reads", legend=false)
hist3a = histogram(abs.(test6.maternal_hapmers_percent - test6.paternal_hapmers_percent), bins=100, title="Brutal Rewrite", xlabel="", ylabel="Num. Reads", legend=false)
hist4a = histogram(abs.(test8.maternal_hapmers_percent - test8.paternal_hapmers_percent), bins=100, title="Peregrine 2021", xlabel="", ylabel="Num. Reads", legend=false)
hist5a = histogram(abs.(test10.maternal_hapmers_percent - test10.paternal_hapmers_percent), bins=100, title="DeChat", xlabel="                                     Absolute value of difference in percentage matching hapmers", ylabel="Num. Reads", legend=false)

combined_plot4 = plot(hist1a, hist2a, hist3a, hist4a, hist5a, layout=(3, 2))

savefig(combined_plot4, "hapmers-histogram-Newsame_reads-full2.svg")
```



Map reads that were analyzed for hapmer analysis using the lr:hqae preset to the entire hg002v1.0.1.fasta genome

`/msc/home/jelber43/WGS_HG002_EZ1_25kb/lrhqae.sbatch`
```bash
#!/bin/bash

#SBATCH -c 200
#SBATCH -p cpu
#SBATCH --mem 300G
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jean.elbers@meduniwien.ac.at

THREADS=200

cd /msc/home/jelber43/WGS_HG002_EZ1_25kb
eval "$(/msc/home/jelber43/.local/bin/micromamba shell hook --shell bash)"
micromamba activate minimap2.28

minimap2 --eqx --secondary=no -t ${THREADS} -Y -c -ax lr:hqae hg002v1.0.1.fasta <(seqtk seq -F "?" WGS_HG002_EZ1_25kb_SUPNew.same_reads.fasta|paste - - - - |fgrep -f Newsame_reads_ids2.txt |tr '\t' '\n') |samtools sort -@${THREADS} > WGS_HG002_EZ1_25kb_SUPNew.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam
samtools index -@ ${THREADS} WGS_HG002_EZ1_25kb_SUPNew.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam

minimap2 --eqx --secondary=no -t ${THREADS} -Y -c -ax lr:hqae hg002v1.0.1.fasta <(seqtk seq -F "?" WGS_HG002_EZ1_25kb.herroNew.same_reads.fasta|paste - - - - |fgrep -f Newsame_reads_ids2.txt |tr '\t' '\n') |samtools sort -@${THREADS} > WGS_HG002_EZ1_25kb.herroNew.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam
samtools index -@ ${THREADS} WGS_HG002_EZ1_25kb.herroNew.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam


minimap2 --eqx --secondary=no -t ${THREADS} -Y -c -ax lr:hqae hg002v1.0.1.fasta <(seqtk seq -F "?" WGS_HG002_EZ1_25kb.brNew.same_reads.fasta |paste - - - - |fgrep -f Newsame_reads_ids2.txt |tr '\t' '\n') |samtools sort -@${THREADS} > WGS_HG002_EZ1_25kb.brNew.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam
samtools index -@ ${THREADS} WGS_HG002_EZ1_25kb.brNew.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam

minimap2 --eqx --secondary=no -t ${THREADS}-Y -c -ax lr:hqae hg002v1.0.1.fasta <(seqtk seq -F "?" WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.same_reads.fasta |paste - - - - |fgrep -f Newsame_reads_ids2.txt |tr '\t' '\n') |samtools sort -@${THREADS} > WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam
samtools index -@ ${THREADS} WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam


minimap2 --eqx --secondary=no -t ${THREADS} -Y -c -ax lr:hqae hg002v1.0.1.fasta <(seqtk seq -F "?" WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.dechat.same_reads.fasta |paste - - - - |fgrep -f Newsame_reads_ids2.txt |tr '\t' '\n') |samtools sort -@${THREADS} > WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.dechat.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam
samtools index -@ ${THREADS} WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.dechat.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam
```


Now run BEST

```bash
~/git/best/target/release/best -t 48 WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.dechat.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam hg002v1.0.1.fasta dechat.best
~/git/best/target/release/best -t 48 WGS_HG002_EZ1_25kb.herro.brk19.pg_asm1x.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam hg002v1.0.1.fasta pg.best
~/git/best/target/release/best -t 48 WGS_HG002_EZ1_25kb.herroNew.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam hg002v1.0.1.fasta herro.best
~/git/best/target/release/best -t 48 WGS_HG002_EZ1_25kb.brNew.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam hg002v1.0.1.fasta br.best
~/git/best/target/release/best -t 48 WGS_HG002_EZ1_25kb_SUPNew.same_reads.fasta_against_hg002v1.0.1.fasta.Newsame_reads_ids2.bam hg002v1.0.1.fasta SUP.best
```

output from *.best.summary_identity_stats.csv

```bash
total_alns  primary_alns  identity  identity_qv  gap_compressed_identity  matches_per_kbp  mismatches_per_kbp  non_hp_ins_per_kbp  non_hp_del_per_kbp  hp_ins_per_kbp  hp_del_per_kbp  data_set
1146030     1145743       0.999741  35.862487    0.999842                 999.803664       0.027125            0.023584            0.022643            0.039365        0.146568        dechat
1146025     1145743       0.999726  35.619356    0.999830                 999.790755       0.031596            0.024639            0.026350            0.040332        0.151300        pg
1156816     1145743       0.999637  34.398948    0.999772                 999.725298       0.051433            0.033581            0.044946            0.054916        0.178323        br
1156746     1145743       0.999637  34.396470    0.999770                 999.725438       0.052275            0.033748            0.044347            0.055095        0.177940        herro
1281942     1145743       0.980335  17.062956    0.984901                 985.021819       7.253553            2.826602            3.413819            1.954716        4.310809        raw
```

Convert WGS_HG002_EZ1_25kb.pod5.bam to cram by aligning to GCF_009914755.1 (T2T-CHM13v2.0_genomic.fna)
to be able to upload to Zenodo in one piece

```bash
micromamba activate minimap2.28
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz

echo "9e6bf6b586bc8954208d1cc1d5f2fc99  ./GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" > checksum
md5sum -c checksum
# ok

pigz -kd GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
minimap2 -t 54 -Y -y --eqx -ax lr:hq GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz <(samtools fastq -@4 -T "*" WGS_HG002_EZ1_25kb.pod5.bam) 2>minimap2.log |samtools sort | samtools view -@8 -h -C -T GCF_009914755.1_T2T-CHM13v2.0_genomic.fna > WGS_HG002_EZ1_25kb.pod5.aligned-to-GCF_009914755.1.cram &

samtools index -@54 WGS_HG002_EZ1_25kb.pod5.aligned-to-GCF_009914755.1.cram
```


Convert CRAM back to BAM if desired to get WGS_HG002_EZ1_25kb.pod5.bam from WGS_HG002_EZ1_25kb.pod5.aligned-to-GCF_009914755.1.cram

```bash
samtools view -@ 4 -h -T GCF_009914755.1_T2T-CHM13v2.0_genomic.fna -O BAM WGS_HG002_EZ1_25kb.pod5.aligned-to-GCF_009914755.1.cram > WGS_HG002_EZ1_25kb.pod5.converted.bam
```
