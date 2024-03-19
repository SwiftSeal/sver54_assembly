process SamtoolsFaidx {
    container 'https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2'
    cpus 1
    memory '4 GB'
    time '1h'
    input:
    path assembly
    output:
    path *.fa
    path *.fai
    script:
    """
    samtools faidx $input
    """
}

process BwaIndex {
    container 'https://depot.galaxyproject.org/singularity/bwa:0.7.16--pl5.22.0_0'
    cpus 1
    memory '4 GB'
    time '1h'
    input:
    path assembly
    output:
    path 'index'
    script:
    """
    mkdir -p index
    bwa index -p index/index $assembly
    """
}

process BwaMem {
    container 'https://depot.galaxyproject.org/singularity/bwa:0.7.16--pl5.22.0_0'
    cpus 32
    memory '64 GB'
    time '24h'
    input:
    path(index)
    path(R1_1)
    path(R1_2)
    path(R2_1)
    path(R2_2)
    output:
    path 'aligned.sam'
    script:
    """
    bwa mem \
      -5SP \
      -t ${task.cpus} \
      $index/index \
      <(zcat $R1_1 $R1_2) <(zcat $R2_1 $R2_2) \
      -o aligned.sam
    """
}

process PairtoolsParse {
    container 'https://depot.galaxyproject.org/singularity/pairtools:1.0.3--py39h9e08559_0'
    cpus 8
    memory '16 GB'
    time '6h'
    input:
    path assembly
    path aligned
    output:
    path 'parsed.pairsam'
    script:
    """
    pairtools parse \
      --min-mapq 40 \
      --walks-policy 5unique \
      --max-inter-align-gap 30 \
      --chroms-path $assembly \
      $aligned \
      > parsed.pairsam
    """
}

process PairtoolsSort {
    container 'https://depot.galaxyproject.org/singularity/pairtools:1.0.3--py39h9e08559_0'
    cpus 8
    memory '16 GB'
    time '6h'
    input:
    path parsed
    output:
    path 'sorted.pairsam'
    script:
    """
    pairtools sort \
      --nproc ${task.cpus} \
      $parsed \
      > sorted.pairsam
    """
}

process PairtoolsDeduplicate {
    publishDir 'scaffolding'
    container 'https://depot.galaxyproject.org/singularity/pairtools:1.0.3--py39h9e08559_0'
    cpus 8
    memory '16 GB'
    time '6h'
    input:
    path sorted
    output:
    path 'dedup.pairsam'
    path 'dedup.stats'
    script:
    """
    pairtools dedup \
      --mark-dups \
      --output-stats dedup.stats \
      --output dedup.pairsam \
      $sorted
    """
}

process PairtoolsSplit {
    container 'https://depot.galaxyproject.org/singularity/pairtools:1.0.3--py39h9e08559_0'
    cpus 8
    memory '16 GB'
    time '6h'
    input:
    path dedup
    output:
    path 'unsorted.sam'
    path 'mapped.pairs'
    script:
    """
    pairtools split \
      --output-pairs mapped.pairs \
      --output-sam unsorted.sam \
      $dedup
    """
}

process SamtoolsSort {
    container 'https://depot.galaxyproject.org/singularity/pairtools:1.0.3--py39h9e08559_0'
    cpus 8
    memory '16 GB'
    time '6h'
    input:
    path unsorted
    output:
    path 'sorted.bam'
    script:
    """
    samtools sort \
      --threads ${task.cpus} \
      -o sorted.bam \
      $unsorted
    """
}

process yahs {
    cpus 1
    memory '4 GB'
    time '1h'
    input:
    path bam
    path assembly
    output:
    path 'yahs_scaffolds_final.fa'
    path 'yahs_scaffolds_final.agp'
    path 'yahs.bin'
    script:
    """
    $APPS/yahs/yahs --no-contig-ec $assembly $bam -o yahs
    """
}

process JuicerPre {
    cpus 1
    memory '4 GB'
    time '1h'
    input:
    path bin
    path agp
    path assembly
    output:
    script:
    """
    $APPS/yahs/juicer pre -a -o juicer $bin $agp $assembly.fai > juicer.log 2>&1
    """
}

process JuicerTools {
    publishDir 'scaffolding'
    cpus 12
    memory '32 GB'
    time '6h'
    input:
    path juicer_txt
    path juicer_log
    output:
    path 'juicer.hic'
    script:
    """
    java -jar -Xmx32G $APPS/juicer_tools.1.9.9_jcuda.0.8.jar pre \
      $juicer_txt \
      juicer.hic \
      <(cat $juicer_log | grep PRE_C_SIZE | awk '{print \$2" "\$3}')
    """
}

workflow {
    index = BwaIndex(Channel.fromPath(params.assembly))
    aligned_reads = BwaMem(
        index,
        Channel.fromPath(params.R1_1),
        Channel.fromPath(params.R1_2),
        Channel.fromPath(params.R2_1),
        Channel.fromPath(params.R2_2)
    )
}
