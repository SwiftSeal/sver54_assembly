params.hifi_reads =
params.nanopore_reads =

process Hifiasm {
    container 'https://depot.galaxyproject.org/singularity/hifiasm:0.19.8--h43eeafb_0'
    cpus 32
    memory '64GB'
    time '12h'
    input:
    path hifi_reads
    output:
    path hifiasm.p_ctg.fa
    script:
    """
    hifiasm --primary -t 32 $hifi_reads -o hifiasm
    awk '/^S/{{print ">"$2;print $3}}' hifiasm.p_ctg.gfa > hifiasm.p_ctg.fa
    """
}

process Flye {
    container
    cpus 16
    memory '32GB'
    time '12h'
    input:
    path nanopore_reads
    output:
    path flye/assembly.fasta
    script:
    """
    flye --nano-hq $nanopore_reads --out-dir flye --threads ${task.cpus}
    """
}

process Quickmerge {

}
