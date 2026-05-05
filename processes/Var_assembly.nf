#!/usr/bin/env nextflow
/*
asssembly of var gene with spades
*/

// 1. Assembly with SPAdes and de-replication with CD-HIT
process AssembleDereplicate {
    tag "SPAdes+CD-HIT - ${sample_name}"
    conda 'bioconda::spades=4.1.0,bioconda::cd-hit=4.8.1'
    
    errorStrategy { task.exitStatus == 137 ? 'ignore' : 'retry' }

    publishDir "${params.out_dir}/Assembly/spades_Cd-hit_out", mode: 'copy', overwrite: true
    
    cpus 12

    input:
    tuple path(read1), path(read2), path(process_counting_size), val(sample_name)

    output:
    tuple path("${sample_name}_spades/contigs.fasta"), val(sample_name), emit: spades_outputs
    tuple path("${sample_name}_spades/${sample_name}_R1.fastq"), 
        path("${sample_name}_spades/${sample_name}_R2.fastq"),
        path("${sample_name}_spades/${sample_name}_unique_contigs_300bp.fasta"), 
        val(sample_name), 
        emit: unique_contigs
    path("${sample_name}_spades/${sample_name}_unique_contigs.fasta")

    script:
    """
    echo "Assembling reads with RNASPAdes for sample: ${sample_name}..."
    rnaspades.py -t ${task.cpus} -o ${sample_name}_spades -k 71 -1 ${read1} -2 ${read2}
    mv ${sample_name}_spades/transcripts.fasta ${sample_name}_spades/contigs.fasta
    
    echo "Removing similar contigs with CD-HIT for sample: ${sample_name}..."
    cd-hit-est -i ${sample_name}_spades/contigs.fasta \
        -o ${sample_name}_spades/${sample_name}_unique_contigs.fasta \
        -c 0.99 -M 0 -T ${task.cpus}
    # count the countigs size and filter out the contigs < 300 bp 
    perl ${process_counting_size} ${sample_name}_spades/${sample_name}_unique_contigs.fasta 300 \
    ${sample_name}_spades/${sample_name}_unique_contigs_300bp.fasta

    echo "Copying input reads to output directory..."
    cp ${read1} ${sample_name}_spades/${sample_name}_R1.fastq
    cp ${read2} ${sample_name}_spades/${sample_name}_R2.fastq
    """
}

// 2. Mapping the polymorphic reads to the full unique contigs with bowtie2 and cclean up the contigs with pilon
process MapPolyReadsContigs {

    tag "Mapping - ${sample_name}"

    conda 'bioconda::bowtie2=2.5.4', 'bioconda::samtools=1.10'
    
    errorStrategy { task.exitStatus == 137 ? 'ignore' : 'retry' }

    publishDir "${params.out_dir}/Assembly/mapped_read_to_unique_contigs", mode: 'copy', overwrite: true
    
    cpus 10

    input:
    tuple path(read1), path(read2), path(unique_contigs), val(sample_name)

    output:
    path("${sample_name}_mapping_stats.txt")
    tuple path("${unique_contigs}"), path("${sample_name}_mapped_to_unique_contigs.bam"), path("${sample_name}_mapped_to_unique_contigs.bam.bai"), val(sample_name), emit: bowtie2_contigs_mapped

    script:
    """

    # Build index
    bowtie2-build ${unique_contigs} ${sample_name}_index

    # Align reads and sort BAM
    bowtie2 -x ${sample_name}_index -1 ${read1} -2 ${read2} --threads ${task.cpus} --sensitive-local | \\
        samtools view -@ ${task.cpus} -S -b - | \\
        samtools sort -@ ${task.cpus} -o ${sample_name}_mapped_to_unique_contigs.bam

    # Compute mapping stats
    samtools flagstat -@ ${task.cpus} ${sample_name}_mapped_to_unique_contigs.bam > ${sample_name}_mapping_stats.txt

    # Index BAM
    samtools index -@ ${task.cpus} ${sample_name}_mapped_to_unique_contigs.bam

    # Clean up
    """
}

// 3. Assembly correction with Pilon
process AssemblyCorrection {
    
    tag "Mapping - ${sample_name}"
    
    conda 'bioconda::pilon=1.24', 'bioconda::assembly-stats=1.0.0'
    
    errorStrategy { task.exitStatus == 137 ? 'ignore' : 'retry' }

    publishDir "${params.out_dir}/Assembly/correction_contigs_pilon", mode: 'copy', overwrite: true

    cpus 2
    memory '64 GB'

    input:
    tuple path(unique_contigs), path(mapped_to_unique_contigs_bam), path(mapped_to_unique_contigs_bam_bai), val(sample_name)

    output:
    tuple path("${sample_name}_genome.fasta"), path(mapped_to_unique_contigs_bam), val(sample_name), emit: pilon_genome
    path("${sample_name}_AssemblyStats.txt")

    script:
    
    """
    export JAVA_OPTS="-Xmx${task.memory.toGiga()}G"
    
    samtools index ${mapped_to_unique_contigs_bam}

    pilon --genome ${unique_contigs} --frags ${mapped_to_unique_contigs_bam} \
        --output ${sample_name}_pilon --chunksize 50000 --fix bases --nostrays --outdir . # add the option "--nostrays" if you have memory issue

    # Clean headers in Pilon output FASTA
    awk '/^>/ {gsub(/_pilon/, "")} {print}' ${sample_name}_pilon.fasta > ${sample_name}_genome.fasta

    # Assembly statistics
    assembly-stats ${sample_name}_genome.fasta > ${sample_name}_AssemblyStats.txt

    # Remove intermediate files
    rm ${sample_name}_pilon.fasta

    """
}


