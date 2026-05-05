#!/usr/bin/env nextflow

//Identify reads derived from non-core genes

/* 
    Assumes that you have reference genomes (human genome, Plasmodium genome, varDB exon 1) 
    fastq files (sample_name_1_trimmed.fastq.gz and sample_name_2_trimmed.fastq.gz) for each sample and samples_list.txt to be process
    1. Remove PfEMP1 and stevor regions from the Plasmodium genome
    (this is done to avoid PfEMP1 and stevor reads in the polymorphic read set)
    3. Extract the polymorphic read 
    3.1. Index genomes (human genome, pf genome, varDB exon 1)
    3.2. Map reads to Human genome and extract unmapped read IDs
    3.3. Map reads to Plasmodium genome and extract unmapped read IDs
    3.4. Map the  reads against var3kb exon 1 database of 2000 field isolate and extract mapped read IDs
    3.5. Combine unmapped and mapped read IDs to get polymorphic read set
    
*/

// 1. Remove PfEMP1 from genome
process RemovePfEMP1fromGenome {
    tag "remove PfEMP1 and stevor from genome"
    conda 'bioconda::bedtools=2.30.0'

    errorStrategy { task.exitStatus == 137 ? 'ignore' : 'retry' }

    publishDir "${params.out_dir}/Reference/Plasmodium_no_pfemp1", mode: 'copy', overwrite: true
    
    input:
    path pf_genome
    path pf_genome_gff_file

    output:
    path "PlasmoDB-68_Pfalciparum3D7_noPfEMP1.fasta", emit: pf_genome_no_pfemp1
    path "pfemp1_regions.bed", emit: pfemp1_bed
    
    script:
    """
    echo "Removing PfEMP1 and stevor regions from the Plasmodium genome..."

    awk -F '\\t' '(\$3 == "protein_coding_gene" || \$3 == "pseudogene") && tolower(\$9) ~ /pfemp1|erythrocyte membrane protein 1|stevor/ {
        start = (\$4 - 1 - 500); 
        end = (\$5 + 500); 
        if (start < 0) start = 0;
        print \$1, start, end
    }' OFS='\\t' ${pf_genome_gff_file} > pfemp1_regions.bed

    bedtools maskfasta -fi ${pf_genome} -bed pfemp1_regions.bed -fo PlasmoDB-68_Pfalciparum3D7_noPfEMP1.fasta
    """
    
}

// 2. Index all genome
process IndexGenomes {
    conda 'bioconda::subread=2.1.1'

    publishDir "${params.out_dir}/Reference/indexes", mode: 'copy', overwrite: true

    errorStrategy { task.exitStatus == 137 ? 'ignore' : 'retry' }
    
    input:
    path human_genome
    path pf_genome_no_pfemp1
    path vardb_exon1

    output:
    path "${human_genome.baseName}.*", emit: human_genome_index
    path "${human_genome}", emit: human_genome_reference
    path "${pf_genome_no_pfemp1.baseName}.*", emit: pf_genome_index
    path "${pf_genome_no_pfemp1}", emit: pf_reference
    path "${vardb_exon1.baseName}.*", emit: vardb_exon1_index
    path "${vardb_exon1}", emit: vardb_exon1_reference
    
    script:
    """
    echo "Indexing human genome, Plasmodium genome (no PfEMP1) and varDB database ..."

    # Index human genome
    echo "1. Indexing human genome..."

    subread-buildindex -o ${human_genome.baseName} ${human_genome}

    # Index the Plasmodium genome (no PfEMP1)
    echo "2. Indexing Plasmodium genome (no PfEMP1)..."

    subread-buildindex -o ${pf_genome_no_pfemp1.baseName} ${pf_genome_no_pfemp1}

    # Index the varDB exon 1 dataset
    echo "3. Indexing varDB exon 1 dataset..."

    subread-buildindex -o ${vardb_exon1.baseName} ${vardb_exon1}
    """
    
}

// 3. Extract polymorphic reads
process ExtractPolymorphicPfemp1Reads {
    tag "$sample_name"
    conda 'bioconda::subread=2.1.1,bioconda::samtools=1.10,bioconda::seqtk=1.3'

    errorStrategy { task.exitStatus == 137 ? 'ignore' : 'retry' }

    publishDir "${params.out_dir}/polymorphic_reads", mode: 'copy', overwrite: true
    
    cpus 4

    input:
    tuple path(read1), path(read2), val(sample_name)
    path human_genome_reference
    path human_genome_index
    path pf_reference
    path pf_genome_index
    path vardb_exon1_reference
    path vardb_exon1_index

    output:
    tuple path("${sample_name}_polymorphic_reads_1.fq"), path("${sample_name}_polymorphic_reads_2.fq"), val(sample_name), emit: polymorphic_reads
    path("${sample_name}_polymorphic_ids.txt"), emit: polymorphic_read_ids
    
    
    script:
    """
    # Step 1: Extract paired reads where BOTH mates are unmapped to human genome
    subread-align -t 0 -i ${human_genome_reference.baseName} -r ${read1} -R ${read2} --SAMoutput -T ${task.cpus} -o ${sample_name}_human_mapped.sam
    samtools view -@ ${task.cpus} -f 4 ${sample_name}_human_mapped.sam | \
        cut -f1 | sort -u > ${sample_name}_human_both_unmapped_read_ids.txt

    [ -s "${sample_name}_human_both_unmapped_read_ids.txt" ] || touch "${sample_name}_human_both_unmapped_read_ids.txt"

    seqtk subseq ${read1} ${sample_name}_human_both_unmapped_read_ids.txt > "${sample_name}_human_unmapped_1.fq"
    seqtk subseq ${read2} ${sample_name}_human_both_unmapped_read_ids.txt > "${sample_name}_human_unmapped_2.fq"

    # Step 2: Extract paired reads where BOTH mates are unmapped to pf genome
    subread-align -t 0 -i ${pf_reference.baseName} -r ${sample_name}_human_unmapped_1.fq -R ${sample_name}_human_unmapped_2.fq --SAMoutput -T ${task.cpus} -o ${sample_name}_mapped_pf_novar.sam
    samtools view -@ ${task.cpus} -f 4 ${sample_name}_mapped_pf_novar.sam | \
        cut -f1 > ${sample_name}_pf_both_unmapped_read_ids.txt

    [ -s "${sample_name}_pf_both_unmapped_read_ids.txt" ] || touch "${sample_name}_pf_both_unmapped_read_ids.txt"

    # Step 3: Extract paired reads where BOTH mates are mapped to exon 1 DB
    subread-align -t 0 -i ${vardb_exon1_reference.baseName} -r ${sample_name}_human_unmapped_1.fq -R ${sample_name}_human_unmapped_2.fq --SAMoutput -T ${task.cpus}  -o ${sample_name}_var3kb_mapped.sam
    samtools view -@ ${task.cpus} -f 4 ${sample_name}_var3kb_mapped.sam > ${sample_name}_vardb_both_mapped_read_ids.txt

    [ -s "${sample_name}_vardb_both_mapped_read_ids.txt" ] || touch "${sample_name}_vardb_both_mapped_read_ids.txt"

    # Step 4: Combine unmapped+mapped read IDs to get polymorphic set
    cat ${sample_name}_pf_both_unmapped_read_ids.txt ${sample_name}_vardb_both_mapped_read_ids.txt | \
        sort -u > ${sample_name}_polymorphic_ids.txt

    if [ ! -s "${sample_name}_polymorphic_ids.txt" ]; then
        echo "WARNING: No polymorphic read IDs found - creating empty outputs" >&2
        touch "${sample_name}_polymorphic_reads_1.fq"
        touch "${sample_name}_polymorphic_reads_2.fq"
    else
        # Step 5: Extract matching reads from FASTQ
        seqtk subseq ${read1} ${sample_name}_polymorphic_ids.txt > ${sample_name}_polymorphic_reads_1.fq
        seqtk subseq ${read2} ${sample_name}_polymorphic_ids.txt > ${sample_name}_polymorphic_reads_2.fq
    fi

    # Cleanup intermediate files
    rm -f ${sample_name}_pf_both_unmapped.sorted.bam ${sample_name}_vardb_both_mapped.sorted.bam \
        ${sample_name}_pf_both_unmapped_read_ids.txt ${sample_name}_vardb_both_mapped_read_ids.txt \
        ${sample_name}_human_unmapped_1.fq ${sample_name}_human_unmapped_2.fq \
        ${sample_name}_human_both_unmapped.sorted.bam ${sample_name}_human_both_unmapped_read_ids.txt \
        ${sample_name}_human_mapped.sam ${sample_name}_mapped_pf_novar.bam \
        ${sample_name}_var3kb_mapped.sam
    """
    
}
