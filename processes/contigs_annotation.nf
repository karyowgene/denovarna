#!/usr/bin/env nextflow
/*
Annotation of the var gene assembly
*/

// 1. Annotate and resolve the assembled contigs
process Annotate_assembled_contigs {
    tag "Annotate transcripts - ${sample_name}"
    conda 'bioconda::hmmer=3.4, bioconda::cath-resolve-hits=4.2.0, bioconda::seqtk=1.3, conda install conda-forge::python'

    publishDir "${params.out_dir}/Annotation/annotated_transcripts", mode: 'copy', overwrite: true

    input:
    tuple path(sspace_out_scaf_seq), path(var_domain), path(process_hmm), val(sample_name)

    output:
    path("Annotated_${sample_name}/${sample_name}_assembled_contigs_hmm"), emit: assembled_contigs_hmm
    path("Annotated_${sample_name}/${sample_name}_assembled_contigs_hmm_cathresolvehits"), emit: assembled_contigs_hmm_cathresolvehits, optional: true
    path("Annotated_${sample_name}/${sample_name}_assembled_id_sig_annotation.txt"), emit: assembled_id_sig_annotation_txt, optional: true
    path("Annotated_${sample_name}/${sample_name}_assembled_id_and_var_annotation.csv"), emit: assembled_id_annotation_csv, optional: true
    path("Annotated_${sample_name}/${sample_name}_assembled_id_sig_annotation.fasta"), emit: assembled_id_sig_annotation_fasta, optional: true

    script:
    """
    mkdir -p Annotated_${sample_name}

    # Annotation of the assembled contigs
    hmmsearch --domtblout Annotated_${sample_name}/${sample_name}_assembled_contigs_hmm -E 1e-5 ${var_domain} ${sspace_out_scaf_seq}

    # Resolve the assembled contigs
    cath-resolve-hits --input-format hmmer_domtblout Annotated_${sample_name}/${sample_name}_assembled_contigs_hmm \
        --hits-text-to-file Annotated_${sample_name}/${sample_name}_assembled_contigs_hmm_cathresolvehits

    # Only continue if cath-resolve-hits output is non-empty
    if [[ -s Annotated_${sample_name}/${sample_name}_assembled_contigs_hmm_cathresolvehits ]]; then
        # Process hmm annotation
        python ${process_hmm} Annotated_${sample_name}/${sample_name}_assembled_contigs_hmm_cathresolvehits ${sample_name}

        # Subset the contigs with significant annotations
        seqtk subseq ${sspace_out_scaf_seq} assembled_id_sig_annotation.txt \
            > ${sample_name}_assembled_id_sig_annotation.fasta

        # Move outputs into final directory
        mv assembled_id_sig_annotation.txt Annotated_${sample_name}/${sample_name}_assembled_id_sig_annotation.txt
        mv ${sample_name}_assembled_id_and_var_annotation.csv Annotated_${sample_name}/${sample_name}_assembled_id_and_var_annotation.csv
        mv ${sample_name}_assembled_id_sig_annotation.fasta Annotated_${sample_name}/${sample_name}_assembled_id_sig_annotation.fasta
    else
        echo "WARNING: No significant HMM hits found for ${sample_name}. Skipping annotation steps." >&2
    fi
    """
}

