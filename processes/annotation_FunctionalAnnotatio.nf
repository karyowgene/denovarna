#!/usr/bin/env nextflow
/*
var annotation

*/
// 1. Augustus gene ad protein annotation
process VargeneAnnotation {
    conda 'bioconda::augustus=3.5.0', 'bioconda::subread=2.1.1', 'bioconda::blast=2.16.0', 'bioconda::seqkit=2.10.0'
    
    errorStrategy { task.exitStatus == 137 ? 'ignore' : 'retry' }

    tag "Vargene annotation - ${sample_name}"

    publishDir "${params.out_dir}/Annotation/augustus_annotation", mode: 'copy', overwrite: true
    
    cpus 6

    input:
    tuple path(genome_fasta), path(assembly_bam), val(sample_name)
    
    output:
    tuple path("${sample_name}.genes.aa.fasta"), path("${sample_name}.exons.na.fasta"), path("${sample_name}.genes.na.fasta"), val(sample_name), emit: augustus_aa_ex_na_outputs
    tuple path("${sample_name}.genes.aa.fasta"), val(sample_name), emit: augustus_aa_outputs 
    path("${sample_name}.out2.augustus.txt")
    tuple path("${sample_name}.annotation.gtf"), path("${sample_name}.annotation.gff"), val(sample_name), emit: annotation_outputs
    tuple path("${sample_name}.counts.txt"), val(sample_name), emit: augustus_stats_outputs
    tuple path("${sample_name}.genome_motif_hits.blastx.txt"), path("${sample_name}.genome_motif_categories.txt"), path("${sample_name}.genome_motif_counts.txt"), val(sample_name), emit: genome_DBLa_motif_count_outputs
    tuple path("${sample_name}.augustus_motif_counts.txt"), path("${sample_name}_LARSFADIG.fasta"), path("${sample_name}_DYVPQYLRW.fasta"), path("${sample_name}_LARSFADIG_DYVPQYLRW.fasta"), val(sample_name), emit: augustus_motif_outputs

    script:
    """
    # Run AUGUSTUS
    augustus --species=pf_VSA ${genome_fasta} > ${sample_name}.out.augustus.txt

    # Add sample name to gene IDs in AUGUSTUS output
    perl -e '
        my \$n = shift;
        while (<>) {
            chomp;
            my @ar = split(/\\t/);
            if (defined \$ar[8]) {
                \$ar[8] =~ s/\\b(g\\d+)/\$n.\$1/g;
                print join("\\t", @ar) . "\\n";
            } else {
                print \$_ . "\\n";
            }
        }
    ' ${sample_name} ${sample_name}.out.augustus.txt > ${sample_name}.out2.augustus.txt

    # Generate FASTA files from AUGUSTUS output
    getAnnoFasta.pl --seqfile=${genome_fasta} ${sample_name}.out2.augustus.txt
    mv ${sample_name}.out2.augustus.aa ${sample_name}.genes.aa.fasta
    mv ${sample_name}.out2.augustus.cdsexons ${sample_name}.exons.na.fasta
    mv ${sample_name}.out2.augustus.codingseq ${sample_name}.genes.na.fasta

    # Create GTF and GFF annotation files for featureCounts
    awk -v OFS='\\t' '\$2 == "AUGUSTUS"' ${sample_name}.out2.augustus.txt \\
        | grep -E "transcript|CDS" \\
        | sed 's/CDS/exon/g' \\
        | perl -e '
            my %CDS;
            while (<STDIN>) {
                if (/exon/) {
                    chomp;
                    my @ar = split(/\\t/);
                    \$CDS{\$ar[0]}++;
                    print \$_, " exon_number \\"", \$CDS{\$ar[0]}, "\\";\\n";
                }
            }' > ${sample_name}.annotation.gtf

    grep -v "^#" ${sample_name}.out2.augustus.txt > ${sample_name}.annotation.gff
        
    # Count features using featureCounts (expects BAM file present with specific name)
    samtools index -@ ${task.cpus} ${assembly_bam}
    featureCounts -p -a ${sample_name}.annotation.gtf -o ${sample_name}.counts.txt ${assembly_bam}


    # Search for LARSFADIG and DYVPQYLRW with up to 2 mismatches
    seqkit grep -s -m 2 -p "LARSFADIG" ${sample_name}.genes.aa.fasta > ${sample_name}_LARSFADIG.fasta
    seqkit grep -s -m 2 -p "DYVPQYLRW" ${sample_name}.genes.aa.fasta > ${sample_name}_DYVPQYLRW.fasta
    seqkit grep -s -m 2 -p "LARSFADIG" ${sample_name}.genes.aa.fasta | \\
        seqkit grep -s -m 2 -p "DYVPQYLRW" > ${sample_name}_LARSFADIG_DYVPQYLRW.fasta

    seqkit stats ${sample_name}_LARSFADIG.fasta ${sample_name}_DYVPQYLRW.fasta \\
        ${sample_name}_LARSFADIG_DYVPQYLRW.fasta > ${sample_name}.augustus_motif_counts.txt

    # Count LARSFADIG and DYVPQYLRW motifs in the genome using blastx

    
    # Create a FASTA file with the motifs
    cat <<EOF > motifs.faa
    >LARSFADIG
    LARSFADIG
    >DYVPQYLRW
    DYVPQYLRW
    EOF

    # Run BLASTX
    blastx -query ${genome_fasta} -subject motifs.faa \\
        -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \\
        -evalue 1e-2 -word_size 2 -max_target_seqs 100000 -gapopen 9 -gapextend 1 \\
        -out ${sample_name}.genome_motif_hits.blastx.txt

    # Count motif matches
    awk '{print \$1, \$2}' ${sample_name}.genome_motif_hits.blastx.txt | sort | uniq > tmp.motif_matches.txt

    # Categorize contigs by motif presence
    awk '
    {
        contig=\$1; motif=\$2
        if (motif == "LARSFADIG") hasLARS[contig]=1
        if (motif == "DYVPQYLRW") hasDYVP[contig]=1
        allContigs[contig]=1
    }
    END {
        for (c in allContigs) {
            if (hasLARS[c] && hasDYVP[c]) {
                print c, "BOTH"
            } else if (hasLARS[c]) {
                print c, "LARSFADIG"
            } else if (hasDYVP[c]) {
                print c, "DYVPQYLRW"
            } else {
                print c, "NONE"
            }
        }
    }' tmp.motif_matches.txt > ${sample_name}.genome_motif_categories.txt

    # Count per category
    awk '{count[\$2]++} END { for (m in count) print m, count[m] }' ${sample_name}.genome_motif_categories.txt \
        > ${sample_name}.genome_motif_counts.txt

    #rm motifs.faa tmp.motif_matches.txt
    """
}

// 2. Domains annotation
process DomainsAnnotation {
    tag "$sample_name"
    conda 'bioconda::hmmer=3.4'
    
    errorStrategy { task.exitStatus == 137 ? 'ignore' : 'retry' }
    
    publishDir "${params.out_dir}/Annotation/domains_annotation", mode: 'copy'
    
    cpus 16

    input:
    tuple path(gene_aa_fasta), path(exons_na_fasta), path(gene_na_fasta), path(domainpathHmmer), path(domainlengh),
        path(vargene_parseDBLa), path(vargene_Domain_hiearchical), path(order_domains), val(sample_name)

    output:
    tuple path(gene_aa_fasta), path(exons_na_fasta), path(gene_na_fasta), path("${sample_name}.Domain.gff"), val(sample_name), emit : domains_gff_outputs
    path "${sample_name}.Seqhits.txt"
    path "${sample_name}.Domhit.domain.txt"
    path "${sample_name}.Domain.association.structured.txt"

    script:
    """
    # Run hmmscan
    # Press the HMM database if needed
    if [ ! -f "${domainpathHmmer}.h3m" ]; then
        hmmpress ${domainpathHmmer}
    fi

    hmmscan --cpu ${task.cpus} --domT 50 -E 1e-6 --tblout ${sample_name}.Seqhits.txt \
        --domtblout ${sample_name}.Domhit.domain.txt \
        ${domainpathHmmer} ${gene_aa_fasta} &> /dev/null

    # Parse HMMER output to GFF
    grep -v "^#" ${sample_name}.Domhit.domain.txt | \
        perl ${vargene_parseDBLa} | \
        sort -nrk 14 | \
        perl ${vargene_Domain_hiearchical} ${domainlengh} | \
        sort -nrk 14 > ${sample_name}.Domain.gff

    # Domain structured summary with correct order
    python ${order_domains} ${sample_name}.Domhit.domain.txt > ${sample_name}.Domain.association.structured.txt

    """
}

// 2. Subdomains annotation
process SubDomainsAnnotation {
    tag "$sample_name"
    conda 'bioconda::hmmer=3.4'
    
    errorStrategy { task.exitStatus == 137 ? 'ignore' : 'retry' }
    
    publishDir "${params.out_dir}/Annotation/subdomains_annotation", mode: 'copy'
    
    cpus 6

    input:
    tuple path(gene_aa_fasta), path(suddomainpathHmmer), path(subdomainlengh),
        path(vargene_Domain_hiearchical), path(order_subdomains), val(sample_name)

    output:
    path "${sample_name}.Seqhits_sub.txt"
    path "${sample_name}.Domhit.suddomain.txt"
    path "${sample_name}.Subdomain.gff"
    path "${sample_name}.Subdomain.association.structured.txt"

    script:
    """
    # Run hmmscan
    # Press the HMM database if needed
    if [ ! -f "${suddomainpathHmmer}.h3m" ]; then
        hmmpress ${suddomainpathHmmer}
    fi

    hmmscan --cpu ${task.cpus} --domT 50 -E 1e-6 --tblout ${sample_name}.Seqhits_sub.txt \
        --domtblout ${sample_name}.Domhit.suddomain.txt \
        ${suddomainpathHmmer} ${gene_aa_fasta} &> /dev/null

    # Parse HMMER output to GFF
    grep -v "^#" ${sample_name}.Domhit.suddomain.txt | \
        sort -nrk 14 | \
        perl ${vargene_Domain_hiearchical} ${subdomainlengh} | \
        sort -nrk 14 | \
        sed 's/Domain\\.//g'> ${sample_name}.Subdomain.gff

    # Domain structured summary with correct order
    sed 's/Domain\\.//g' ${sample_name}.Domhit.suddomain.txt | python ${order_subdomains} > ${sample_name}.Subdomain.association.structured.txt

    """
}

// 3. Blast annotation

process BlastAnnotation {
    tag "$sample_name"
    conda 'bioconda::blast=2.16.0'
    
    errorStrategy { task.exitStatus == 137 ? 'ignore' : 'retry' }
    
    publishDir "${params.out_dir}/Annotation/blast_3D7ref_annotation", mode: 'copy'
    
    cpus 6

    input:
    tuple path(gene_aa_fasta), path(gene_na_fasta), path(exons_na_fasta), path(domain_gff), path(Pf3D7_aa_fa),
        path(List_3D7_Annotation), path(vargene_blastAnntationtransfer_pl), val(sample_name)

    output:
    path "${sample_name}.comp.aa.3d7.blast"
    path "${sample_name}.3D7ref.genes.na.fasta"
    path "${sample_name}.3D7ref.exons.na.fasta"
    path "${sample_name}.3D7ref.genes.aa.fasta"
    path "${sample_name}.VARgenes.na.fasta"
    path "${sample_name}.VARgenes.aa.fasta"

    script:
    """
    # Prepare BLAST database from Pf3D7_aa_fa if not already formatted
    makeblastdb -in ${Pf3D7_aa_fa} -dbtype prot -out Pf3D7db

    # Run BLASTP
    blastp -db Pf3D7db -query ${gene_aa_fasta} -outfmt 6 -num_threads ${task.cpus} \
        -max_target_seqs 1 > ${sample_name}.comp.aa.3d7.blast

    # Transfer annotations to gene, exon, and aa fasta
    cat ${gene_na_fasta} | sed '/^>/ s/\\.cds[0-9]*//g' | perl ${vargene_blastAnntationtransfer_pl} \
        ${List_3D7_Annotation} ${sample_name}.comp.aa.3d7.blast | sed 's/VAR1CSA/VAR/g' > ${sample_name}.3D7ref.genes.na.fasta

    cat ${exons_na_fasta} | sed '/^>/ s/\\.cds[0-9]*//g' | perl ${vargene_blastAnntationtransfer_pl} \
        ${List_3D7_Annotation} ${sample_name}.comp.aa.3d7.blast | sed 's/VAR1CSA/VAR/g' > ${sample_name}.3D7ref.exons.na.fasta

    cat ${gene_aa_fasta} | perl ${vargene_blastAnntationtransfer_pl} \
        ${List_3D7_Annotation} ${sample_name}.comp.aa.3d7.blast | sed 's/VAR1CSA/VAR/g' > ${sample_name}.3D7ref.genes.aa.fasta

    # Extract var genes
    n=\$(cut -f 1 ${domain_gff} | sort -u | paste -sd ' ' -)

    samtools faidx ${sample_name}.3D7ref.genes.na.fasta \$n | perl ${vargene_blastAnntationtransfer_pl} \
        ${List_3D7_Annotation} ${sample_name}.comp.aa.3d7.blast | sed 's/VAR1CSA/VAR/g' > ${sample_name}.VARgenes.na.fasta

    samtools faidx ${sample_name}.3D7ref.genes.aa.fasta \$n | perl ${vargene_blastAnntationtransfer_pl} \
        ${List_3D7_Annotation} ${sample_name}.comp.aa.3d7.blast | sed 's/VAR1CSA/VAR/g' > ${sample_name}.VARgenes.aa.fasta
    """
}
