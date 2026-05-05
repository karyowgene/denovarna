#!/usr/bin/env nextflow
/*
 * This is the main Nextflow pipeline for identifying non-core reads in Plasmodium genome sequencing data.
 * It includes processes for indexing the human genome and extracting unmapped reads.
 * This pipeline is designed for DSL2 and uses Nextflow's built-in features for process management and data handling.
 */

nextflow.enable.dsl=2

// Include modules
include { RemovePfEMP1fromGenome as RPE } from './processes/identify_non_core_reads.nf'
include { IndexGenomes as IG } from './processes/identify_non_core_reads.nf'
include { ExtractPolymorphicPfemp1Reads as Poly_reads } from './processes/identify_non_core_reads.nf'
include { AssembleDereplicate as VAR_ASEM_SPAdes_CDHIT } from './processes/Var_assembly.nf'
include { MapPolyReadsContigs as VAR_MAP_bowtie2 } from './processes/Var_assembly.nf'
include { AssemblyCorrection as VAR_CORR_Pilon} from './processes/Var_assembly.nf'
include { VargeneAnnotation as VAR_AA_Augustus} from './processes/annotation_FunctionalAnnotatio.nf'
include { DomainsAnnotation as VAR_AA_Domains} from './processes/annotation_FunctionalAnnotatio.nf'
include { SubDomainsAnnotation as VAR_AA_Subdomains} from './processes/annotation_FunctionalAnnotatio.nf'
include { BlastAnnotation as VAR_AA_Blast_3D7} from './processes/annotation_FunctionalAnnotatio.nf'

workflow {
    // Define input channels
    Channel.fromPath(params.human_genome).set { human_genome_ch }
    Channel.fromPath(params.pf_genome).set { pf_genome_ch }
    Channel.fromPath(params.pf_genome_gff).set { pf_genome_gff_ch }
    Channel.fromPath(params.vardb_exon1).set { vardb_exon1_ch }


    // Load sample list
    def sample_list = file(params.sample_list).readLines()*.trim()*.replaceAll(/[^a-zA-Z0-9_-]/, '')

    Channel
        .fromFilePairs("${params.trimmed_reads}/*_{R1,R2}_001.fastq.gz", size: 2)
        .filter { sample_name, reads -> sample_list.contains(sample_name) }
        .map { sample_name, reads -> tuple(reads[0], reads[1], sample_name) }
        .set { trimmed_reads_ch }

    // === Genome-Level Processes ===
    RPE(pf_genome_ch, pf_genome_gff_ch)
    IG(human_genome_ch, RPE.out.pf_genome_no_pfemp1, vardb_exon1_ch)

    // define the outputs directly to chanels
        // Get reference sequences (as value channels)
    def human_genome_ref_ch = IG.out.human_genome_reference.first()
    def pf_reference_ch = IG.out.pf_reference.first()
    def vardb_exon1_ref_ch = IG.out.vardb_exon1_reference.first()

    // Get complete index sets
    def human_genome_index_ch = IG.out.human_genome_index.first()
    def pf_genome_index_ch = IG.out.pf_genome_index.first()
    def vardb_exon1_index_ch = IG.out.vardb_exon1_index.first()

    // === Sample-Level Processes: extract polymorphic reads ===
    Poly_reads(trimmed_reads_ch, human_genome_ref_ch, human_genome_index_ch, pf_reference_ch, pf_genome_index_ch, vardb_exon1_ref_ch, vardb_exon1_index_ch)

    // === Sample-Level Processes: Assembly ( SPAdes + CD-HIT + countigs size and filter out the contigs < 300 bp pipeline) ===
            // add the option "--nostrays" if you have memory issue with pilon process at the " 3. Assembly correction with Pilon"
    def assembly_ch = Poly_reads.out.polymorphic_reads
    assembly_ch
    .map { read1, read2, sample_name -> 
        def process_counting_size = file(params.process_counting_contigs_size)
        tuple(read1, read2, process_counting_size, sample_name)
    }
    .set { assembly_dereplicate_ch }

    VAR_ASEM_SPAdes_CDHIT(assembly_dereplicate_ch)

    //=== Sample-Level Processes: Map the polymorphic reads to the full unique contigs with bowtie2 and cclean up the contigs with pilon ===
    VAR_MAP_bowtie2(VAR_ASEM_SPAdes_CDHIT.out.unique_contigs)
    VAR_CORR_Pilon(VAR_MAP_bowtie2.out.bowtie2_contigs_mapped)

    //=== Sample-Level Processes: Vargene annotation ===
            // Augustus gene and protein annotation
    VAR_AA_Augustus(VAR_CORR_Pilon.out.pilon_genome)

            // Domains annotation
    def gene_aa_ch = VAR_AA_Augustus.out.augustus_aa_ex_na_outputs
    gene_aa_ch
        .map { gene_aa, exon_na, gene_na, sample_name ->
            def vargene_parseDBLa = file(params.vargene_MP_parseDBLa_pl)
            def vargene_Domain_hiearchical = file(params.vargene_Domain_hmm2gff_hiearchical_pl)
            def order_domains = file(params.order_domains)
            def domainpathHmmer = file(params.domain_Path_Hmmer)
            def domainlengh = file(params.domain_Path_lengh)
            tuple(gene_aa, exon_na, gene_na, domainpathHmmer, domainlengh, vargene_parseDBLa, vargene_Domain_hiearchical, order_domains, sample_name)
        }
        .set {domains_annotation_ch}
    VAR_AA_Domains(domains_annotation_ch)

            // Subdomains annotation
    def gene_aa_s_ch = VAR_AA_Augustus.out.augustus_aa_outputs
    gene_aa_s_ch
        .map { gene_aa, sample_name ->
            def vargene_Domain_hiearchical = file(params.vargene_Domain_hmm2gff_hiearchical_pl)
            def sudomainpathHmmer = file(params.subdomain_Path_Hmmer)
            def subdomainlengh = file(params.subdomain_Path_lengh)
            def order_subdomains = file(params.order_subdomains)
            tuple(gene_aa, sudomainpathHmmer, subdomainlengh, vargene_Domain_hiearchical, order_subdomains, sample_name)
        }
        .set {subdomains_annotation_ch}
    VAR_AA_Subdomains(subdomains_annotation_ch)

        // Blast 3D7 annotation
    def gene_exon_an_ch = VAR_AA_Domains.out.domains_gff_outputs
    gene_exon_an_ch
        .map { gene_aa, exon_na, gene_na, domain_gff, sample_name ->
            def Pf3D7_aa_fa = file(params.pf_proteins_aa)
            def List_3D7_Annotation = file(params.list_3D7_annotation)
            def vargene_blastAnntationtransfer_pl = file(params.vargene_blastAnntationtransfer_pl)
            tuple(gene_aa, exon_na, gene_na, domain_gff, Pf3D7_aa_fa, List_3D7_Annotation, vargene_blastAnntationtransfer_pl, sample_name)
        }
        .set {blast_3d7_annotation}
    VAR_AA_Blast_3D7(blast_3d7_annotation)
}
