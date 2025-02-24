version 1.0

import "AlignmentPipeline.wdl" as AlignAndMarkDuplicates

workflow AlignAndCall {
  meta {
    description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }

  input {
    File unmapped_bam
    String base_name

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File mt_amb
    File mt_ann
    File mt_bwt
    File mt_pac
    File mt_sa
    File blacklisted_sites
    File blacklisted_sites_index

    #Shifted reference is used for calling the control region (edge of mitochondria reference).
    #This solves the problem that BWA doesn't support alignment to circular contigs.
    File mt_shifted_dict
    File mt_shifted_fasta
    File mt_shifted_fasta_index
    File mt_shifted_amb
    File mt_shifted_ann
    File mt_shifted_bwt
    File mt_shifted_pac
    File mt_shifted_sa

    File shift_back_chain

    File? gatk_override
    String? m2_extra_args
    String? m2_filter_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Boolean compress_output_vcf
    Float? verifyBamID
  }

  parameter_meta {
    unmapped_bam: "Unmapped and subset bam, optionally with original alignment (OA) tag"
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToMt {
    input:
      input_bam = unmapped_bam,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_amb = mt_amb,
      mt_ann = mt_ann,
      mt_bwt = mt_bwt,
      mt_pac = mt_pac,
      mt_sa = mt_sa
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToShiftedMt {
    input:
      input_bam = unmapped_bam,
      mt_dict = mt_shifted_dict,
      mt_fasta = mt_shifted_fasta,
      mt_fasta_index = mt_shifted_fasta_index,
      mt_amb = mt_shifted_amb,
      mt_ann = mt_shifted_ann,
      mt_bwt = mt_shifted_bwt,
      mt_pac = mt_shifted_pac,
      mt_sa = mt_shifted_sa
  }

  call M2 as CallMt {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bai = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      # Everything is called except the control region.
      m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:576-16024 "
  }

  call M2 as CallShiftedMt {
    input:
      input_bam = AlignToShiftedMt.mt_aligned_bam,
      input_bai = AlignToShiftedMt.mt_aligned_bai,
      ref_fasta = mt_shifted_fasta,
      ref_fai = mt_shifted_fasta_index,
      ref_dict = mt_shifted_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      # Everything is called except the control region.
      m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:8025-9144 "
  }

  call LiftoverAndCombineVcfs {
    input:
      shifted_vcf = CallShiftedMt.raw_vcf,
      vcf = CallMt.raw_vcf,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      shift_back_chain = shift_back_chain
  }

  call MergeStats {
    input:
      shifted_stats = CallShiftedMt.stats,
      non_shifted_stats = CallMt.stats,
      gatk_override = gatk_override
  }

  call Filter as InitialFilter {
    input:
      raw_vcf = LiftoverAndCombineVcfs.merged_vcf,
      raw_vcf_index = LiftoverAndCombineVcfs.merged_vcf_index,
      raw_vcf_stats = MergeStats.stats,
      base_name = base_name,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      vaf_filter_threshold = 0.01,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta,
      run_contamination = false
  }

 
  call SplitMultiAllelicsAndRemoveNonPassSites {
    input:
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      filtered_vcf = InitialFilter.filtered_vcf,
      gatk_override = gatk_override
  }

  call GetContamination {
    input:
      input_vcf = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker
  }

  call Filter as FilterContamination {
    input:
      raw_vcf = InitialFilter.filtered_vcf,
      raw_vcf_index = InitialFilter.filtered_vcf_idx,
      raw_vcf_stats = MergeStats.stats,
      run_contamination = true,
      hasContamination = GetContamination.hasContamination,
      contamination_major = GetContamination.major_level,
      contamination_minor = GetContamination.minor_level,
      verifyBamID = verifyBamID,
      base_name = base_name,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      vaf_filter_threshold = vaf_filter_threshold,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta
  }

  output {
    File mt_aligned_bam = AlignToMt.mt_aligned_bam
    File mt_aligned_bai = AlignToMt.mt_aligned_bai
    File mt_aligned_shifted_bam = AlignToShiftedMt.mt_aligned_bam
    File mt_aligned_shifted_bai = AlignToShiftedMt.mt_aligned_bai
    File out_vcf = FilterContamination.filtered_vcf
    File out_vcf_index = FilterContamination.filtered_vcf_idx
    File input_vcf_for_haplochecker = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker
    File duplicate_metrics = AlignToMt.duplicate_metrics
    File contamination_metrics = GetContamination.contamination_file
    String major_haplogroup = GetContamination.major_hg
    Float contamination = FilterContamination.contamination
  }
}


task GetContamination {
  input {
    File input_vcf
  }

  meta {
    description: "Uses new Haplochecker to estimate levels of contamination in mitochondria"
  }
  parameter_meta {
    input_vcf: "Filtered and split multi-allelic sites VCF for mitochondria"
  }
  command <<<
  set -e
  PARENT_DIR="$(dirname "~{input_vcf}")"
  java -jar /path_to/haplocheckCLI.jar "${PARENT_DIR}"

  sed 's/\"//g' output >  output-noquotes

  grep "SampleID" output-noquotes > headers
  FORMAT_ERROR="Bad contamination file format"
  if [ `awk '{print $2}' headers` != "Contamination" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $6}' headers` != "HgMajor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $8}' headers` != "HgMinor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $14}' headers` != "MeanHetLevelMajor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $15}' headers` != "MeanHetLevelMinor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi

  grep -v "SampleID" output-noquotes > output-data
  awk -F "\t" '{print $2}' output-data > contamination.txt
  awk -F "\t" '{print $6}' output-data > major_hg.txt
  awk -F "\t" '{print $8}' output-data > minor_hg.txt
  awk -F "\t" '{print $14}' output-data > mean_het_major.txt
  awk -F "\t" '{print $15}' output-data > mean_het_minor.txt
  >>>
  output {
    File contamination_file = "output-noquotes"
    String hasContamination = read_string("contamination.txt") 
    String major_hg = read_string("major_hg.txt")
    String minor_hg = read_string("minor_hg.txt")
    Float major_level = read_float("mean_het_major.txt")
    Float minor_level = read_float("mean_het_minor.txt")
  }
}

task LiftoverAndCombineVcfs {
  input {
    File shifted_vcf
    File vcf
    String basename = basename(shifted_vcf, ".vcf")

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File shift_back_chain
  }

  meta {
    description: "Lifts over shifted vcf of control region and combines it with the rest of the chrM calls."
  }
  parameter_meta {
    shifted_vcf: "VCF of control region on shifted reference"
    vcf: "VCF of the rest of chrM on original reference"
    ref_fasta: "Original (not shifted) chrM reference"
    shift_back_chain: "Chain file to lift over from shifted reference to original chrM"
  }
  command<<<
    set -e

    java -XX:ActiveProcessorCount=1 -XX:ConcGCThreads=1 -jar /path_to/picard.jar LiftoverVcf \
      I=~{shifted_vcf} \
      O=~{basename}.shifted_back.vcf \
      R=~{ref_fasta} \
      CHAIN=~{shift_back_chain} \
      REJECT=~{basename}.rejected.vcf

    java -XX:ActiveProcessorCount=1 -XX:ConcGCThreads=1 -jar /path_to/picard.jar MergeVcfs \
      I=~{basename}.shifted_back.vcf \
      I=~{vcf} \
      O=~{basename}.merged.vcf
    >>>
    output{
        # rejected_vcf should always be empty
        File rejected_vcf = "~{basename}.rejected.vcf"
        File merged_vcf = "~{basename}.merged.vcf"
        File merged_vcf_index = "~{basename}.merged.vcf.idx"
    }
}

task M2 {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File input_bam
    File input_bai
    Int max_reads_per_alignment_start = 75
    String? m2_extra_args
    Boolean? make_bamout
    Boolean compress
    File? gatk_override
    Int? mem
  }

  String output_vcf = "raw" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 7000
  Int command_mem = machine_mem - 500

  meta {
    description: "Mutect2 for calling Snps and Indels"
  }
  parameter_meta {
    input_bam: "Aligned Bam"
  }
  command <<<
      set -e

      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

      gatk --java-options "-Xmx~{command_mem}m -XX:ActiveProcessorCount=1 -XX:ConcGCThreads=1" Mutect2 \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O ~{output_vcf} \
        ~{true='--bam-output bamout.bam' false='' make_bamout} \
        ~{m2_extra_args} \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --max-reads-per-alignment-start ~{max_reads_per_alignment_start} \
        --max-mnp-distance 0 \
        --dont-use-soft-clipped-bases
  >>>
  output {
      File raw_vcf = "~{output_vcf}"
      File raw_vcf_idx = "~{output_vcf_index}"
      File stats = "~{output_vcf}.stats"
      File output_bamOut = "bamout.bam"
  }
}

task Filter {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File raw_vcf
    File raw_vcf_index
    File raw_vcf_stats
    Boolean compress
    Float? vaf_cutoff
    String base_name

    String? m2_extra_filtering_args
    Int max_alt_allele_count
    Int? min_reads_per_strand
    Float? vaf_filter_threshold
    Float? f_score_beta

    Boolean run_contamination
    String? hasContamination
    Float? contamination_major
    Float? contamination_minor
    Float? verifyBamID
     
    File blacklisted_sites
    File blacklisted_sites_index

    File? gatk_override
  }

  String output_vcf = base_name + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  Float hc_contamination = if run_contamination && hasContamination == "YES" then (if contamination_major == 0.0 then contamination_minor else 1.0 - contamination_major) else 0.0
  Float max_contamination = if defined(verifyBamID) && verifyBamID > hc_contamination then verifyBamID else hc_contamination

  meta {
    description: "Mutect2 Filtering for calling Snps and Indels"
  }
  parameter_meta {
      vaf_filter_threshold: "Hard cutoff for minimum allele fraction. All sites with VAF less than this cutoff will be filtered."
      f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
  }
  command <<<
      set -e
      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

      gatk --java-options "-Xmx2500m -XX:ActiveProcessorCount=1 -XX:ConcGCThreads=1" FilterMutectCalls \
        -V ~{raw_vcf} \
        -R ~{ref_fasta} \
        -O filtered.vcf \
        --stats ~{raw_vcf_stats} \
        ~{m2_extra_filtering_args} \
        --max-alt-allele-count ~{max_alt_allele_count} \
        --mitochondria-mode \
        ~{"--min-reads-per-strand " + min_reads_per_strand} \
        ~{"--min-allele-fraction " + vaf_filter_threshold} \
        ~{"--f-score-beta " + f_score_beta} \
        ~{"--contamination-estimate " + max_contamination}

      gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ConcGCThreads=1" VariantFiltration \
        -V filtered.vcf \
        -O ~{output_vcf} \
        --apply-allele-specific-filters \
        --mask ~{blacklisted_sites} \
        --mask-name "blacklisted_site"

  >>>
  output {
      File filtered_vcf = "~{output_vcf}"
      File filtered_vcf_idx = "~{output_vcf_index}"
      Float contamination = "~{hc_contamination}"
  }
}

task MergeStats {
  input {
    File shifted_stats
    File non_shifted_stats
    File? gatk_override
  }

  command{
    set -e
    gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ConcGCThreads=1" MergeMutectStats \
      --stats ~{shifted_stats} \
      --stats ~{non_shifted_stats} \
      -O raw.combined.stats
  }
  output {
    File stats = "raw.combined.stats"
  }
}

task SplitMultiAllelicsAndRemoveNonPassSites {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File filtered_vcf
    File? gatk_override
  }

  command {
    set -e
    gatk --java-options "-XX:ActiveProcessorCount=1 -XX:ConcGCThreads=1" LeftAlignAndTrimVariants \
      -R ~{ref_fasta} \
      -V ~{filtered_vcf} \
      -O split.vcf \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac

      gatk SelectVariants \
        -V split.vcf \
        -O splitAndPassOnly.vcf \
        --exclude-filtered
  
  }
  output {
    File vcf_for_haplochecker = "splitAndPassOnly.vcf"
  }
}