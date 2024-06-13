version 1.0

workflow GenerateConsensus {
  meta {
    description: "Generate personal reference genome for realignment and recalling."
  }

  input {
    File vcf_for_haplochecker
    File mt_fasta
  }

  call GenerateFasta {
    input:
      input_vcf = vcf_for_haplochecker,
      ref_fasta = mt_fasta
  }

  call GenerateShiftedFasta {
    input:
      consensus_fasta = GenerateFasta.consensus_fasta
  }

  output {
    File consensus_fasta = GenerateFasta.consensus_fasta
    File consensus_fasta_index = GenerateFasta.consensus_fasta_index
    File consensus_dict = GenerateFasta.consensus_dict
    File consensus_amb = GenerateFasta.consensus_amb
    File consensus_ann = GenerateFasta.consensus_ann
    File consensus_bwt = GenerateFasta.consensus_bwt
    File consensus_pac = GenerateFasta.consensus_pac
    File consensus_sa = GenerateFasta.consensus_sa
    File consensus_shifted_fasta = GenerateShiftedFasta.consensus_shifted_fasta
    File consensus_shifted_fasta_index = GenerateShiftedFasta.consensus_shifted_fasta_index
    File consensus_shifted_dict = GenerateShiftedFasta.consensus_shifted_dict
    File consensus_shifted_amb = GenerateShiftedFasta.consensus_shifted_amb
    File consensus_shifted_ann = GenerateShiftedFasta.consensus_shifted_ann
    File consensus_shifted_bwt = GenerateShiftedFasta.consensus_shifted_bwt
    File consensus_shifted_pac = GenerateShiftedFasta.consensus_shifted_pac
    File consensus_shifted_sa = GenerateShiftedFasta.consensus_shifted_sa
  }
}

task GenerateFasta {
  input {
    File input_vcf
    File ref_fasta
  }

  meta {
    description: "Use bcftools to generate consensus fasta"
  }

  command <<<
    set -e

    bgzip -c ~{input_vcf} > PASS.vcf.gz
    tabix -p vcf PASS.vcf.gz
    bcftools consensus \
      -i 'TYPE="snp" && AD[0:0]<AD[0:1]' \
      -f ~{ref_fasta} \
      PASS.vcf.gz \
      > consensus.fasta
    
    samtools faidx consensus.fasta
    samtools dict consensus.fasta -o consensus.dict
    bwa index consensus.fasta
  >>>
  output {
    File consensus_fasta = "consensus.fasta"
    File consensus_fasta_index = "consensus.fasta.fai"
    File consensus_dict = "consensus.dict"
    File consensus_amb = "consensus.fasta.amb"
    File consensus_ann = "consensus.fasta.ann"
    File consensus_bwt = "consensus.fasta.bwt"
    File consensus_pac = "consensus.fasta.pac"
    File consensus_sa = "consensus.fasta.sa"
  }

}

task GenerateShiftedFasta {
  input {
    File consensus_fasta
  }

  meta {
    description: "Shift consensus fasta by 8000 bases"
  }

  command <<<
    set -e 

    python <<CODE
    with open("~{consensus_fasta}") as f:
        fasta_in = f.readlines()

    header = fasta_in[0:1]
    old_seq = fasta_in[1:]
    old_seq = "".join(old_seq).replace("\n","")
    new_seq = old_seq[8000:] + old_seq[:8000]
    new_seq = [new_seq[i:i+100] + "\n" for i in range(0, len(new_seq), 100)]
    new_seq = header + new_seq

    with open('consensus_shifted.fasta','w') as f2:
        for items in new_seq:
            f2.write(items)
    CODE

    samtools faidx consensus_shifted.fasta
    samtools dict consensus_shifted.fasta -o consensus_shifted.dict
    bwa index consensus_shifted.fasta
  >>>
  output {
    File consensus_shifted_fasta = "consensus_shifted.fasta"
    File consensus_shifted_fasta_index = "consensus_shifted.fasta.fai"
    File consensus_shifted_dict = "consensus_shifted.dict"
    File consensus_shifted_amb = "consensus_shifted.fasta.amb"
    File consensus_shifted_ann = "consensus_shifted.fasta.ann"
    File consensus_shifted_bwt = "consensus_shifted.fasta.bwt"
    File consensus_shifted_pac = "consensus_shifted.fasta.pac"
    File consensus_shifted_sa = "consensus_shifted.fasta.sa"
  }
}