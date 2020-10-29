#!/usr/bin/env nextflow

process lenstools_make_peptides {

  conda 'bioconda::pyvcf bioconda::biopython'

  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(vcf)
  path tx_aa

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.aa.fa"), emit: peptide_fastas

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py make-peptides -v ${vcf} -t ${tx_aa} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.aa.fa
  """
}

process lenstools_filter_expressed_variants {

  conda 'bioconda::pyvcf bioconda::biopython'

  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(vcf), path(quant)
  val parstr

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.vcf"), emit: peptide_fastas

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py expressed-variants -v ${vcf} -a ${quant} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.annot.expfilt.vcf
  """
}
