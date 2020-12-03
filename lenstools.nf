#!/usr/bin/env nextflow

process lenstools_make_peptides {

  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy'

  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(vcf)
  path tx_aa

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.aa.fa"), emit: peptide_fastas

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py make-peptides -v ${vcf} -t ${tx_aa} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.snvs.aa.fa
  """
}

process lenstools_make_indel_peptides {

  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy'

  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(vcf)
  path tx_aa

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.aa.fa"), emit: peptide_fastas

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py make-indel-peptides -v ${vcf} -t ${tx_aa} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.indels.aa.fa
  """
}

process lenstools_filter_expressed_variants {

  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy'

  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(vcf), path(quant)
  val parstr

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.vcf"), emit: expressed_vcf

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py expressed-variants -v ${vcf} -a ${quant} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.annot.expfilt.vcf
  """
}

process lenstools_filter_isolated_variants {

  conda 'bioconda::pyvcf bioconda::biopython'

  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(cand_vcf), path(somatic_vcf), path(germline_vcf)
  val parstr

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.vcf"), emit: isolated_vcfs

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py expressed-variants -c ${cand_vcf} -s ${somatic_vcf} -g ${germlline_vcf}  #-o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.annot.isofilt.vcf
  """
}
