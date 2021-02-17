#!/usr/bin/env nextflow

process lenstools_make_peptides {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'

  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(vcf)
  path tx_aa

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.mt_aa.fa"), emit: mutant_fastas
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.wt_aa.fa"), emit: wildtype_fastas
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.wt_aa.fa"), path("*.mt_aa.fa"), emit: peptide_fastas

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py make-peptides -v ${vcf} -t ${tx_aa} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.snvs.mt_aa.fa -w ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.snvs.wt_aa.fa
  """
}


process lenstools_make_indel_peptides {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  cache false

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

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy'

  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(vcf), path(quant)
  val parstr

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.vcf"), emit: expressed_vcfs

  script:
  """
  AVCF=`echo ${vcf}`
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py expressed-variants ${parstr} -v ${vcf} -a ${quant} -o \${AVCF%.vcf*}.expfilt.vcf --abundance-threshold 32.9
  #python ${params.project_dir}/workflow/lenstools/bin/lenstools.py expressed-variants ${parstr} -v ${vcf} -a ${quant} -o \${AVCF%.vcf*}.expfilt.vcf
  sed -i 's/^#\\t/1\\t/g' \${AVCF%.vcf*}.expfilt.vcf
  """
}


process lenstools_filter_isolated_variants {

  label "lenstools"
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'

  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(somatic_vcf), path(germline_vcf)
  val parstr

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.vcf"), emit: isolated_vcfs

  script:
  """
  AVCF=`echo ${somatic_vcf}`
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py isolated-variants -s ${somatic_vcf} -g ${germline_vcf}  -o \${AVCF%.vcf*}.isofilt.vcf
  #Post-processing to deal with weird pyvcf bug where first entry is turned into header.
  sed -i 's/^#hr/#CHROM\tPOS\tID\tREF\tALT\\\nchr/'  \${AVCF%.vcf*}.isofilt.vcf
  """
}


process lenstools_calculate_agretopicity {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'

  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(wt_nmp), path(mt_nmp)

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*agreto.netmhcpan.txt"), emit: agreto_netmhcpans

  script:
  """
  MTNMP=`echo ${mt_nmp}`
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py calculate-agretopicity -w ${wt_nmp} -m ${mt_nmp} -o  \${MTNMP%.netmhcpan.txt}.agreto.netmhcpan.txt
  """
}


process lenstools_filter_peptides {

  label "lenstools"
  
  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(binding_affinities)

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*filt_peps*"), emit: filtered_peptides

  script:
  """
  BA_TMP=`echo ${binding_affinities}`
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py filter-peptides -i ${binding_affinities} -o  \${BA_TMP%.netmhcpan.txt}.agreto.netmhcpan.txt
  """
}

process lenstools_rna_covered_variants {

  label "lenstools"
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  cache false
  
  input:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(isec_vcf), path(rna_bam), path(rna_bai)

  output:
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*rna_cov.vcf*"), emit: rna_cov_vcfs

  script:
  """
  AVCF=\$(echo ${isec_vcf})
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py rna-covered-variants -v ${isec_vcf} -b ${rna_bam} -o  \${AVCF%.vcf}.rna_cov.vcf
  """
}

process lenstools_make_pyclonevi_inputs {                                                           
                                                                                                    
  label "lenstools"                                                                                 
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
                                                                                                    
  input:                                                                                            
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path(candidate_vcf), path(mutect_vcf), path(sequenza_segments)
  val parstr                                                                                        
                                                                                                    
  output:                                                                                           
  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*pcvi_input"), emit: pcvi_inputs
                                                                                                    
  script:                                                                                           
  """                                                                                               
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py make-pyclone-vi-inputs -c ${candidate_vcf} -m ${mutect_vcf} -s ${sequenza_segments} --samp-id ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.pcvi_input
  """                                                                                               
} 
