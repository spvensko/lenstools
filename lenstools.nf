#!/usr/bin/env nextflow

process lenstools_make_peptides {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}"
  cache false

  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path(vcf)
  path tx_aa

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("*.mt_aa.fa"), emit: mutant_fastas
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("*.wt_aa.fa"), emit: wildtype_fastas
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("*.wt_aa.fa"), path("*.mt_aa.fa"), emit: peptide_fastas
//  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.mt_aa.fa"), emit: mutant_fastas
//  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.wt_aa.fa"), emit: wildtype_fastas
//  tuple val(pat_name), val(dataset), val(norm_prefix), val(tumor_prefix), path("*.wt_aa.fa"), path("*.mt_aa.fa"), emit: peptide_fastas

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py make-peptides -v ${vcf} -t ${tx_aa} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.snvs.mt_aa.fa -w ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.snvs.wt_aa.fa
  """
}


process lenstools_make_indel_peptides {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}"

  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path(vcf)
  path tx_aa

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("*.aa.fa"), emit: peptide_fastas

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py make-indel-peptides -v ${vcf} -t ${tx_aa} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.indels.aa.fa
  """
}


process lenstools_filter_expressed_variants {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}"

  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path(vcf), path(quant)
  val parstr

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("*.vcf"), emit: expressed_vcfs

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
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}"

  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path(somatic_vcf), path(germline_vcf)
  val parstr

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("*.vcf"), emit: isolated_vcfs

  script:
  """
  AVCF=`echo ${somatic_vcf}`
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py isolated-variants -s ${somatic_vcf} -g ${germline_vcf}  --proximity 42 -o \${AVCF%.vcf*}.isofilt.vcf
  #Post-processing to deal with weird pyvcf bug where first entry is turned into header.
  sed -i 's/^#hr/#CHROM\tPOS\tID\tREF\tALT\\\nchr/'  \${AVCF%.vcf*}.isofilt.vcf
  """
}


process lenstools_calculate_agretopicity {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}"
  cache false

  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path(wt_nmp), path(mt_nmp)

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("*agreto.netmhcpan.txt"), emit: agreto_netmhcpans

  script:
  """
  MTNMP=`echo ${mt_nmp}`
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py calculate-agretopicity -w ${wt_nmp} -m ${mt_nmp} -o  \${MTNMP%.netmhcpan.txt}.agreto.netmhcpan.txt
  """
}


process lenstools_filter_peptides {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy'
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}"

  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path(binding_affinities)

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("*filt_peps*"), emit: filtered_peptides

  script:
  """
  BA_TMP=`echo ${binding_affinities}`
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py filter-peptides -i ${binding_affinities} -o  \${BA_TMP%.netmhcpan.txt}.agreto.netmhcpan.txt
  """
}

process lenstools_rna_covered_variants {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}"

  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path(isec_vcf), path(rna_bam), path(rna_bai)

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("*rna_cov.vcf*"), emit: rna_cov_vcfs

  script:
  """
  AVCF=\$(echo ${isec_vcf})
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py rna-covered-variants -v ${isec_vcf} -b ${rna_bam} -o  \${AVCF%.vcf}.rna_cov.vcf
  """
}

process lenstools_make_pyclonevi_inputs {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}"
  cache false

  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path(candidate_vcf), path(mutect_vcf), path(sequenza_segments), path(sequenza_solutions)
  val parstr

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("*pcvi_input"), emit: pcvi_inputs

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py make-pyclone-vi-inputs -c ${candidate_vcf} -m ${mutect_vcf} -s ${sequenza_segments} --samp-id ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix} --sequenza-solutions ${sequenza_solutions} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.pcvi_input
  """
}

process lenstools_add_snv_metadata {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}_${rna_prefix}"
  cache false

  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(rna_prefix), val(dataset), path(netmhcpan_input), path(mutant_fasta), path(pcvi), path(mutant_nucs), path(quants)
  path gtf
  val parstr

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(rna_prefix), val(dataset), path("*metadata.txt"), emit: metadatas

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py add-snv-metadata -n ${netmhcpan_input} -m ${mutant_fasta} -q ${quants} -u ${mutant_nucs} -g ${gtf} -c ${pcvi}  -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}_${rna_prefix}.snv.metadata.txt
  """
}

process lenstools_add_indel_metadata {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}_${rna_prefix}"
  cache false

  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(rna_prefix), val(dataset), path(netmhcpan_input), path(mutant_fasta), path(pcvi), path(quants)
  path gtf
  val parstr

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(rna_prefix), val(dataset), path("*metadata.txt"), emit: metadatas

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py add-indel-metadata -n ${netmhcpan_input} -m ${mutant_fasta} -q ${quants} -g ${gtf} -c ${pcvi}  -o "${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}_${rna_prefix}.indel.metadata.txt"
  """
}

process lenstools_make_genomic_context {
  
  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${norm_prefix}_${tumor_prefix}"
  cache false
  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path(vcf)
  path tx_cds

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("*.nuc.fa"), emit: mutant_nucs

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py make-genomic-context -v ${vcf} -c ${tx_cds} -o ${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}.snvs.nuc.fa
  """
}

process lenstools_filter_expressed_hervs {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${prefix}"
  cache false

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(quant)
  val parstr

  output:
  tuple val(pat_name), val(prefix), val(dataset), path('*.expressed_hervs.txt'), emit: expressed_hervs

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py expressed-hervs ${parstr} -a ${quant} -o ${dataset}-${pat_name}-${prefix}.expressed_hervs.txt
  """
}

process lenstools_make_herv_peptides {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${prefix}"
  cache false

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(expressed_hervs)
  path herv_ref
  val parstr

  output:
  tuple val(pat_name), val(prefix), val(dataset), path('*.hervs.peptides.fa'), emit: herv_peptides

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py make-herv-peptides ${parstr} -e ${expressed_hervs} -r ${herv_ref} -o ${dataset}-${pat_name}-${prefix}.hervs.peptides.fa
  """
}

process lenstools_get_herv_metadata {

  label "lenstools"
  conda 'bioconda::pyvcf bioconda::biopython anaconda::numpy anaconda::scipy bioconda::pysam'
  tag "${dataset}/${pat_name}/${prefix}"
  cache false

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(netmhcpan), path(abundances)
  val parstr

  output:
  tuple val(pat_name), val(prefix), val(dataset), path('*.hervs.metadata.txt'), emit: herv_metadata

  script:
  """
  python ${params.project_dir}/workflow/lenstools/bin/lenstools.py get-herv-metadata ${parstr} -a ${abundances} -n ${netmhcpan} -o ${dataset}-${pat_name}-${prefix}.hervs.metadata.txt
  """
}
