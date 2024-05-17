 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

process split_tss_samples{

  //Docker Image
  container = 'prc992/gc_calc_par:v1.0'
  label 'default_mem'

  publishDir "$params.output/$sampleId", mode : 'copy', pattern : '*.csv'

  input:
  tuple val(sampleId),path(faFile),path(bedFileIn)
  
  output:
  path("*.csv")

  script:
  """
  python $params.python_prog -i $bedFileIn -o $params.output -c $params.cores
  """
}

workflow {

    chSampleInfo = Channel.fromPath(params.samples) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId,row.faFile,row.bedFileIn) }

    split_tss_samples(chSampleInfo)
}

