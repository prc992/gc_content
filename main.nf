 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

process split_tss_samples{

  tag "Sample - $sampleId - $cores" 

  //Docker Image
  container = 'prc992/gc_calc_par:v1.1'
  label 'default_mem'
  cpus {cores}

  publishDir "$params.output/$sampleId", mode : 'copy', pattern : '*.csv'

  input:
  tuple val(sampleId),path(faFile),path(bedFileIn),val(cores)
  
  output:
  path("*.csv")

  script:
  """
  python $params.python_prog -i $bedFileIn -o $params.output -c $cores
  """
}

workflow {

    chSampleInfo = Channel.fromPath(params.samples) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId,row.faFile,row.bedFileIn,row.cores) }

    split_tss_samples(chSampleInfo)
}

