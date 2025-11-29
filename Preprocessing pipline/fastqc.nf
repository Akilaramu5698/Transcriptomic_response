params.fastq="/home/bic/akila/fastq/meropenem-gse190441/9hrs/treated/*.fastq"
params.qc_report="/home/bic/akila/fastq/meropenem-gse190441/9hrs/treated/fastqc_report"
process QC {

publishDir("${params.qc_report}", mode:'copy')

input:
 path fastq

output:
 path "*"

script:
"""

fastqc $fastq

"""

} 

workflow {
fastq_ch=Channel.fromPath(params.fastq)
QC(fastq_ch)
QC.out.view()

}
