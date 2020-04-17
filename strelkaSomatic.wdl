version 1.0

workflow strelkaSomatic {

    input {
	File tumorBam
	File normalBam
	File refFasta
	File? regionsBed
	String outputFileNamePrefix = "strelkaSomatic"
    }
    
    parameter_meta {
	tumorBam: "Input BAM file with tumor data"
	normalBam: "Input BAM file with normal data"
	refFasta: "Reference FASTA file"
	regionsBed: "BED file with call regions"
	outputFileNamePrefix: "Prefix for output files"
    }

    call configureAndRun {
	input:
	tumorBam = tumorBam,
	normalBam = normalBam,
	refFasta = refFasta,
	regionsBed = regionsBed,
	outputFileNamePrefix = outputFileNamePrefix
    }

    output {
	File snvsVcf = configureAndRun.snvsVcf
	File indelsVcf = configureAndRun.indelsVcf
    }

    meta {
	author: "Iain Bancarz"
	email: "ibancarz@oicr.on.ca"
	description: "Strelka variant caller in somatic mode"
	dependencies: [
	{
	    name: "samtools/1.9",
	    url: "https://github.com/samtools/samtools"
	},
	{
	    name: "strelka/2.9.10",
	    url: "https://github.com/Illumina/strelka/releases/tag/v2.9.10"
	},
	{
	    name: "python/2.7",
	    url: "https://www.python.org/downloads/release/python-2716/"
	}  
	]
    }

}

task configureAndRun {

    input {
	File tumorBam
	File normalBam
	File refFasta
	File? regionsBed
	String outputFileNamePrefix
	String modules = "python/2.7 samtools/1.9 strelka/2.9.10"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    String regionsBedArg = if defined(regionsBed) then "--callRegions ~{regionsBed}" else ""
    
    command <<<
	set -eo pipefail

	samtools faidx ~{refFasta}
	samtools index ~{normalBam}
	samtools index ~{tumorBam}

	configureStrelkaSomaticWorkflow.py \
	--normalBam ~{normalBam} \
	--tumorBam ~{tumorBam} \
	--referenceFasta ~{refFasta} \
	~{regionsBedArg} \
	--runDir .

	./runWorkflow.py -m local -j ~{threads}
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File snvsVcf = "./results/variants/somatic.snvs.vcf.gz"
	File indelsVcf = "./results/variants/somatic.indels.vcf.gz"
    }
}

