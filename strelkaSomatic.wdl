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

    call configure {
	input:
	tumorBam = tumorBam,
	normalBam = normalBam,
	refFasta = refFasta,
	regionsBed = regionsBed,
	outputFileNamePrefix = outputFileNamePrefix
    }

    call run {
	input:
	script = configure.script,
	outputFileNamePrefix = outputFileNamePrefix
    }

    output {
	File snvsVcf = run.snvsVcf
	File indelsVcf = run.indelsVcf
    }

    meta {
	author: "Iain Bancarz"
	email: "ibancarz@oicr.on.ca"
	description: "Strelka variant caller in somatic mode"
	dependencies: [
	{
	    name: "strelka/2.9.10",
	    url: "https://github.com/samtools/samtools"
	},
	{
	    name: "python/2.7",
	    url: "https://www.python.org/downloads/release/python-2716/"
	}  
	]
    }

}

task configure {

    input {
	File tumorBam
	File normalBam
	File refFasta
	File? regionsBed
	String outputFileNamePrefix
	String modules = "python/2.7 strelka/2.9.10"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    String regionsBedArg = if defined(regionsBed) then "--callRegions ~{regionsBed}" else ""
    
    command <<<
	configureStrelkaSomaticWorkflow.py \
	--normalBam ~{normalBam} \
	--tumorBam ~{tumorBam} \
	--referenceFasta ~{refFasta} \
	~{regionsBedArg} \
	--runDir .
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File script = "runWorkflow.py"
    }
}

task run {

    input {
	File script
	String outputFileNamePrefix
	String modules = "python/2.7 strelka/2.9.10"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    # TODO will the "-m sge" option work in conjunction with Cromwell?
    
    command <<<
	~{script} -m local -j ~{threads}
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File snvsVcf = "somatic.snvs.vcf.gz"
	File indelsVcf = "somatic.indels.vcf.gz"
    }
    
}
