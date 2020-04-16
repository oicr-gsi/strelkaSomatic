version 1.0

workflow strelkaSomatic {

    input {
	File tumorBam
	File normalBam
	File refFasta
	File? regionsBed
	String tmpDir = "/scratch2/users/ibancarz/strelka/workdir"
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
	tmpDir = tmpDir,
	outputFileNamePrefix = outputFileNamePrefix
    }

    call run {
	input:
	script = configure.script,
	tmpDir = tmpDir,
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

task configure {

    input {
	File tumorBam
	File normalBam
	File refFasta
	File? regionsBed
	String tmpDir
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

	rm -rf ~{tmpDir}/*

	configureStrelkaSomaticWorkflow.py \
	--normalBam ~{normalBam} \
	--tumorBam ~{tumorBam} \
	--referenceFasta ~{refFasta} \
	~{regionsBedArg} \
	--runDir ~{tmpDir}
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File script = "~{tmpDir}/runWorkflow.py"
    }
}

task run {

    input {
	File script
	String outputFileNamePrefix
	String tmpDir
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
	File snvsVcf = "~{tmpDir}/results/variants/somatic.snvs.vcf.gz"
	File indelsVcf = "~{tmpDir}/results/variants/somatic.indels.vcf.gz"
    }
    
}
