version 1.0

workflow strelkaGermline {

    input {
	Array[File] bamInput
	File fastaRef
	String outputFileNamePrefix = "strelkaGermline"
    }
    
    parameter_meta {
	bamInput: "Array of 2 or more input BAM files"
	fastaRef: "Reference FASTA file"
	outputFileNamePrefix: "Prefix for output files"
    }

    call configure {
	input:
	bamInput = bamInput,
	fastaRef = fastaRef,
	outputFileNamePrefix = outputFileNamePrefix
    }

    call run {
	input:
	script = configure.script,
	outputFileNamePrefix = outputFileNamePrefix
    }

    output {
	File variantsVcf = run.variantsVcf
	Array[File] genomeVcf = run.genomeVcf
    }

    meta {
	author: "Iain Bancarz"
	email: "ibancarz@oicr.on.ca"
	description: "Strelka variant caller in germline mode"
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
	Array[File] bamInput,
	File fastaRef,
	String outputFileNamePrefix,
	String modules = "python/2.7 strelka/2.9.10"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    command <<<
	
    >>>

    
}

