version 1.0

workflow strelkaSomatic {

    input {
	File tumorBam
	File normalBam
	File refFasta
	File? bedFile
	Int? numChunk
	String outputFileNamePrefix = "strelkaSomatic"
    }
    
    parameter_meta {
	tumorBam: "Input BAM file with tumor data"
	normalBam: "Input BAM file with normal data"
	refFasta: "Reference FASTA file"
	bedFile: "BED file designating regions to process"
	numChunk: "If BED file given, number of chunks in which to split each chromosome"
	outputFileNamePrefix: "Prefix for output files"
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
	},
	{
	    name: "gatk/4.1.1.0",
	    url: "https://software.broadinstitute.org/gatk/download/index"
	}
	]
    }

    # Interval file provided, perform scatter/gather
    if(defined(bedFile)) {

	call splitIntervals {
	    input:
            intervals = bedFile,
            scatterCount = numChunk
	}

	Array[File] intervals = splitIntervals.intervalFiles

	scatter(interval in intervals) {

	    call configureAndRun as configureAndRunParallel {
		input:
		tumorBam = tumorBam,
		normalBam = normalBam,
		refFasta = refFasta,
		regionsBed = interval,
		outputFileNamePrefix = outputFileNamePrefix
	    }
	}

	# TODO should we concatenate and output the TSV and XML stats files?

	call vcfGather as snvsVcfGather {
	    input:
	    vcfs = configureAndRunParallel.snvsVcf,
	    variantType = "snvs"
	}
	call vcfGather as indelsVcfGather {
	    input:
	    vcfs = configureAndRunParallel.indelsVcf,
	    variantType = "indels"
	}
    }

    if(!defined(bedFile)) {

	call configureAndRun as configureAndRunSingle {
	    input:
	    tumorBam = tumorBam,
	    normalBam = normalBam,
	    refFasta = refFasta,
	    outputFileNamePrefix = outputFileNamePrefix
	}
    }

    output {
	File snvsVcf = select_first([snvsVcfGather.result, configureAndRunSingle.snvsVcf])
	File indelsVcf = select_first([indelsVcfGather.result, configureAndRunSingle.indelsVcf])
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

    parameter_meta {
	tumorBam: "BAM file with aligned reads from tumor sample"
	normalBam: "BAM file with aligned reads from normal sample"
	refFasta: "FASTA reference file"
	regionsBed: "BED file designating regions to process"
	modules: "Environment module names and version to load (space separated) before command execution"
	jobMemory: "Memory allocated for job"
	timeout: "Hours before task timeout"
	threads: "Number of threads for processing"
    }

    meta {
	output_meta: {
	    snvsVcf: "VCF file with SNVs, .gz compressed",
	    indelsVcf: "VCF file with indels, .gz compressed"
	}
    }

    String bedGrep = if defined(regionsBed) then "grep -v \"^@\" ~{regionsBed} | bgzip > regions.bed.gz" else ""
    #String regionsBedArg = if defined(regionsBed) then "--callRegions regions.bed" else ""
    String regionsBedArg = if defined(regionsBed) then "--callRegions regions.bed.gz" else ""
    #String indexFeatureFile = if defined(regionsBed) then "gatk --java-options \"-Xmx~{jobMemory}g\" IndexFeatureFile -F regions.bed" else ""
    String indexFeatureFile = if defined(regionsBed) then "tabix regions.bed.gz" else ""
    
    command <<<
	set -eo pipefail

	samtools faidx ~{refFasta}
	samtools index ~{normalBam}
	samtools index ~{tumorBam}
	~{bedGrep}
	~{indexFeatureFile}

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

task splitIntervals {

    input {
	String modules = "gatk/4.1.2.0 hg19/p13"
	String refFasta = "$HG19_ROOT/hg19_random.fa"
	String refFai = "$HG19_ROOT/hg19_random.fa.fai"
	String refDict = "$HG19_ROOT/hg19_random.dict"
	File? intervals
	Int? scatterCount
	String? splitIntervalsExtraArgs
	Int memory = 16
	Int timeout = 72
    }

    parameter_meta {
	modules: "Environment module names and version to load (space separated) before command execution"
	refFasta: "Path to the reference fasta"
	intervals: "Interval file to split for scattering"
	scatterCount: "Number of files to split the interval file into"
	memory: "Memory allocated for job"
	timeout: "Hours before task timeout"
    }

    meta {
	output_meta: {
	    intervals: "Interval files that were split to be used for scattering."
	}
    }

    String scatterArg = if defined(scatterCount) then "--scatter-count ~{scatterCount}" else ""

    command <<<
	set -eo pipefail

	mkdir interval-files
	ln -s ~{refFai}
	ln -s ~{refDict}
	gatk --java-options "-Xmx~{memory}g" SplitIntervals \
	-R ~{refFasta} \
	~{"-L " + intervals} \
	~{scatterArg} \
	-O interval-files \
	~{splitIntervalsExtraArgs}

	cp interval-files/*.interval_list .
    >>>

    output {
	Array[File] intervalFiles = glob("*.interval_list")
    }
}

task vcfGather {
    # TODO do not use hard-coded genome reference module name
    input {
	String modules = "gatk/4.1.2.0 hg19/p13"
	String gatk = "$GATK_ROOT/bin/gatk"
	String refFasta = "$HG19_ROOT/hg19_random.fa"
	Array[File] vcfs
	String variantType
	Int memory = 32
	Int timeout = 72
    }

    parameter_meta {
	modules: "Environment module names and version to load (space separated) before command execution"
	gatk: "GATK to use"
	vcfs: "VCFs from scatter to merge together"
	memory: "Memory allocated for job"
	timeout: "Hours before task timeout"
    }

    meta {
	output_meta: {
	    unfilteredVcf: "Merged vcf, unfiltered.",
	    unfilteredVcfIndex: "Merged vcf index, unfiltered.",
	    filteredVcf: "Merged vcf, processed through gatk FilterMutectCalls.",
	    filteredVcfIndex: "Merged vcf index, processed through gatk FilterMutectCalls."
	}
    }

    #String outputPrefix = basename(select_first([vcfs[0], ""]), ".vcf")
    String outputPrefix = basename(vcfs[0], ".vcf.gz")
    String outputName = "~{outputPrefix}.vcf"

    command <<<
	set -eo pipefail

	~{gatk} GatherVcfs \
	-I ~{sep=" -I " vcfs} \
	-O ~{outputName}
    >>>

    runtime {
	memory:  "~{memory} GB"
	modules: "~{modules}"
	timeout: "~{timeout}"
    }

    output {
	File result = "~{outputName}"
    }
}
