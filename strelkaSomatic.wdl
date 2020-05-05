version 1.0

workflow strelkaSomatic {

    input {
	File tumorBam
	File tumorBai
	File normalBam
	File normalBai
	File refFasta
	File refIndex
	File refDict
	String refModule
	File? bedFile
	Int? numChunk
	String outputFileNamePrefix = "strelkaSomatic"
    }
    
    parameter_meta {
	tumorBam: "Input BAM file with tumor data"
	tumorBai: "BAM index file for tumor data"
	normalBam: "Input BAM file with normal data"
	normalBai: "BAM index file for normal data"
	refFasta: "Reference FASTA file"
	refIndex: "Reference FAI index"
	refModule: "Name of genome reference environment module"
	refDict: "Reference DICT file"
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
	    name: "gatk/4.1.6.0",
	    url: "https://software.broadinstitute.org/gatk/download/index"
	}
	]
    }

    # Interval file provided, perform scatter/gather
    if(defined(bedFile)) {
	call splitIntervals {
	    input:
	    refFasta = refFasta,
	    refFai = refIndex,
	    refDict = refDict,
	    refModule = refModule,
            intervals = bedFile,
            scatterCount = numChunk
	}
	Array[File] intervals = splitIntervals.intervalFiles
	scatter(interval in intervals) {
	    call configureAndRun as configureAndRunParallel {
		input:
		tumorBam = tumorBam,
		tumorBai = tumorBai,
		normalBam = normalBam,
		normalBai = normalBai,
		refFasta = refFasta,
		refIndex = refIndex,
		refModule = refModule,
		regionsBed = interval,
		outputFileNamePrefix = outputFileNamePrefix
	    }
	}
	call vcfGather as snvsVcfGather {
	    input:
	    vcfs = configureAndRunParallel.snvsVcf,
	}
	call vcfGather as indelsVcfGather {
	    input:
	    vcfs = configureAndRunParallel.indelsVcf,
	}
    }
    # No interval file, run as single process
    if(!defined(bedFile)) {
	call configureAndRun as configureAndRunSingle {
	    input:
	    tumorBam = tumorBam,
	    tumorBai = tumorBai,
	    normalBam = normalBam,
	    normalBai = normalBai,
	    refFasta = refFasta,
	    refIndex = refIndex,
	    refModule = refModule,
	    outputFileNamePrefix = outputFileNamePrefix
	}
    }

    # Do not output the TSV and XML stats files; contents not needed
    output {
	File snvsVcf = select_first([snvsVcfGather.result, configureAndRunSingle.snvsVcf])
	File indelsVcf = select_first([indelsVcfGather.result, configureAndRunSingle.indelsVcf])
    }
}

task configureAndRun {

    input {
	File tumorBam
	File tumorBai
	File normalBam
	File normalBai
	File refFasta
	File refIndex
	File refModule
	File? regionsBed
	String outputFileNamePrefix
	String nonRefModules = "python/2.7 samtools/1.9 strelka/2.9.10"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	tumorBam: "BAM file with aligned reads from tumor sample"
	tumorBai: "BAM index file for tumor data"
	normalBam: "BAM file with aligned reads from normal sample"
	normalBai: "BAM index file for normal data"
	refFasta: "FASTA reference file"
	refIndex: "FAI reference index file"
	refModule: "Genome reference module name"
	regionsBed: "BED file designating regions to process"
	nonRefModules: "Environment module names other than genome reference"
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

    String modules = "~{refModule} ~{nonRefModules}"
    String bedName = "regions.bed.gz"
    # remove header from the .bed file and compress, for compatibility with Strelka
    String writeBed = if defined(regionsBed) then "grep -v \"^@\" ~{regionsBed} | bgzip > ~{bedName}" else ""
    String indexFeatureFile = if defined(regionsBed) then "tabix ~{bedName}" else ""
    String regionsBedArg = if defined(regionsBed) then "--callRegions ~{bedName}" else ""

    # index files not explicitly used, but needed by configureStrelkaSomaticWorkflow.py
    # ie. tumorBai, normalBai, refIndex; tabix output on bedfile (if any)

    command <<<
	set -eo pipefail

	~{writeBed}
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
	String refFasta
	String refFai
	String refDict
	String refModule
	String gatk = "$GATK_ROOT/bin/gatk"
	File? intervals
	Int? scatterCount
	String? splitIntervalsExtraArgs
	String nonRefModules = "gatk/4.1.2.0"
	Int memory = 32
	Int timeout = 72
    }

    parameter_meta {
	refFasta: "Path to the reference fasta"
	refFai: "Path to the reference .fai index"
	refDict: "Path to the reference .dict dictionary"
	refModule: "Genome reference module"
	gatk: "GATK executable path"
	intervals: "Interval file to split for scattering"
	scatterCount: "Number of files to split the interval file into"
	splitIntervalsExtraArgs: "Additional arguments for the 'gatk SplitIntervals' command"
	memory: "Memory allocated for job"
	timeout: "Hours before task timeout"
	nonRefModules: "Environment modules other than the genome refence"
    }

    meta {
	output_meta: {
	    intervals: "Interval files that were split to be used for scattering."
	}
    }

    String modules = "~{refModule} ~{nonRefModules}"
    String intervalsArg = if defined(intervals) then "-L ${intervals }" else ""
    String scatterArg = if defined(scatterCount) then "--scatter-count ~{scatterCount}" else ""

    command <<<
	set -eo pipefail

	mkdir interval-files
	ln -s ~{refFai}
	ln -s ~{refDict}
	~{gatk} --java-options "-Xmx~{memory-8}g" SplitIntervals \
	-R ~{refFasta} \
	~{intervalsArg} \
	~{scatterArg} \
	--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
	-O interval-files \
	~{splitIntervalsExtraArgs}

	cp interval-files/*.interval_list .
    >>>

    runtime {
	memory:  "~{memory} GB"
	modules: "~{modules}"
	timeout: "~{timeout}"
    }

    output {
	Array[File] intervalFiles = glob("*.interval_list")
    }
}

task vcfGather {

    input {
	String modules = "gatk/4.1.2.0"
	String gatk = "$GATK_ROOT/bin/gatk"
	Array[File] vcfs
	Int memory = 16
	Int timeout = 12
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

    String outputPrefix = basename(vcfs[0], ".vcf.gz")
    String outputName = "~{outputPrefix}.vcf"

    command <<<
	set -eo pipefail

	~{gatk} GatherVcfs \
	-I ~{sep=" -I " vcfs} \
	-O ~{outputName}

	gzip ~{outputName}
    >>>

    runtime {
	memory:  "~{memory} GB"
	modules: "~{modules}"
	timeout: "~{timeout}"
    }

    output {
	File result = "~{outputName}.gz"
    }
}
