version 1.0

struct referenceResources {
  String refDict
  String refIndex 
  String indelsVcfGather_refIndex 
  String snvsVcfGather_refIndex 
  String refFasta 
  String refModule
  String bedFile
}

workflow strelkaSomatic {

   input {
      File tumorBam
      File tumorBai
      File normalBam
      File normalBai
      String reference
      String? bedFile
      Int? numChunk
      String outputFileNamePrefix = "strelkaSomatic"
    }

    # store genome references as String, not File
    # allows correct interpolation of environment variables
    # Eg. for "strelkaSomatic.refDict": "$HG38_ROOT/hg38_random.dict"
    Map [String,referenceResources] resources = {
      "hg19": {
        "refDict": "$HG19_ROOT/hg19_random.dict",
        "refIndex": "$HG19_ROOT/hg19_random.fa.fai",
        "indelsVcfGather_refIndex": "$HG19_ROOT/hg19_random.fa.fai",
        "snvsVcfGather_refIndex": "$HG19_ROOT/hg19_random.fa.fai",
        "refFasta": "$HG19_ROOT/hg19_random.fa",
        "refModule": "hg19/p13 samtools/1.9",
        "bedFile": "$HG19_ROOT/hg19.chrom.sizes.bed"
      },
      "hg38": {
        "refDict": "$HG38_ROOT/hg38_random.dict",
        "refIndex": "$HG38_ROOT/hg38_random.fa.fai",
        "indelsVcfGather_refIndex": "$HG38_ROOT/hg38_random.fa.fai",
        "snvsVcfGather_refIndex": "$HG38_ROOT/hg38_random.fa.fai",
        "refFasta": "$HG38_ROOT/hg38_random.fa",
        "refModule": "hg38/p12 samtools/1.9",
        "bedFile": "$HG38_ROOT/hg38.chrom.sizes.bed"
      },
      "mm10": {
        "refDict": "$MM10_ROOT/mm10.dict",
        "refIndex": "$MM10_ROOT/mm10.fa.fai",
        "indelsVcfGather_refIndex": "$MM10_ROOT/mm10.fa.fai",
        "snvsVcfGather_refIndex": "$MM10_ROOT/mm10.fa.fai",
        "refFasta": "$MM10_ROOT/mm10.fa",
        "refModule": "mm10/p6 samtools/1.9",
        "bedFile": "$MM10_ROOT/mm10.chrom.sizes.bed"
      }
    }   
 
    parameter_meta {
	tumorBam: "Input BAM file with tumor data"
	tumorBai: "BAM index file for tumor data"
	normalBam: "Input BAM file with normal data"
	normalBai: "BAM index file for normal data"
	reference: "Reference assembly id"
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
	    name: "python/3.7",
	    url: "https://www.python.org/downloads/release/python-378/"
	},
	{
	    name: "gatk/4.1.6.0",
	    url: "https://software.broadinstitute.org/gatk/download/index"
	}
	]
    output_meta: {
    snvsVcf: {
        description: "VCF file with SNVs, .gz compressed",
        vidarr_label: "snvsVcf"
    },
    indelsVcf: {
        description: "VCF file with indels, .gz compressed",
        vidarr_label: "indelsVcf"
    }
}
    }

	call splitIntervals {
	    input:
	    refFasta = resources[reference].refFasta,
	    refFai = resources[reference].refIndex,
	    refDict = resources[reference].refDict,
	    refModule = resources[reference].refModule,
	    intervals = select_first([bedFile, resources[reference].bedFile]),
            scatterCount = numChunk
	}
	Array[File] intervals = splitIntervals.intervalFiles
	call convertIntervalsToBed {
	    input:
	    intervalFiles = intervals
	}
	Array[File] bedIntervals = convertIntervalsToBed.bedFiles
	scatter(interval in bedIntervals) {
	    call configureAndRun as configureAndRunParallel {
		input:
		tumorBam = tumorBam,
		tumorBai = tumorBai,
		normalBam = normalBam,
		normalBai = normalBai,
		refFasta = resources[reference].refFasta,
		refIndex = resources[reference].refIndex,
		refModule = resources[reference].refModule,
		regionsBed = interval,
		outputFileNamePrefix = outputFileNamePrefix
	    }
	}
	call vcfGather as snvsVcfGather {
	    input:
	    vcfs = configureAndRunParallel.snvsVcf,
            refIndex = resources[reference].snvsVcfGather_refIndex
	}
	call vcfGather as indelsVcfGather {
	    input:
	    vcfs = configureAndRunParallel.indelsVcf,
            refIndex = resources[reference].indelsVcfGather_refIndex
	}
    
    # Do not output the TSV and XML stats files; contents not needed
    output {
	File snvsVcf = snvsVcfGather.result
	File indelsVcf = indelsVcfGather.result
    }
}

task configureAndRun {

    input {
	File tumorBam
	File tumorBai
	File normalBam
	File normalBai
	String refFasta
	String refIndex
	String refModule
	String? regionsBed
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

task convertIntervalsToBed {

    # convert GATK .interval_files to .bed by subtracting 1 from each interval start point

    input {
	Array[File] intervalFiles
	String modules = "python/3.7"
	Int memory = 16
	Int timeout = 4
    }

    parameter_meta {
	intervalFiles: "Files in GATK .interval_file format"
	modules: "Environment modules"
	memory: "Memory allocated for job"
	timeout: "Hours before task timeout"
    }

    meta {
	output_meta: {
	    bedFiles: "Output files converted to .bed format"
	}
    }

    command <<<
        python3 <<CODE
        import os, re
        intervalFiles = re.split(",", "~{sep=',' intervalFiles}")
        for intervalFile in intervalFiles:
            items = re.split("\.", os.path.basename(intervalFile))
            items.pop() # remove .interval_list suffix
            bedName = ".".join(items)+".bed"
            with open(intervalFile, 'r') as inFile, open(bedName, 'w') as outFile:
                for line in inFile:
                    if not re.match("@", line): # omit the GATK header
                        fields = re.split("\t", line.strip())
                        fields[1] = str(int(fields[1]) - 1)
                        outFile.write("\t".join(fields)+"\n")
        CODE
    >>>

    runtime {
	memory:  "~{memory} GB"
	modules: "~{modules}"
	timeout: "~{timeout}"
    }

    output {
	Array[File] bedFiles = glob("*.bed")
    }

}

task splitIntervals {

    input {
	String refFasta
	String refFai
	String refDict
	String refModule
	String gatk = "$GATK_ROOT/bin/gatk"
	String? intervals
	Int? scatterCount
	String? splitIntervalsExtraArgs
	String nonRefModules = "gatk/4.1.2.0"
	Int memory = 32
        Int overhead = 8
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
        overhead: "Memory overhead for running on a node"
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
	~{gatk} --java-options "-Xmx~{memory-overhead}g" SplitIntervals \
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
        String refIndex
	Array[File] vcfs
	Int memory = 16
	Int timeout = 12
    }

    parameter_meta {
	modules: "Environment module names and version to load (space separated) before command execution"
	gatk: "GATK to use"
        refIndex: "fai of the reference assembly, determines the ordering by chromosome"
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
	-R ~{refIndex} \
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
