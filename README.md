# strelkaSomatic

Strelka variant caller in somatic mode

## Overview

## Dependencies

* [samtools 1.9](https://github.com/samtools/samtools)
* [strelka 2.9.10](https://github.com/Illumina/strelka/releases/tag/v2.9.10)
* [python 2.7](https://www.python.org/downloads/release/python-2716/)
* [python 3.7](https://www.python.org/downloads/release/python-378/)
* [gatk 4.1.6.0](https://software.broadinstitute.org/gatk/download/index)


## Usage

### Cromwell
```
java -jar cromwell.jar run strelkaSomatic.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`tumorBam`|File|Input BAM file with tumor data
`tumorBai`|File|BAM index file for tumor data
`normalBam`|File|Input BAM file with normal data
`normalBai`|File|BAM index file for normal data
`reference`|String|Reference assembly id


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`bedFile`|String?|None|BED file designating regions to process
`numChunk`|Int?|None|If BED file given, number of chunks in which to split each chromosome
`outputFileNamePrefix`|String|"strelkaSomatic"|Prefix for output files


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`splitIntervals.gatk`|String|"$GATK_ROOT/bin/gatk"|GATK executable path
`splitIntervals.splitIntervalsExtraArgs`|String?|None|Additional arguments for the 'gatk SplitIntervals' command
`splitIntervals.nonRefModules`|String|"gatk/4.1.2.0"|Environment modules other than the genome refence
`splitIntervals.memory`|Int|32|Memory allocated for job
`splitIntervals.overhead`|Int|8|Memory overhead for running on a node
`splitIntervals.timeout`|Int|72|Hours before task timeout
`convertIntervalsToBed.modules`|String|"python/3.7"|Environment modules
`convertIntervalsToBed.memory`|Int|16|Memory allocated for job
`convertIntervalsToBed.timeout`|Int|4|Hours before task timeout
`configureAndRunParallel.nonRefModules`|String|"python/2.7 samtools/1.9 strelka/2.9.10"|Environment module names other than genome reference
`configureAndRunParallel.jobMemory`|Int|16|Memory allocated for job
`configureAndRunParallel.threads`|Int|4|Number of threads for processing
`configureAndRunParallel.timeout`|Int|4|Hours before task timeout
`snvsVcfGather.modules`|String|"gatk/4.1.2.0"|Environment module names and version to load (space separated) before command execution
`snvsVcfGather.gatk`|String|"$GATK_ROOT/bin/gatk"|GATK to use
`snvsVcfGather.memory`|Int|16|Memory allocated for job
`snvsVcfGather.timeout`|Int|12|Hours before task timeout
`indelsVcfGather.modules`|String|"gatk/4.1.2.0"|Environment module names and version to load (space separated) before command execution
`indelsVcfGather.gatk`|String|"$GATK_ROOT/bin/gatk"|GATK to use
`indelsVcfGather.memory`|Int|16|Memory allocated for job
`indelsVcfGather.timeout`|Int|12|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`snvsVcf`|File|VCF file with SNVs, .gz compressed|vidarr_label: snvsVcf
`indelsVcf`|File|VCF file with indels, .gz compressed|vidarr_label: indelsVcf


## Commands
 
This section lists command(s) run by strelkaSomatic workflow
 
* Running strelkaSomatic
 
### Main task, configuring and running Strelka2
 
```
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
```
 
### Converting interval files to .bed
 
.bed files use 0-based index, we need to do a proper conversion of interval files 
 
```
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
```
 
### Splitting intervals
 
SplitIntervals is used to produce a requested number of files
to make the analysis parallel, decreasing demand for resources per chunk
and improving speed

```
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
```
 
### Combining vcf files after scatter
 
This task uses GatherVcfs tools from GATK suite which needs the files to be ordered. 
No overlap between covered features allowed.
 
```
    set -eo pipefail
 
    ~{gatk} GatherVcfs \
    -I ~{sep=" -I " vcfs} \
    -R ~{refFasta}
    -O ~{outputName}
 
    gzip ~{outputName}
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
