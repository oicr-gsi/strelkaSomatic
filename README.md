# strelkaSomatic

Strelka variant caller in somatic mode

## Overview

## Dependencies

* [samtools 1.9](https://github.com/samtools/samtools)
* [strelka 2.9.10](https://github.com/Illumina/strelka/releases/tag/v2.9.10)
* [python 2.7](https://www.python.org/downloads/release/python-2716/)
* [python 3.7](https://www.python.org/downloads/release/python-378/)
* [gatk 4.1.2.0](https://software.broadinstitute.org/gatk/download/index)
* [bcftools 1.9](https://github.com/samtools/bcftools)


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
`outputFileNamePrefix`|String|Prefix for output files


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`bedFile`|String?|None|BED file designating regions to process
`numChunk`|Int?|None|If BED file given, number of chunks in which to split each chromosome
`mode`|String|"WG"|WG (default), exome or targeted


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
`snvsVcfGather.modules`|String|"gatk/4.1.2.0 tabix/1.9"|Environment module names and version to load (space separated) before command execution
`snvsVcfGather.gatk`|String|"$GATK_ROOT/bin/gatk"|GATK to use
`snvsVcfGather.memory`|Int|16|Memory allocated for job
`snvsVcfGather.timeout`|Int|12|Hours before task timeout
`indelsVcfGather.modules`|String|"gatk/4.1.2.0 tabix/1.9"|Environment module names and version to load (space separated) before command execution
`indelsVcfGather.gatk`|String|"$GATK_ROOT/bin/gatk"|GATK to use
`indelsVcfGather.memory`|Int|16|Memory allocated for job
`indelsVcfGather.timeout`|Int|12|Hours before task timeout
`vcfCombine.modules`|String|"bcftools/1.9 tabix/1.9"|environment modules
`vcfCombine.jobMemory`|Int|16|Memory allocated for job
`vcfCombine.threads`|Int|4|Number of threads for processing
`vcfCombine.timeout`|Int|4|Hours before task timeout
`injectFields.modules`|String|"strelkasomatic-scripts/1.0 tabix/1.9"|environment modules
`injectFields.jobMemory`|Int|16|Memory allocated for job
`injectFields.threads`|Int|4|Number of threads for processing
`injectFields.timeout`|Int|4|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`vcfSnvs`|File|VCF file with SNVs, .gz compressed|vidarr_label: vcfSnvs
`vcfIndels`|File|VCF file with indels, .gz compressed|vidarr_label: vcfIndels
`vcfAll`|File|VCF file with SNVs and indels, .gz compressed|vidarr_label: vcfAll
`vcfAllIndex`|File|tabix index for VCF file with SNVs and indels|vidarr_label: vcfAllIndex
`vcfExtended`|File|extended VCF file with SNVs and indels, with added gt + ad fields .gz compressed|vidarr_label: vcfExtended
`vcfExtendedIndex`|File|tabix index for extended VCF|vidarr_label: vcfExtendedIndex


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
 
 
 ### Merging the SNVs and Indels into a single file
 
 Strelka will produce separate files for SNVs and Indels.  
 Downstream processing sometimes needs these in the same file
 The vcfCombine tasks combines the data and indexes the bgzipped output
 
 ```
 set -eo pipefail
 	
 	bcftools concat -a -o ~{outputFileNamePrefix}.strelka2_all.vcf ~{vcfSnvs} ~{vcfIndels}
 	bgzip ~{outputFileNamePrefix}.strelka2_all.vcf
 	tabix ~{outputFileNamePrefix}.strelka2_all.vcf.gz
 ```
 
 ### Add in GT and AD fields
 
 Strelka does not utilize GT or AD as part of the sample fields.
 These are useful to have and sometimes required for downstream processing.
 Appropriate values can be generated in a variety of ways from existing information in the INFO and Sample FIELDS
 A helper python script, scripts/strelka_add_gt_ad.py, is used to create and extended.vcf file in the injectFields task
 
 ```
 	set -eo pipefail
 	python3 strelka_add_gt_ad.py -i ~{vcfIn} -o ~{outputFileNamePrefix}.strelka2_all.extended.vcf
 	bgzip ~{outputFileNamePrefix}.strelka2_all.extended.vcf
 	tabix ~{outputFileNamePrefix}.strelka2_all.extended.vcf.gz
 ```
 
 
 
 
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
