# strelkaSomatic

Strelka variant caller in somatic mode

## Overview

## Dependencies

* [samtools 1.9](https://github.com/samtools/samtools)
* [strelka 2.9.10](https://github.com/Illumina/strelka/releases/tag/v2.9.10)
* [python 2.7](https://www.python.org/downloads/release/python-2716/)
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
`normalBam`|File|Input BAM file with normal data
`refFasta`|File|Reference FASTA file
`refIndex`|File|Reference FAI index


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`bedFile`|File?|None|BED file designating regions to process
`numChunk`|Int?|None|If BED file given, number of chunks in which to split each chromosome
`outputFileNamePrefix`|String|"strelkaSomatic"|Prefix for output files


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`splitIntervals.modules`|String|"gatk/4.1.2.0 hg19/p13"|Environment module names and version to load (space separated) before command execution
`splitIntervals.refFasta`|String|"$HG19_ROOT/hg19_random.fa"|Path to the reference fasta
`splitIntervals.refFai`|String|"$HG19_ROOT/hg19_random.fa.fai"|Path to the reference .fai index
`splitIntervals.refDict`|String|"$HG19_ROOT/hg19_random.dict"|Path to the reference .dict dictionary
`splitIntervals.splitIntervalsExtraArgs`|String?|None|Additional arguments for the 'gatk SplitIntervals' command
`splitIntervals.memory`|Int|16|Memory allocated for job
`splitIntervals.timeout`|Int|4|Hours before task timeout
`configureAndRunParallel.modules`|String|"python/2.7 samtools/1.9 strelka/2.9.10"|Environment module names and version to load (space separated) before command execution
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
`configureAndRunSingle.regionsBed`|File?|None|BED file designating regions to process
`configureAndRunSingle.modules`|String|"python/2.7 samtools/1.9 strelka/2.9.10"|Environment module names and version to load (space separated) before command execution
`configureAndRunSingle.jobMemory`|Int|16|Memory allocated for job
`configureAndRunSingle.threads`|Int|4|Number of threads for processing
`configureAndRunSingle.timeout`|Int|4|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`snvsVcf`|File|None
`indelsVcf`|File|None


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
