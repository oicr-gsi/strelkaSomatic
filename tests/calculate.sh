#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

zcat regressionTest.strelka2_snvs.vcf.gz | egrep -v "^##(cmdline|fileDate|reference|startTime)" | md5sum
zcat regressionTest.strelka2_indels.vcf.gz | egrep -v "^##(cmdline|fileDate|reference|startTime)" | md5sum
zcat regressionTest.strelka2_all.vcf.gz | egrep -v "^##(cmdline|fileDate|reference|startTime)" | md5sum
zcat regressionTest.strelka2_all.extended.vcf.gz | egrep -v "^##(cmdline|fileDate|reference|startTime)" | md5sum


