#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

zcat somatic.indels.vcf.gz | egrep -v "^##(cmdline|fileDate|reference|startTime)" | md5sum
zcat somatic.snvs.vcf.gz | egrep -v "^##(cmdline|fileDate|reference|startTime)" | md5sum
