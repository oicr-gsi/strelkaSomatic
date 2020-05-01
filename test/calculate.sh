#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

md5sum somatic.indels.vcf.gz somatic.snvs.vcf.gz

