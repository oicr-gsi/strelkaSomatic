# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.3.0] - 2025-02-28
### Added
- [GRD-869](https://jira.oicr.on.ca/browse/GRD-869)
- modifications to output filenames, indicating strelka2
- modifications to output filenames to include the outputFileNamePrefix, this was not working properly
- added a mode parameter, that gives ability to set --exome and --targeted flags
- added task vcfCombine, which combines the snvs and indels into a single all output
- added a python script (strelka_add_gt_ad.py) for generating missing gt and ad fields; this is based on prior work in the variantMerging workflow.  the python script is available through the strelkaSomatic-scripts module
- added an injectFields task that runs the python script
- added two additonal outputs, the _all.vcf.gz (combined snvs + indels) and _all.extended.vcf.gz (with gt + ad informaition)
- extended regression tests to account for new outputs

## [1.2.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add [idarr labels to outputs (changes to medata only)

## [1.1.0] 2022-06-28
### Added
- assembly-specific choices move into wdl

## [1.0.0] 2020-10-24
### Changed
- bed file as String, minor changes. Upgrade for Vidarr

## [0.2.0] 2020-09-02
### Added
- GP-2466 Convert .interval_list files to .bed format before input to Strelka

## [0.1.0] 2020-07-15
### Added
- Initial production release of Strelka 2 somatic mode workflow
