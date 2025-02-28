import argparse
import vcf
import re


parser = argparse.ArgumentParser(description='output file')
parser.add_argument('-i', '--input', help='input vcf file path', required=True)
parser.add_argument('-o', '--output', help='output vcf file path', required=True)
#parser.add_argument('-r', '--reference',help="reference assembly, this will be added to the strelka header", required=True)


args = parser.parse_args()


def strelka_genotypes(ref, alt, info):
    """

    Retrieve standard 0/0, 0/1, 1/1 style genotypes from INFO field.
    Normal -- NT field (ref, het, hom, conflict)
    Tumor -- SGT field
      - for SNPs specified as GG->TT for the normal and tumor diploid alleles. These
        can also represent more complex alleles in which case we set at heterozygotes
        pending longer term inclusion of genotypes in Strelka2 directly
        (https://github.com/Illumina/strelka/issues/16)
      - For indels, uses the ref, het, hom convention
    ref: The REF allele from a VCF line
    alt: A list of potentially multiple ALT alleles (rec.ALT.split(";"))
    info: The VCF INFO field
    fname, coords: not currently used, for debugging purposes

    """
    known_names = set(["het", "hom", "ref", "conflict"])

    def name_to_gt(val):
        if val.lower() == "het":
            return "0/1"
        elif val.lower() == "hom":
            return "1/1"
        elif val.lower() in set(["ref", "conflict"]):
            return "0/0"
        else:
            """ Non-standard representations, het is our best imperfect representation """
            return "0/1"

    def alleles_to_gt(val):
        gt_indices = {}
        for j, gT in enumerate([ref] + alt):
            gt_indices.update({str(gT).upper(): j})
        tumor_gts = [gt_indices[x.upper()] for x in val if x in gt_indices]
        if tumor_gts and val not in known_names:
            if max(tumor_gts) == 0:
                tumor_gt = "0/0"
            elif 0 in tumor_gts:
                tumor_gt = "0/%s" % min([x for x in tumor_gts if x > 0])
            else:
                tumor_gt = "%s/%s" % (min(tumor_gts), max(tumor_gts))
        else:
            tumor_gt = name_to_gt(val)
        return tumor_gt

    nt_val = info['NT']
    normal_gt = name_to_gt(nt_val)
    sgt_val = info['SGT']
    if not sgt_val:
        tumor_gt = "0/0"
    else:
        sgt_val = sgt_val.split("->")[-1]
        tumor_gt = alleles_to_gt(sgt_val)
    
    gt={"TUMOR":tumor_gt,"NORMAL":normal_gt}
    return gt


def strelka_allele_depth(ref,alt,sample_info):
    """
      this will vary in Strelka files for snvs vc indels
      
      Retrieve allele reads and generate AD in form of REF_count,ALT1_count,ALT2_count etc.
      ref: The REF allele from a VCF line
      alt: A list of potentially multiple ALT alleles (rec.ALT.split(";"))
    """
    ## check if this is a snv record by the presence of specific fields
    if all( hasattr(sample_info.data,x) for x in ['AU','CU','GU','TU']):
       depth_ref = ','.join([str(sample_info[x + 'U'][0]) for x in ref])
       depth_alt = "0"
       if alt[0] is not None:
          depth_alt = ','.join([str(sample_info[str(x) + 'U'][0]) for x in alt])

    elif all( hasattr(sample_info.data,x) for x in ['TAR','TIR','TOR']):
       ### see notes on https://groups.google.com/g/strelka-discuss/c/Pz9EYKnE6u8?pli=1   TAR = reference + reads not supporting INDEL, TIR = reads supporting INDEL
       ### using only the tier- value [0]
       depth_ref = str(sample_info['TAR'][0])
       depth_alt = str(sample_info['TIR'][0])
    else:
       depth_ref="."
       depth_alt="."

    return [depth_ref,depth_alt]


def generate_header(vcf_reader: vcf.Reader, samples: list):
    ### file format
    header_lines = ["##fileformat=" + vcf_reader.metadata['fileformat'] + "\n"]
    

    ### FILTERS
    for flt in vcf_reader.filters:
        header_lines.append("##FILTER=<ID=" + vcf_reader.filters[flt].id +
                            ",Description=\"" + vcf_reader.filters[flt].desc + "\">\n")
    ### FORMATS
    for f in vcf_reader.formats:
        num = vcf_reader.formats[f].num
        if vcf_reader.formats[f].num is None:
            num = '.'
        elif isinstance(vcf_reader.formats[f].num, int) and vcf_reader.formats[f].num < 0:
            num = _special_character(vcf_reader.formats[f].num)

        header_lines.append("##FORMAT=<ID=" + vcf_reader.formats[f].id +
                            ",Number=" + str(num) + ",Type=" + vcf_reader.formats[f].type +
                            ",Description=\"" + vcf_reader.formats[f].desc + "\">\n")

    ### add in the GT and AD formats
    header_lines.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\""
                            "Genotype, constructed from Strelka SGT INFO\">\n")
    header_lines.append("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\""
                            "Allelic depths for the ref and alt alleles in the order listed, computed from Strelka Sample information\">\n")

    ### INFO
    for i in vcf_reader.infos:
        num = vcf_reader.infos[i].num
        if vcf_reader.infos[i].num is None:
            num = '.'
        elif isinstance(vcf_reader.infos[i].num, int) and vcf_reader.infos[i].num < 0:
            num = _special_character(vcf_reader.infos[i].num)
        header_lines.append("##INFO=<ID=" + vcf_reader.infos[i].id +
                            ",Number=" + str(num) + ",Type=" + vcf_reader.infos[i].type +
                            ",Description=\"" + vcf_reader.infos[i].desc + "\">\n")

    ### contigs
    for c in vcf_reader.contigs:
        header_lines.append("##contig=<ID=" + vcf_reader.contigs[c].id + ",length=" + str(vcf_reader.contigs[c].length) + ">\n")
    
    ### everything else
    for m in vcf_reader.metadata:
        if m == "fileformat": continue
        v=",".join(map(str, vcf_reader.metadata[m])) if isinstance(vcf_reader.metadata[m], list) else vcf_reader.metadata[m]
        line = "##" + m + "=" + v + "\n"
        header_lines.append(line)

    header_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    header_fields.extend(samples)
    header_lines.append("#" + "\t".join(header_fields) + "\n")
    return header_lines



vcf_reader = vcf.Reader(filename=args.input, compressed=True)

revised_records=[]
### capture the samples in this file, to inject into the reformed header
for rec in vcf_reader:
    
  sampleIds=[rec.samples[0].sample,rec.samples[1].sample]

  ### start with CHROM and POS
  revised_fields = [rec.CHROM, str(rec.POS)]
  
  ### ID
  value = "." if rec.ID is None else record.ID
  revised_fields.append(value)

  ### REF  
  revised_fields.append(rec.REF)

  ### ALT  
  value = ",".join(map(str, rec.ALT))
  revised_fields.append(value)

  ### QUAL  
  value = "." if rec.QUAL is None else rec.QUAL
  revised_fields.append(value)

  ### FILTER
  value = "."
  if rec.FILTER is None or len(rec.FILTER) == 0:
      value = "." if len(vcf_reader.filters) == 0 else 'PASS'
  else: 
      value = ";".join(rec.FILTER)
  revised_fields.append(value)

  ### INFO
  value = ";".join(list(map(lambda k: k + "=" + str(rec.INFO[k]),rec.INFO)))
  revised_fields.append(value)
  
  ### FORMAT
  ### Add the GT and AD Fields to the start of indel records
  value = "GT:AD:" + rec.FORMAT
  revised_fields.append(value)
  
  gts=strelka_genotypes(rec.REF,rec.ALT,rec.INFO)
  
  
  ### SAMPLES
  ### for the GT field we need to know which sample is the Tumour and which is the Normal
  for sample_info in rec.samples:
      gt=gts[sample_info.sample]
      ad=",".join(strelka_allele_depth(rec.REF,rec.ALT,sample_info))
      sample_data=[gt,ad]
      for k in rec.FORMAT.split(":"):
          ## form the value for scalar vs list differently, creating a comma separated string of list values
          v=",".join(map(str, sample_info[k])) if isinstance(sample_info[k], list) else sample_info[k]
          sample_data.append(v)

      value=":".join(map(str,sample_data))
      revised_fields.append(value)

  ### create the revised tab separated record from the new fields
  
  
  revised_record="\t".join(revised_fields) + "\n"
  revised_records.append(revised_record)


## generate the header lines, this requires 
header_lines = generate_header(vcf_reader, sampleIds)

with open(args.output, mode='+w') as out:
    out.writelines(header_lines)
    out.writelines(revised_records)





 
   



