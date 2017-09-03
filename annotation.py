#!/Users/Ryan/anaconda2/bin/python

import vcf
import json
import urllib2
import os
import getopt
import sys
from time import localtime, strftime

# vcf_file = 'task/Challenge_data.vcf'
# anno_file = 'task/anno_data.vcf'

def print_help():
    print "Annotate a given vcf file with follow information. Output an annotated vcf file."
    print "1. Type of variation. The most deleterious possibility from SnpEff output. [Effect_Type] in INFO field."
    print "2. Depth of sequence coverage at the site of variation. [DP] in INFO field"
    print "3. Number of reads supporting the variant. [AO] in INFO field"
    print "4. Percentage of reads supporting the variant. [Percent_Var] in INFO field"
    print "5. Percentage of reads supporting reference reads. [Percent_Ref] in INFO field"
    print "6. Allele frequency of variant from Broad Institute ExAC Project API. [ExAC_AF] in INFO field"
    print "7. Allele count of variant from Broad Institute ExAC Project API. [ExAC_AC] in INFO field"
    print "8. Total number of alleles from Broad Institute ExAC Project API. [ExAC_AN] in INFO field"
#    print "9. SnpEff annotations of highest possible impact. [ANN] in INFO field"

    print "\n"
    print "USAGE: python annotation.py -i Challenge_data.vcf -o Challenge_data.ann.vcf"
    print "OPTIONS: "
    print "-i, --input     str     input vcf file"
    print "-o, --output    str     output vcf file after annotation"
    print "-h, --help"

    
def run_SnpEff(vcf_file, snpEff_dir, genome, out_file):
    command = "java -jar -Xmx4g " + snpEff_dir + "/snpEff.jar " + genome + " " + vcf_file + " > " + out_file
    os.system(command)
    return


def fetch_ExAC_info(varID):
    url = "http://exac.hms.harvard.edu/" + "rest/variant/" + varID
    data = json.load(urllib2.urlopen(url))
    var = data['variant']
    ExAC_info = {}

    if 'allele_freq' in var.keys():
        ExAC_info['ExAC_AF'] = var['allele_freq']
    else:
        ExAC_info['ExAC_AF'] = ''
        
    if 'allele_count' in var.keys():
        ExAC_info['ExAC_AC'] = var['allele_count']
    else:
        ExAC_info['ExAC_AC'] = ''
    
    if 'allele_num' in var.keys():
        ExAC_info['ExAC_AN'] = var['allele_num']
    else:
        ExAC_info['ExAC_AN'] = ''
        
    return ExAC_info


def max_impact(ANN):
    
    if len(ANN) == 0:
        max_ANN = []
        effect_type = ''
    else:
        impact_grade = {
            'HIGH': 4,
            'MODERATE': 3,
            'LOW': 2,
            'MODIFIER': 1
        }
        impact = [impact_grade[item.split('|')[2]] for item in ANN]
        max_impact = max(impact)
        max_ANN = [item for item in ANN if impact_grade[item.split('|')[2]] == max_impact]
        effect_type = [item.split('|')[1] for item in max_ANN][0]
        
    return max_ANN, effect_type


def annotate(vcf_file, anno_file, SNPEFF_DIR, GENOME_VER, INFO_TO_KEEP):
    
    filename = os.path.splitext(vcf_file)[0]
    
    ## perform snpEff annotaion of variants' deleteriousness
    snpEff_ann_vcf_file = filename + '.snpeff.ann.vcf'
    
    now = strftime("%Y%m%d-%H%M%S", localtime())
    sys.stdout.write(now + ": Start running SnpEff ...\n")
    
    run_SnpEff(vcf_file, SNPEFF_DIR, GENOME_VER, snpEff_ann_vcf_file)
    
    now = strftime("%Y%m%d-%H%M%S", localtime())
    sys.stdout.write(now + ": SnpEff annotation finished successfully ...\n")
    
    
    
    ## reformat the output vcf file with desired output items
    vcf_read = vcf.Reader(open(snpEff_ann_vcf_file, 'r'))
    
    ## modify INFO section in vcf header
    vcf_read.infos = {key:value for key, value in vcf_read.infos.items() if key in INFO_TO_KEEP}
    vcf_read.infos['Percent_Var'] = vcf.parser._Info(
        id = 'Percent_Var',
        num = 'A',
        type = 'Float',
        desc = 'Percentage of reads supporting the variant',
        source = None,
        version = None        
    )
    vcf_read.infos['Percent_Ref'] = vcf.parser._Info(
        id = 'Percent_Ref',
        num = 1,
        type = 'Float',
        desc = 'Percentage of reads supporting the reference',
        source = None,
        version = None        
    )
    vcf_read.infos['ExAC_AF'] = vcf.parser._Info(
        id = 'ExAC_AF',
        num = 'A',
        type = 'Float',
        desc = 'ExAC allele frequency',
        source = None,
        version = None        
    )
    vcf_read.infos['ExAC_AC'] = vcf.parser._Info(
        id = 'ExAC_AC',
        num = 'A',
        type = 'Integer',
        desc = 'ExAC allele count',
        source = None,
        version = None        
    )
    vcf_read.infos['ExAC_AN'] = vcf.parser._Info(
        id = 'ExAC_AN',
        num = 'A',
        type = 'Integer',
        desc = 'ExAC total number of alleles',
        source = None,
        version = None        
    )    
    vcf_read.infos['Effect_Type'] = vcf.parser._Info(
        id = 'Effect_Type',
        num = 'A',
        type = 'String',
        desc = 'Effect type of variant',
        source = None,
        version = None        
    )    
    
    vcf_write = vcf.Writer(open(anno_file, 'w'), vcf_read)
    
    
    ## extracting desired information and fetching from ExAC database through API
    now = strftime("%Y%m%d-%H%M%S", localtime())
    sys.stdout.write(now + ": Start annotating with ExAC information ...\n")    
    
    for record in vcf_read:
        
        new_record = vcf.model._Record(
            CHROM = record.CHROM,
            POS = record.POS,
            ID = record.ID,
            REF = record.REF,
            ALT = record.ALT,
            QUAL = record.QUAL,
            FILTER = record.FILTER,
            INFO = {},
            FORMAT = record.FORMAT,
            sample_indexes = '',
            samples = record.samples
        )

        new_record.INFO['AO'] = record.INFO['AO']
        new_record.INFO['RO'] = record.INFO['RO']
        new_record.INFO['DP'] = record.INFO['DP']
        new_record.INFO['Percent_Ref'] = float(record.INFO['RO'])/record.INFO['DP']*100
        
        new_record.INFO['Percent_Var'] = []
        new_record.INFO['ExAC_AF'] = []
        new_record.INFO['ExAC_AC'] = []
        new_record.INFO['ExAC_AN'] = []       
        new_record.INFO['Effect_Type'] = []
        
        max_ANN = []
        for i in range(len(record.ALT)):
            var = str(record.ALT[i])
            ANN = [item for item in record.INFO['ANN'] if item.split('|')[0] == var]
            max_impact_ANN = max_impact(ANN)
            max_ANN.extend(max_impact_ANN[0])
            e_type = max_impact_ANN[1]
            
            varID = '-'.join([
                str(record.CHROM),
                str(record.POS),
                str(record.REF),
                var
            ])
            
            ExAC_info = fetch_ExAC_info(varID)
            new_record.INFO['Percent_Var'].append(float(record.INFO['AO'][i])/record.INFO['DP']*100)
            new_record.INFO['ExAC_AF'].append(ExAC_info['ExAC_AF'])
            new_record.INFO['ExAC_AC'].append(ExAC_info['ExAC_AC'])
            new_record.INFO['ExAC_AN'].append(ExAC_info['ExAC_AN'])
            new_record.INFO['Effect_Type'].append(e_type)
            
#        new_record.INFO['ANN'] = max_ANN

        vcf_write.write_record(new_record)
        
    now = strftime("%Y%m%d-%H%M%S", localtime())
    sys.stdout.write(now + ": Finished successfully!\n")
    
    return
   

def main():
    
    ## a few configuration variables
    
    # please replace SNPEFF_DIR with the SnpEff directory
    SNPEFF_DIR = "~/BioSoftwares/source/snpEff/"
    
    # version of the genome for SnpEff annotation
    GENOME_VER = "GRCh37.75"
    
    # items to keep in the INFO field
    INFO_TO_KEEP = ['AO', 'RO', 'DP', 'Percent_Var', 'Percent_Ref', 'ExAC_AF', 'ExAC_AC', 'ExAC_AN', 'Effect_Type']
    
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "i:o:h", ["input=", "output=", "help"])
    except getopt.GetoptError as err:
        print str(err)
        print_help()
        sys.exit()
    
    for o,a in opts:
        if o in ("-i", "--input"):
            vcf_file = a
        elif o in ("-o", "--output"):
            anno_file = a        
        elif o in ("-h", "--help"):
            print_help()
            sys.exit()

    if len(sys.argv[1:]) > 0:
        annotate(vcf_file, anno_file, SNPEFF_DIR, GENOME_VER, INFO_TO_KEEP)
    else:
        print_help()       
        
if __name__=="__main__":
    main()
