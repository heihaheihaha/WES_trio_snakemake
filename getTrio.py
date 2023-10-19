#Rewrite from getTrio.pl
# 2023/10/14
import argparse
import os
import csv
import subprocess

parser = argparse.ArgumentParser(description="Get DNMs\nUsage:\n--input <preDNMs file> <initial.extreme> <initial.Trios.extreme>")
parser.add_argument("-i","--input", required=True, nargs="+", help="input file")
parser.add_argument("-o","--output", required=True, help="output file")

args = parser.parse_args()

if len(args.input) != 3:
	print("Usage:\n--input <preDNMs file> <initial.extreme> <initial.Trios.extreme>, 3 files are needed")
	exit(1)

output = args.output
output_path = os.path.abspath(output)
output_dir = os.path.dirname(output_path)
tmp_full_output = output_dir+"/DNMs.xls"

print("OUTPUT will be saved at", tmp_full_output)

OUTPUT = args.output


# Check if the input files exist
def test_file(file):
	if not os.path.exists(file):
		print(f"{file} does not exist")
		exit(1)

for i in args.input:
	test_file(i)

"""
.preDNMs has columns like this:
#Chr    start   end     ref     alt     DeNovoCNN probability   Child coverage  Father coverage Mother coverage
which is the filered output of DeNovoCNN

FIlter is done by rule Filter_DenovoCNN, which is also the `lib/FilterDNMs.pl`

Two parameters are used:
- DeNovoCNN probability
- Average coverage of child, father and mother

FIlter_DenovoCNN also strandized the variants,here is an example:
in predictions.csv:
# reference    variant
T       TGG

it should be -    A  in the preDNMs file

predictions.csv have more columns:
Chromosome      Start position  Reference       Variant DeNovoCNN probability   Child coverage  Father coverage Mother coverage Child BAM       Father BAM      Mother BAM
"""
#Read in the pos of denovo variants
tp = open(args.input[0], "r") # preDNMs file
DNM_POS = []
for line in tp:
	if line.startswith("#"):
		continue
	line = line.strip().split("\t")
	DNM_POS.append("\t".join(line[0:5]))
tp.close()
"""
Trio.initial.tmp.xls is the output of ExtremeVar2
It has columns like this:
Extreme Ratio_of_tools (ReVe,gt,0.7:VEST3,gt,0.6:REVEL,gt,0.4:CADD,gt,20:GERP++,gt,2)   Rare_or_Common  Chr     Start   End     Ref     Alt     Func.refGene    Gene.refGene        GeneDetail.refGene      ExonicFunc.refGene      AAChange.refGene        GeneFullName.refGene    GeneFunction.refGene    GeneExpressionTissue.refGeneGeneDisease.refGene     OMIM.refGene    RVIS_percentile.refGene RVIS.refGene    GDI-Phred.refGene       GDI.refGene     Tissue_specificity(Uniprot).refGene     P(rec).refGene      P(HI).refGene   pLI.refGene     pRec.refGene    pNull.refGene   Expression(GNF/Atlas).refGene   cytoBand        rmsk    avsnp150        SIFT    PolyPhen2-HDIV      PolyPhen2-HVAR  LRT     MutationTaster  MutationAssessor        FATHMM  PROVEAN VEST3   MetaSVM MetaLR  M-CAP   CADD    DANN    fathmm-MKL      Eigen       GenoCanyon      fitCons GERP++  phyloP  phastCons       SiPhy   REVEL   ReVe    Interpro_domain gene4denovo(phenotype|Pumbed_ID)        ChinaMap        gnomAD_exome_ALL    gnomAD_exome_AFR        gnomAD_exome_SAS        gnomAD_exome_AMR        gnomAD_exome_EAS        gnomAD_exome_NFE        gnomAD_exome_FIN        gnomAD_exome_ASJ    gnomAD_exome_OTH        gnomad312_AF    gnomad312_AF_raw        gnomad312_AF_XX gnomad312_AF_XY gnomad312_AF_popmax     gnomad312_AF_faf95_popmax  gnomad312_AF_afr gnomad312_AF_ami        gnomad312_AF_amr        gnomad312_AF_asj        gnomad312_AF_eas        gnomad312_AF_fin        gnomad312_AF_mid        gnomad312_AF_nfe    gnomad312_AF_oth        gnomad312_AF_sas        ExAC_ALL        ExAC_AFR        ExAC_AMR        ExAC_EAS        ExAC_FIN        ExAC_NFE        ExAC_OTH    ExAC_SAS        esp6500siv2_all 1000g2015aug_all        1000g2015aug_afr        1000g2015aug_amr        1000g2015aug_eas        1000g2015aug_eur        1000g2015aug_sas    CLNID   CLNDN   CLNDISDB        CLNREVSTAT      CLNSIG  CLNSIGCONF      HGMD_based_tag  HGMD_based_disease      HGMD_based_pub  InterVar_automated dbscSNV_ADA_SCORE        dbscSNV_RF_SCORE        Otherinfo
N       0.00    Common  chr1    69511   69511   A       G       exonic  OR4F5   -       nonsynonymous SNV       OR4F5:NM_001005484:exon1:c.421A>G:p.T141A       olfactory receptor family 4 subfamily F member 5    FUNCTION: Odorant receptor. {ECO:0000305}.;     -       -       -       -       -       2.31572 113.16669       -  0.07159  0.06092 0.176329298172162       0.644086264144236       0.179584437683602       -       1p36.33 -       rs2691305       0.652:T 0.0:B   0.0:B   0.000:N 1.000:P     -0.855:T        1.26:T  1.54:N  0.012:T -0.997:T        0.000:T -:-     0.047:T 0.454:T 0.003:N -1.341:T        0.000:T 0.487:T 1.15:N  -0.684:N        0.000:N     4.198:N 0.053:T 0.041:T GPCR: rhodopsin-like: 7TM       -       0.998295        0.9497  0.6075  0.9854  0.9514  0.9994  0.9726  0.9916  0.9767  0.9506  0.8460      0.8614  0.8427  0.8496  0.9998  0.9752  0.5948  0.9980  0.8951  0.9784  0.9998  0.9907  0.9     0.9674  0.8624  0.9772  0.9394  0.5942  0.9507  0.9994  0.9907      0.9716  0.9597  0.9832  0.7598  -       -       -       -       -       -       -       -       -       -       -       -       -       -       -       Benign(PM1|BA1|BS1|BP4)     -       -       1.000000        10520.73        Tsample=3,Csample=3,Hetsample=0,Homsample=3     CHILD=Hom1,Dpth=131,GQ=99,HetR=1;FU0=Hom1,Dpth=104,GQ=99,HetR=1;MU0=Hom1,Dpth=121,GQ=99,HetR=1
"""

writer = csv.writer(open(tmp_full_output, "w"), delimiter="\t")
# Read in the header of Trio.initial.tmp.xls
with open(args.input[1], "r") as infile:
	reader = csv.reader(infile, delimiter="\t")
	title = next(reader)
	writer.writerow(["Inherited_type"]+title) #Write the title row
	for row in reader:
		# Check if the variant is denovo
		# If it is denovo, marked is as DNMs in column 1
		if "\t".join(row[2:7]) in DNM_POS:
			writer.writerow(["DNM"]+row) # Only write the denovo variants

with open(args.input[2], "r") as infile: #initial.Trios.extreme
	reader = csv.reader(infile, delimiter="\t")
	title = next(reader)
	for row in reader:
		# Check if the variant is denovo
		# If it is denovo, marked is as DNMs in column 1
		if "Inherited_type" not in row:
			writer.writerow(row)
		elif "DNMs" not in row:
			writer.writerow(row)

bash_command = "cut -f 1-11,13-27,31-34,38,42,46,52-54,56-61,65,70,80,86,89,94-95,98,101-112,116 "+tmp_full_output+" > "+OUTPUT
subprocess.run(bash_command, shell=True)