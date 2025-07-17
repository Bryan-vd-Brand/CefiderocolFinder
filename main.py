import argparse
import logging
import os
import sys
import subprocess
import csv
import pandas as pd
import glob



def setup_logger(log_file):
    """Setup logger to log messages to a file and stdout."""
    logger = logging.getLogger("CefiderocolFinder")
    logger.setLevel(logging.DEBUG)

    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')
    console_handler.setFormatter(console_formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger

def run_subprocess(command):
    """Helper function to run subprocess commands and handle errors."""
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        logger.info(f"Command succeeded: {command}")
        return result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {e.cmd}\nError: {e.stderr}")
        sys.exit(1)


def read_config(config_file):
    """Read the configuration file."""
    if not os.path.exists(config_file):
        logger.error(f"Configuration file not found: {config_file}")
        sys.exit(1)

    config = {}
    try:
        with open(config_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    key, value = line.strip().split('\t', 1)
                    config[key.strip()] = value.strip()
    except Exception as e:
        logger.error(f"Error reading configuration file: {e}")
        sys.exit(1)

    return config

def delete_intermediates(output_dir,name):
    """Delete intermediate files in the output directory."""
    try:
        for root, dirs, files in os.walk(output_dir):
            for file in files:
                if name in file:
                    if file.endswith(".sam") or file.endswith(".bam"):
                        os.remove(os.path.join(root, file))
                        logger.info(f"Deleted file: {file}")
    except Exception as e:
        logger.error(f"Error deleting intermediate files: {e}")
        sys.exit(1)

def read_adaptations(tsv_file, species):
    """Read a TSV file containing species adaptations and return all adaptations for the specified species."""
    #Species	Gene	Adaptation	Nomenclature
    if not os.path.exists(tsv_file):
        logger.error(f"Adaptations file not found: {tsv_file}")
        sys.exit(1)
    
    adaptations_df = pd.read_table(tsv_file,header=0,sep='\t')

    if not isNotEmptyDF(adaptations_df):
        logger.warning(f"No adaptations found for species: {species}")
    
    adaptations = adaptations_df[adaptations_df["Species"] == species]

    return adaptations

def isNotEmptyDF (dataframe):
    return dataframe.shape[0] != 0

def check_gene(df,geneName,query):
    """Check for gene in ANN column"""
    df_gene = df[df['ANN'].str.contains(geneName)]
    df_gene = df_gene[df_gene['ANN'].str.contains(query)]
 
    if isNotEmptyDF(df_gene):
        return df_gene
    else:
        return pd.DataFrame()

def AddFilterDataframe(longDF,newDF,Name,Gene):
    #helper method takes DF's creates new row
    #remove extra columns of VCF table
    newDF = newDF[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","AC","AF","AN","BaseQRankSum","DP","ExcessHet","FS","InbreedingCoeff","MLEAC","MLEAF","MQ","MQRankSum","QD","ReadPosRankSum","SOR","ANN","LOF","NMD"]]
    newDF.insert(0,'Name',Name)
    newDF.insert(1,'Gene',Gene)
    longDF = pd.concat([longDF,newDF],ignore_index=True)
    return longDF

def parse_adaptations(variants_table, adaptations, species, output_dir, name):
    """Match variants from the tabularized VCF file to adaptations."""
    logger.info("Parsing table for known adaptations against cefiderocol")
    foundVariants_tsv = os.path.join(output_dir, f"{name}_foundVariants.tsv")
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	AC	AF	AN	BaseQRankSum	DP	ExcessHet	FS	InbreedingCoeff	MLEAC	MLEAF	MQ	MQRankSum	QD	ReadPosRankSum	SOR	ANN	LOF	NMD	TestKP.GT	TestKP.AD	TestKP.DP	TestKP.GQ	TestKP.PL
    variants_df = pd.read_table(variants_table,sep='\t',header=0)
    foundVariant_df = pd.DataFrame(columns=["Name","Gene","CHROM","POS","ID","REF","ALT","QUAL","FILTER","AC","AF","AN","BaseQRankSum","DP","ExcessHet","FS","InbreedingCoeff","MLEAC","MLEAF","MQ","MQRankSum","QD","ReadPosRankSum","SOR","ANN","LOF","NMD"])

    for gene in adaptations["Gene"]:
        adaptation = adaptations.loc[adaptations["Gene"] == gene, "Adaptation"].iloc[0]
        gene_df = check_gene(variants_df,gene,adaptation)
        if isNotEmptyDF(gene_df):
            foundVariant_df = AddFilterDataframe(foundVariant_df,gene_df,name,gene)
    foundVariant_df.to_csv(foundVariants_tsv,sep="\t")
    return foundVariant_df

def check_adaptations(variants_df,adaptations,species,output_dir, name, minDP):
    """Check variants that matched adapatations for reinstating frameshifts"""
    logger.info("Parsing found adaptations for reinstating frameshifts and other subvariants")
    foundVariants_tsv = os.path.join(output_dir, f"{name}_frameshiftVariants.tsv")
    otherVariant_tsv = os.path.join(output_dir, f"{name}_otherVariants.tsv")
   #	Name	Gene	CHROM	POS	ID	REF	ALT	QUAL	FILTER	AC	AF	AN	BaseQRankSum	DP	ExcessHet	FS	InbreedingCoeff	MLEAC	MLEAF	MQ	MQRankSum	QD	ReadPosRankSum	SOR	ANN	LOF	NMD
    foundVariant_df = pd.DataFrame(columns=["Name","Gene","CHROM","POS","ID","REF","ALT","QUAL","FILTER","AC","AF","AN","BaseQRankSum","DP","ExcessHet","FS","InbreedingCoeff","MLEAC","MLEAF","MQ","MQRankSum","QD","ReadPosRankSum","SOR","ANN","LOF","NMD"])
    otherVariant_df = pd.DataFrame(columns=["Name","Gene","CHROM","POS","ID","REF","ALT","QUAL","FILTER","AC","AF","AN","BaseQRankSum","DP","ExcessHet","FS","InbreedingCoeff","MLEAC","MLEAF","MQ","MQRankSum","QD","ReadPosRankSum","SOR","ANN","LOF","NMD"])
    #loop all variants in DB
    for gene in adaptations["Gene"]:
        #all results of current gene
        gene_df = variants_df[variants_df["Gene"] == gene]
        #filter out results that do not meet minDP
        gene_df = gene_df[gene_df.DP > minDP]
        #adaptation type of current gene
        adaptation = adaptations.loc[adaptations["Gene"] == gene, "Adaptation"].iloc[0]
        if isNotEmptyDF(gene_df):
            #frameshift, count ref/alt of all framehifts in gene for reinstating frameshifts
            difference = 0
            if "frameshift_variant" in adaptation:
                for index, row in gene_df.iterrows():
                    len_ref = len(row["REF"])
                    len_alt = len(row["ALT"])
                    if row["ALT"] == ".":
                        len_alt = 0
                    difference += (len_ref - len_alt)
                #frameshift?
                if difference % 3 != 0:
                    foundVariant_df = pd.concat([foundVariant_df,gene_df],ignore_index=True)
                else:
                    logger.info(f"Found reinstating frameshifts for {species} {name} {gene}; removed")


            #missense_variant, pass through
            if "missense_variant" in adaptation or "conservative_inframe_insertion" in adaptation:
                otherVariant_df = pd.concat([otherVariant_df,gene_df],ignore_index=True)
    foundVariant_df.to_csv(foundVariants_tsv,sep="\t")
    otherVariant_df.to_csv(otherVariant_tsv,sep="\t")
    return foundVariant_df


def quality_control(fastq_files, output_dir):
    """Perform quality control using FastQC."""
    logger.info("Running FastQC for quality control")
    run_subprocess(f"fastqc {' '.join(fastq_files)} -o {output_dir}")


def align_reads(fastq_files, ref_genome, output_dir, threads, name):
    """Align reads using BWA-MEM and convert to BAM."""
    logger.info("Aligning reads using BWA MEM")
    sam_file = os.path.join(output_dir, f"{name}_aligned.sam")
    bam_file = os.path.join(output_dir, f"{name}_aligned.bam")
    sorted_bam = os.path.join(output_dir, f"{name}_sorted.bam")
    run_subprocess(f"bwa mem -t {threads} -R '@RG\\tID:{name}\\tSM:{name}\\tLB:01\\tPL:ILLUMINA' {ref_genome} {' '.join(fastq_files)} -o {sam_file}")
    run_subprocess(f"samtools view -@ {threads} -b -o {bam_file} {sam_file}")
    run_subprocess(f"samtools sort -@ {threads} {bam_file} -o {sorted_bam}")
    run_subprocess(f"samtools index -@ {threads} -b {sorted_bam}")
    return sorted_bam

def call_variants(bam_file, ref_genome, output_dir, threads, name, memory):
    """Call variants using GATK HaplotypeCaller."""
    logger.info("Calling variants with GATK HaplotypeCaller")
    vcf_file = os.path.join(output_dir, f"{name}_variants.vcf")
    #testing
    gatk_bamout = os.path.join(output_dir,f"{name}_GATKbamout.bam")
    run_subprocess(f"./gatk --java-options {memory} HaplotypeCaller -ploidy 1 --native-pair-hmm-threads {threads} -R {ref_genome} -I {bam_file} -bamout {gatk_bamout} -O {vcf_file}")
    return vcf_file


def annotate_variants(vcf_file, output_dir, ref_genome, name):
    """Annotate variants using SnpEff."""
    logger.info("Annotating variants with SnpEff")
    annotated_vcf = os.path.join(output_dir, f"{name}_snpEff.vcf")
    run_subprocess(f"snpEff -ud 100 -noStats {ref_genome} {vcf_file} > {annotated_vcf}")
    return annotated_vcf

def convert_vcf_to_table(vcf_file, output_dir, name):
    """Convert VCF file to a tabular format using GATK VariantsToTable."""
    logger.info("convert vcf to table")
    output_table = os.path.join(output_dir, f"{name}_ann_norm_vcf.tsv")
    run_subprocess(f"./gatk VariantsToTable -V {vcf_file} -O {output_table}")
    return output_table

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="CefiderocolFinder: Variant Discovery Tool")
    parser.add_argument("--config", type=str, help="Path to the configuration file.")
    parser.add_argument("--output", type=str, required=True, help="Directory to store output files.")
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate files.")
    parser.add_argument("--species", type=str, required=True, help="Species to filter adaptations.")
    parser.add_argument("--adaptations", type=str, help="Path to adaptations file.")
    parser.add_argument("--genome-ref", type=str, help="Path to reference genome file.")
    parser.add_argument("--reads", type=str, nargs='+', required=True, help="List of input read files (FASTQ format).")
    parser.add_argument("--threads", type=str, help="# of threads/cores to use when processing; default 8")
    parser.add_argument("--name", type=str, required=True, help="Identifier of isolate/file")
    parser.add_argument("--memory", type=str, help="Working memory size for Haplotypecaller (Java: -Xmx24g)")
    parser.add_argument("--fastqc", type=bool, help="Toggle FASTQC")
    parser.add_argument("--minDP", type = str, help="Set the minimum read depth for a variant (DP); default 30")
    return parser.parse_args()

def main():
    args = parse_arguments()

    global logger
    # Setup logging
    log_file = os.path.join(args.output, "CefiderocolFinder.log")
    os.makedirs(args.output, exist_ok=True)
    logger = setup_logger(log_file)

    logger.info("Starting CefiderocolFinder")

    # Read configuration file
    if not args.config:
        #default config
        args.config = "config.yml"
    config = read_config(args.config)
    logger.info("Configuration loaded successfully.")

    #threads
    if not args.threads:
        args.threads = "8"

    #Read genome reference, if not set use config default for supplied spieces
    if not args.genome_ref:
        print(config[args.species])
        args.genome_ref = config[args.species]
        print(f"set genome_ref to {args.genome_ref}")

    # Read adaptations
    if not args.adaptations:
        adaptations_file = os.path.join("referenceDB", "adaptations.tsv") # Default adaptations file
    else:
        adaptations_file = args.adaptations
    adaptations = read_adaptations(adaptations_file, args.species)
    logger.info(f"Found {len(adaptations)} adaptations for species: {args.species}")

    #set working memory of Haplotypecaller
    if not args.memory:
        args.memory = "-Xmx24g"

    #set fastqc
    if not args.fastqc:
        args.fastqc = True

    #set DP
    if not args.minDP:
        args.minDP = 10

    #Perform quality control
    if args.fastqc:
        quality_control(args.reads, args.output)

    #Align reads
    aligned_bam = align_reads(args.reads, args.genome_ref, args.output, args.threads, args.name)
    logger.info(f"Aligned {args.reads} to {args.genome_ref}")

    #Call variants
    variants_vcf = call_variants(aligned_bam, args.genome_ref, args.output, args.threads, args.name, args.memory)

    #Annotate variants
    genome_ref_accession = config[f"{args.species}_accession"]
    annotated_vcf = annotate_variants(variants_vcf, args.output, genome_ref_accession, args.name)

    #Convert to Table
    variants_table = convert_vcf_to_table(annotated_vcf, args.output, args.name)
    logger.info(f"Converted VCF to table: {variants_table}")

    #Parse variants and match to adaptations
    found_adaptations = parse_adaptations(variants_table, adaptations, args.species, args.output, args.name)
    logger.info(f"Matched adaptations to variants")
    print(found_adaptations)

    #check for minDP, reinstating frameshifts, assuming substitutions caused by this are not phenotypically relevant
    filtered_adaptations = check_adaptations(found_adaptations,adaptations,args.species,args.output,args.name, args.minDP)
    logger.info(f"Filtered adaptations for reinstating frameshifts")

    #Report back sucesses and generate success file
    logger.info("Generating finish file")
    finish_file = os.path.join(args.output, f"{args.name}_finish.touch")
    run_subprocess(f"touch {finish_file}")

    # Cleanup intermediates if not keeping them
    if not args.keep_temp:
        delete_intermediates(args.output,args.name)

    logger.info("CefiderocolFinder completed successfully.")

if __name__ == "__main__":
    main()
