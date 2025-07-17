def get_Escherichia_coli(wildcards):
	return config["Escherichia_coli"][wildcards.sample]

def get_Klebsiella_pneumoniae(wildcards):
	return config["Klebsiella_pneumoniae"][wildcards.sample]

def get_Acinetobacter_baumannii(wildcards):
	return config["Acinetobacter_baumannii"][wildcards.sample]

def get_Pseudomonas_aeruginosa(wildcards):
	return config["Pseudomonas_aeruginosa"][wildcards.sample]

rule all:
  input: 
    expand("results/Escherichia_coli/{sample}/{sample}_finish.touch", sample = config['Escherichia_coli']),
    expand("results/Klebsiella_pneumoniae/{sample}/{sample}_finish.touch", sample = config['Klebsiella_pneumoniae']),
    expand("results/Pseudomonas_aeruginosa/{sample}/{sample}_finish.touch", sample = config['Pseudomonas_aeruginosa']),
    expand("results/Acinetobacter_baumannii/{sample}/{sample}_finish.touch", sample = config['Acinetobacter_baumannii'])

rule CefiderocolFinder_Escherchia_coli:
	input:
		fastq = get_Escherichia_coli
	output:
		outputDir = directory("results/Escherichia_coli/{sample}/"),
		tsv = "results/Escherichia_coli/{sample}/{sample}_frameshiftVariants.tsv",
		finish = "results/Escherichia_coli/{sample}/{sample}_finish.touch"
	resources:
		mem_mb=config['mem_mb']['generic'],
		runtime_min=config['runtime_min']['generic'],
		queue=config['queue']['generic']
	threads:
		config['threads']['generic']
	shell:
		'''
		python main.py --config config_cefiderocolFinder.yml --keep-temp --output {output.outputDir} --species Escherichia_coli --reads {input.fastq} --name {wildcards.sample}
		'''
rule CefiderocolFinder_Pseudomonas_aeruginosa:
	input:
		fastq = get_Pseudomonas_aeruginosa
	output:
		outputDir = directory("results/Pseudomonas_aeruginosa/{sample}/"),
		tsv = "results/Pseudomonas_aeruginosa/{sample}/{sample}_frameshiftVariants.tsv",
		finish = "results/Pseudomonas_aeruginosa/{sample}/{sample}_finish.touch"
	resources:
		mem_mb=config['mem_mb']['generic'],
		runtime_min=config['runtime_min']['generic'],
		queue=config['queue']['generic']
	threads:
		config['threads']['generic']
	shell:
		'''
		python main.py --config config_cefiderocolFinder.yml --keep-temp --output {output.outputDir} --species Pseudomonas_aeruginosa --reads {input.fastq} --name {wildcards.sample}
		'''

rule CefiderocolFinder_Klebsiella_pneumoniae:
	input:
		fastq = get_Klebsiella_pneumoniae
	output:
		outputDir = directory("results/Klebsiella_pneumoniae/{sample}/"),
		tsv = "results/Klebsiella_pneumoniae/{sample}/{sample}_frameshiftVariants.tsv",
		finish = "results/Klebsiella_pneumoniae/{sample}/{sample}_finish.touch"
	resources:
		mem_mb=config['mem_mb']['generic'],
		runtime_min=config['runtime_min']['generic'],
		queue=config['queue']['generic']
	threads:
		config['threads']['generic']
	shell:
		'''
		python main.py --config config_cefiderocolFinder.yml --keep-temp --output {output.outputDir} --species Klebsiella_pneumoniae --reads {input.fastq} --name {wildcards.sample}
		'''


rule CefiderocolFinder_Acinetobacter_baumannii:
	input:
		fastq = get_Acinetobacter_baumannii
	output:
		outputDir = directory("results/Acinetobacter_baumannii/{sample}/"),
		tsv = "results/Acinetobacter_baumannii/{sample}/{sample}_frameshiftVariants.tsv",
		finish = "results/Acinetobacter_baumannii/{sample}/{sample}_finish.touch"
	resources:
		mem_mb=config['mem_mb']['generic'],
		runtime_min=config['runtime_min']['generic'],
		queue=config['queue']['generic']
	threads:
		config['threads']['generic']
	shell:
		'''
		python main.py --config config_cefiderocolFinder.yml --keep-temp --output {output.outputDir} --species Acinetobacter_baumannii --reads {input.fastq} --name {wildcards.sample}
		'''
  
