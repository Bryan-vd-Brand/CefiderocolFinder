## <a name="singularity">Singularity Container</a>
There is a singularity container build recipe (.def) with the below guidance on the conda environment and reference genomes integrated. 
* Its recommended to bind the folder with read data and your output directory
* When using singularity exec the working directory must be set to the internal wd in the container

*   singularity exec --cleanenv \
      --bind /path/to/your/readfolder:/reads \
      --bind /path/to/your/outputfolder:/CefiderocolFinderOutput \
      --pwd /opt/CefiderocolFinder \
      CefiderocolFinder.sif \
      CefiderocolFinder \
      --reads /reads/sample_R1.fastq.gz /reads/sample_R2.fastq.gz \
      --output /CefiderocolFinderOutput \
      --species Escherichia_coli \
      --name sample_identifier \
      --config config_cefiderocolFinder.yml

or alternatively through shell when testing

*   singularity shell --cleanenv \
      --bind /path/to/your/readfolder:/reads \
      --bind /path/to/your/outputfolder:/CefiderocolFinderOutput \
      --pwd /opt/CefiderocolFinder \
      CefiderocolFinder.sif

or test the container

* singularity run CefiderocolFinder.sif --help

## <a name="singularity build">Building the Singularity Container</a>
* Clone the git
* Please download the GATK variant caller: https://github.com/broadinstitute/gatk/releases
* Place the gatk-package-4.6.1.0-local.jar in the CefiderocolFinder folder.
* Place the gatkPythonPackageArchive.zip in the folder.
* Place the GATK wrapper file "gatk" in the folder.
* Build the singularity container using the following build command;
* singularity build --fakeroot CefiderocolFinder.sif CefiderocolFinder.def


## <a name="conda">Conda Environment</a>
CefiderocolFinder relies on a conda environment for analysis:
* Please download the GATK variant caller: https://github.com/broadinstitute/gatk/releases
* Place the gatk-package-4.6.1.0-local.jar in the CefiderocolFinder folder.
* Place the gatkPythonPackageArchive.zip in the folder.
* Place the GATK wrapper file "gatk" in the folder.
* Build your environment using the CefiderocolFinder_conda.yml file.

## <a name="RefGenome">Reference Genomes</a>
CefiderocolFinder uses NCBI reference genomes, these have to be added to snpEff.
Find your snpEff folder of your CefiderocolFinder conda environment (/conda/envs/CefiderocolFinder/share/snpeff*) and run the following commands to download and configure reference genomes:
* Acinetobacter_baumannii: ./scripts/buildDbNcbi.sh NZ_CP045110
* Escherichia_coli: ./scripts/biuldDbNcbi.sh NC_000913
* Klebsiella_pneumoniae: ./scripts/buildDbNcbi.sh NZ_CP117227
* Pseudomonas_aeruginosa: ./scripts/buildDbNcbi.sh NC_002516
* The Refseq Identifiers for CRAB, CPEC, CPKP and CPPA were GCF_009035845, GCF_000005845, GCF_028596205 and GCF_000006765 respectively.
* When changing reference genomes aditionally update the chromosome accession in the config_cefiderocolFinder.yml file.

## <a name="cefiderocolFinder">Running CefiderocolFinder</a>
With your conda environment set up and the reference genomes added the the snpEff database (see above) you can now run CefiderocolFinder:
* In your own application through the command line: python main.py --help
* By default FastQC is enabled for read QC; Disable by --fastqc False
* By default intermediary alignment files are removed, keep them by using the --keep-temp parameter
* The found genetic adaptations are output as a table in frameshiftVariants.tsv and missenseVariants.tsv.

## <a name="cefiderocolFinderSM">Running CefiderocolFinder with snakemake</a>
Alternatively using the provided snakemake file and the associated run file: run_pipeline.sh
* When using the run_pipeline.sh, you have to set the PATH_TO_PIPE variable in the file and add snakemake to the environment.
* add your files to species_config_example.yaml.
* set the config file as file to use in run_pipeline.sh.
* in the profile subfolder config.yaml file the HPC command is located, by default bsub.
* Some variables for the bsub request are set in config_snakemake.yaml.

## <a name="species">Supported Species</a>
CefiderocolFinder supports the following species:
* Acinetobacter_baumannii
* Escherichia_coli
* Klebsiella_pneumoniae
* Pseudomonas_aeruginosa

## <a name="adaptations">Adaptations</a>
CefiderocolFinder supports overwriting the adaptations file.
Species indentifier must match with the config file.
Expected format is the following:
* Species Gene    Adaptation Nomenclature 
* TestSpecies TestGene    TestAdaption    TestNomenclature

## <a name="faq">FAQ</a>
* Building the environment failed, the gatkPythonPackageArchive.zip was not found:
* Depending on your conda .bashrc settings your env might be build at a different location, replace the gatkPythonPackageArchive.zip in CefiderocolFinder_conda.yaml with a full path to the file.


## <a name="versions">CefiderocolFinder versions</a>
* V1.0.1: Fixed an issue with multiple adaptations in a single species & gene combination. Added singularity container support. Added stdOUT passthrough for container shell usage. 