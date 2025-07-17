PATH_TO_PIPE="/home/your_path/to_checkout_of_CefiderocolFinder" #location of cefiderocolFinder
SNAKEFILE_PATH="$PATH_TO_PIPE/Snakefile" #where to find the Snakefile
ENV_PREFIX="$PATH_TO_PIPE/env_db" #where to store environment after built
ISOLATES_FILE="$PATH_TO_PIPE/species_config_example.yaml" #config file listing the isolates and their ID  to be analysed

snakemake -s $SNAKEFILE_PATH --rerun-incomplete \
    --configfiles \
    $PATH_TO_PIPE/config_snakemake.yaml \
    $ISOLATES_FILE \
    --conda-prefix $ENV_PREFIX \
    --singularity-prefix $ENV_PREFIX \
    --singularity-args "--bind $PATH_TO_PIPE" \
    --profile $PATH_TO_PIPE/profile/ \
