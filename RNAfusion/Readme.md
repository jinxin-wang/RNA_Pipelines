## Install nf-core/rnafusion Pipeline

#### Download RNAfusion package

```
[j_wang@n12 ~]$ mkdir -p lib/nfcore/rnafusion/
[j_wang@n12 ~]$ cd lib/nfcore/rnafusion/
[j_wang@n12 rnafusion]$ wget https://github.com/nf-core/rnafusion/archive/refs/tags/2.1.0.tar.gz
[j_wang@n12 rnafusion]$ tar -xf 2.1.0.tar.gz
[j_wang@n12 rnafusion]$ mv rnafusion-2.1.0 2.1.0
```

#### Modify [base.conf](https://github.com/jinxin-wang/INSERM_U981_Pipelines/blob/main/RNAfusion/base.config)

```
[j_wang@n12 rnafusion]$ cd 2.1.0/conf
[j_wang@n12 conf]$ cp base.conf base.conf.bak
[j_wang@n12 conf]$ emacs base.conf

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnafusion Nextflow base config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A 'blank slate' config file, appropriate for general use on most high performance
    compute environments. Assumes that all software is installed and available on
    the PATH. Runs in `local` mode - all jobs will be run on the logged in environment.
----------------------------------------------------------------------------------------
*/

process {

    cpus   = { check_max( 16, 'cpus'   ) }
    memory = { check_max( 128.GB, 'memory' ) }
    time   = { check_max( 168.h, 'time'   ) }
    queue = 'longq'

    errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_low {
        cpus   = { check_max( 8, 'cpus'   ) }
    	memory = { check_max( 32.GB, 'memory' ) }
        time   = { check_max( 6.h, 'time') }
	queue = 'shortq'
    }
    withLabel:process_medium {
        cpus   = { check_max( 16, 'cpus'   ) }
	memory = { check_max( 64.GB, 'memory' ) }
        time   = { check_max( 24.h, 'time') }
        queue = 'mediumq'	
    }
    withLabel:process_high {
        cpus   = { check_max( 16, 'cpus'   ) }
	memory = { check_max( 128.GB, 'memory' ) }
        time   = { check_max( 168.h, 'time') }
        queue = 'longq'	
    }
    withLabel:process_long {
        cpus   = { check_max( 16, 'cpus'   ) }
	memory = { check_max( 128.GB, 'memory' ) }
        time   = { check_max( 168.h, 'time') }
	queue = 'longq'
    }
    withLabel:process_verylong {
        cpus   = { check_max( 16, 'cpus'   ) }
	memory = { check_max( 128.GB, 'memory' ) }
        time   = { check_max( 1440.h, 'time') }
	queue = 'verylongq'
    }
}

```

#### Set Bash Variable and Export to Global Shell Environment

add the following lines to ~/.bashrc
```
# Our main temporary directory.
BIGR_DEFAULT_TMP="/mnt/beegfs/userdata/${USER}/tmp"

# Create default tmp directory if missing
if [ ! -d "${BIGR_DEFAULT_TMP}" ]; then
    mkdir --parents --verbose "${BIGR_DEFAULT_TMP}";
fi

# Used in many bash / Python scripts
if [ -z ${TMP} ]; then
  declare -x TMP
  TMP="${BIGR_DEFAULT_TMP}"
  export TMP
fi

# Used in some bash / R / perl / Python scripts
if [ -z ${TEMP} ]; then
  declare -x TEMP
  TEMP="${BIGR_DEFAULT_TMP}"
  export TEMP
fi

# Used in some bash / R / perl / Python scripts
if [ -z ${TMPDIR} ]; then
  declare -x TMPDIR
  TMPDIR="${BIGR_DEFAULT_TMP}"
  export TMPDIR
fi

# Used in some bash / R / perl scripts
if [ -z ${TEMPDIR} ]; then
  declare -x TEMPDIR
  TEMPDIR="${BIGR_DEFAULT_TMP}"
  export TEMPDIR
fi

## Added by `nf-core download` v2.7.2 ##
export NXF_SINGULARITY_CACHEDIR="/mnt/beegfs/scratch/${USER}/.singularity_cache"
export NXF_TEMP="/mnt/beegfs/scratch/${USER}/.tmp_dir"
```

#### Build the reference

1. create an account in the [site](https://cancer.sanger.ac.uk/cosmic)
2. put your username and password in the [nextflow.config](https://github.com/jinxin-wang/INSERM_U981_Pipelines/blob/main/RNAfusion/nextflow.config)
 - cosmic_username = your username
 - cosmic_passwd   = your password
3. modify [nextflow.config]
```
[j_wang@n12 rnafusion]$ cd ..
[j_wang@n12 2.1.0]$ mv nextflow.config nextflow.config.bak
[j_wang@n12 2.1.0]$ cp /home/j_wang@intra.igr.fr/RNA_Pipelines/RNAfusion/nextflow.config nextflow.config
```

4. start the pipeline

```
module load java/17.0.4.1
module load singularity/3.6.3
module load nextflow/21.10.6

nextflow run /home/j_wang@intra.igr.fr/lib/nfcore/rnafusion/2.1.0/main.nf --starindex --build_references --genome GRCh38 -profile singularity --outdir $PWD -resume 
```

#### Start RNAfusion pipeline 

```
module load java/17.0.4.1
module load singularity/3.6.3
module load nextflow/21.10.6

nextflow run /home/j_wang@intra.igr.fr/lib/nfcore/rnafusion/2.1.0/main.nf --starindex --genome GRCh38 -profile singularity --outdir $PWD -resume --arriba --star_fusion
```
