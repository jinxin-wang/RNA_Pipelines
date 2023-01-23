## Install nf-core/rnafusion Pipeline

#### Download RNAfusion package

```
[j_wang@n12 ~]$ cd /mnt/beegfs/scratch/j_wang/
[j_wang@n12 j_wang]$ mkdir -p lib/nfcore/rnafusion/
[j_wang@n12 j_wang]$ cd lib/nfcore/rnafusion/
[j_wang@n12 rnafusion]$ wget https://github.com/nf-core/rnafusion/archive/refs/tags/2.1.0.tar.gz
[j_wang@n12 rnafusion]$ tar -xf 2.1.0.tar.gz
[j_wang@n12 rnafusion]$ mv rnafusion-2.1.0 2.1.0
```


#### Build Refereces 
1. create an account in the [site](https://cancer.sanger.ac.uk/cosmic)
2. put your username and password in the [nextflow.config](https://github.com/jinxin-wang/INSERM_U981_Pipelines/blob/main/RNAfusion/nextflow.config)
 - cosmic_username = your username
 - cosmic_passwd   = your password

#### Modify [nextflow.config](https://github.com/jinxin-wang/INSERM_U981_Pipelines/blob/main/RNAfusion/nextflow.config)

```
[j_wang@n12 rnafusion]$ emacs nextflow.config

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

    withLabel:process_low {
        time  = { check_max( 6.h, 'time') }
        queue = 'shortq'	
    }
    withLabel:process_medium {
        time  = { check_max( 24.h, 'time') }
        queue = 'mediumq'
    }
    withLabel:process_long {
        time  = { check_max( 168.h, 'time') }
        queue = 'longq'
    }
    withLabel:process_verylong {
        time   = { check_max( 1440.h, 'time') }
        queue  = 'verylongq'
    }

}
```

#### Set Bash Variable and Export to Global Shell Environment

add the following lines to ~/.bashrc
```
## Added by `nf-core download` v2.7.2 ##
export NXF_SINGULARITY_CACHEDIR="/mnt/beegfs/scratch/j_wang/.singularity_cache"
export NXF_TEMP="/mnt/beegfs/scratch/j_wang/.tmp_dir"
export TMPDIR="/mnt/beegfs/scratch/j_wang/.tmp_dir"
export TEMP="/mnt/beegfs/scratch/j_wang/.tmp_dir"
export TMP="/mnt/beegfs/scratch/j_wang/.tmp_dir"
export TEMPDIR="/mnt/beegfs/scratch/j_wang/.tmp_dir"
```

#### Start the pipeline 
```
cat rnafusion.sh
#!/bin/bash

module load java/17.0.4.1
module load singularity/3.6.3
module load nextflow/21.10.6

nextflow run /mnt/beegfs/scratch/j_wang/lib/nfcore/rnafusion/2.1.0/main.nf --starindex --build_references --all --genome GRCh38 -profile singularity -resume
```
