# NuMetWG

## Who is this package for?

This package is for scientists experienced with command-line tools who want to analyze bisulfite- and oxidative-bisulfite-treated DNA libraries from NuGEN (Ovation® Ultralow Methyl-Seq Library Systems) or CEGX (TrueMethyl® Whole Genome) Whole Genome kits.

## What does this package do?

  - Map to bisulfite converted genome.
  - Detect and quantify PCR duplicates.
  - If applicable - perform bisulfite conversion controls analysis and reporting.

This package runs [`bismark`](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) and [`cegx_bsexpress`](https://bitbucket.org/cegx-bfx/cegx_bsexpress) to map bisulfite reads to a reference genome and perform spike-in controls analysis (if applicable), respectively. PCR duplication-detection is accomplished with a custom function: `analysisUtils.countPCRdups`, contained herein.

## Installation instructions

Note that installing docker requires administrative privileges.

### Obtain and install docker using system-specific instructions:

  - Docker is a system to deploy software with a reliable run-time environment. Read more here: [What is Docker?](https://www.docker.com/what-docker)
  - Installation instructions can be found here: [Docker install page](https://docs.docker.com/engine/installation/)
  
### Obtain the docker image

  The docker image for `Bridge` analysis is located in the container repository at the address:
    
    gcr.io/nugen-production/github-nugentechnologies-bridge:hg19
    
  Development version of the image will be:
  
    gcr.io/nugen-production/github-nugentechnologies-bridge:hg19
    
### How does it work?

  The container looks for `R1*.fq` files mounted into the `/input` path inside the
  container and runs the NuGEN analysis pipeline and CEGX controls analysis on each file in serial.
  
### Usage

  1) Make a directory containing the files you wish to analyze.  
  2) Run the following command from INSIDE the directory with fastq files:
  
    docker run -it --rm -v $PWD:/input gcr.io/nugen-production/github-nugentechnologies-bridge:hg19

  Please note:
  
  - This only works if rddata is mounted on the docker host. (`nt-srv-virtual-01` or `nugen-gcp-virtual`)
  - Make sure that the fastq files have an "R1" prefix and ".fq" suffix, that is how the script finds files.

## Special cases

To analyze files with a prefix other than `R1`, add `-e PREFIX=<somethingelse>` to the docker call, for example `-e PREFIX=R2`
