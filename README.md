# ArchR_docker

__THIS IS NOT MEANT FOR PUBLIC USE YET. This repository is pre-release and is not intended for external use but has been made public for various internal reasons. Thanks for your patience!__

To build the singularity container:
Set your working directory to this github repository and run:

`sudo singularity build archr_test.img Singularity`

Where "Singularity" is the file contained in the github repository and "archr_test.img" is the output signularity image file
On Pelayo, normal (non-sudo) users will not be able to do this.

To run the built container:

`singularity run archr_test.img`
