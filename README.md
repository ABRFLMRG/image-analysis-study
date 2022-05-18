# image-analysis-study

Configuration files for the `Cytopacq` software, used to generate the original images used in the [LMRG Image Analysis Study](https://sites.google.com/view/lmrg-image-analysis-study). 

You can find all the original image files, label (ground truth) files, and PSF files used in the study here: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6560910.svg)](https://doi.org/10.5281/zenodo.6560910)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6560759.svg)](https://doi.org/10.5281/zenodo.6560759)


## Nuclei Images
Nuclei images can be produced using a containerized version of `Cytopacq` here: [https://github.com/tp81/cytopacq](https://github.com/tp81/cytopacq). The page also contains brief instructions. Use the code under "Using the Software," replacing the ini file name with the `hl60` ini file that you wish to use. The nuclei ini files call a specific PSF file that can be downloaded using the link, above. 

Note that the original images used in the study had a distorted PSF due to faulty interpolation of the .tif version of the PSF in Cytopacq. The problem is completely solved by using the .ics version of the PSF. We provide both file versions, above. You can learn more about the issue in [this image.sc thead](https://forum.image.sc/t/3d-image-analysis-tools-and-reproducibility-event-april-27th-via-zoom/64671/8). 

The different ini files refer to different object clustering and dynamic range (~signal-to-noise). In the filenames, `c` indicates clustering (00 = low clustering) and `dr` indicates dynamic range (90 = high dynamic range, high SNR).

## FISH in *C. elegans Images* (punctate, sub-diffraction staining)
FISH images can be produced using [this singularity container](https://s3.msi.umn.edu/umii-tpengo-containers/cytogen-ubuntu-2020-05-11-c72890b5b480.sif). We are currently working to create a codebase version that can then be containerized, but in the meantime, this should work. You can replace the ini file with any of the `celegans` ini files in this repository. The celegans files call a specific PSF files that can be downloaded using the link, above.

Example code
```
singularity exec -B /data/tpengo /home/umii/tpengo/SW/cytogen-ubuntu-2020-05-11-c72890b5b480.sif cytopacq -c celegans_dyn-10_ceff-0.ini -l ./20201105/celegans_dyn-10_ceff-0/celegans_dyn-10_ceff-0_label.ics -f ./20201105/celegans_dyn-10_ceff-0/celegans_dyn-10_ceff-0_final.ics -e ./20201105/celegans_dyn-10_ceff-0/celegans_dyn-10_ceff-0_error.log -w ./20201105/celegans_dyn-10_ceff-0/celegans_dyn-10_ceff-0_warning.log -d /usr/local/plugins
```

The different ini files refer to different object clustering and dynamic range (~signal-to-noise). In the filenames, `ceff` indicates clustering (00 = low clustering) and `dyn` indicates dynamic range (90 = high dynamic range, high SNR).
