## PhenoCAT image preparation

Script 'crop_images' generates RGB images (256x256 pix) of nuclei of HeLa cells.

#### Input: 
Raw images, nuclear segmentation, and illumnination correction file from 1 plate of MOTC screen (Berchtold D. et al, Mol Cell 2018). Original (image) data is not published (~3.7 TB).

#### Output:
10 batch folders of each about 62,000 images. Every batch folder contains random images from plate.

#### Functions:
cropNuclei_saveRG_batches.m
IllumCorrect.m
fixNonNumericalValueInImage.m

#### cropNuclei_saveRG_batches.m
Nuclei are selected according to metaData file (contains only IDs of filtered nuclei; nuclei from apoptotic, mitotic and mis-segmented cells are removed, see publication). Nuclei are croppped according to nuclear segmentation (mask), red and green channel are illumination corrected, B/C adjusted, and finally merged and saved as RGB. Final images are 256x256 pixels, background around nuclear segmentation is black.









 