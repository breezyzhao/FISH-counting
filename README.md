# Automated Fish counting for KaryoCreate
### The method is for KaryoCreate <br />
Sarah Keegan,Institute for Systems Genetics, NYUMC <br />
Xin Zhao,Institute for Systems Genetics, NYUMC <br />		
## Parameters:<br />
|Parameters| Note|
|----------|-----|
|spot_distance_cutoff| max distance between for 2 spots to be considered 1|
|**NUCLEI DETECTION PARAMETERS:**|
|nucl_id_contrast_enh_type| method to use for contrast enhancement of DAPI image - 2 options: none (no contrast-enhance) or rescale_percentile|
|nucl_id_ce_percentile| the top/bottom percentile of pixels to set to the image maximum/minimum when rescaling pixel values.  Similar to % “saturated pixels” in Imagej’s Contrast Enhance.  *NOTE* if this is set to 0, then it is the same as setting none for nucl_id_contrast_enh_type (no contrast-enhance)|
|nucl_id_med_filter_size| size of median filter to apply prior to nucleus detection|
|nucl_id_watershed 0/1| use watershed to separate touching nuclei|
|nucl_id_ws_gauss_sigma| stdev for  gaussian kernel, gauss filter is applied to distance transform prior to watershed|
|nucli_id_ws_min_dist| peaks in the distance map must be at least this distance away (if closer, only larger peak will be kept)|
|nucl_id_th| thresholding algorithm to define nucleus objects (e.g. “otsu” or “li”)|
|nucl_id_min_solidity| minimum solidity to be included as a nuclei (otherwise object is discarded for the analysis)|
|nucl_id_min_area| minimum area of nuclei, objects smaller will be discarded|
|nucl_id_max_area| maximum area of nuclei, objects larger will be discarded|
|**BLOB DETECTION PARAMETERS:**|
|blob_min_sigma| min.stdev for gaussian kernel.  keep low to detect smaller blobs|
|blob_max_sigma| max.stdev for gaussian kernel|
|blob_num_sigma| number of intermediate values of stdev to consider|
|blob_th_GFP| intensity threshold for GFP image, local maxima smaller than thresh are ignored. Reduce this to detect blobs with less intensities.|
|blob_th_RFP| intensity threshold for RFP image, local maxima smaller than thresh are ignored. Reduce this to detect blobs with less intensities.|
|blob_overlap (0-1)| if area of 2 blobs overlaps by a fraction greater than this threshold, the smaller blob is eliminated|
|white_tophat 0/1| apply white tophat - returns the bright spots of an image that are smaller than the structuring element|
|tophat_disk_size| structuring element is a disk, specify size here|
|GFP_ce_percentile_ll| percentile of pixels saturated from the bottom (low intensity)|
|GFP_ce_percentile_ul| percentile of pixels saturated from the top (high intensity)|
|RFP_ce_percentile_ll| same as GFP|
|GFP_ce_percentile_ul| same as GFP|
|count_0| whether 0 spot will be include, 1 for yes, 0 for no|

## How to run it: <br />
1. Generate a spot_counting_parameters.txt (you may see the spot_counting_parameters.txt as a reference) <br />
2. For each FISH imaging sildes, make sure they are named end with ..GFP ..RFP and ..DAPI (this script can only accept TIF file) <br />
<img width="446" alt="image" src="https://user-images.githubusercontent.com/50238955/117856223-5b634d00-b259-11eb-89c4-656e51e3a63f.png">
3. Put the spot_counting_parameters.txt into the slides directory <br />
<img width="850" alt="image" src="https://user-images.githubusercontent.com/50238955/117856460-a1b8ac00-b259-11eb-8fee-301b95aaede8.png">
4. Change the PATH in FISH-counting.py
<img width="1076" alt="image" src="https://user-images.githubusercontent.com/50238955/117856578-c90f7900-b259-11eb-86be-be74b77514cf.png">
5. Run the script: <br />
Python Fis-counting.py

## Results examples: <br />

The out put will contain count_summary.txt and percentage_summary.txt, which should looks like table below: <br />
| 0(spot) | 1(spot) | 2(spot)  | 3(spot)| 4(spot)| 5(spot) | 6(spot) | 7(spot)| 8(spot) | type | sample  |
|---------|---------|----------|--------|--------|---------|---------|--------|---------|------|---------|
| 0 |10.62|88.5|0.88| 0 | 0 | 0 | 0 | 0 | GFP  | sampleA |
| 0 |5.20 |22.0|72.5| 0 | 0 | 0 | 0 | 0 | RFP  | sampleA |

There will also have some barplots show the percentage for different spots(for example, 10.62% of the imaging only had 1 spot and 88.5% of the imaging had 2 spots) in each sub-directory
<img width="1084" alt="image" src="https://user-images.githubusercontent.com/50238955/117858265-b0a05e00-b25b-11eb-8d27-553098d317e0.png">


