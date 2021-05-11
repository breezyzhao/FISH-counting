# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 11:50:49 2016

@author: sarahkeegan
"""

import numpy as np
from scipy import ndimage
import mahotas as mh
from skimage.feature import peak_local_max
from skimage import io, filters, exposure, measure, segmentation, img_as_ubyte, img_as_uint
import matplotlib.pyplot as plt 
import warnings
from skimage.morphology import disk, binary_erosion
from skimage.filters import rank
from skimage.exposure import equalize_adapthist, rescale_intensity
import ij_thresholds as ij_th

warnings.filterwarnings('ignore', message='.+ is a low contrast image', category=UserWarning)

class identify_nuclei_class:
    """given an image (greyscale), of nuclei, identifies nuclei, separates with watershed,
    returns nuceli count, nuclei total area, and labeled mask of nuceli"""
           
    def __init__(self):

        self.save_dir=''
        self.labeling_structure = [[1,1,1],[1,1,1],[1,1,1]] #[[0,1,0],[1,1,1],[0,1,0]] #[[1,1,1],[1,1,1],[1,1,1]]
        self.filter_type='gauss' #no filter if blank, else 'gauss' or 'median'
        self.filter_radius=2
        self.contrast_enh_type='' #no contrast enhance if blank, else 'equalize_adapthist' or 'rescale_percentile'
        self.contrast_enh_kernel_size=50
        self.contrast_enh_clip_limit=0.02
        self.contrast_enh_rescale_perc=0.1
        self.min_area=0
        self.max_area=-1
        self.min_solidity=0
        self.watershed=True
        self.ws_use_erosion=False
        self.ws_filter_mask = False
        self.ws_filter_mask_size = 2
        self.ws_gauss_filter_size = 10
        self.ws_local_min_distance = 10
        self.remove_edge=False
        self.threshold_type='otsu' #'li' or 'otsu' or 'yen' or 'local_otsu'
        self.local_th_disk_size=80
        self.two_pass_th=False
        self.two_pass_ratio_lim=0.75

        self.reset_proc_params()


    def reset_proc_params(self):
        self.nuclei_image = []
        self.nuclei_mask = []
        self.nuclei_count = 0
        self.nuclei_total_area = 0
        self.labeled_clusters = []
        self.nuclei_props = []
        self.outline_img = []
        self.proc_nuclei_image = []
        self.full_background_mask = []
        self.orig_image_mask = []


    def do_rescale(self, img, perc):
        if(perc > 0):
            ll = perc
            ul = 100 - perc
            pmin, pmax = np.percentile(img, (ll, ul))
            return exposure.rescale_intensity(img, in_range=(int(pmin), int(pmax)))
        else:
            return img

    def filter_min_area(self, props, mask, min_area):
        for prop in props:
            if(prop.area < min_area):
                for c in prop.coords:
                    #set coordinates of this segment to False in the image mask
                    mask[c[0],c[1]] = False

    def filter_max_area(self, props, mask, max_area):
        for prop in props:
            if(prop.area > max_area):
                for c in prop.coords:
                    #set coordinates of this segment to False in the image mask
                    mask[c[0],c[1]] = False

    def filter_min_solidity(self, props, mask, min_solidity):
        for prop in props:
            if(prop.solidity < min_solidity):
                for c in prop.coords:
                    #set coordinates of this segment to False in the image mask
                    mask[c[0],c[1]] = False

    def save_maxima(self, distances_, labeled_clusters_, min_dist, save_file):
        maxima_ = peak_local_max(distances_, min_distance=min_dist, indices=False, exclude_border=False)
        maxima_ = maxima_.astype('uint8') * 255
        labeled_mask_ = labeled_clusters_ > 0
        labeled_display_ = labeled_clusters_ * labeled_mask_
        maxima_marked_ = segmentation.mark_boundaries(maxima_, labeled_display_, color=[1, 0, 0], mode='inner')
        io.imsave(save_file, img_as_ubyte(maxima_marked_))

    def rem_matching_objects(self, mask1, mask2, props, ratio_lim=0.75):
        for prop in props:
            in_mask=0
            coords = prop.coords
            for c in coords:
                if(mask1[c[0], c[1]]):
                    in_mask+=1
                    #in_mask=True
                    #break
            if(in_mask/prop.area > ratio_lim):
                for c in coords:
                    mask2[c[0], c[1]] = False

    def enhance_and_smooth(self, img, img_type, save_dir=''):
        # local histogram equalization
        if (self.contrast_enh_type == 'equalize_adapthist'):

            img = equalize_adapthist(img, kernel_size=self.contrast_enh_kernel_size, clip_limit=self.contrast_enh_clip_limit)

            if (img_type == 'uint16'): img = img_as_uint(img)
            else: img = img_as_ubyte(img)

            if (save_dir): io.imsave(save_dir + '/nuclei_enh1.tif', img)

            displ_img = img.copy()
        elif (self.contrast_enh_type == 'rescale_percentile'):

            img = self.do_rescale(img, self.contrast_enh_rescale_perc)

            if (save_dir): io.imsave(save_dir + '/nuclei_enh2.tif', img)

            displ_img = img.copy()
        else:  # just enhance the display image...
            pmin, pmax = np.percentile(img, (1, 99))
            displ_img = rescale_intensity(img, in_range=(pmin, pmax))

        # image blur
        if (self.filter_type == 'gauss'):
            img = ndimage.gaussian_filter(img, self.filter_radius)
            if (save_dir): io.imsave(save_dir + '/nuclei_filter1.tif', img)
        elif (self.filter_type == 'median'):
            img = filters.rank.median(img, disk(self.filter_radius))
            if (save_dir): io.imsave(save_dir + '/nuclei_filter2.tif', img)

        return img, displ_img

    def process_image(self, nuclei_image, save_dir=''):
        #segment the nuclei

        img_type = nuclei_image.dtype

        self.reset_proc_params()

        self.nuclei_image = nuclei_image.copy()

        if (save_dir):
            io.imsave(save_dir + '/nuclei_orig.tif', self.nuclei_image)

        nuclei_image, displ_img = self.enhance_and_smooth(nuclei_image, img_type, save_dir)

        if (self.threshold_type == 'local_otsu'):
            local_otsu = filters.rank.otsu(nuclei_image, disk(self.local_th_disk_size))
            image_mask = nuclei_image > local_otsu
        elif(self.threshold_type == 'li'):
            image_mask = nuclei_image > ij_th.ImageJ_Li_threshold(nuclei_image) #filters.threshold_li(nuclei_image,tolerance=0.5)
            #print(nuclei_image)
            #print(ij_th.ImageJ_Li_threshold(nuclei_image))
        elif(self.threshold_type == 'yen'):
            image_mask = nuclei_image > filters.threshold_yen(nuclei_image)
        else: #otsu is default
            image_mask = nuclei_image > filters.threshold_otsu(nuclei_image)
        # from sklearn.cluster import KMeans
        # est = KMeans(n_clusters=3)
        # arr = nuclei_image.flatten()
        # arr = np.reshape(arr, (-1,1))
        # est.fit(arr)
        # #est.labels_
        # #est.cluster_centers_
        # classes=[[],[],[]]
        # for i,j in zip(nuclei_image.flatten(), est.labels_):
        #     classes[j].append(i)

        if(self.two_pass_th):

            image_mask = ndimage.binary_fill_holes(image_mask)

            labeled_clusters, num_clusters = ndimage.label(image_mask, structure=self.labeling_structure)
            nuclei_props = measure.regionprops(labeled_clusters)
            if (self.min_area > 0): self.filter_min_area(nuclei_props, image_mask, self.min_area) # ** important to do this **

            if (save_dir):
                io.imsave(save_dir + '/nuclei_thresh1pass.tif', 255 * image_mask.astype('uint8'))

            # remove the detected objects and threshold again
            # use this second threshold instead
            rev_mask = (image_mask == False)

            nuclei_image_rem = nuclei_image * rev_mask.astype(nuclei_image.dtype)

            if (save_dir): io.imsave(save_dir + '/nuclei_image_rem.tif', nuclei_image_rem)

            if (self.threshold_type == 'local_otsu'):
                local_otsu = filters.rank.otsu(nuclei_image_rem, disk(self.local_th_disk_size))
                nuclei_image_rem_mask = nuclei_image > local_otsu
            elif (self.threshold_type == 'li'):
                li = ij_th.ImageJ_Li_threshold(nuclei_image_rem) #filters.threshold_li(nuclei_image_rem,tolerance=0.5)
                nuclei_image_rem_mask = nuclei_image > li
            elif (self.threshold_type == 'yen'):
                yen=filters.threshold_yen(nuclei_image_rem)
                nuclei_image_rem_mask = nuclei_image > yen
            else:  # otsu is default
                global_otsu = filters.threshold_otsu(nuclei_image_rem)
                nuclei_image_rem_mask = nuclei_image > global_otsu

            nuclei_image_rem_mask = ndimage.binary_fill_holes(nuclei_image_rem_mask)

            if (save_dir):
                io.imsave(save_dir + '/nuclei_image_rem_mask_pre.tif', 255 * nuclei_image_rem_mask.astype('uint8'))

            #get rid of objects already found with thresh1
            lc, n = ndimage.label(nuclei_image_rem_mask, structure=self.labeling_structure)
            p = measure.regionprops(lc)
            self.rem_matching_objects(image_mask, nuclei_image_rem_mask, p, self.two_pass_ratio_lim)

            if (save_dir):
                io.imsave(save_dir + '/nuclei_image_rem_mask_post.tif', 255 * nuclei_image_rem_mask.astype('uint8'))

            image_mask = image_mask | nuclei_image_rem_mask

            if (save_dir):
                io.imsave(save_dir + '/nuclei_thresh2pass.tif', 255 * image_mask.astype('uint8'))

        orig_image_mask = image_mask.copy()
        image_mask = ndimage.binary_fill_holes(image_mask)
        
        labeled_clusters, num_clusters = ndimage.label(image_mask, structure=self.labeling_structure)
        nuclei_props = measure.regionprops(labeled_clusters)

        #the full background mask before removing objects < min_area etc.
        self.full_background_mask = (image_mask == False)

        
        # filter by min area (pre-watershed), and if no watershed, then max_area (else will do after ws)
        if(self.min_area > 0 or (self.max_area > 0 and self.watershed==False)):
            if(self.min_area > 0): self.filter_min_area(nuclei_props, image_mask, self.min_area)
            if(self.max_area > 0 and self.watershed==False): self.filter_max_area(nuclei_props, image_mask, self.max_area)

            labeled_clusters, num_clusters = ndimage.label(image_mask, structure=self.labeling_structure)
            nuclei_props = measure.regionprops(labeled_clusters)
            
        if(save_dir):
            io.imsave(save_dir + '/nuclei_thresh.tif', 255*image_mask.astype('uint8'))
            
        if(self.watershed):
            #WATERSHED to SEPARATE TOUCHING NUCLEI

            if (self.ws_use_erosion):
                new_image_mask = binary_erosion(image_mask, )
            else:
                new_image_mask = image_mask

            if (self.ws_filter_mask):
                new_image_mask = filters.rank.median(new_image_mask.astype('uint8'), disk(self.ws_filter_mask_size))
                if (save_dir):
                    io.imsave(save_dir + '/med_filt_mask.tif', 255 * new_image_mask)
                new_image_mask = new_image_mask > 0

            distances = mh.distance(new_image_mask)
            distances = distances/float(distances.ptp()) * 255 #stretch values to full range, don't know if this matters, but helps for visualization

            if(save_dir):
                fig1 = plt.figure()
                DefaultSize = fig1.get_size_inches()
                fig1.set_size_inches( (DefaultSize[0]*2, DefaultSize[1]*2) )
                DPI = fig1.get_dpi()
                fig1.set_dpi(DPI*2)

                self.save_maxima(distances, labeled_clusters, 1, save_dir + '/maxima.tif')

                plt.imshow(distances, cmap='jet')
                fig1.savefig(save_dir + '/distances_raw.pdf')
                fig1.clf()

            #smooth distance map to avoid over-segmentation
            distances = ndimage.gaussian_filter(distances, self.ws_gauss_filter_size)
            if (save_dir):
                plt.imshow(distances, cmap='jet')
                fig1.savefig(save_dir + '/distances_gauss.pdf')
                fig1.clf()

            if (save_dir):
                self.save_maxima(distances, labeled_clusters, 1, save_dir + '/maxima_gauss.tif')
                self.save_maxima(distances, labeled_clusters, self.ws_local_min_distance, save_dir + '/maxima_guass_dist.tif')

            maxima = peak_local_max(distances, min_distance=self.ws_local_min_distance, indices=False, exclude_border=False)
            spots, n = mh.label(maxima)
            surface = distances.max() - distances
            areas = mh.cwatershed(surface, spots)
            
            joined_labels = segmentation.join_segmentations(areas, labeled_clusters)
            labeled_clusters = joined_labels * image_mask
            nuclei_props = measure.regionprops(labeled_clusters)
            
            #re-apply the min area filter in case watershed segmented more small spots
            #also apply max area filter.
            if(self.min_area > 0):
                self.filter_min_area(nuclei_props, labeled_clusters, self.min_area)

            if(self.max_area > 0):
                self.filter_max_area(nuclei_props, labeled_clusters, self.max_area)

        # don't examine nuclei touching the edge
        if (self.remove_edge):
            labeled_clusters = segmentation.clear_border(labeled_clusters)

        #filter for solidity
        if(self.min_solidity>0):
            self.filter_min_solidity(nuclei_props, labeled_clusters, self.min_solidity)

        nuclei_props = measure.regionprops(labeled_clusters, self.nuclei_image)
        
        labeled_mask = labeled_clusters > 0
        labeled_display = labeled_clusters * labeled_mask

        marked = segmentation.mark_boundaries(displ_img, labeled_display, color=[1, 0, 0], mode='inner')
        marked = img_as_ubyte(marked)
            
        if(save_dir):
            io.imsave(save_dir + '/nuclei_outlines.tif', marked)
            #io.imsave(save_dir + '/test.tif', displ_img)
        
        #get nuclei area, count
        self.nuclei_count = len(nuclei_props)
        self.nuclei_total_area = np.sum(labeled_mask.astype('uint8')) #changes True/False to 1/0
        self.labeled_clusters = labeled_clusters.astype('uint16') # 16 bit is enough, i don't think there will be > 65535 nuclei
        self.nuclei_props=nuclei_props
        self.outline_img = marked
        self.display_img = displ_img
        self.proc_nuclei_image = nuclei_image
        self.nuclei_mask = labeled_mask
        self.orig_image_mask = orig_image_mask


#################
run_test=False
if(run_test):
    from skimage import io, exposure, img_as_ubyte, filters

    contrast_enh='rescale_percentile'
    ce_clip_limit=0.01
    ce_rescale_perc=0.1 # 1,2
    filter='median' #gauss
    filter_r=2
    use_ws=True
    ws_gf_size=2 # +1 or 2?
    ws_min_d=4 # 4,6
    th_type='otsu'
    local_th_disk_size=200
    min_area=50

    # filter='median'
    # ce_rescale_perc=0.1
    # ws_min_d = 4
    # ws_gf_size = 3

    nucl_id = identify_nuclei_class()

    nucl_id.contrast_enh_type = contrast_enh
    nucl_id.contrast_enh_rescale_perc = ce_rescale_perc
    nucl_id.contrast_enh_clip_limit = ce_clip_limit  # 0.01 #0.02

    nucl_id.filter_type = filter
    nucl_id.filter_radius = filter_r

    nucl_id.watershed = use_ws
    nucl_id.ws_gauss_filter_size = ws_gf_size
    nucl_id.ws_local_min_distance = ws_min_d

    nucl_id.threshold_type = th_type  # 'otsu' #'local_otsu'
    nucl_id.local_th_disk_size = local_th_disk_size  # (if using local_otsu)

    nucl_id.min_area = min_area

    img = io.imread('/Users/sarahkeegan/Downloads/images/4_KRas_3W0171.tif')
    cur_img = exposure.rescale_intensity(img)

    #####
    nucl_id.ws_use_erosion=True
    nucl_id.process_image(cur_img, '/Users/sarahkeegan/Downloads/images/test/')
    outline_img = nucl_id.outline_img
    io.imsave('/Users/sarahkeegan/Downloads/images/OUTLINE_er.tif', outline_img) # add more smoothing (+2) show with and without erosion
    #####
    nucl_id.ws_use_erosion = False
    nucl_id.process_image(cur_img, '/Users/sarahkeegan/Downloads/images/test/')
    outline_img = nucl_id.outline_img
    io.imsave('/Users/sarahkeegan/Downloads/images/OUTLINE_no-er.tif',outline_img)  # add more smoothing (+2) show with and without erosion







                
                
        