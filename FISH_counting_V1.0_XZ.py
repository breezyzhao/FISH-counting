"""
This script will count bright spots that are inside nuclei
For each image in folder:
(1) Identify nuclei - blue channels
(2) Identify spots (by blob detection) - green and red channels
(3) Get count of spots/nuclei
(4) Output summary file
(5) Output each image with nuclei and spots outlined - so these can be examined for correctness

@Author: Sarah Keegan,Institute for Systems Genetics,NYUMC
        Xin Zhao,Institute for Systems Genetics,NYUMC
@version 0.1
"""

import glob
import os
from skimage import io, exposure, feature, img_as_ubyte, segmentation, draw, filters, morphology
import identify_nuclei as idn
import matplotlib.pyplot as plt
import numpy as np
from cv2 import cv2
import pandas as pd
from matplotlib.ticker import PercentFormatter ### new added
np.warnings.filterwarnings('ignore')

def get_files_from_dirs(start_dir):
    # gets list of folders in start_dir
    # opens each, gets files corresponding to green, red and blue channels (in that order)
    file_list=[]
    files_and_folders = os.listdir(start_dir)
    for ff in files_and_folders:
        if(os.path.isdir(start_dir+'/'+ff)):
            cur_file_list=[]
            for i,ext in enumerate(['d1','d2','d3']):
                flist = glob.glob(start_dir + '/' + ff + '/*_D_*'+ext+'.TIF')
                if(len(flist)==0):
                    flist = glob.glob(start_dir + '/' + ff + '/*_D_*'+ext+'.tif')
                if(len(flist)==0):
                    if(ext == 'd1'):
                        old_ext=ext
                        ext='d0'
                        flist = glob.glob(start_dir + '/' + ff + '/*_D_*' + ext + '.TIF')
                        if (len(flist) == 0):
                            flist = glob.glob(start_dir + '/' + ff + '/*_D_*' + ext + '.tif')
                        if (len(flist) == 0):
                            print("Error: Did not find tif/TIF file ending with " + old_ext + "/" + ext + " in ", start_dir + '/' + ff)
                        else:
                            cur_file_list.append(flist[0])
                    else:
                        print("Error: Did not find tif/TIF file ending with "+ext+" in ", start_dir + '/' + ff)
                else:
                    cur_file_list.append(flist[0])
            if(len(cur_file_list) == 3):
                cur_file_list.append(ff)
                file_list.append(cur_file_list)

    return file_list

def do_rescale(img, ll_perc, ul_perc):
    ll = ll_perc
    ul = 100 - ul_perc
    pmin, pmax = np.percentile(img, (ll, ul))
    return exposure.rescale_intensity(img, in_range=(int(pmin), int(pmax)))

def read_parameters_from_file(file_name):
    params = {}

    df = pd.read_table(file_name, sep='\t')
    for row in df.iterrows():
        row = row[1]### change [num] to which parameters you want
        #subfolder = row['Folder_name'] + "_para_" +str(row[1])
        subfolder = row['Folder_name']
        params[subfolder]={}
        for col in df.columns:
            if(col != subfolder):
                params[subfolder][col]=row[col]
    return params

def get_props_and_label(img, props, text_c=[0, 255, 0]):
    areas=[]
    eccentricities=[]
    solidities=[]
    ax_lens=[]
    text_size=0.6
    text_lw=2
    label_w = 260  # note w/h depend on text size/lw, and length of string used to label (this is hard-coded here)
    label_h = 20

    # label each object in outline image with its area, eccentricity and solidity
    for prop in props:
        r, c = prop.centroid
        areas.append(prop.area)
        eccentricities.append(prop.eccentricity)
        solidities.append(prop.solidity)
        ax_lens.append(prop.major_axis_length)
        label = (str(prop.label)+': a=' + str(round(prop.area, 0)) +
                 ' e=' + str(round(prop.eccentricity, 2)) +
                 ' s=' + str(round(prop.solidity, 2)))
        if ((r + label_h) >= len(outline_img)): r_pos = r - ((r + label_h) - len(outline_img))
        else: r_pos = r
        if ((c + label_w) >= len(outline_img[0])): c_pos = c - ((c + label_w) - len(outline_img[0]))
        else: c_pos = c
        cv2.putText(img,label,(int(c_pos),int(r_pos)),cv2.FONT_HERSHEY_SIMPLEX,
                    text_size,text_c,text_lw,cv2.LINE_AA)

    return (img, areas, eccentricities, solidities, ax_lens)

def draw_dot_on_img(draw_img, row, col, color, thickness=0):
    row = int(row)
    col = int(col)
    if(thickness==0):
        draw_img[row][col] = color
    else:
        for cur_r in range(row-thickness, row+thickness+1, 1):
            for cur_c in range(col-thickness, col+thickness+1,1):
                if(cur_r >=0 and cur_c >= 0 and cur_r < len(draw_img) and cur_c < len(draw_img[0])):
                    draw_img[cur_r][cur_c] = color

def reshape_to_rgb(grey_img):
    #makes single color channel image into rgb
    ret_img = np.zeros(shape=[grey_img.shape[0],grey_img.shape[1],3], dtype='uint8')
    grey_img_=img_as_ubyte(grey_img)

    ret_img[:,:,0]=grey_img_
    ret_img[:, :, 1] = grey_img_
    ret_img[:, :, 2] = grey_img_
    return ret_img

save_extra_images=False # this is for spot detection, set to False unless testing
temp_dir='' # this is for the nuclei detection, leave blank unless testing
do_finding=True #set to True to find spots
plot_results_hists=True #set to True to make plots of the spot-finding results

new_folder_stucture=False  #set to True if analyzing images directly from the microscope folders

#base_dir=input()
#work_dir=input()

base_dir= "/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/15.FISH_spot_finding/raw_data/Roni/Abbot Probes FISH 7-7-2021/Slides/"
work_dir= "/Users/zhaox12/Dropbox (NYU Langone Health)/Xin_backup/Teresa_lab/project/15.FISH_spot_finding/test/test"

params = read_parameters_from_file(base_dir + '/spot_counting_parameters.txt')
folders=list(params.keys())
if(do_finding):
    spot_df=pd.DataFrame()

    for folder_i, folder in enumerate(folders):
        if(not os.path.isdir(base_dir + '/' + folder)):
            print(folder, "not found, skipping...")
            continue
        print(folder)

        #set up output directory
        if (not os.path.exists(work_dir + '/' + folder)):
            os.mkdir(work_dir + '/' + folder)
        cur_work_dir = work_dir + '/' + folder

        blob_th_GFP = float(params[folder]['blob_th_GFP'])
        blob_th_RFP = float(params[folder]['blob_th_RFP'])
        spot_dist = int(params[folder]['spot_distance_cutoff'])
        if(params[folder]['nucl_id_contrast_enh_type'] == 'none'):
            DAPI_ce_perc = 0
        else:
            DAPI_ce_perc = float(params[folder]['nucl_id_ce_percentile'])

        rescale_intensity_perc_GFP = [params[folder]['GFP_ce_percentile_ll'],params[folder]['GFP_ce_percentile_ul']]
        rescale_intensity_perc_RFP = [params[folder]['RFP_ce_percentile_ll'],params[folder]['RFP_ce_percentile_ul']]
        # new added -211006
        if('count_from_0' in params[folder]):
            count_from_0 = int(params[folder]['count_from_0'])
        else:
            count_from_0 = 0

        #open 3 files, one for each channel
        if(new_folder_stucture):
            file_list  = get_files_from_dirs(base_dir + '/' + folder)  # [ [g,r,b], [g,r,b], ...]
            if (len(file_list) == 0):
                print("No files found in folder:",folder)
                continue
        else:
            file_list = glob.glob(base_dir + '/' + folder + '/*_DAPI.tif')
            if(len(file_list) == 0):
                print("No files found in folder.")
                continue
        all_nucl_areas=[]
        all_nucl_eccentricities=[]
        all_nucl_solidities=[]
        all_nucl_maj_ax_lens=[]
        for file_ in file_list:
            if(new_folder_stucture):
                DAPI_file = os.path.split(file_[2])[1]
                GFP_file = os.path.split(file_[0])[1]
                RFP_file = os.path.split(file_[1])[1]
                nucl_outline_file = DAPI_file[:-4] + '_nucl.tif'
                fnames = [GFP_file[:-4], RFP_file[:-4]]
                in_folder = file_[3]
            else:
                DAPI_file = os.path.split(file_)[1]
                base_file_name = DAPI_file[:-9]
                GFP_file = base_file_name+'_GFP.tif'
                RFP_file = base_file_name + '_RFP.tif'
                nucl_outline_file = base_file_name + '_nucl.tif'
                fnames = [base_file_name, base_file_name]
                in_folder=''

            #load files - each file has 3 channels (2 are blank)
            try:

                DAPI_img = io.imread(base_dir + '/' + folder + '/' + in_folder + '/' + DAPI_file)[:,:,2]
                GFP_img = io.imread(base_dir + '/' + folder + '/' + in_folder + '/' + GFP_file)[:,:,1]
                RFP_img = io.imread(base_dir + '/' + folder + '/' + in_folder + '/' + RFP_file)[:,:,0]
            except FileNotFoundError as e:
                print("GFP/RFP files not found for '",DAPI_file, "' in folder '", folder, "'")
                continue

            if (save_extra_images):
                if (not os.path.exists(cur_work_dir + '/nuclei')):
                    os.mkdir(cur_work_dir + '/nuclei')

            DAPI_img = exposure.rescale_intensity(DAPI_img)
            #GFP_img = exposure.rescale_intensity(GFP_img)
            #RFP_img = exposure.rescale_intensity(RFP_img)

            #identify nuclei in DAPI - setup params
            nucl_id = idn.identify_nuclei_class()
            if(params[folder]['nucl_id_contrast_enh_type'] == 'none'):
                nucl_id.contrast_enh_type = ''
            else:
                nucl_id.contrast_enh_type = params[folder]['nucl_id_contrast_enh_type']

            nucl_id.contrast_enh_rescale_perc = params[folder]['nucl_id_ce_percentile']
            nucl_id.filter_type = 'median'
            nucl_id.filter_radius = params[folder]['nucl_id_med_filter_size']
            nucl_id.ws_use_erosion = False

            nucl_id.watershed = bool(params[folder]['nucl_id_watershed'])
            nucl_id.ws_gauss_filter_size = float(params[folder]['nucl_id_ws_gauss_sigma'])
            nucl_id.ws_local_min_distance = float(params[folder]['nucl_id_ws_min_dist'])
            nucl_id.min_solidity = float(params[folder]['nucl_id_min_solidity'])
            nucl_id.min_area = float(params[folder]['nucl_id_min_area'])
            nucl_id.max_area = float(params[folder]['nucl_id_max_area'])
            nucl_id.threshold_type = params[folder]['nucl_id_th']
            nucl_id.remove_edge = True

            nucl_id.process_image(DAPI_img, temp_dir)
            outline_img = np.copy(nucl_id.outline_img)

            io.imsave(cur_work_dir + '/' + nucl_outline_file[:-4]+'_pre-filtered_mask_'+str(DAPI_ce_perc)+'.tif',
                      nucl_id.orig_image_mask.astype('uint8')*255)

            # label each object in outline image with its area, eccentricity and solidity
            (outline_img, areas, eccentricities, solidities, maj_ax_lens) = get_props_and_label(outline_img, nucl_id.nuclei_props)
            all_nucl_areas.extend(areas)
            all_nucl_eccentricities.extend(eccentricities)
            all_nucl_solidities.extend(solidities)
            all_nucl_maj_ax_lens.extend(maj_ax_lens)

            io.imsave(cur_work_dir + '/' + nucl_outline_file[:-4]+'_'+str(DAPI_ce_perc)+'.tif', outline_img)

            #identify spots
            #for each identified nuclei, crop to bounding box
            labels=['GFP','RFP']
            spot_list = []
            for meas_img_i,meas_img in enumerate([GFP_img,RFP_img]):
                meas_img_displ = segmentation.mark_boundaries(meas_img, (nucl_id.nuclei_mask*nucl_id.labeled_clusters),
                                                              color=[1, 0, 0], mode='inner')
                meas_img_displ = img_as_ubyte(meas_img_displ)
                meas_img_displ2 = img_as_ubyte(meas_img_displ).copy()

                #blob th can be different for GFP/RFP
                if(meas_img_i == 0):
                    blob_th=blob_th_GFP
                    rescale_intensity_perc = rescale_intensity_perc_GFP
                else:
                    blob_th=blob_th_RFP
                    rescale_intensity_perc = rescale_intensity_perc_RFP

                for prop in nucl_id.nuclei_props:
                    nucl_spot_count = 0
                    (min_row, min_col, max_row, max_col) = prop.bbox
                    obj_mask = prop.image
                    cur_img = meas_img[min_row:max_row,min_col:max_col]
                    nucl_y,nucl_x = prop.centroid

                    cur_img_displ = reshape_to_rgb(cur_img)
                    if(save_extra_images):
                        io.imsave(cur_work_dir + '/nuclei/orig_' + nucl_outline_file[:-4] + "_" + str(prop.label) +
                                  '_' + labels[meas_img_i] + '.tif', cur_img_displ)

                    cur_img = do_rescale(cur_img, rescale_intensity_perc[0],rescale_intensity_perc[1])
                    cur_img_displ = reshape_to_rgb(cur_img)
                    if (save_extra_images):
                        io.imsave(cur_work_dir + '/nuclei/rescale_int_' + nucl_outline_file[:-4] + "_" + str(prop.label) +
                            '_' + labels[meas_img_i] + '.tif', cur_img_displ)

                    if(params[folder]['white_tophat']):
                        # apply white top hat to subtract background/spots below a minimum size
                        cur_img = morphology.white_tophat(cur_img, selem=morphology.disk(params[folder]['tophat_disk_size']), )
                        cur_img = exposure.rescale_intensity(cur_img)
                        cur_img_displ = reshape_to_rgb(cur_img)
                        if (save_extra_images):
                            io.imsave(cur_work_dir + '/nuclei/white_th_' + nucl_outline_file[:-4] + "_" + str(prop.label) +
                                      '_' + labels[meas_img_i] + '.tif', cur_img_displ)

                    # blob detection
                    blobs=feature.blob_log(cur_img, min_sigma=float(params[folder]['blob_min_sigma']),
                                           max_sigma=float(params[folder]['blob_max_sigma']),
                                           num_sigma=int(params[folder]['blob_num_sigma']),
                                           threshold=blob_th,
                                           overlap=float(params[folder]['blob_overlap']),
                                           exclude_border=0) ## setting exclude_border to > 0 finds NO blobs ???
                    # The radius of each blob is approximately sq-root of 2 * sigma
                    nucl_spots=[]
                    for blob_i,blob in enumerate(blobs):
                        radius=(2*blob[2])**0.5
                        if(obj_mask[int(blob[0])][int(blob[1])]):
                            nucl_spot_count += 1

                            nucl_spots.append([fnames[meas_img_i],labels[meas_img_i],prop.label,nucl_x,nucl_y,blob_i,blob[0],blob[1],radius,-1,0])

                            #draw dot on entire image (easier to see)
                            #draw_dot_on_img(meas_img_displ, blob[0]+min_row, blob[1]+min_col, [0, 0, 255], 0)
                            rr, cc = draw.circle_perimeter(int(blob[0]+min_row), int(blob[1]+min_col), int(radius*3),
                                                           method='bresenham', shape=meas_img_displ.shape)
                            meas_img_displ[rr, cc] = [0, 0, 255]

                            #draw circle on nuclei image
                            rr, cc = draw.circle_perimeter(int(blob[0]), int(blob[1]), int(radius),
                                                           method='bresenham', shape=cur_img_displ.shape)
                            cur_img_displ[rr, cc] = [0, 0, 255]

                    #calculate nearest distance each blob to other blob within nuclei
                    drawn_labels=[]
                    for spot_i,spot in enumerate(nucl_spots):
                        min_d = -1
                        min_d_id = -1
                        for spot_j,spot_ in enumerate(nucl_spots):
                            if(spot_i!=spot_j):
                                d=((spot[6]-spot_[6])**2 + (spot[7]-spot_[7])**2)**0.5
                                if((min_d==-1) or (d < min_d)):
                                    min_d=d
                                    min_d_id=spot_[5]
                        nucl_spots[spot_i][9]=min_d_id
                        nucl_spots[spot_i][10] = min_d

                        #draw spot dist to nearest as label
                        if(not ({spot[5],min_d_id} in drawn_labels)):
                            drawn_labels.append({spot[5],min_d_id})
                            cv2.putText(meas_img_displ, str(np.round(min_d,1)), (int(spot[7])+min_col+5, int(spot[6])+min_row+5), cv2.FONT_HERSHEY_SIMPLEX,
                                    0.4, [0,255,0], 1, cv2.LINE_AA)

                    # add to running spot list
                    if(count_from_0 and len(nucl_spots)==0): # add nuclei that have zero spots, if desired #add 211006
                        nucl_spots.append([fnames[meas_img_i],labels[meas_img_i],prop.label,nucl_x,nucl_y,-1,0,0,0,-1,0])

                    spot_list.extend(nucl_spots)
                    if (save_extra_images):
                        io.imsave(cur_work_dir + '/nuclei/final_' + nucl_outline_file[:-4] + "_" + str(prop.label) +
                                  '_' + labels[meas_img_i] + '.tif', cur_img_displ)

                params_txt=str(spot_dist)+'_' +str(blob_th)+'_'+str(rescale_intensity_perc[0])+'_'+\
                           str(rescale_intensity_perc[1])+'_'+str(DAPI_ce_perc)
                if(meas_img_i==0):
                    #if(save_extra_images):
                    if(new_folder_stucture):
                        io.imsave(cur_work_dir + '/' + GFP_file[:-4]+'_GFP_'+params_txt+'_blob_marked.tif', meas_img_displ)
                        io.imsave(cur_work_dir + '/' + GFP_file[:-4]+'_GFP_'+params_txt+'_counts_marked.tif', meas_img_displ2)
                    else:
                        io.imsave(cur_work_dir + '/' + GFP_file[:-4]+'_'+params_txt+'_blob_marked.tif', meas_img_displ)
                        io.imsave(cur_work_dir + '/' + GFP_file[:-4]+'_'+params_txt+'_counts_marked.tif', meas_img_displ2)
                        print(GFP_file[:-4]) #### we need to remove labels[meas_img_i], if GFP_file[:-4] contains GFP
                else:
                    #if (save_extra_images):
                    if(new_folder_stucture):
                        io.imsave(cur_work_dir + '/' + RFP_file[:-4]+'_RFP_'+params_txt+'_blob_marked.tif', meas_img_displ)
                        io.imsave(cur_work_dir + '/' + RFP_file[:-4]+'_RFP_'+params_txt+'_counts_marked.tif', meas_img_displ2)
                    else:
                        io.imsave(cur_work_dir + '/' + RFP_file[:-4] +'_'+params_txt+'_blob_marked.tif', meas_img_displ)
                        io.imsave(cur_work_dir + '/' + RFP_file[:-4] +'_'+params_txt+'_counts_marked.tif', meas_img_displ2)
                        print(RFP_file[:-4]) #### we need to remove labels[meas_img_i], if GFP_file[:-4] contains GFP

            #save spot distances
            spot_cols=['file_name','type','nuclei_label','nucl_x','nucl_y','spot_id','spot_x','spot_y','spot_r','dist_nearest_id','dist_nearest']
            cur_spot_df = pd.DataFrame(spot_list, columns=spot_cols)
            cur_spot_df['folder']=folder
            ext_cols=['folder',]
            ext_cols.extend(spot_cols)
            cur_spot_df = cur_spot_df[ext_cols]
            spot_df = spot_df.append(cur_spot_df)

        #make histogram of nuclei properties:
        plt.hist(all_nucl_areas, bins=np.arange(np.min(all_nucl_areas), np.max(all_nucl_areas), 50))
        plt.savefig(cur_work_dir + '/' + str(DAPI_ce_perc) + '_nuclei_areas_hist.pdf')
        plt.clf()

        plt.hist(all_nucl_eccentricities, bins=np.arange(np.min(all_nucl_eccentricities), np.max(all_nucl_eccentricities), 0.1))
        plt.savefig(cur_work_dir + '/' + str(DAPI_ce_perc)+ '_nuclei_ecc_hist.pdf')
        plt.clf()

        plt.hist(all_nucl_solidities, bins=np.arange(np.min(all_nucl_solidities), np.max(all_nucl_solidities), 0.01))
        plt.savefig(cur_work_dir + '/' + str(DAPI_ce_perc) + '_nuclei_sol_hist.pdf')
        plt.clf()

        plt.hist(all_nucl_maj_ax_lens, bins=np.arange(np.min(all_nucl_maj_ax_lens), np.max(all_nucl_maj_ax_lens), 5))
        plt.savefig(cur_work_dir + '/' + str(DAPI_ce_perc) + '_nuclei_maj_ax_len_hist.pdf')
        plt.clf()

    # ** NOTE: spot_dist and blob_th_GFP/RFP are set separately for each folder
    # if there is > 1 folder listed in spot_counting_paramters.txt, the name of the file below
    # will contain the spot_dist and blob_th_GFP/RFP for the LAST folder that is processed. **
    spot_df.index=range(len(spot_df))
    spot_df.to_csv(work_dir + '/'+str(spot_dist)+'_GFP_' +str(blob_th_GFP)+'_RFP_'+str(blob_th_RFP)+'_'+
                   str(DAPI_ce_perc)+'_all_data.txt', sep='\t')

if(plot_results_hists):
    print("Processing results...")
    #read all data and find valid spots (counting close-by spots as a single spot), get spot counts per nuclei
    folder=folders[len(folders)-1] #file was saved with parameters for last folder in list
    last_blob_th_GFP = float(params[folder]['blob_th_GFP'])
    last_blob_th_RFP = float(params[folder]['blob_th_RFP'])
    last_spot_dist = int(params[folder]['spot_distance_cutoff'])
    if (params[folder]['nucl_id_contrast_enh_type'] == 'none'): last_DAPI_ce_perc = 0
    else: last_DAPI_ce_perc = float(params[folder]['nucl_id_ce_percentile'])
    df = pd.read_csv(work_dir + '/'+str(last_spot_dist)+'_GFP_' +str(last_blob_th_GFP)+'_RFP_'+str(last_blob_th_RFP)+
                     '_'+str(last_DAPI_ce_perc)+'_all_data.txt', sep='\t', index_col=0)

    df['valid_spot']=0
    df['spot_count']=0

    max_spots={}
    max_spots['RFP']=8
    max_spots['GFP']=8

    #get spot counts per nuclei, using spot distance cutoff to combine spots
    spot_counts={}
    for folder_i, folder in enumerate(folders):
        print(folder)
        spot_counts[folder]={}
        cur_df = df[df['folder']==folder]
        spot_dist=int(params[folder]['spot_distance_cutoff'])
        blob_th_GFP = float(params[folder]['blob_th_GFP'])
        blob_th_RFP = float(params[folder]['blob_th_RFP'])
        blob_rescl_perc_GFP = [params[folder]['GFP_ce_percentile_ll'],params[folder]['GFP_ce_percentile_ul']]
        blob_rescl_perc_RFP = [params[folder]['RFP_ce_percentile_ll'],params[folder]['RFP_ce_percentile_ul']]
        if (params[folder]['nucl_id_contrast_enh_type'] == 'none'): DAPI_ce_perc=0
        else: DAPI_ce_perc = float(params[folder]['nucl_id_ce_percentile'])
        if ('count_from_0' in params[folder]):
            count_from_0 = int(params[folder]['count_from_0'])
        else:
            count_from_0 = 0

        nucl_per_spot_count={}
        for type in df['type'].unique(): # type: GFP or RFP
            if(type == 'GFP'):
                blob_th = blob_th_GFP
                rescale_intensity_perc = blob_rescl_perc_GFP
            else:
                blob_th = blob_th_RFP
                rescale_intensity_perc = blob_rescl_perc_RFP
            nucl_per_spot_count[type] = []
            spot_counts[folder][type]=[]
            type_df = cur_df[cur_df['type'] == type]
            for fn in type_df['file_name'].unique():
                fn_df = type_df[type_df['file_name']==fn]
                hist_arr=[0,]*(max_spots[type]+1)

                #open image to label with number of spots on nuclei
                img = io.imread(work_dir + '/' + folder + '/' + str(fn) + '_'+type +'_'+str(spot_dist)+'_' +
                                str(blob_th)+'_'+str(rescale_intensity_perc[0])+'_'+str(rescale_intensity_perc[1])+'_'+
                                str(DAPI_ce_perc)+'_counts_marked.tif')
                for nucl in fn_df['nuclei_label'].unique():
                    cur_spots = fn_df[fn_df['nuclei_label']==nucl]

                    if(len(cur_spots)==1 and cur_spots.iloc[0]['spot_id']==-1):
                        num_valid_spots = 0 # no spots found for this nuclei
                    else:
                        valid_spot_list = []
                        for row_i,row in enumerate(cur_spots.iterrows()):
                            data = row[1]
                            if(row_i > 0):
                                drop_spot=False
                                for valid_spot in valid_spot_list:
                                    d=((data['spot_x']-valid_spot[1])**2+(data['spot_y']-valid_spot[2])**2)**0.5
                                    if(d <= spot_dist):
                                        drop_spot=True
                                        break
                                if(not drop_spot):
                                    valid_spot_list.append((row[0], data['spot_x'], data['spot_y']))
                            else:
                                valid_spot_list.append((row[0], data['spot_x'], data['spot_y']))

                        # set those valid to 1
                        for valid_spot in valid_spot_list:
                            df.at[valid_spot[0], 'valid_spot'] = 1

                        num_valid_spots = len(valid_spot_list)

                    for row_i, row in enumerate(cur_spots.iterrows()):
                        df.at[row[0], 'spot_count'] = num_valid_spots

                    spot_counts[folder][type].append(num_valid_spots)

                    if(num_valid_spots > max_spots[type]):
                        hist_arr[max_spots[type]]+=1
                    else:
                        hist_arr[num_valid_spots]+=1

                    text_size = 0.8
                    text_c=[0,255,0]
                    text_lw = 1
                    cv2.putText(img, str(len(valid_spot_list)), (int(cur_spots.iloc[0]['nucl_x']),
                                int(cur_spots.iloc[0]['nucl_y'])), cv2.FONT_HERSHEY_SIMPLEX,
                                text_size, text_c, text_lw, cv2.LINE_AA)
                hist_arr.insert(0,fn)
                nucl_per_spot_count[type].append(hist_arr)
                io.imsave(work_dir + '/' + folder + '/' + str(fn) + '_' + type +'_'+str(spot_dist)+'_' +str(blob_th)+
                          '_'+str(rescale_intensity_perc[0])+'_'+str(rescale_intensity_perc[1])+ '_'+str(DAPI_ce_perc)+
                          '_counts_marked.tif', img)

        #save hist to csv
        for type in df['type'].unique():
            if (type == 'GFP'):
                blob_th = blob_th_GFP
                rescale_intensity_perc = blob_rescl_perc_GFP
            else:
                blob_th = blob_th_RFP
                rescale_intensity_perc = blob_rescl_perc_RFP
            cols = ['file_name']
            cols.extend(range(0, max_spots[type] + 1))
            df_hist = pd.DataFrame(nucl_per_spot_count[type], columns=cols)
            if(not count_from_0): # drop the zero column
                df_hist.drop(labels=0,axis=1,inplace=True)
            df_hist=df_hist.sort_values(by='file_name')
            df_hist.to_csv(work_dir + '/' + folder + '/' + type+'_'+str(spot_dist)+'_' +str(blob_th)+'_'+
                           str(rescale_intensity_perc[0]) + '_'+str(rescale_intensity_perc[1])+'_'+
                           str(DAPI_ce_perc)+'_counts_by_filename.txt', sep='\t')
    df.to_csv(work_dir + '/'+str(last_spot_dist)+'_GFP_' + str(last_blob_th_GFP) + '_RFP_' + str(last_blob_th_RFP)+'_'+
              str(last_DAPI_ce_perc)+'_all_data_with_counts.txt', sep='\t')

    #plot spot counts per nuclei histogram, for RFP and GFP
    spot_counts_all={}
    for type in df['type'].unique():
        spot_counts_all[type]=[]
    for folder_i, folder in enumerate(folders):
        spot_dist = int(params[folder]['spot_distance_cutoff'])
        if (params[folder]['nucl_id_contrast_enh_type'] == 'none'): DAPI_ce_perc=0
        else: DAPI_ce_perc = float(params[folder]['nucl_id_ce_percentile'])

        blob_th_GFP = float(params[folder]['blob_th_GFP'])
        blob_th_RFP = float(params[folder]['blob_th_RFP'])

        blob_rescl_perc_GFP = [params[folder]['GFP_ce_percentile_ll'], params[folder]['GFP_ce_percentile_ul']]
        blob_rescl_perc_RFP = [params[folder]['RFP_ce_percentile_ll'], params[folder]['RFP_ce_percentile_ul']]
        if ('count_from_0' in params[folder]):
            count_from_0 = int(params[folder]['count_from_0'])
        else:
            count_from_0 = 0

        for type in spot_counts[folder].keys():
            if (type == 'GFP'):
                blob_th = blob_th_GFP
                rescale_intensity_perc=blob_rescl_perc_GFP
            else:
                blob_th = blob_th_RFP
                rescale_intensity_perc = blob_rescl_perc_RFP
            if(len(spot_counts[folder][type])>0):
                if(count_from_0):
                    min_=0
                else:
                    min_=1
                ret=plt.hist(spot_counts[folder][type], np.arange(min_, 10, .25), histtype='step')
                x_freq,y_bins,z_patches=plt.hist(spot_counts[folder][type], np.arange(0, 10, 1), histtype='step')
                plt.savefig(work_dir + '/' + folder + '/spot_count_hist_' + type + '_' + folder.replace('/','-') +
                            '_'+str(spot_dist)+ '_' +str(blob_th)+'_'+str(rescale_intensity_perc[0]) + '_'+
                            str(rescale_intensity_perc[1])+'_'+str(DAPI_ce_perc)+'.pdf')
                plt.clf()
                #percent=plt.hist(spot_counts[folder][type], np.arange(min_, 8, .25),
                #                 weights=np.ones(len(spot_counts[folder][type])) / len(spot_counts[folder][type]))
                #plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
                #plt.ylim(0,1)
                #plt.savefig(work_dir + '/' + folder + '/spot_perc_hist_' + type + '_' + folder.replace('/','-') +
                #            '_'+str(spot_dist)+ '_' +str(blob_th)+'_'+str(rescale_intensity_perc[0]) + '_'+
                #            str(rescale_intensity_perc[1])+'_'+str(DAPI_ce_perc)+'.pdf')
                #plt.clf()
                #spot_counts_all[type].extend(spot_counts[folder][type])
                print(x_freq)
                spot1=pd.DataFrame([x_freq])
                spot1["type"]=type
                spot1["sample"]=folder
                spot1.drop(spot1.columns[[0]], axis=1) 
                spot1.to_csv(work_dir + '/count_summary.txt', sep='\t',mode='a', header=True,index=False)
                freq, bins, patches=plt.hist(spot_counts[folder][type], np.arange(0,10,1),
                                 weights=np.ones(len(spot_counts[folder][type])) / len(spot_counts[folder][type]))
                plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
                bin_centers = np.diff(bins)*0.5 + bins[:-1]
                n = 0
                percent = []
                df_list=[]
                #f = open(work_dir + "/perc.txt", "a")  
                for fr, x, patch in zip(freq, bin_centers, patches):
                    height = "{:.2f}".format(float(freq[n])*100)
                    perc = float(height)
                    percent.append(perc)
                    plt.annotate(height,
                                xy = (x, float(height)/100),             # top left corner of the histogram bar
                                xytext = (0,0.2),             # offsetting label position above its bar
                                textcoords = "offset points", # Offset (in points) from the *xy* value
                                ha = 'center', va = 'bottom'
                                )
                    n=n+1
                perc1=pd.DataFrame([percent])
                perc1["type"]=type
                perc1["sample"]=folder
                perc1.drop(perc1.columns[[0]], axis=1) 
                print(str(percent))
                #spot1.to_csv(work_dir + '/count_summary.txt', sep='\t',mode='a', header=True,index=False)
                perc1.to_csv(work_dir + '/percentage_summary.txt', sep='\t',mode='a', header=True,index=False)
                #f.write(str([perc1]))
                #f.close()
                    #n = n+1
                #print(perc1)
                plt.ylim(0,1)
                plt.savefig(work_dir + '/' + folder + '/spot_perc_hist_' + type + '_' + folder.replace('/','-') +
                            '_'+str(spot_dist)+ '_' +str(blob_th)+'_'+str(rescale_intensity_perc[0]) + '_'+
                            str(rescale_intensity_perc[1])+'_'+str(DAPI_ce_perc)+'.pdf')
                plt.clf()
                spot_counts_all[type].extend(spot_counts[folder][type])
            else:
                print("No count data for ", folder, " ", type)

    #plot combined histogram for all folders
    # for type in spot_counts_all.keys():
    #     ret=plt.hist(spot_counts_all[type], np.arange(0,8,0.25), histtype='step')
    #     plt.savefig(base_dir + '/spot_count_hist_' + type + '_all_' + str(spot_dist) + '.pdf')
    #     plt.clf()

    #plot for all spots, (not only 'valid' spots), the distance-between-spots histogram
    for folder_i, folder in enumerate(folders):
        spot_dist = int(params[folder]['spot_distance_cutoff'])
        if (params[folder]['nucl_id_contrast_enh_type'] == 'none'): DAPI_ce_perc=0
        else: DAPI_ce_perc = float(params[folder]['nucl_id_ce_percentile'])

        blob_th_GFP = float(params[folder]['blob_th_GFP'])
        blob_th_RFP = float(params[folder]['blob_th_RFP'])
        blob_rescl_perc_GFP = [params[folder]['GFP_ce_percentile_ll'], params[folder]['GFP_ce_percentile_ul']]
        blob_rescl_perc_RFP = [params[folder]['RFP_ce_percentile_ll'], params[folder]['RFP_ce_percentile_ul']]

        cur_df = df[df['folder']==folder]
        for type in df['type'].unique():
            if (type == 'GFP'):
                blob_th = blob_th_GFP
                rescale_intensity_perc = blob_rescl_perc_GFP
            else:
                blob_th = blob_th_RFP
                rescale_intensity_perc = blob_rescl_perc_RFP
            type_df = cur_df[cur_df['type'] == type]
            type_df = type_df[type_df['dist_nearest']>0]
            if(len(type_df)>0):
                ret=plt.hist(type_df['dist_nearest'], bins=np.arange(np.min(type_df['dist_nearest']), np.max(type_df['dist_nearest']), 1))
                plt.savefig(work_dir + '/' + folder + '/spot_dist_hist_'+type+'_'+folder.replace('/', '-') +
                            '_'+str(spot_dist)+'_' +str(blob_th)+'_'+str(rescale_intensity_perc[0]) + '_'+
                            str(rescale_intensity_perc[1])+'_'+str(DAPI_ce_perc)+'.pdf')
                plt.clf()
            else:
                print("No distance data for ", folder, " ", type)



