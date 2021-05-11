import numpy as np
import math

def bimodalTest(input_hist):
    hist_len=len(input_hist)
    b=False
    modes=0
    for k in range(1, hist_len-1, 1):
        if(input_hist[k-1] < input_hist[k] and input_hist[k+1] < input_hist[k]):
            modes += 1
            if(modes > 2):
                return False;

    if (modes == 2):
        b=True
    return b

def ImageJ_Intermodes_threshold(img):
    # get image histogram, bins are from 0 to 255
    img_px = img.flatten()
    img_hist, bin_edges = np.histogram(img_px, bins=range(0, 256 + 1, 1), )
    hist_len = len(img_hist)
    minbin = -1
    maxbin = -1
    for i in range(hist_len):
        if(img_hist[i] > 0): maxbin = i
    for i in range(hist_len-1,-1, -1):
        if (img_hist[i] > 0): minbin = i

    length = (maxbin - minbin) + 1

    hist = np.zeros(length, dtype='float')
    for i in range(minbin,maxbin+1,1):
        hist[i-minbin] = img_hist[i]

    iter = 0
    threshold = -1
    while (not bimodalTest(hist)):
        # smooth with a 3 point running mean filter
        previous = 0
        current = 0
        next = hist[0]
        for i in range(length-1):
            previous = current
            current = next
            next = hist[i+1]
            hist[i] = (previous+current+next) / 3
        hist[length-1] = (current + next) / 3
        iter+=1
        if(iter > 10000):
            threshold = -1
            print("Intermodes Threshold not found after 10000 iterations.")
            return threshold

    # The threshold is the mean between the 2 peaks
    tt = 0
    for i in range(1,length-1): #(int i=1; i < length - 1; i++)
        if(hist[i-1] < hist[i] and hist[i+1] < hist[i]):
            tt += i
            # print("mode:",i)

    threshold = int(math.floor(tt / 2.0))
    return threshold + minbin

def ImageJ_Li_threshold(img, tolerance=0.5, initial_guess=-1):
    # Implements Li's Minimum Cross Entropy thresholding method
    # This implementation is based on the iterative version (Ref. 2) of the algorithm.
    # 1) Li C.H. and Lee C.K. (1993) "Minimum Cross Entropy Thresholding"
    # Pattern Recognition, 26(4): 617-625
    # 2) Li C.H. and Tam P.K.S. (1998) "An Iterative Algorithm for Minimum
    # Cross Entropy Thresholding"Pattern Recognition Letters, 18(8): 771-776
    # 3) Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding
    # Techniques and Quantitative Performance Evaluation" Journal of
    # Electronic Imaging, 13(1): 146-165
    # http://citeseer.ist.psu.edu/sezgin04survey.html
    # Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
    # Ported to Python by Sarah Keegan

    # get image histogram, bins are from 0 to 255
    img_px=img.flatten()
    img_hist,bin_edges=np.histogram(img_px, bins=range(0,256+1,1), )
    num_pixels = len(img_px)
    mean=np.mean(img_px)
    #tolerance=0.5
    num_pixels=0

    if(initial_guess < 0):
        new_thresh = mean
    else:
        new_thresh=initial_guess

    while(True):
        old_thresh = new_thresh
        threshold = int(old_thresh + 0.5)  # range, rounding to nearest int

        # Calculate the means of background and object pixels
        # Background
        img_px_back = img_px[img_px <= threshold]
        if(len(img_px_back) > 0):
            mean_back = np.mean(img_px_back)
        else:
            mean_back = 0

        # Object
        img_px_obj = img_px[img_px > threshold]
        if (len(img_px_obj) > 0):
            mean_obj = np.mean(img_px_obj)
        else:
            mean_obj = 0
        #print(mean_obj)
        #print(mean_back)
        temp = (mean_back - mean_obj) / (math.log(mean_back) - math.log(mean_obj))
        #print(temp)

        #DBL_EPSILON = 2.220446049250313E-16
        if (temp < -2.220446049250313E-16): # check this!
            new_thresh = int(temp - 0.5)
        else:
            new_thresh = int(temp + 0.5)

        # Stop the iterations when the difference between the
        # new and old threshold values is less than the tolerance
        if(abs(new_thresh - old_thresh) <= tolerance):
            break
    #print(threshold)
    return threshold
