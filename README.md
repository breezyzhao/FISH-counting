# FISH-counting
Author:Sarah Keegan,Institute for Systems Genetics,NYUMC <br />
        Xin Zhao,Institute for Systems Genetics,NYUMC <br />		
## Parameters:<br />
<img width="856" alt="image" src="https://user-images.githubusercontent.com/50238955/117855191-3e7a4a00-b258-11eb-842f-04cd6b0e7472.png">
bold-fonted parameters mean potential important parameters for algorithm <br />

## How to run it: <br />
1. Generate a spot_counting_parameters.txt (you may see the spot_counting_parameters.txt as a reference) <br />
2. For each FISH imaging sildes, make sure they are named end with ..GFP ..RFP and ..DAPI (this script can only accept TIF file) <br />
<img width="446" alt="image" src="https://user-images.githubusercontent.com/50238955/117856223-5b634d00-b259-11eb-89c4-656e51e3a63f.png">
3. Put the spot_counting_parameters.txt into the slides directory <br />
<img width="850" alt="image" src="https://user-images.githubusercontent.com/50238955/117856460-a1b8ac00-b259-11eb-8fee-301b95aaede8.png">
4. Change the PATH in FISH-counting.py
<img width="1076" alt="image" src="https://user-images.githubusercontent.com/50238955/117856578-c90f7900-b259-11eb-86be-be74b77514cf.png">
5. Run the script

## results examples:

barplots show the percentage for different spots(for example, 10% imaing slides with 1 spot)
<img width="1084" alt="image" src="https://user-images.githubusercontent.com/50238955/117858265-b0a05e00-b25b-11eb-8d27-553098d317e0.png">
<img width="747" alt="image" src="https://user-images.githubusercontent.com/50238955/117858302-bf871080-b25b-11eb-87db-21076f7975b6.png">

