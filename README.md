# Automated Fish counting for KaryoCreate
### The method is for KaryoCreate <br />
Sarah Keegan,Institute for Systems Genetics, NYUMC <br />
Xin Zhao,Institute for Systems Genetics, NYUMC <br />		
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


