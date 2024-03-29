# Model based B1 correction of MTsat Maps

**This MATLAB code is designed to generate a correction curve to be able to correct MTsat for B1 errors.**
## Requirements:
- MTsat data (MT-weighted image, low flip angle and high flip angle image)
- relative B1 map (where the prescribed B1 value == 1)
- Load/Save MRI data: https://github.com/ulrikls/niak (If you use other code, you will just need to change the load and save functions)
- For datafitting: polyfitn toolbox https://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn
- and get equation from https://www.mathworks.com/matlabcentral/fileexchange/9577-symbolic-polynomial-manipulation

**Optional code to improve image quality:**
- Denoising Matlab code: https://github.com/sunenj/MP-PCA-Denoising
- Unring code: https://github.com/josephdviviano/unring)
- Note the sample data was already processed with the denoising and unringing code. 

## General Steps:
1. Simulate  the sequence with known timing and saturation parameters. 
    - as in simSeq_M0b_R1obs.m
2. Generate MTsat and R1 values for 1 subject, making sure all flip angles used in the calculation are modified for relative B1. (flip angle achieved = flip angle prescribed at scanner multiplied by the B1map value). 
    - as in sampleCode_calc_M0bappVsR1.m
3. Fit the simulation results to the MTsat dataset to solve for M0B,app
    - as in sampleCode_calc_M0bappVsR1.m
4. Correlate M0B,app with the R1 values used in the fit. 
    - as in sampleCode_calc_M0bappVsR1.m
5. Be sure that the results of step 2 and 4 stored to a fitResults structure
    - as in sampleCode_calc_M0bappVsR1.m
6. Now you are able to correct future MTsat data using this imaging protocol.
    - as in sample_code_correct_MTsat.m

Note that many of the calculations in step 6 are repeated from the previous script. 
This is intended as steps 2-5 should only need to be completed on one dataset
for a given imaging protocol. Then the fit results can be reused for step 6
in subsequent datasets. 

* Note the parts of the code will need to be modified for use.
- Matlab code files are broken up into sections. Generally, code that needs 
to be adjusted is placed towards the top of the section. 
 
**Lines that require user intervention are labeled with:**
"-> USER DEFINED"
Commands that are optional (mainly extra processing step) are labelled with:
"-> USER DEFINED OPTION"

Code included in PCA folder is used to recreate figures from paper, but is not needed for B1 correction of data. 

**If used, please reference the following publication:**
**For additional help with this code, contact christopher.rowley@mcgill.ca**
