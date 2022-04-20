# README #

Model equation BLUP (meBLUP) is a new method to estimate genetic breeding values (GEBVs) combining relationships form pedigree and genotype information. In comparison to traditional methods, it is not necessary to set up any relationship matrix or its inverse. For further information check the publications listed below.

The current repository includes the scripts developed and used in the Master's thesis of Himmelbauer, 2022 [3].

## Publications

[1] E. **Groeneveld** and A. **Neumaier**. "BLUP without (inverse) relationsship matrix". In: *World Congress on Genetics Applied to Livestock Production.* Electronic Poster Session-Theory to Application 3. Auckland, New Zealand, 2018, p.21. URL: http://www.wcgalp.org/system/files/proceedings/2018/blup-without-inverse-relationship-matrix.pdf 

[2] A. **Neumaier** and E. **Groeneveld**. "Model equations for prediction with pedigree (unpublished work)". 2021

[3] J. **Himmelbauer**. "Parameter estimation in model equations for animal breeding". Master's thesis, University of Vienna, 2022. URL: 

## Description

### meCODE
The folder *meCODE* contains all R-scripts to run meBLUP as it is described in the the Thesis Himmelbauer, 2022 [3]. 
The file 

    meSTART.r

can be used to sstart a simple evaluation, estimating GEBVs, SNP-effects and varaince parameters (e.g. $h^2$) simultaneously. 

The file 

    meSTART_EBV.r

contains the whole appraoch developed and tested in [3]. This appraoch is based on several steps and it can be used to estimate GEBVs, SNP-effects and some varaince parameters (e.g. polygenetic effect), if the heritability ($h^2$) or an estimator for $h^2$ is available for the respective dataset.

In these two start-Scripts several functions are used, wich are defined in the scripts provided in the folder *meMML*. For detailed information please read [3] and the description in the header of the respective R-script. 

To run the start scripts, it is necessary to download the whole folder *meCODE* and save it localy. The start scripts have to be executed in the folder *meCode*. Otherwise the Variable *PATH_meBLUP* in the scripts has to be replaced with the full path to the folder *meMML*.

#### Toy-Example
In the folder *meCODE/toy_data* a very small toy dataset is provided. The two scripts can be tested on this toy dataset with the following code lines, assuming a $h^2$ of 0.3 and using all 4 SNPs. For all other parameters the default values are used. To run a single meBLUP evaluations this script can be used:
    
    R CMD BATCH '--args "toy_data/toy" "0.3" "all" '  meSTART.r

To applay the developed approach to estimate EBVs the following script is used:

    R CMD BATCH '--args "toy_data/toy" "0.3" "all" '  meSTART_EBV.r

### meSIM
The folder *meSIM* contains all R-scripts to simulate a dataset, taht can be used for meBLUP. The simulation was used in [3] to get an independent dataset for test purpose. The script *start.sh* can be used to start the simulation:

    $ sh start.sh

All parameters for the simulation can be changed in *scripts/Sim.r*. The parameters to change are explained in this file.

**NOTE:** The calculation of EBVs is done with the mix99 program suite (MiX99 Development Team. "Solving Large Mixed Model Equations". Jokionene, Finnland, 2019. URL: https://www.luke.fi/en/services/mix99-solving-large-mixed-model-equations). If these programs are not available, the scripts to estimate (G)EBVs for the simulation have to be adapted.  
