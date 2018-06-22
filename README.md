# Well decline curve analysis  
This program reads well header data and production logs (e.g. exported from Drilling Info as .csv files) and
walks the user through the genreation of decline curves for each well provided in the input data. Decine curves
are fit with a the hyperbolic curve that is estimated using an iterative least squares method.  

### ENVIRONMENT  
**Python 2.7**  
*Required Libraries*  
* sys  
* os  
* numpy  
* pandas  
* matplotlib version 1.5.3  
* geopy  
* mpl_toolkits.basemap  
  
## 0.0 Provide .csv of Well production data and .csv of well header data
Examples in `./data` 
prepare data by deleting all the commas from the values in the 'Well production data' and 'well header data' files
This can be done easily in excel with find-and-replace.  
  
## 1.0 Select wells based on critieria available in well data  
Follow the prompts after starting program  
  
## 2.0 Compute type curve  
Based on https://www.uky.edu/KGS/emsweb/devsh/production/decline_obj.py  
  
## 3.0 plot type curve  
Results are saved to `./results`  

## 4.0 Plot map of wells used in type curve  
NOTE: this function will only work if shapefiles are provided in a folder located at '../well_decline_curve_shapefiles/\*'  
The shapefiles used in this function are:  
1. map.readshapefile('../well_decline_curve_shapefiles/PB_County/PB_County', 'PB_County', drawbounds = False)  
2. map.readshapefile('../well_decline_curve_shapefiles/PB_County/PB_County', 'PB_State', drawbounds = False)  
3. map.readshapefile('../well_decline_curve_shapefiles/PB_County/PB_County', 'PB_survey', drawbounds = False)  

Mapping results are saved to `./results`  


Open Stories:
3. batch-process multiple GPS locations