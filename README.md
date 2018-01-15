# Well decline curve analysis  
### ENVIRONMENT  
**Python 2.7**  
*Required Libraries*  
*sys  
*os  
*numpy  
*pandas  
*matplotlib version 1.5.3
*geopy  
*mpl_toolkits.basemap  
  
## 0.0 Provide .csv of Well production data and .csv of well header data
Examples in `./data` 
prepare data by:
1 - going into excel and deleting all commas from the cvs file

## 1.0 Select wells based on critieria available in well data  
Follow the prompts after starting program  
  
## 2.0 Compute type curve  
Based on https://www.uky.edu/KGS/emsweb/devsh/production/decline_obj.py  
  
## 3.0 plot type curve  
Results are saved to `./results`  

## 4.0 Plot map of wells used in type curve  
Results are saved to `./results`  


Updates:
1. make new column 
BOE per month (barrels of oil equivalent per day= ((Liquid (bbl) + (Gas (mcf)/6)) / Perforated Interval Length) * 10000 * (1/30.4)
change y-axis to be this new column units

2. import date columns from header data as datetime format

3. batch-process multiple GPS locations