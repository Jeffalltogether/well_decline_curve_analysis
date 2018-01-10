### drilling info analysis

### Boiler-plate imports and code
import sys
sys.path.append('C:/Users/JeffPC2/Documents/Sync/PythonScripts/forJessica/utils/')
import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import pylab
from geopy.distance import vincenty
# import python module for type curve analysis
from decline import DeclineObj

# import tools and custom code
from tools import load_merge_header_and_production_csv, swap_production_dates_for_time_delta, current_selection, decline_curve
from tools import handle_numerical_variables, handle_dateTime_variables, plot_map, decline_curve_for_fitting, fit_decline_curve
# set plot text size
matplotlib.rcParams.update({'font.size': 12})

class Quick_TypeCurve_Analysis(object):
	'''
	Type curve analysis based on Jessica's work.

	Decline curve estimates from a python module available at:
	     http://www.uky.edu/KGS/emsweb/devsh/production/decline_obj.py
	'''

	def __init__(self, headerCSV, productionCSV):
		self.wellDF = load_merge_header_and_production_csv(headerCSV, productionCSV)
		self.userLocation = []

	def subset_wells_by_distance(self):
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		# obtain longitude and latitudes from user
		while  len(self.userLocation) != 2:
			while True:
				try:
					self.userLocation = raw_input('\nDefine the center of your radius in Latitude (WGS84), and Longitude (WGS84) (separate by comma): ')
					self.userLocation = [x.strip() for x in self.userLocation.split(',')]
					self.userLocation = [float(x) for x in self.userLocation]
				except ValueError:
					print 'Please enter numbers'
					continue
				else:
					break

		# obtain the selection radius from user
		while True:
			try:
				userRadius = float(raw_input('\nDefine the radius within which you will keep all nearby wells (in miles): '))
			except ValueError:
				print 'Please enter numbers'
				continue
			else:
				break	

		# add vicintiy column to data set
		dist = np.zeros(len(self.wellDF['API/UWI']))
		for i,(lat,lon) in enumerate(zip(self.wellDF['Surface Latitude (WGS84)'], self.wellDF['Surface Longitude (WGS84)'])):
			dist[i] = vincenty([lat, lon], self.userLocation).miles

		self.wellDF['vicinity'] = dist

		# keep only wells withing the user selected radius
		self.wellDF = self.wellDF.loc[self.wellDF['vicinity'] <= userRadius]

		# notify user of changes to current selection
		print '%i wells selected' %(len(set(self.wellDF['API/UWI'])))

	def subset_by_well_name(self):
		allWells = list(set(self.wellDF['API/UWI']))

		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '\nSelect one or more of the followig wells by API/UWI number\n'
		print 'all wells available...'
		for i,well in enumerate(allWells):
			print '%i -- %s' %(i, well)


		selection = raw_input('well selection [separate by commas]:\n')
		selectionList = [x.strip() for x in selection.split(',')]

		self.wellDF = self.wellDF[self.wellDF['API/UWI'].isin(selectionList)]
		current_selection(self.wellDF)

		# notify user of changes to current selection
		print '%i wells selected' %(len(set(self.wellDF['API/UWI'])))

	def subset_well_by_variable(self):
		allVariables = self.wellDF.columns.values

		print '\nSelect one or more of the followig variables\n'
		print 'all variables available...'
		
		# generate dictionary of variables
		variableDict = dict()
		for i,var in enumerate(allVariables):
			print '%i -- %s' %(i, var)
			variableDict.update({i:var})

		selectedVars = []
		while  len(selectedVars) == 0:
			try:
				selection = raw_input('Select the variables by their number [separate multiple selections by commas]:\n')
				selectionList = [x.strip() for x in selection.split(',')]
				selectedVars = [variableDict.get(int(key)) for key in selectionList]
			except ValueError:
				print 'Please enter variables by their number'
				continue
			else:
				break

		print 'you selected the following variables: '
		print selectedVars

		for colName in selectedVars:
			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
			print '\nthe variable \"%s\" is of type \"%s\"' %(colName, self.wellDF[colName].dtypes)

			if str(self.wellDF[colName].dtypes) in ['float64', 'int64']:
				self.wellDF = handle_numerical_variables(self.wellDF, colName)

			elif str(self.wellDF[colName].dtypes) in ['object']:
				print 'I do not know how to handle objects yet'

			elif str(self.wellDF[colName].dtypes) in ['datetime64', 'timedelta[ns]','datetime64[ns]']:
				self.wellDF = handle_dateTime_variables(self.wellDF, colName)

			else:
				print 'data type not recognized, skipping variable'
				continue

		# notify user of changes to current selection
		print '%i wells selected' %(len(set(self.wellDF['API/UWI'])))

	def plot_decline_curves(self):
		# get time dela column from seleccted wells
		self.wellDF = swap_production_dates_for_time_delta(self.wellDF)
		
		plotDF = self.wellDF[['Time Delta', 'Liquid (bbl)']]

		# plot well data
		fig, ax = plt.subplots(figsize = (12,8))
		for API in set(self.wellDF['API/UWI']):
			plotData = self.wellDF.loc[self.wellDF['API/UWI'] == API, ['Time Delta', 'Liquid (bbl)']]
			days = plotData['Time Delta'].dt.days
			liquid = np.array(plotData['Liquid (bbl)'])
			ax.semilogy(days, liquid, 'o-', label = API)

		# add titles and legend
		ax.set_xlabel('Time [Days]')
		ax.set_ylabel('Rate\n[Liquid (bbl)]')
		plt.legend()

		# save and display plot
		plt.savefig('./results/production_data.png')
		plt.close()

	def generate_type_curve(self):
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		# get time dela column from seleccted wells
		self.wellDF = swap_production_dates_for_time_delta(self.wellDF)

		# decline curve estiamged parameters
		qi, b, di, r2 = fit_decline_curve(self.wellDF)

		# times to estimate for the plot in int(days)
		time_0 = 0
		time_n = np.timedelta64(self.wellDF['Time Delta'].max())
		decline_t = np.arange(time_0, time_n, np.timedelta64(10,'D'))
		decline_t = (decline_t / np.timedelta64(1, 'D')).astype(int)

		# estimated decline curve
		decline_y = decline_curve(decline_t, qi, b, di)

		# plot well data
		fig, ax = plt.subplots(figsize = (15,8))
		for API in set(self.wellDF['API/UWI']):
			plotData = self.wellDF.loc[self.wellDF['API/UWI'] == API, ['Time Delta', 'Liquid (bbl)']]
			days = plotData['Time Delta'].dt.days
			liquid = np.array(plotData['Liquid (bbl)'])
			ax.semilogy(days, liquid, 'o-', label = API)

		# add decline estimate
		ax.plot(decline_t, decline_y, '-', color='black', linewidth=5.0, label = 'Estimated Decline')
		
		# set axis limits
		xmin = (self.wellDF['Time Delta'].min() / np.timedelta64(1, 'D')).astype(int)
		xmin = xmin*0.15
		xmax = (self.wellDF['Time Delta'].max() / np.timedelta64(1, 'D')).astype(int)
		xmax = xmax*1.06
		ax.set_xlim([xmin, xmax])

		# add titles and legend
		ax.set_xlabel('Time [Days]')
		ax.set_ylabel('Rate\n[Liquid (bbl)]')
		ax.set_title('Decline Curve Parameters: qi=%.2f, b=%.4f, di=%.4f, r2=%.3f' %(qi, b, di, r2))
		ax.legend(bbox_to_anchor=(1.28, 1.05))
		
		# Customize the major grid
		ax.grid(which='major', linestyle='-', linewidth='0.5', color='grey')

		# Customize the minor grid
		ax.grid(which='minor', linestyle=':', linewidth='0.5', color='grey')

		# eliminate unnecessary white space
		plt.subplots_adjust(left=0.07, right=0.8, top=0.9, bottom=0.1)

		# save and display plot
		plt.savefig('./results/decline_estimate.png')
		plt.close()

	def map_selected_wells(self):
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print 'generating map, this may take a minute...'
		# send data to mapping function
		if not(self.userLocation):
			plot_map(self.wellDF)
		else:
			plot_map(self.wellDF, self.userLocation)

	def save_selected_data(self):
		# # select by well number
		# wellByName = raw_input ('would you like to save the sorted data as a CSV? [y/n]: ')
		# # check user input
		# while wellByName not in ('y', 'n', 'Y', 'N'):
		# 	wellByName = raw_input('please try again [y/n]?  ')
		# # save wells
		# if wellByName == 'y' or wellByName == 'Y':
		# 	self.wellDF.to_csv('./results/selected_wells.csv')
		self.wellDF.to_csv('./results/selected_wells.csv')



if __name__ == '__main__':
	### well data files
	headerCSV = 'C:/Users/JeffPC2/Documents/Sync/PythonScripts/forJessica/data/Well_header_data.csv'
	productionCSV = 'C:/Users/JeffPC2/Documents/Sync/PythonScripts/forJessica/data/Production_Time_Series.CSV'

	analysis = Quick_TypeCurve_Analysis(headerCSV, productionCSV)

	# select nearby wells with a circular radius
	wellByName = raw_input ('would you like to select wells near a GPS location? [y/n]: ')
	# check user input
	while wellByName not in ('y', 'n', 'Y', 'N'):
		wellByName = raw_input('please try again [y/n]?  ')

	if wellByName == 'y' or wellByName == 'Y':
		analysis.subset_wells_by_distance()


	# # select by well number
	# wellByName = raw_input ('would you like to select individual wells by API-UWI number? [y/n]: ')
	# # check user input
	# while wellByName not in ('y', 'n', 'Y', 'N'):
	# 	wellByName = raw_input('please try again [y/n]?  ')

	# if wellByName == 'y' or wellByName == 'Y':
	# 	analysis.subset_by_well_name()


	# # select by variable ranges
	# wellByVariable = raw_input ('would you like to subset wells by column values? [y/n]: ')
	# # check user input
	# while wellByVariable not in ('y', 'n', 'Y', 'N'):
	# 	wellByVariable = raw_input('please try again [y/n]?  ')

	# if wellByVariable == 'y' or wellByVariable == 'Y':
	# 	analysis.subset_well_by_variable()

	# plot raw decline curves
	# analysis.plot_decline_curves()

	# plot type curve and raw decline curves
	analysis.generate_type_curve()

	# # plot map
	# analysis.map_selected_wells()

	# # save csv
	# analysis.save_selected_data()	

	# analysis.fit_decline_curve()