### drilling info analysis

### Boiler-plate imports and code
import sys
sys.path.append('./utils/')
import os, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from geopy.distance import vincenty

# import tools and custom code
from tools import load_merge_header_and_production_csv, swap_production_dates_for_time_delta
from tools import current_selection, decline_curve, handle_numerical_variables, handle_dateTime_variables
from tools import handle_object_variables, plot_map, fit_decline_curve, add_BOE_per_day_column, nominal_decline

def main(headerCSV, productionCSV):
	analysis = Quick_TypeCurve_Analysis(headerCSV, productionCSV)
	print '\n********************************************************************************'
	print '*                                                                              *'	
	print '*                       Well Type Curve Analysis                               *'
	print '*                                                                              *'
	print '* Quit this program anytime by pressing `ctrl+C`                               *\n'
	print 'reading well header data from: %s' %headerCSV
	print 'reading production data from: %s' %productionCSV

	# select by well number
	print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
	wellByName = raw_input ('would you like to select individual wells by API-UWI number? [y/n]: ')
	# check user input
	while wellByName not in ('y', 'n', 'Y', 'N'):
		wellByName = raw_input('please try again [y/n]?  ')

	if wellByName == 'y' or wellByName == 'Y':
		analysis.subset_by_well_name()


	# select nearby wells with a circular radius
	print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
	wellByName = raw_input ('would you like to select wells near a GPS location? [y/n]: ')
	# check user input
	while wellByName not in ('y', 'n', 'Y', 'N'):
		wellByName = raw_input('please try again [y/n]?  ')

	if wellByName == 'y' or wellByName == 'Y':
		analysis.subset_wells_by_distance()


	# select by variable ranges
	print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
	wellByVariable = raw_input ('would you like to subset wells by column values? [y/n]: ')
	# check user input
	while wellByVariable not in ('y', 'n', 'Y', 'N'):
		wellByVariable = raw_input('please try again [y/n]?  ')

	if wellByVariable == 'y' or wellByVariable == 'Y':
		analysis.subset_well_by_variable()


	# plot type curve for all selected wells
	print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
	b_value = None
	# determine if user wants to pre-specify any of the decline curve aprameters
	fixed_b = raw_input ('would you like to pre-specify the decline curve b-factor? [y/n]: ')
	# check user input
	while fixed_b not in ('y', 'n', 'Y', 'N'):
		fixed_b = raw_input('please try again [y/n]?  ')

	if fixed_b.upper() == 'Y':
		while True:
			try:
				b_value = float(raw_input('Enter value for b-factor: '))
			except ValueError:
				print 'Please enter a number'
				continue
			else:
				break	

	analysis.generate_type_curve(b_value)

	# plot map
	print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
	analysis.map_selected_wells()

	# save csv
	print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
	analysis.save_selected_data()	

	# plot wells individually
	print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
	analysis.plot_individual_wells_and_type_curves()

	return

class Quick_TypeCurve_Analysis(object):
	'''
	Type curve analysis based on Jessica's work.

	Decline curve estimates from a python module available at:
	     http://www.uky.edu/KGS/emsweb/devsh/production/decline_obj.py
	'''

	def __init__(self, headerCSV, productionCSV):
		self.wellDF = load_merge_header_and_production_csv(headerCSV, productionCSV)
		self.wellDF = add_BOE_per_day_column(self.wellDF)
		self.userLocation = []

	def subset_wells_by_distance(self):
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

		return

	def subset_by_well_name(self):
		allWells = list(set(self.wellDF['API/UWI']))

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

		return

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
				self.wellDF = handle_object_variables(self.wellDF, colName)

			elif str(self.wellDF[colName].dtypes) in ['datetime64', 'timedelta[ns]','datetime64[ns]']:
				self.wellDF = handle_dateTime_variables(self.wellDF, colName)

			else:
				print 'data type not recognized, skipping variable'
				continue

		# notify user of changes to current selection
		print '%i wells selected' %(len(set(self.wellDF['API/UWI'])))

		return

	def generate_type_curve(self, b_value = None):
		# get time dela column from seleccted wells
		self.wellDF = swap_production_dates_for_time_delta(self.wellDF)

		# decline curve estiamged parameters
		qi, b, di, r2 = fit_decline_curve(self.wellDF, fixed_b_factor = b_value)

		d_nominal = nominal_decline(qi, b, di)

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
			plotData = self.wellDF.loc[self.wellDF['API/UWI'] == API, ['Time Delta', 'BOE per day']]
			days = plotData['Time Delta'].dt.days
			liquid = np.array(plotData['BOE per day'])
			ax.semilogy(days, liquid, '-', label = API)

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
		ax.set_ylabel('BOE per Day\n[Barrels of Oil Equivalent per Day]')
		ax.set_title('Decline Curve Parameters: qi=%.2f, b=%.4f, nominal decline rate=%.1f, r2=%.3f' %(qi, b, d_nominal, r2))
		num_col = math.ceil(len(set(self.wellDF['API/UWI']))/40.0) # number of columns to put in legend
		num_col = int(num_col)
		ax.legend(bbox_to_anchor=(1.26, 0.9), ncol = num_col, fontsize = 9-num_col, labelspacing=0.2)
		
		# Customize the major grid
		ax.grid(which='major', linestyle='-', linewidth='0.5', color='grey')

		# Customize the minor grid
		ax.grid(which='minor', linestyle=':', linewidth='0.5', color='grey')

		# eliminate unnecessary white space
		plt.subplots_adjust(left=0.07, right=0.8, top=0.9, bottom=0.1)

		# save and display plot
		plt.savefig('./results/Average_decline_estimate.png')
		plt.close()

		return

	def map_selected_wells(self):
		print 'generating map, this may take a minute...'
		# send data to mapping function
		if not(self.userLocation):
			plot_map(self.wellDF)
		else:
			plot_map(self.wellDF, self.userLocation)

		return

	def save_selected_data(self):
		print 'saving selected wells to .csv'		
		self.wellDF.to_csv('./results/selected_wells.csv')
		return

	def plot_individual_wells_and_type_curves(self):
		print 'generating plots for all selected wells'
		# get time dela column from seleccted wells
		self.wellDF = swap_production_dates_for_time_delta(self.wellDF)
		declineFit = []

		for well in np.unique(self.wellDF['API/UWI']):
			print 'fitting well # %s' %(str(well))
			wellData = self.wellDF[self.wellDF['API/UWI'] == well]

			# decline curve estiamged parameters
			qi, b, di, r2 = fit_decline_curve(wellData)

			# compute Nominal decline
			d_nominal = nominal_decline(qi, b, di)

			# add data to list for saving to excel
			declineFit.append([wellData, qi, b, d_nominal, di, r2])

			# times to estimate for the plot in int(days)
			time_0 = 0
			time_n = np.timedelta64(wellData['Time Delta'].max())
			decline_t = np.arange(time_0, time_n, np.timedelta64(10,'D'))
			decline_t = (decline_t / np.timedelta64(1, 'D')).astype(int)

			# estimated decline curve
			decline_y = decline_curve(decline_t, qi, b, di)

			# plot well data
			fig, ax = plt.subplots(figsize = (15,8))
			days = wellData['Time Delta'].dt.days
			liquid = np.array(wellData['BOE per day'])
			ax.semilogy(days, liquid, 'o-', label = well)

			# add decline estimate
			ax.plot(decline_t, decline_y, '-', color='black', linewidth=5.0, label = 'Estimated Decline')
			
			# set axis limits
			xmin = (wellData['Time Delta'].min() / np.timedelta64(1, 'D')).astype(int)
			xmin = xmin*0.15
			xmax = (wellData['Time Delta'].max() / np.timedelta64(1, 'D')).astype(int)
			xmax = xmax*1.06
			ax.set_xlim([xmin, xmax])

			# add titles and legend
			ax.set_xlabel('Time [Days]')
			ax.set_ylabel('BOE per Day\n[Barrels of Oil Equivalent per Day]')
			ax.set_title('Decline Curve Parameters: qi=%.2f, b=%.4f, nominal decline rate=%.1f, r2=%.3f' %(qi, b, d_nominal, r2))
			ax.legend(bbox_to_anchor=(1.28, 1.05))
			
			# Customize the major grid
			ax.grid(which='major', linestyle='-', linewidth='0.5', color='grey')

			# Customize the minor grid
			ax.grid(which='minor', linestyle=':', linewidth='0.5', color='grey')

			# eliminate unnecessary white space
			plt.subplots_adjust(left=0.07, right=0.8, top=0.9, bottom=0.1)

			# save and display plot
			plt.savefig('./results/' + str(well) + '_decline_estimate.png')
			plt.close()	

		declineFitDF = pd.DataFrame(declineFit, columns = ['API/UWI', 'qi', 'b', 'nominal decline rate', 'effective decline rate[di]', 'r2'])
		declineFitDF.to_csv('./results/individual_well_decline_curves.csv')

		return

if __name__ == '__main__':
	### well data files
	headerCSV = './data/Well_header_data.csv'
	productionCSV = './data/Production_Time_Series.CSV'

	main(headerCSV, productionCSV)
