## helper functions for well analysis

def load_merge_header_and_production_csv(headerCSV, productionCSV):
	import pandas as pd
	import datetime
	# load header data
	headerDF = pd.read_csv(headerCSV, dtype={'API/UWI': 'O'})

	# convert all date column to datetime
	for column in headerDF.columns.values:
		if any(word in column for word in ['date','Date']):
			headerDF[column] = pd.to_datetime(headerDF[column], errors='coerce')

	# force certain columns to be float
	floatColumns = ['BE Oil EUR', 'BE Oil EUR NORM 10kft', 'BE Oil Delta EUR', 'BE Oil Delta EUR Percent', 'BE Gas EUR',
					'BE Oil GAS NORM 10kft', 'BE Gas Delta EUR', 'BE Gas Delta EUR Percent', 'BE B-Factor', 'BE Initial Rate',
					'BE Final Rate', 'BE Initial Decline', 'BE Oil RRR', 'BE Gas RRR', 'Perforated Interval Length']

	for column in headerDF.columns.values:
		if column in floatColumns:
				headerDF[column] = pd.to_numeric(headerDF[column], errors = 'coerce')
	# load time-series data
	timeSeriesDF = pd.read_csv(productionCSV, dtype={'API/UWI': 'O'})

	# drop unneeded columns
	timeSeriesDF = timeSeriesDF.drop(['Entity ID', 'API/UWI List', 'Days'], axis = 1)

	# convert production date to dtype date time and API to string
	#################### loop through all columns with the word "date" and convert to datetime #######################
	timeSeriesDF['Production Date'] = pd.to_datetime(timeSeriesDF['Production Date'])

	# merge header data and timeseries data
	wellDF = pd.merge(timeSeriesDF, headerDF, on='API/UWI', how = 'right')
	wellDF['Production Date'] = pd.to_datetime(wellDF['Production Date'])
	wellDF['API/UWI'] = wellDF['API/UWI'].astype(object)

	print '%i wells available' %(len(set(wellDF['API/UWI'])))

	return wellDF

def add_BOE_per_day_column(wellDF):
	wellDF['BOE per day'] = ( (wellDF['Liquid (bbl)']+(wellDF['Gas (mcf)']/6.0) ) / wellDF['Perforated Interval Length']) * 10000 * (1/30.4)
	return wellDF

def swap_production_dates_for_time_delta(wellDF):
	import pandas as pd

	# generate new column with time delta rather than production date
	wellDF['Time Delta'] = '69 days 00:00:00'
	wellDF['Time Delta'] = pd.to_timedelta(wellDF['Time Delta'])

	for API in set(wellDF['API/UWI']):
		# get peak production value
		peakProduction = wellDF.loc[wellDF['API/UWI'] == API, 'Liquid (bbl)'].max()
		#find date of peak production
		startDate = wellDF.loc[(wellDF['API/UWI'] == API) & (wellDF['Liquid (bbl)'] == peakProduction), 'Production Date'].min()
		#subtract this date from all production dates
		delta = wellDF.loc[wellDF['API/UWI'] == API, 'Production Date'] - startDate
		# update the Time Delta column with this nuber
		wellDF.loc[wellDF['API/UWI'] == API, 'Time Delta'] = delta

	return wellDF

def current_selection(dataFrame):
	import pandas as pd

	print 'You currently have these wells selected:'
	# list current well selection
	allWells = list(set(dataFrame['API/UWI']))
	for i,well in enumerate(allWells):
		print '%i -- %s' %(i, well)

	print '%i wells selected' %(len(allWells))
	return

def handle_numerical_variables(dataFrame, colName):
	# describe the data to the user
	print '\nThe reamining wells in your selection have the following detailes for the column \"%s\" \n(note: mising data is excluded): ' %(colName)
	print dataFrame[colName].dropna().describe()

	# provide subsetting options to user and get their desired option
	print '\nYou have the following options for this variable \"%s\":' %colName
	print "\n  1 -- greater than or equal too [>=]\n  2 -- less than or equal too [<=] \n  3 -- equal too [==]\n  4 -- between x0 and x1 [x0 < variable < x1]"
	selection = raw_input('\nSelect how you would like to subset the data based on the variable \"%s\":' %colName)

	# check user input
	while selection not in ('1', '2', '3', '4'):
		selection = raw_input('please select 1, 2, 3, or 4: ')

	# execute subsetting request 
	# greater than
	if selection == '1':
		while True:
			try:
				criteria = float(raw_input('\nSelect all wells with \"%s\" >= ' %colName))
			except ValueError:
				print 'Please enter a number'
				continue
			else:
				break
		# subset wells				
		dataFrame = dataFrame.loc[dataFrame[colName] >= criteria]

	# less than
	if selection == '2':
		# get user input for threshold value
		while True:
			try:
				criteria = float(raw_input('\nSelect all wells with \"%s\" <= ' %colName))
			except ValueError:
				print 'Please enter a number'
				continue
			else:
				break						
		# subset wells
		dataFrame= dataFrame.loc[dataFrame[colName] <= criteria]

	# equal too
	if selection == '3':
		# get user input for threshold value
		while True:
			try:
				criteria = float(raw_input('\nSelect all wells with \"%s\" == ' %colName))
			except ValueError:
				print 'Please enter a number'
				continue
			else:
				break						
		# subset wells
		criteria = raw_input('\nSelect all wells with \"%s\" == ' %colName)
		dataFrame= dataFrame.loc[dataFrame[colName] == criteria]

	# between x0 and x1
	if selection == '4':
		# get user input for lower bound threshold value
		Limits = []
		while  len(Limits) != 2:
			while True:
				try:
					Limits = raw_input('\nSelect all wells with \"%s\" between the following values: <lower limit, upper limit> ' %colName)
					Limits = [x.strip() for x in Limits.split(',')]
					Limits = [float(x) for x in Limits]
				except ValueError:
					print 'Please enter numbers'
					continue
				else:
					break		
		# subset wells
		dataFrame= dataFrame.loc[(dataFrame[colName] >= float(Limits[0])) & (dataFrame[colName] >= float(Limits[1]))]

	return dataFrame

def handle_dateTime_variables(dataFrame, colName):
	import numpy as np

	# describe the data to the user
	print '\nThe reamining wells in your selection have the following detailes for the column \"%s\" \n(note: mising data is excluded): ' %(colName)
	print dataFrame[colName].dropna().describe()

	# provide subsetting options to user and get their desired option
	print '\nYou have the following options for this variable \"%s\":' %colName
	print "\n  1 -- greater than or equal too [>=]\n  2 -- less than or equal too [<=] \n  3 -- equal too [==]\n  4 -- between x0 and x1 [x0 < variable < x1]"
	selection = raw_input('\nSelect how you would like to subset the data based on the variable \"%s\":' %colName)

	# check user input
	while selection not in ('1', '2', '3', '4'):
		selection = raw_input('please select 1, 2, 3, or 4: ')

	# execute subsetting request 
	# greater than
	if selection == '1':
		while True:
			try:
				criteria = np.datetime64(raw_input('\nSelect all wells with \"%s\" >= (yyyy-mm-dd)' %colName) + 'T00:00:00')
			except ValueError:
				print 'Please enter a date in the correct format'
				continue
			else:
				break
		# subset wells			
		dataFrame = dataFrame.loc[dataFrame[colName] >= criteria]

	# less than
	if selection == '2':
		# get user input for threshold value
		while True:
			try:
				criteria = np.datetime64(raw_input('\nSelect all wells with \"%s\" <= (yyyy-mm-dd)' %colName) + 'T00:00:00')
			except ValueError:
				print 'Please enter a date in the correct format'
				continue
			else:
				break						
		# subset wells
		dataFrame= dataFrame.loc[dataFrame[colName] <= criteria]

	# equal too
	if selection == '3':
		# get user input for threshold value
		while True:
			try:
				criteria = np.datetime64(raw_input('\nSelect all wells with \"%s\" started on date (yyyy-mm-dd) \"%s\" == ' %colName) + 'T00:00:00')
			except ValueError:
				print 'Please enter a date in the correct format'
				continue
			else:
				break						
		# subset wells
		dataFrame= dataFrame.loc[dataFrame[colName] == criteria]

	# between x0 and x1
	if selection == '4':
		# get user input for lower bound threshold value
		Limits = []
		while  len(Limits) != 2:
			while True:
				try:
					Limits = raw_input('\nSelect all wells with \"%s\" between the following dates: <lower limit, upper limit> ' %colName)
					Limits = [x.strip() for x in Limits.split(',')]
					Limits = [np.datetime64(x + 'T00:00:00') for x in Limits]
				except ValueError:
					print 'Please enter a date in the correct format'
					continue
				else:
					break		
		# subset wells
		dataFrame= dataFrame.loc[(dataFrame[colName] >= np.datetime64(Limits[0] + 'T00:00:00')) & (dataFrame[colName] >= np.datetime64(Limits[1] + 'T00:00:00'))]

	return dataFrame

def handle_object_variables(dataFrame, colName):
	# describe the data to the user
	print '\nThe reamining wells in your selection have the following detailes for the column \"%s\" \n(note: mising data is excluded): ' %(colName)
	print dataFrame[colName].dropna().describe()
	print 'All unique values for this variable, \"%s\", are listed below: \n' %(colName)
	uniqueValues = dataFrame[colName].dropna().unique()
	print uniqueValues

	# provide subsetting options to user and get their desired option
	print '\nYou have the following options for this variable \"%s\":' %colName
	print '\n  1 -- inclue your selection \n  2 -- exclude your selection'
	selection = raw_input('\nSelect how you would like to subset the data based on the variable \"%s\":' %colName)

	# check user input
	while selection not in ('1', '2'):
		selection = raw_input('please select 1 or 2: ')

	# execute subsetting request 
	# Include
	if selection == '1':
		criteria = raw_input('\nInput all levels of the variable \"%s\" that you would like to INCLUDE (separated by commas; names are case sensitive; do not include quotations): ' %colName)
		criteria = [x.strip() for x in criteria.split(',')]
		while set(criteria).issubset(uniqueValues) != True:
			print 'Check spelling and case'
			criteria = raw_input('\nInput all levels of the variable \"%s\" that you would like to INCLUDE (separated by commas; names are case sensitive; do not include quotations): ' %colName)
			criteria = [x.strip() for x in criteria.split(',')]

		# subset wells			
		dataFrame = dataFrame[dataFrame[colName].isin(criteria)]

	# Exclude
	if selection == '2':
		criteria = raw_input('\nInput all levels of the variable \"%s\" that you would like to EXCLUDE (separated by commas; names are case sensitive; do not include quotations): ' %colName)
		criteria = [x.strip() for x in criteria.split(',')]
		while set(criteria).issubset(uniqueValues) != True:
			print 'Check spelling and case'
			criteria = raw_input('\nInput all levels of the variable \"%s\" that you would like to EXCLUDE (separated by commas; names are case sensitive; do not include quotations): ' %colName)
			criteria = [x.strip() for x in criteria.split(',')]

		# subset wells			
		dataFrame = dataFrame[~dataFrame[colName].isin(criteria)] #`~` specifies 'not in'

	return dataFrame

def plot_map(latLongDF, userLocation = None):
	# https://peak5390.wordpress.com/2012/12/08/matplotlib-basemap-tutorial-plotting-points-on-a-simple-map/
	# https://stackoverflow.com/questions/23751635/good-python-toolkit-for-plotting-points-on-a-city-map
	from mpl_toolkits.basemap import Basemap
	import matplotlib.pyplot as plt
	import matplotlib.patches as mpatches
	import numpy as np

	labels = latLongDF['API/UWI']
	latitudes = latLongDF['Surface Latitude (WGS84)']
	longitudes = latLongDF['Surface Longitude (WGS84)']
	production = (latLongDF['Cum Oil'] / float(latLongDF['Cum Oil'].max()))*50

	# select map area
	upper_right_lon = np.min(longitudes) - 0.05
	lower_left_lon = np.max(longitudes) + 0.05
	lower_left_lat = np.min(latitudes) - 0.05
	upper_right_lat = np.max(latitudes) + 0.05
	 
	# draw map
	map = Basemap(projection='merc', lat_0 = 57, lon_0 = -135,
		resolution = 'h', area_thresh = 0.1,
		llcrnrlon=lower_left_lon, llcrnrlat=lower_left_lat,
		urcrnrlon=upper_right_lon, urcrnrlat=upper_right_lat)
	 
	map.drawcoastlines()
	map.drawcountries()
	map.fillcontinents(color = 'coral')
	map.drawmapboundary()

	# add roads
	# map.readshapefile('../well_decline_curve_shapefiles/tl_2010_48_prisecroads/tl_2010_48_prisecroads', 'Streets', drawbounds = False)

	# for shape in map.Streets:
	# 	xx, yy, = zip(*shape)
	# 	map.plot(xx, yy, linewidth = 1.5, color='green', alpha=.75) 

	# add PB_block
	map.readshapefile('../well_decline_curve_shapefiles/PB_County/PB_block', 'PB_block', drawbounds = False)

	for shape in map.PB_block:
		xx, yy, = zip(*shape)
		map.plot(xx, yy, linewidth = 1.5, color='black', alpha=.75) 

	# add PB_County
	map.readshapefile('../well_decline_curve_shapefiles/PB_County/PB_County', 'PB_County', drawbounds = False)

	for shape in map.PB_County:
		xx, yy, = zip(*shape)
		map.plot(xx, yy, linewidth = 1.5, color='gray', alpha=.75) 

	# add PB_State
	map.readshapefile('../well_decline_curve_shapefiles/PB_County/PB_County', 'PB_State', drawbounds = False)

	for shape in map.PB_State:
		xx, yy, = zip(*shape)
		map.plot(xx, yy, linewidth = 1.5, color='orange', alpha=.75) 

	# add PB_survey
	map.readshapefile('../well_decline_curve_shapefiles/PB_County/PB_County', 'PB_survey', drawbounds = False)

	for shape in map.PB_survey:
		xx, yy, = zip(*shape)
		map.plot(xx, yy, linewidth = 1.5, color='green', alpha=.75) 


	# add wells to map
	for lab,lat,lon,prod in zip(labels, latitudes, longitudes, production):
		x,y = map(lon, lat)
		map.plot(x, y, 'o', markerfacecolor = 'none', markeredgecolor = 'blue', markersize=prod)
		map.plot(x, y, 'bx', markersize=5)

	# manually enter legend 
	selectedWells = mpatches.Patch(color='blue', label='selected wells')
	block = mpatches.Patch(color='black', label='block')
	survey = mpatches.Patch(color='green', label='survey')
	county = mpatches.Patch(color='gray', label='county')
	
	patches = [selectedWells, block, survey, county]

	if userLocation:
		x,y = map(userLocation[1], userLocation[0])
		map.plot(x, y, 'ro', markersize=20, alpha=0.75)
		map.plot(x, y, 'x', markersize=10, label = 'Selected Location')
		myLocation = mpatches.Patch(color='red', label='my location\n lon: %.5f\n lat: %.5f' %(userLocation[1], userLocation[0]))
		patches.append(myLocation)


	# organize legend and plot locations
	plt.subplots_adjust(left=0.01, right=0.7, top=0.99, bottom=0.01)
	plt.legend(handles=patches, bbox_to_anchor=(1.42, 1.0))

	plt.savefig('./results/map.png')
	plt.close()

	return

def decline_curve(t, qi, b, di):
	# unit of t is in days
	return qi*(1.0-b*di*t)**(-1.0/b)

def fit_decline_curve(wellDF):
	import numpy as np
	from sklearn.metrics import r2_score
	from decline import DeclineObj

	# select only the decline portion of the well's production for analysis
	declineData = wellDF[wellDF['Time Delta'] >= '0 days']
	declineData = declineData.sort_values(['Time Delta'])

	days = declineData['Time Delta'].dt.days
	x_train = np.array(days)
	y_train = np.array(declineData['BOE per day'])

	#start params
	print 'fitting decline curve parameters'
	
	startParams = np.array([5000.0, 0.1, 0.04]) #[qi, b, di]

	# eliminate outliners from curves to fix convergence issues
	smooth_t = []
	smooth_y = []
	for API in set(declineData['API/UWI']):
		# Eliminate data that is below 3 standard deviations of the mean
		stdev = np.std(declineData.loc[(declineData['API/UWI'] == API), 'BOE per day'])
		mean = np.std(declineData.loc[(declineData['API/UWI'] == API), 'BOE per day'])
		goodData = declineData.loc[((declineData['API/UWI'] == API) & (declineData['BOE per day'] >= mean - 3.0*stdev))]
		
		# get decline data for individual well
		y_train = np.array(goodData['BOE per day'])
		t_train = np.array(goodData['Time Delta'].dt.days)

		# append to lists
		smooth_t.append(t_train)
		smooth_y.append(y_train)			

	# flatten list
	smooth_t = np.array([j for i in smooth_t for j in i])
	smooth_y = np.array([j for i in smooth_y for j in i])

	# sort from smallest times to largest times
	smooth = zip(smooth_t, smooth_y)
	smooth = sorted(smooth)

	# separate into x and y
	x = np.array([t for t,y in smooth])
	y = np.array([y for t,y in smooth])

	# compute decline curve
	model = DeclineObj(x,y,[y.max(),1.4,-0.75], model="HYP")
	results = model(limits=[(0,y.max()*10.0),(1.0,4.0),(-0.99,-0.001)])
	# results = model(limits = [[0,y.max()*10.0],[1.4,1.41],[-2.0,1.0]])
	qi, b, di = model.parameters
	
	# compute goodness-of-fit parameter
	coefficient_of_dermination = r2_score(smooth_y, decline_curve(smooth_t, qi, b, di))

	# display accuracy of fit to user
	print '\nR2 value of the fit is %.2f \n' %coefficient_of_dermination
	print '	R-squared is a statistical measure of how close the data are to the fitted regression line. \n\
	It is also known as the coefficient of determination, or the coefficient of multiple determination for multiple regression.\n\
	It is the percentage of the response variable variation that is explained by a linear model.\n\
		R-squared = Explained variation / Total variation\n\
		R-squared is always between 0 and 100%:\n\
	0.0 indicates that the model explains none of the variability of the response data around its mean.\n\
	1.00 indicates that the model explains all the variability of the response data around its mean.\n\
	In general, the higher the R-squared, the better the model fits your data.\n\n\
	If your R-squared value is low it can mean your selected wells have very different decline curve trends,\n\
	or that your selected wells have an unusual decline curve shape.'

	# return solution
	return qi, b, di, coefficient_of_dermination

def nominal_decline(qi, b, di):
	import numpy as np

	x = 365 #days
	q_oneYear = qi*(1.0-b*di*x)**(-1.0/b)
	return 100 * ((qi/q_oneYear)**b - 1.0) / (b * 12)