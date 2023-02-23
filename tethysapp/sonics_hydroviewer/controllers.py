from tethys_sdk.gizmos import *
from django.shortcuts import render
from tethys_sdk.gizmos import PlotlyView
from django.http import HttpResponse, JsonResponse
from tethys_sdk.routing import controller
import os
import sys
import json
import math
import requests
import geoglows
import pandas as pd
import xarray as xr
import netCDF4 as nc
import datetime as dt
from glob import glob
from lmoments3 import distr
import plotly.graph_objs as go
from .app import SonicsHydroviewer as app


@controller(name='home',url='sonics-hydroviewer')
def home(request):
	"""
	Controller for the app home page.
	"""

	folder = app.get_custom_setting('folder')
	print(folder)
	forecast_nc_list = sorted(glob(os.path.join(folder, "*.nc")))
	print(forecast_nc_list)

	dates_array = []

	for file in forecast_nc_list:
		dates_array.append(file[len(folder) + 1 + 23:-3])

	dates = []

	for date in dates_array:
		date_f = dt.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8])).strftime('%Y-%m-%d')
		dates.append([date_f, date])

	dates.append(['Select Date', dates[-1][1]])
	dates.reverse()

	date_picker = DatePicker(name='datesSelect',
							 display_text='Date',
							 autoclose=True,
							 format='yyyy-mm-dd',
							 start_date=dates[-1][0],
							 end_date=dates[1][0],
							 start_view='month',
							 today_button=True,
							 initial='')

	region_index = json.load(open(os.path.join(os.path.dirname(__file__), 'public', 'geojson', 'index.json')))
	regions = SelectInput(
		display_text='Zoom to a Region:',
		name='regions',
		multiple=False,
		#original=True,
		options=[(region_index[opt]['name'], opt) for opt in region_index],
        initial='',
        select2_options={'placeholder': 'Select a Region', 'allowClear': False}
	)

	# Select Province
	province_index = json.load(open(os.path.join(os.path.dirname(__file__), 'public', 'geojson2', 'index2.json')))
	provinces = SelectInput(
		display_text='Zoom to a Province:',
		name='provinces',
		multiple=False,
		# original=True,
		options=[(province_index[opt]['name'], opt) for opt in province_index],
		initial='',
		select2_options={'placeholder': 'Select a Province', 'allowClear': False}
	)

	# Select Basins
	basin_index = json.load(open(os.path.join(os.path.dirname(__file__), 'public', 'geojson3', 'index3.json')))
	basins = SelectInput(
		display_text='Zoom to a Basin:',
		name='basins',
		multiple=False,
		# original=True,
		options=[(basin_index[opt]['name'], opt) for opt in basin_index],
		initial='',
		select2_options={'placeholder': 'Select a Basin', 'allowClear': False}
	)

	# Select SubBasins
	subbasin_index = json.load(open(os.path.join(os.path.dirname(__file__), 'public', 'geojson4', 'index4.json')))
	subbasins = SelectInput(
		display_text='Zoom to a Subbasin:',
		name='subbasins',
		multiple=False,
		# original=True,
		options=[(subbasin_index[opt]['name'], opt) for opt in subbasin_index],
		initial='',
		select2_options={'placeholder': 'Select a Subbasin', 'allowClear': False}
	)

	context = {
		"date_picker": date_picker,
		"regions": regions,
		"provinces": provinces,
		"basins": basins,
		"subbasins": subbasins,
	}

	return render(request, 'sonics_hydroviewer/home.html', context)


def gve_1(loc: float, scale: float, shape: float, rp: int or float) -> float:
	"""
	Solves the Gumbel Type I probability distribution function (pdf) = exp(-exp(-b)) where b is the covariate. Provide
	the standard deviation and mean of the list of annual maximum flows. Compare scipy.stats.gumbel_r
	Args:
	  std (float): the standard deviation of the series
	  xbar (float): the mean of the series
	  skew (float): the skewness of the series
	  rp (int or float): the return period in years
	Returns:
	  float, the flow corresponding to the return period specified
	"""

	return ((scale / shape) * (1 - math.exp(shape * (math.log(-math.log(1 - (1 / rp))))))) + loc

@controller(name='get_hydrographs',url='sonics-hydroviewer/get-hydrographs')
def get_hydrographs(request):
	try:
		get_data = request.GET
		# get stream attributes
		comid = get_data['comid']
		region = get_data['region']
		subbasin = get_data['subbasin']
		watershed = get_data['watershed']

		folder = app.get_custom_setting('folder')
		forecast_nc_list = sorted(glob(os.path.join(folder, "*.nc")), reverse=True)
		nc_file = forecast_nc_list[0]

		qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_hist
		time_dataset = qout_datasets.time_hist

		historical_simulation_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values,
												columns=['Streamflow (m3/s)'])
		historical_simulation_df.index.name = 'Datetime'

		'''Getting Return Periods'''
		return_periods_values = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).threshold
		return_periods_values = return_periods_values.values
		return_periods = [10, 5, 2.33]

		d = {'rivid': [comid], 'return_period_10': [return_periods_values[2]],
			 'return_period_5': [return_periods_values[1]], 'return_period_2_33': [return_periods_values[0]]}

		rperiods = pd.DataFrame(data=d)
		rperiods.set_index('rivid', inplace=True)

		'''Plotting hydrograph'''
		hydroviewer_figure = geoglows.plots.historic_simulation(historical_simulation_df)

		x_vals = (
			historical_simulation_df.index[0], historical_simulation_df.index[len(historical_simulation_df.index) - 1],
			historical_simulation_df.index[len(historical_simulation_df.index) - 1], historical_simulation_df.index[0])
		max_visible = max(historical_simulation_df.max())

		'''Getting Return Periods'''
		r2_33 = int(rperiods.iloc[0]['return_period_2_33'])

		colors = {
			'2.33 Year': 'rgba(243, 255, 0, .4)',
			'5 Year': 'rgba(255, 165, 0, .4)',
			'10 Year': 'rgba(255, 0, 0, .4)',
		}

		if max_visible > r2_33:
			visible = True
			hydroviewer_figure.for_each_trace(
				lambda trace: trace.update(visible=True) if trace.name == "Maximum & Minimum Flow" else (), )
		else:
			visible = 'legendonly'
			hydroviewer_figure.for_each_trace(
				lambda trace: trace.update(visible=True) if trace.name == "Maximum & Minimum Flow" else (), )

		def template(name, y, color, fill='toself'):
			return go.Scatter(
				name=name,
				x=x_vals,
				y=y,
				legendgroup='returnperiods',
				fill=fill,
				visible=visible,
				line=dict(color=color, width=0))

		r5 = int(rperiods.iloc[0]['return_period_5'])
		r10 = int(rperiods.iloc[0]['return_period_10'])

		hydroviewer_figure.add_trace(
			template('Return Periods', (r10 * 0.05, r10 * 0.05, r10 * 0.05, r10 * 0.05), 'rgba(0,0,0,0)', fill='none'))
		hydroviewer_figure.add_trace(template(f'2.33 Year: {r2_33}', (r2_33, r2_33, r5, r5), colors['2.33 Year']))
		hydroviewer_figure.add_trace(template(f'5 Year: {r5}', (r5, r5, r10, r10), colors['5 Year']))
		hydroviewer_figure.add_trace(template(f'10 Year: {r10}', (
		r10, r10, max(r10 + r10 * 0.05, max_visible), max(r10 + r10 * 0.05, max_visible)), colors['10 Year']))

		hydroviewer_figure['layout']['xaxis'].update(autorange=True)

		chart_obj = PlotlyView(hydroviewer_figure)

		context = {
			'gizmo_object': chart_obj,
		}

		return render(request, 'sonics_hydroviewer/gizmo_ajax.html', context)

	except Exception as e:
		exc_type, exc_obj, exc_tb = sys.exc_info()
		print("error: " + str(e))
		print("line: " + str(exc_tb.tb_lineno))

		return JsonResponse({
			'error': f'{"error: " + str(e), "line: " + str(exc_tb.tb_lineno)}',
		})

@controller(name='get_simulated_discharge_csv',url='sonics-hydroviewer/get-simulated-discharge-csv')
def get_simulated_discharge_csv(request):
	"""
	Get historic simulations from ERA Interim
	"""

	try:
		get_data = request.GET
		# get stream attributes
		comid = get_data['comid']
		region = get_data['region']
		subbasin = get_data['subbasin']
		watershed = get_data['watershed']

		folder = app.get_custom_setting('folder')
		forecast_nc_list = sorted(glob(os.path.join(folder, "*.nc")), reverse=True)
		nc_file = forecast_nc_list[0]

		qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_hist
		time_dataset = qout_datasets.time_hist

		historical_simulation_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
		historical_simulation_df.index.name = 'Datetime'

		pairs = [list(a) for a in zip(historical_simulation_df.index, historical_simulation_df.iloc[:, 0])]

		response = HttpResponse(content_type='text/csv')
		response['Content-Disposition'] = 'attachment; filename=simulated_discharge_{0}.csv'.format(comid)

		historical_simulation_df.to_csv(encoding='utf-8', header=True, path_or_buf=response)

		return response

	except Exception as e:
		exc_type, exc_obj, exc_tb = sys.exc_info()
		print("error: " + str(e))
		print("line: " + str(exc_tb.tb_lineno))

		return JsonResponse({
			'error': f'{"error: " + str(e), "line: " + str(exc_tb.tb_lineno)}',
		})

@controller(name='get-time-series',url='sonics-hydroviewer/get-time-series')
def get_time_series(request):
	try:
		get_data = request.GET
		# get stream attributes
		comid = get_data['comid']
		region = get_data['region']
		subbasin = get_data['subbasin']
		watershed = get_data['watershed']
		startdate = get_data['startdate']

		folder = app.get_custom_setting('folder')

		'''Getting Forecast Stats'''
		if startdate != '':
			nc_file = folder + '/PISCO_HyD_ARNOVIC_v1.0_' + startdate + '.nc'
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_hist
			time_dataset = qout_datasets.time_hist
			historical_simulation_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			initial_condition = historical_simulation_df.loc[historical_simulation_df.index == pd.to_datetime(historical_simulation_df.index[-1])]

			'''ETA eqm Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_eta_eqm
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_eta_eqm_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_eta_eqm_df.index.name = 'Datetime'
			forecast_eta_eqm_df = forecast_eta_eqm_df.append(initial_condition)
			forecast_eta_eqm_df.sort_index(inplace=True)

			'''ETA scal Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_eta_scal
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_eta_scal_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_eta_scal_df.index.name = 'Datetime'
			forecast_eta_scal_df = forecast_eta_scal_df.append(initial_condition)
			forecast_eta_scal_df.sort_index(inplace=True)

			'''GFS Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_gfs
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_gfs_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_gfs_df.index.name = 'Datetime'
			forecast_gfs_df = forecast_gfs_df.append(initial_condition)
			forecast_gfs_df.sort_index(inplace=True)

			'''WRF Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_wrf
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_wrf_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_wrf_df.index.name = 'Datetime'
			forecast_wrf_df = forecast_wrf_df.append(initial_condition)
			forecast_wrf_df.sort_index(inplace=True)

		else:

			forecast_nc_list = sorted(glob(os.path.join(folder, "*.nc")), reverse=True)
			nc_file = forecast_nc_list[0]
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_hist
			time_dataset = qout_datasets.time_hist
			historical_simulation_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			initial_condition = historical_simulation_df.loc[historical_simulation_df.index == pd.to_datetime(historical_simulation_df.index[-1])]

			'''ETA eqm Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_eta_eqm
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_eta_eqm_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_eta_eqm_df.index.name = 'Datetime'
			forecast_eta_eqm_df = forecast_eta_eqm_df.append(initial_condition)
			forecast_eta_eqm_df.sort_index(inplace=True)

			'''ETA scal Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_eta_scal
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_eta_scal_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_eta_scal_df.index.name = 'Datetime'
			forecast_eta_scal_df = forecast_eta_scal_df.append(initial_condition)
			forecast_eta_scal_df.sort_index(inplace=True)

			'''GFS Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_gfs
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_gfs_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_gfs_df.index.name = 'Datetime'
			forecast_gfs_df = forecast_gfs_df.append(initial_condition)
			forecast_gfs_df.sort_index(inplace=True)

			'''WRF Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_wrf
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_wrf_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_wrf_df.index.name = 'Datetime'
			forecast_wrf_df = forecast_wrf_df.append(initial_condition)
			forecast_wrf_df.sort_index(inplace=True)

		'''Getting Return Periods'''
		return_periods_values = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).threshold
		return_periods_values = return_periods_values.values
		return_periods = [10, 5, 2.33]

		d = {'rivid': [comid], 'return_period_10': [return_periods_values[2]],
			 'return_period_5': [return_periods_values[1]], 'return_period_2_33': [return_periods_values[0]]}

		rperiods = pd.DataFrame(data=d)
		rperiods.set_index('rivid', inplace=True)

		'''Plotting Forecast'''

		dates = forecast_gfs_df.index.tolist()
		startdate = dates[0]
		enddate = dates[-1]

		gfs_forecast = go.Scatter(name='GFS Forecast',
								  x=forecast_gfs_df.index,
								  y=forecast_gfs_df['Streamflow (m3/s)'],
								  showlegend=True,
								  line=dict(color='black', dash='dash'))

		eta_eqm_forecast = go.Scatter(name='ETA eqm Forecast',
									  x=forecast_eta_eqm_df.index,
									  y=forecast_eta_eqm_df['Streamflow (m3/s)'],
									  showlegend=True,
									  line=dict(color='blue', dash='dash'))

		eta_scal_forecast = go.Scatter(name='ETA scal Forecast',
									   x=forecast_eta_scal_df.index,
									   y=forecast_eta_scal_df['Streamflow (m3/s)'],
									   showlegend=True,
									   line=dict(color='green', dash='dash'))

		wrf_forecast = go.Scatter(name='WRF Forecast',
								  x=forecast_wrf_df.index,
								  y=forecast_wrf_df['Streamflow (m3/s)'],
								  showlegend=True,
								  line=dict(color='brown', dash='dash'))

		layout = go.Layout(title='SONICS Forecast at {0}'.format(comid),
						   xaxis=dict(title='Dates', ),
						   yaxis=dict(title='Streamflow (m<sup>3</sup>/s)', autorange=True),
						   showlegend=True)

		hydroviewer_figure = go.Figure(data=[gfs_forecast, eta_eqm_forecast, eta_scal_forecast, wrf_forecast], layout=layout)

		x_vals = (forecast_gfs_df.index[0], forecast_gfs_df.index[len(forecast_gfs_df.index) - 1],
				  forecast_gfs_df.index[len(forecast_gfs_df.index) - 1], forecast_gfs_df.index[0])
		max_visible = max(forecast_gfs_df.max())
		max_visible = max(max(forecast_eta_eqm_df.max()), max_visible)
		max_visible = max(max(forecast_eta_scal_df.max()), max_visible)
		max_visible = max(max(forecast_wrf_df.max()), max_visible)

		'''Adding Recent Days'''
		records_df = historical_simulation_df.loc[
			historical_simulation_df.index >= pd.to_datetime(startdate - dt.timedelta(days=4))]
		records_df = records_df.loc[records_df.index <= pd.to_datetime(startdate + dt.timedelta(days=2))]

		if len(records_df.index) > 0:
			hydroviewer_figure.add_trace(go.Scatter(
				name='1st days forecasts SONICS',
				x=records_df.index,
				y=records_df.iloc[:, 0].values,
				line=dict(color='#FFA15A', )
			))

			x_vals = (records_df.index[0], forecast_gfs_df.index[len(forecast_gfs_df.index) - 1],
					  forecast_gfs_df.index[len(forecast_gfs_df.index) - 1], records_df.index[0])
			max_visible = max(max(records_df.max()), max_visible)

		'''Getting Return Periods'''
		r2_33 = int(rperiods.iloc[0]['return_period_2_33'])

		colors = {
			'2.33 Year': 'rgba(243, 255, 0, .4)',
			'5 Year': 'rgba(255, 165, 0, .4)',
			'10 Year': 'rgba(255, 0, 0, .4)',
		}

		if max_visible > r2_33:
			visible = True
			hydroviewer_figure.for_each_trace(
				lambda trace: trace.update(visible=True) if trace.name == "Maximum & Minimum Flow" else (), )
		else:
			visible = 'legendonly'
			hydroviewer_figure.for_each_trace(
				lambda trace: trace.update(visible=True) if trace.name == "Maximum & Minimum Flow" else (), )

		def template(name, y, color, fill='toself'):
			return go.Scatter(
				name=name,
				x=x_vals,
				y=y,
				legendgroup='returnperiods',
				fill=fill,
				visible=visible,
				line=dict(color=color, width=0))

		r5 = int(rperiods.iloc[0]['return_period_5'])
		r10 = int(rperiods.iloc[0]['return_period_10'])

		hydroviewer_figure.add_trace(template('Return Periods', (r10 * 0.05, r10 * 0.05, r10 * 0.05, r10 * 0.05), 'rgba(0,0,0,0)', fill='none'))
		hydroviewer_figure.add_trace(template(f'2.33 Year: {r2_33}', (r2_33, r2_33, r5, r5), colors['2.33 Year']))
		hydroviewer_figure.add_trace(template(f'5 Year: {r5}', (r5, r5, r10, r10), colors['5 Year']))
		hydroviewer_figure.add_trace(template(f'10 Year: {r10}', (r10, r10, max(r10 + r10 * 0.05, max_visible), max(r10 + r10 * 0.05, max_visible)), colors['10 Year']))

		hydroviewer_figure['layout']['xaxis'].update(autorange=True)

		chart_obj = PlotlyView(hydroviewer_figure)

		context = {
			'gizmo_object': chart_obj,
		}

		return render(request, 'sonics_hydroviewer/gizmo_ajax.html', context)


	except Exception as e:
		exc_type, exc_obj, exc_tb = sys.exc_info()
		print("error: " + str(e))
		print("line: " + str(exc_tb.tb_lineno))

		return JsonResponse({
			'error': f'{"error: " + str(e), "line: " + str(exc_tb.tb_lineno)}',
		})

@controller(name='get_forecast_data_csv',url='sonics-hydroviewer/get-forecast-data-csv')
def get_forecast_data_csv(request):
	"""""
	Returns Forecast data as csv
	"""""

	try:
		get_data = request.GET
		# get stream attributes
		comid = get_data['comid']
		region = get_data['region']
		subbasin = get_data['subbasin']
		watershed = get_data['watershed']
		startdate = get_data['startdate']

		folder = app.get_custom_setting('folder')

		'''Getting Forecast Stats'''
		if startdate != '':
			nc_file = folder + '/PISCO_HyD_ARNOVIC_v1.0_' + startdate + '.nc'
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_hist
			time_dataset = qout_datasets.time_hist
			historical_simulation_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			initial_condition = historical_simulation_df.loc[historical_simulation_df.index == pd.to_datetime(historical_simulation_df.index[-1])]

			'''ETA eqm Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_eta_eqm
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_eta_eqm_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_eta_eqm_df.index.name = 'Datetime'
			forecast_eta_eqm_df = forecast_eta_eqm_df.append(initial_condition)
			forecast_eta_eqm_df.sort_index(inplace=True)
			forecast_eta_eqm_df.rename(columns={"Streamflow (m3/s)": "ETA eqm Streamflow (m3/s)"}, inplace=True)

			'''ETA scal Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_eta_scal
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_eta_scal_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_eta_scal_df.index.name = 'Datetime'
			forecast_eta_scal_df = forecast_eta_scal_df.append(initial_condition)
			forecast_eta_scal_df.sort_index(inplace=True)
			forecast_eta_scal_df.rename(columns={"Streamflow (m3/s)": "ETA scal Streamflow (m3/s)"}, inplace=True)

			'''GFS Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_gfs
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_gfs_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_gfs_df.index.name = 'Datetime'
			forecast_gfs_df = forecast_gfs_df.append(initial_condition)
			forecast_gfs_df.sort_index(inplace=True)
			forecast_gfs_df.rename(columns={"Streamflow (m3/s)": "GFS Streamflow (m3/s)"}, inplace=True)

			'''WRF Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_wrf
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_wrf_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values,
										   columns=['Streamflow (m3/s)'])
			forecast_wrf_df.index.name = 'Datetime'
			forecast_wrf_df = forecast_wrf_df.append(initial_condition)
			forecast_wrf_df.sort_index(inplace=True)
			forecast_wrf_df.rename(columns={"Streamflow (m3/s)": "WRF Streamflow (m3/s)"}, inplace=True)

		else:

			forecast_nc_list = sorted(glob(os.path.join(folder, "*.nc")), reverse=True)
			nc_file = forecast_nc_list[0]
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_hist
			time_dataset = qout_datasets.time_hist
			historical_simulation_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values,columns=['Streamflow (m3/s)'])
			initial_condition = historical_simulation_df.loc[historical_simulation_df.index == pd.to_datetime(historical_simulation_df.index[-1])]

			'''ETA eqm Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_eta_eqm
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_eta_eqm_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_eta_eqm_df.index.name = 'Datetime'
			forecast_eta_eqm_df = forecast_eta_eqm_df.append(initial_condition)
			forecast_eta_eqm_df.sort_index(inplace=True)
			forecast_eta_eqm_df.rename(columns={"Streamflow (m3/s)": "ETA eqm Streamflow (m3/s)"}, inplace=True)

			'''ETA scal Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_eta_scal
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_eta_scal_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_eta_scal_df.index.name = 'Datetime'
			forecast_eta_scal_df = forecast_eta_scal_df.append(initial_condition)
			forecast_eta_scal_df.sort_index(inplace=True)
			forecast_eta_scal_df.rename(columns={"Streamflow (m3/s)": "ETA scal Streamflow (m3/s)"}, inplace=True)

			'''GFS Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_gfs
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_gfs_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_gfs_df.index.name = 'Datetime'
			forecast_gfs_df = forecast_gfs_df.append(initial_condition)
			forecast_gfs_df.sort_index(inplace=True)
			forecast_gfs_df.rename(columns={"Streamflow (m3/s)": "GFS Streamflow (m3/s)"}, inplace=True)

			'''WRF Forecast'''
			qout_datasets = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).qr_wrf
			time_dataset = xr.open_dataset(nc_file, autoclose=True).sel(comid=comid).time_frst
			forecast_wrf_df = pd.DataFrame(qout_datasets.values, index=time_dataset.values, columns=['Streamflow (m3/s)'])
			forecast_wrf_df.index.name = 'Datetime'
			forecast_wrf_df = forecast_wrf_df.append(initial_condition)
			forecast_wrf_df.sort_index(inplace=True)
			forecast_wrf_df.rename(columns={"Streamflow (m3/s)": "WRF Streamflow (m3/s)"}, inplace=True)

		forecast_df = pd.concat([forecast_eta_eqm_df, forecast_eta_scal_df, forecast_gfs_df, forecast_wrf_df], axis=1)

		# Writing CSV
		response = HttpResponse(content_type='text/csv')
		response['Content-Disposition'] = 'attachment; filename=streamflow_forecast_{0}_{1}.csv'.format(comid, startdate)

		forecast_df.to_csv(encoding='utf-8', header=True, path_or_buf=response)

		return response

	except Exception as e:
		exc_type, exc_obj, exc_tb = sys.exc_info()
		print("error: " + str(e))
		print("line: " + str(exc_tb.tb_lineno))

		return JsonResponse({
				'error': f'{"error: " + str(e), "line: " + str(exc_tb.tb_lineno)}',
		})