import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.colors as colors
from matplotlib import cm

mpl.rcParams['hatch.linewidth'] = 0.5

# ================================================================
# Read and interpolate data
# ================================================================

def smooth_spikes(input_data):

	nx = np.shape(input_data)[1]
	ny = np.shape(input_data)[0]

	output_data = np.copy(input_data)

	made_replacements = True
	num_iter = 0
	max_iter = 10

	deviation = 2.0

	print('====Starting Smoothing====')

	while made_replacements == True and num_iter < max_iter:

		made_replacements = False
		num_replacements = 0
		num_iter += 1

		for k in range(1, ny-1):
			for j in range(1, nx-1):

				neighbors = np.array([input_data[k-1, j],
									input_data[k+1, j],
									input_data[k, j-1],
									input_data[k, j+1],
									input_data[k-1, j-1],
									input_data[k+1, j-1],
									input_data[k-1, j+1],
									input_data[k+1, j+1]])

				mean = np.mean(neighbors)
				std = np.std(neighbors)


				if output_data[k, j] < mean-deviation*std or mean+deviation*std < output_data[k, j]:
					output_data[k, j] = mean
					made_replacements = True
					num_replacements += 1

		print('Finished iteration %02d with %02d replacements.' % (num_iter, num_replacements))

		input_data = np.copy(output_data)

	print('====Finished Smoothing====')

	return output_data


# path_to_files = 'columnwise_runs'
# path_to_files = 'fine_col'
# path_to_files = 'fine_col_new'
# path_to_files = 'fine_col_iso_5'
path_to_files = 'fine_col_iso_5_financials'

# path_to_files = 'new_runs/rowwise'
num_procs = 72

stages = [1, 2, 3, 4, 5]
nn = len(stages)

for k, stage_num in enumerate(stages):

	raw_data = np.genfromtxt('%s/raw_data_%dstage.csv' % (path_to_files, stage_num), skip_header=1, delimiter=',')

	if k == 0:
		x_mesh = np.reshape(raw_data[:, 0], (num_procs, num_procs)).T
		y_mesh = np.reshape(raw_data[:, 1], (num_procs, num_procs)).T

		raw_ec = np.zeros((nn, num_procs, num_procs))
		raw_lcow = np.zeros((nn, num_procs, num_procs))

		interp_ec = np.zeros((nn, num_procs, num_procs))
		interp_lcow = np.zeros((nn, num_procs, num_procs))

		clean_ec = np.zeros((nn, num_procs, num_procs))
		clean_lcow = np.zeros((nn, num_procs, num_procs))

	# Mask which excluded NaN values
	mask = np.isfinite(raw_data[:, 2])

	# Known valid points
	x0 = raw_data[mask, 0:2]

	# Known data points
	y0_ec = raw_data[mask, 2]
	y0_lcow = raw_data[mask, 3]

	raw_ec[k, :, :] = np.reshape(raw_data[:, 2], (num_procs, num_procs)).T
	raw_lcow[k, :, :] = np.reshape(raw_data[:, 3], (num_procs, num_procs)).T

	interp_ec[k, :, :] = griddata(x0, y0_ec, (x_mesh, y_mesh), method='linear')
	interp_lcow[k, :, :] = griddata(x0, y0_lcow, (x_mesh, y_mesh), method='linear')

	clean_ec[k, :, :] = smooth_spikes(interp_ec[k, :, :])
	clean_lcow[k, :, :] = smooth_spikes(interp_lcow[k, :, :])


fig, ax = plt.subplots(3, nn, figsize=(18, 12), dpi=100)
for k in range(nn):
	ax[0, k].contourf(x_mesh, y_mesh, raw_ec[k, :, :])
	ax[0, k].set_title('Raw Data, %d stage' % (k+2))
	ax[1, k].contourf(x_mesh, y_mesh, interp_ec[k, :, :])
	ax[1, k].set_title('Holes Filled, %d stage' % (k+2))
	ax[2, k].contourf(x_mesh, y_mesh, clean_ec[k, :, :])
	ax[2, k].set_title('Outliers Smoothed, %d stage' % (k+2))
plt.savefig('data_cleaning_ec_3.pdf')

fig, ax = plt.subplots(3, nn, figsize=(18, 12), dpi=100)
for k in range(nn):
	ax[0, k].contourf(x_mesh, y_mesh, raw_lcow[k, :, :])
	ax[0, k].set_title('Raw Data, %d stage' % (k+2))
	ax[1, k].contourf(x_mesh, y_mesh, interp_lcow[k, :, :])
	ax[1, k].set_title('Holes Filled, %d stage' % (k+2))
	ax[2, k].contourf(x_mesh, y_mesh, clean_lcow[k, :, :])
	ax[2, k].set_title('Outliers Smoothed, %d stage' % (k+2))
plt.savefig('data_cleaning_lcow_3.pdf')


# Merge all stages into 1 plot
compiled_ec = np.zeros((num_procs, num_procs))
compiled_lcow = np.zeros((num_procs, num_procs))
compiled_stage = np.zeros((num_procs, num_procs))

nx = np.shape(compiled_ec)[1]
ny = np.shape(compiled_ec)[0]

for k in range(ny):
	for j in range(nx):
		ec_options = clean_ec[:, k, j]
		lcow_options = clean_lcow[:, k, j]

		# Find the minimum value from all possible stages
		try:
			idx = np.nanargmin(lcow_options)

		except:
			compiled_ec[k, j] = 1e6
			compiled_lcow[k, j] = 1e6
			compiled_stage[k, j] = 6

		else:
			compiled_ec[k, j] = ec_options[idx]
			compiled_lcow[k, j] = lcow_options[idx]
			compiled_stage[k, j] = stages[idx]



final_mask = np.zeros((num_procs, num_procs))

#==== Mask option 3, use EC above a certain point to mask
final_mask[compiled_ec < 300] = 1

for k in range(1, ny-1):
	for j in range(1, nx-1):
		if final_mask[k, j-1] == 1 and final_mask[k, j+1] == 1 and final_mask[k, j] == 0:
			final_mask[k, j] = 1
			compiled_ec[k, j] = 0.5*(compiled_ec[k, j-1] + compiled_ec[k, j+1])
			compiled_lcow[k, j] = 0.5*(compiled_lcow[k, j-1] + compiled_lcow[k, j+1])


#==== Mask option 1, straight line cutoff
# left_top_pt = 66.0/1000.0
# right_bottom_pt = 135.0/1000.0

# m = (np.amax(y_mesh) - np.amin(y_mesh))/(left_top_pt - right_bottom_pt)
# b = np.amax(y_mesh) - m*left_top_pt

# y_calc = m*x_mesh + b

# final_mask[y_mesh < y_calc] = 1

#==== Mask option 2, try to figure out slopes of lcow
# for j in range(nx):
# 	inflection_pt = False

# 	for k in range(1, ny):
		# slope = compiled_lcow[k, j] - compiled_lcow[k-1, j]

		# if k == 1:
		# 	prev_slope = slope

		# if slope > 3.0*np.abs(prev_slope):
		# 	inflection_pt = True

		# if inflection_pt == True:
		# 	final_mask[k, j] = 1.0

		# prev_slope = slope




compiled_ec[final_mask == 0] = np.nan
compiled_lcow[final_mask == 0] = np.nan
compiled_stage[final_mask == 0] = 6

# ================================================================
# Generate and format plots
# ================================================================
ylgnbu = cm.get_cmap('YlGnBu')
ylgnbu_shift = ylgnbu(np.linspace(0, 1, 256)**(1.0/1.5))
ylgnbu_shift = colors.ListedColormap(ylgnbu_shift)

fig, ax = plt.subplots(1, 2, figsize=(12, 6), dpi=100)

# Bins for contourf, 
cont_levels = [1.5, 2.5, 3.5, 4.5, 5.01]

cbar = ax[0].contourf(x_mesh, y_mesh, compiled_ec, [0, 2, 5, 10, 20, 30, 40, 50], extend='max', cmap=ylgnbu_shift)
# ax[0].plot(points[0], points[1], '.', color='k', markersize=.5, alpha=0.5)
ax[0].contourf(x_mesh, y_mesh, compiled_stage, cont_levels, colors='none', hatches=['/', '\\', '//', '\\\\'])
ax[0].contour(x_mesh, y_mesh, compiled_stage, cont_levels, colors='k', linewidths=0.75)
fig.colorbar(cbar, ax=ax[0])

xticks = np.linspace(0.005, 0.155, 6)
xlabels = ['%.0f' % (1000.0*k) for k in xticks]
ax[0].set_xticks(xticks)
ax[0].set_xticklabels(xlabels)
ax[0].set_xlabel('Feed Concentration (parts per thousand)')

yticks = np.linspace(0.3, 0.7, 9)
ylabels = ['%.0f' % (100.0*k) for k in yticks]
ax[0].set_yticks(yticks)
ax[0].set_yticklabels(ylabels)
ax[0].set_ylabel('Water Recovery (%)')

ax[0].set_title('Energy Consumption (kWh/m$^3$)')
ax[0].grid(alpha=0.2, zorder=10)


cbar = ax[1].contourf(x_mesh, y_mesh, compiled_lcow, [0, 0.5, 1, 2, 3, 4, 5, 6], extend='max', cmap=ylgnbu_shift)
# ax[1].plot(points[0], points[1], '.', color='k', markersize=.5, alpha=0.5)
ax[1].contourf(x_mesh, y_mesh, compiled_stage, cont_levels, colors='none', hatches=['/', '\\', '//', '\\\\'])
ax[1].contour(x_mesh, y_mesh, compiled_stage, cont_levels, colors='k', linewidths=0.75)
fig.colorbar(cbar, ax=ax[1])

ax[1].set_xticks(xticks)
ax[1].set_xticklabels(xlabels)
ax[1].set_xlabel('Feed Concentration (parts per thousand)')

ax[1].set_yticks(yticks)
ax[1].set_yticklabels(ylabels)

ax[1].set_title('Levelized Cost of Water (\$/m$^3$)')
ax[1].grid(alpha=0.2, zorder=10)


patch_handles = []
patch_handles.append(mpatches.Patch(edgecolor='k', facecolor='none',hatch='',label='1-Stage'))
patch_handles.append(mpatches.Patch(edgecolor='k', facecolor='none',hatch='//',label='2-Stage'))
patch_handles.append(mpatches.Patch(edgecolor='k', facecolor='none',hatch='\\\\',label='3-Stage'))
patch_handles.append(mpatches.Patch(edgecolor='k', facecolor='none',hatch='////',label='4-Stage'))
patch_handles.append(mpatches.Patch(edgecolor='k', facecolor='none',hatch='\\\\\\\\',label='5-Stage'))

ax[0].legend(handles = patch_handles, loc=8, ncol=5, bbox_to_anchor=(1.25, -.23))

plt.savefig('contour_ec_and_lcow_extended_col_plus_3.png', bbox_inches='tight')
plt.savefig('contour_ec_and_lcow_extended_col_plus_3.pdf', bbox_inches='tight')

