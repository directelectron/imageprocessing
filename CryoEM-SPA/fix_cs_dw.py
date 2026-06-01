#!/usr/bin/env python3
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.ndimage import uniform_filter1d
import os
import shutil
import sys

def find_linear_negative_region(
	arr: np.ndarray,
	window_size: int = None,
	min_window: int = 5,
	smoothing_sigma: int = 3,
	r2_threshold: float = 0.9,
) -> dict:
	"""
	Find the region of a 1D numpy array that is most closely linear
	with a negative slope, then fit a linear equation to it.
	
	Strategy:
		1. Smooth the array to reduce noise for region detection.
		2. Slide a window across the array; for each window, fit a line
			and score it by R² (goodness of fit) penalised if slope is
			non-negative.
		3. Select the window with the highest score (best linear fit AND
			negative slope).
		4. Fit the final linear model to the RAW (unsmoothed) data in
			that window.
	
	Parameters
	----------
	arr : np.ndarray
		1D input array (noisy).
	window_size : int, optional
		Number of points in the sliding window.  If None, defaults to
		len(arr) // 4 (25% of the array), with a minimum of min_window.
	min_window : int
		Minimum allowed window size.  Default 5.
	smoothing_sigma : int
		Width of the uniform smoothing filter applied before region
		detection (does NOT affect the final fit, which uses raw data).
	r2_threshold : float
		Expansion stops when R² of raw fit falls below
		r2_threshold * best_R² from the sliding window step.
		Default 0.9 (i.e. 90% of best R²).
	
	Returns
	-------
	dict with keys:
		'slope'		 : float  — slope of the fitted line
		'intercept'	: float  — intercept of the fitted line
		'r_value'	  : float  — Pearson r of fit in the selected region
		'r_squared'	: float  — R² of fit in the selected region
		'idx_start'	: int	 — start index of the best region
		'idx_end'	  : int	 — end index (inclusive) of the best region
		'fitted_full' : ndarray — line evaluated over ALL indices
		'fitted_region': ndarray — line evaluated over the best region only
		'x_region'	 : ndarray — indices of the best region
		'y_region'	 : ndarray — raw values in the best region
	"""
	n = len(arr)
	
	if window_size is None:
		window_size = max(min_window, n // 4)
	window_size = min(window_size, n)
	
	# ── Step 1: smooth for region detection only ──────────────────────────
	smoothed = uniform_filter1d(arr.astype(float), size=smoothing_sigma)
	
	# ── Step 2: slide window and score each position ──────────────────────
	best_score      = -np.inf
	best_start      = 0
	best_r2_initial = 0.0
	
	for i in range(n - window_size + 1):
		x = np.arange(window_size, dtype=float)
		y = smoothed[i : i + window_size]
		
		slope, intercept, r_value, _, _ = stats.linregress(x, y)
		
		# Only consider windows with a negative slope
		if slope >= 0:
			continue
		
		r_squared = r_value ** 2
		
		# Score: R² weighted by how negative the normalised slope is.
		# Normalising the slope by the data range avoids units issues.
		data_range = np.ptp(smoothed) if np.ptp(smoothed) > 0 else 1.0
		norm_slope = abs(slope) / (data_range / window_size)
		score = r_squared + 0.1 * np.tanh(norm_slope)
		
		if score > best_score:
				best_score = score
				best_start = i
				best_r2_initial = r_squared
	
	best_end = best_start + window_size - 1
	
	# ── Step 3: expand window, stopping when R² < threshold * best_R² ────
	# The threshold is computed from the R² of the raw data in the
	# initial best window (not from the smoothed fit used for scoring).
	x_init = np.arange(best_start, best_end + 1, dtype=float)
	y_init = arr[best_start : best_end + 1].astype(float)
	_, _, r_init, _, _ = stats.linregress(x_init, y_init)
	r2_raw_initial = r_init ** 2

	# The R² value that expansion must stay above
	r2_limit = r2_threshold * r2_raw_initial

	exp_start = best_start
	exp_end   = best_end

	def r2_for_region(start: int, end: int) -> float:
		"""Compute R² of a linear fit to raw data between start and end."""
		if end - start < 1:
			return 0.0
		x = np.arange(start, end + 1, dtype=float)
		y = arr[start : end + 1].astype(float)
		_, _, r, _, _ = stats.linregress(x, y)
		return r ** 2

	# Expand one index at a time in the better direction first, then
	# alternate left/right, stopping whichever side drops below the limit.
	can_expand_left = exp_start > 0
	can_expand_right = exp_end < n - 1

	while can_expand_left or can_expand_right:
		# Try expanding left
		r2_left = (
			r2_for_region(exp_start - 1, exp_end)
			if can_expand_left else -np.inf
		)
		# Try expanding right
		r2_right = (
			r2_for_region(exp_start, exp_end + 1)
			if can_expand_right else -np.inf
		)

		# Choose whichever expansion gives the higher R²
		if r2_left >= r2_right:
			if r2_left >= r2_limit:
				exp_start -= 1
				can_expand_left = exp_start > 0
			else:
				can_expand_left = False
			# Also check right side with updated window
			if can_expand_right:
				if r2_for_region(exp_start, exp_end + 1) >= r2_limit:
					exp_end += 1
					can_expand_right = exp_end < n - 1
				else:
					can_expand_right = False
		else:
			if r2_right >= r2_limit:
				exp_end += 1
				can_expand_right = exp_end < n - 1
			else:
				can_expand_right = False
			# Also check left side with updated window
			if can_expand_left:
				if r2_for_region(exp_start - 1, exp_end) >= r2_limit:
					exp_start -= 1
					can_expand_left = exp_start > 0
				else:
					can_expand_left = False
	
	# ── Step 4: final linear fit on raw data in expanded region ───────────
	x_region = np.arange(best_start, best_end + 1, dtype=float)
	y_region = arr[best_start : best_end + 1].astype(float)
	
	slope, intercept, r_value, p_value, std_err = stats.linregress(x_region, y_region)
	
	# ── Step 5: evaluate over full array ──────────────────────────────────
	x_full		 = np.arange(n, dtype=float)
	fitted_full  = slope * x_full + intercept
	fitted_region = slope * x_region + intercept
	
	return {
		'slope':         slope,
		'intercept':     intercept,
		'r_value':       r_value,
		'r_squared':     r_value ** 2,
		'p_value':       p_value,
		'std_err':       std_err,
		'idx_start':     best_start,
		'idx_end':       best_end,
		'fitted_full':   fitted_full,
		'fitted_region': fitted_region,
		'x_region':      x_region,
		'y_region':      y_region,
		
	}


def plot_dose_weights(dw: np.ndarray, frame_doses: np.array, pixel_size: float, output_path: str) -> None:
	
	MIN_BIN_FACTOR = 32
	
	n_frames, n_bins = dw.shape
	n_bins_to_ignore = int(n_bins / MIN_BIN_FACTOR)
	n_bins_to_plot = n_bins - n_bins_to_ignore
	
	fig, ax = plt.subplots(figsize=(10, 6))
	ax.imshow(dw[:,n_bins_to_ignore:], vmin=0.0, vmax=1.0, cmap='jet', aspect='auto', interpolation='nearest')
	
	ax.set_title('Empirical dose weights', fontsize=13, fontweight='bold', pad=10)
	ax.set_xlabel('Resolution (Å)', fontsize=11, labelpad=6)
	ax.set_ylabel('Cumulative Exposure (e/Å²)', fontsize=11, labelpad=6)
	
	f_min = 0.5 / (MIN_BIN_FACTOR * pixel_size)
	f_max = 0.5 / pixel_size
	freqs = np.linspace(f_min, f_max, n_bins_to_plot)
	xtick_pos = range(0, len(freqs), int(n_bins_to_plot / 8))
	xtick_res  = [1.0 / freqs[i] for i in xtick_pos]
	ax.set_xticks(np.array(xtick_pos))
	ax.set_xticklabels([f'{r:.2f}' for r in xtick_res])
	
	ytick_pos = range(0, n_frames, int(n_frames / 8))
	ytick_dose  = [frame_doses[i] for i in ytick_pos]
	ax.set_yticks(np.array(ytick_pos))
	ax.set_yticklabels([f'{r:.2f}' for r in ytick_dose])
	
	plt.tight_layout()
	fig.savefig(output_path, dpi=150, bbox_inches='tight')


if __name__ == "__main__":
	
	print('')
	print('fix_cs_dw.py')
	print('Direct Electron LP')
	print('')
	
	parser = argparse.ArgumentParser()
	parser.add_argument('kV', type=int, help='Accelerating voltage')
	parser.add_argument('epA2', type=float, help='Total exposure in e/A2')
	parser.add_argument("use", choices=['O', 'M', 'I', 'G', 'C', 'o', 'm', 'i', 'g', 'c'], help="O = original, M = modified, I = No initial downweighting, G = Grigorieff, C = constant")
	parser.add_argument('--norm', action='store_true', help='Normalize each resolution bin to have max = 1')
	args = parser.parse_args()
	
	OUTPUT_NAME_ORIGINAL = 'original'
	OUTPUT_NAME_MODIFIED = 'modified'
	OUTPUT_NAME_NOINITIALDOWNWEIGHTING = 'noinitialdownweighting'
	OUTPUT_NAME_GRIGORIEFF = 'grigorieff'
	OUTPUT_NAME_CONSTANT = 'constant'
	
	kV = args.kV
	total_dose = args.epA2
	output_name = OUTPUT_NAME_ORIGINAL
	if (args.use == 'M') or (args.use == 'm'):
		output_name = OUTPUT_NAME_MODIFIED
	elif (args.use == 'I') or (args.use == 'i'):
		output_name = OUTPUT_NAME_NOINITIALDOWNWEIGHTING
	elif (args.use == 'G') or (args.use == 'g'):
		output_name = OUTPUT_NAME_GRIGORIEFF
	elif (args.use == 'C') or (args.use == 'c'):
		output_name = OUTPUT_NAME_CONSTANT
	print(f'Accel. voltage (kV):    {kV}')
	print(f'Total dose (e/A2):      {total_dose}')
	print(f'Output type:            {output_name}')
	
	if len(sys.argv) < 4:
		print('Usage: fix_cs_dw.py <kV> <e/A2> <use>')
		exit()
	kV = int(sys.argv[1])
	total_dose = float(sys.argv[2])
	
	grigorieff_a = 0.490
	grigorieff_b = 1.665
	if kV == 100:
		grigorieff_a = 0.312
	elif kV == 120:
		grigorieff_a = 0.328
	elif kV == 200:
		grigorieff_a = 0.368
	elif kV == 300:
		grigorieff_a = 0.490
	else:
		print('ERROR: Invalid kV!')
		exit()
	
	if total_dose < 0.1:
		print('ERROR: Invalid e/A2!')
		exit()
	
	MIN_BIN_FACTOR = 32
	PLOT_FREQUENCY = 16
	MIN_LOG_DW = -20
	
	CS_LOG_FILE_IN = 'job.log'
	CS_HYPERPARAMS_FILE_IN = 'hyperparams.cs'
	CS_DW_FILE_IN = 'refm_empirical_dw.npy'
	
	input_filename = f'{CS_DW_FILE_IN[:-4]}_{OUTPUT_NAME_ORIGINAL}.npy'
	if not os.path.isfile(input_filename):
		input_filename = CS_DW_FILE_IN
	
	print(f'Job log file:           {CS_LOG_FILE_IN}')
	print(f'Hyperparams input file: {CS_HYPERPARAMS_FILE_IN}')
	print(f'DW input file:          {input_filename}')
	
	if not os.path.isfile(CS_LOG_FILE_IN):
		print(f'ERROR: {CS_LOG_FILE_IN} does not exist in the current directory!')
		exit()
	if not os.path.isfile(CS_HYPERPARAMS_FILE_IN):
		print(f'ERROR: {CS_HYPERPARAMS_FILE_IN} does not exist in the current directory!')
		exit()
	if not os.path.isfile(input_filename):
		print(f'ERROR: {input_filename} does not exist in the current directory!')
		exit()
	
	hyperparams_metadata = np.load(CS_HYPERPARAMS_FILE_IN)
	ang_per_pixel = float(hyperparams_metadata[0][3])
	print(f'Angstroms per pixel:    {ang_per_pixel:.2f}')
	
	box_size = 0
	with open(CS_LOG_FILE_IN, 'r') as file:
		for line in file:
			if ('dims' in line) and('ndims' not in line):
				box_size = int(line.strip().split(' ')[-2])
	print(f'Box size (pixels):      {box_size}')
	
	dose_weights = np.load(input_filename)
	n_frames, n_res_bins = dose_weights.shape
	print(f'Frames per movie:       {n_frames}')
	print(f'Resolution bins:        {n_res_bins}')
	dw_file_dtype = dose_weights.dtype
	frame_numbers = np.arange(n_frames)
	res_bin_numbers = np.arange(n_res_bins)
	
	if (ang_per_pixel <= 0.0) or (box_size < 1) or (n_frames < 1) or (n_res_bins < 1):
		print('ERROR: Invalid parameters!')
		exit()
	
	min_log_dw_sharpness = 32.0 / float(n_frames)
	
	print('')
	
	if OUTPUT_NAME_ORIGINAL not in input_filename:
		print(f'Backing up the original {CS_DW_FILE_IN} file...')
		shutil.copy(CS_DW_FILE_IN, f'{CS_DW_FILE_IN[:-4]}_original.npy')
	
	print('Setting up processing...')
	
	frame_doses = (frame_numbers.astype(float) + 1.0) * total_dose / float(n_frames)
	
	dose_weights_modified = np.zeros_like(dose_weights)
	dose_weights_nodownweightinitial = np.zeros_like(dose_weights)
	dose_weights_grigorieff = np.zeros_like(dose_weights)
	dose_weights_constant = np.ones_like(dose_weights)
	colors = cm.viridis(np.linspace(0.1, 0.9, 5))
	
	n_bins_to_ignore = int(n_res_bins / MIN_BIN_FACTOR)
	n_bins_to_plot = n_res_bins - n_bins_to_ignore
	
	f_min = 0.5 / (MIN_BIN_FACTOR * ang_per_pixel)
	f_max = 0.5 / ang_per_pixel
	freqs = np.linspace(f_min, f_max, n_bins_to_plot)
	tick_pos = range(0, len(freqs), int(n_res_bins / 8))
	tick_res  = [1.0 / freqs[i] for i in tick_pos]
	
	print('Processing each resolution bin...')
	for res_bin in range(n_res_bins):
		
		res_bin_angstrom = (2.0 * box_size * ang_per_pixel) / (res_bin + 1.0)
		print(f'  {res_bin}: {res_bin_angstrom:.2f} Å')
		
		dw_original = np.log(dose_weights[:, res_bin])
		dw_original_peak = np.max(uniform_filter1d(dw_original[:int(n_frames / 2)].astype(float), size=2))
		if args.norm:
			dw_original = np.clip(dw_original - dw_original_peak, a_min=None, a_max=0.0)
		
		result = find_linear_negative_region(dw_original)
		dw_fitted = result['fitted_full']
		
		if result['r_squared'] < 0.85:
			print(f'    ERROR: R² = {result["r_squared"]:.4f} is too low!')
			exit()
		
		dw_merged = np.zeros(n_frames)
		for i in range(n_frames):
			if i < result['idx_start']:
				dw_merged[i] = dw_original[i]
			elif i < result['idx_end']:
				dw_merged[i] = (dw_original[i] * (result['idx_end'] - i) + dw_fitted[i] * (i - result['idx_start'])) / (result['idx_end'] - result['idx_start'])
			else:
				dw_merged[i] = dw_fitted[i]
		dw_merged_peakloc = np.argmax(uniform_filter1d(dw_merged.astype(float), size=2))
		if args.norm:
			dw_merged_peak = np.max(uniform_filter1d(dw_merged.astype(float), size=2))
			dw_merged -= dw_merged_peak
		z = min_log_dw_sharpness * (dw_merged - MIN_LOG_DW)
		softplus = np.log1p(np.exp(np.clip(z, -500, 500))) / min_log_dw_sharpness
		dw_merged = np.clip(MIN_LOG_DW + softplus, a_min=None, a_max=0.0)
		dose_weights_modified[:, res_bin] = np.exp(dw_merged)
		
		dw_nodownweightinitial = np.zeros(n_frames)
		dw_nodownweightinitial[dw_merged_peakloc:] = dw_merged[dw_merged_peakloc:]
		dw_nodownweightinitial = np.clip(dw_nodownweightinitial, a_min=None, a_max=0.0)
		dose_weights_nodownweightinitial[:, res_bin] = np.exp(dw_nodownweightinitial)
		
		dw_grigorieff = (-1.0 * frame_doses * ((1.0 / res_bin_angstrom) ** grigorieff_b)) / grigorieff_a
		dw_grigorieff -= dw_grigorieff[dw_merged_peakloc]
		dw_grigorieff[:dw_merged_peakloc] = dw_merged[:dw_merged_peakloc] / 4.0
		z = min_log_dw_sharpness * (dw_grigorieff - MIN_LOG_DW)
		softplus = np.log1p(np.exp(np.clip(z, -500, 500))) / min_log_dw_sharpness
		dw_grigorieff = np.clip(MIN_LOG_DW + softplus, a_min=None, a_max=0.0)
		dose_weights_grigorieff[:, res_bin] = np.exp(dw_grigorieff)
		
		if (res_bin % int(n_res_bins / PLOT_FREQUENCY)) == 0:
			fig, ax = plt.subplots(figsize=(8, 5))
			ax.plot(frame_doses, dw_original, label=OUTPUT_NAME_ORIGINAL, color=colors[0], marker='o', markersize=3)
			ax.plot(frame_doses, dw_fitted, label='fitted', color=colors[1], linestyle=':', linewidth=2)
			ax.plot(frame_doses, dw_merged, label=OUTPUT_NAME_MODIFIED, color=colors[2], linestyle='--', linewidth=2)
			ax.plot(frame_doses, dw_nodownweightinitial, label=OUTPUT_NAME_NOINITIALDOWNWEIGHTING, color=colors[3], linestyle='--', linewidth=2)
			ax.plot(frame_doses, dw_grigorieff, label=OUTPUT_NAME_GRIGORIEFF, color=colors[4], linestyle='-.', linewidth=2)
			ax.set_title(f'Resolution bin {res_bin} ({res_bin_angstrom:.2f} Å)', fontsize=13, fontweight='bold', pad=10)
			ax.set_xlabel('Cumulative Exposure (e/Å²)', fontsize=11, labelpad=6)
			ax.set_ylabel('Empirical Dose Weight (log)', fontsize=11, labelpad=6)
			ax.grid(True, alpha=0.3)
			plt.legend(loc='upper right')
			ax.set_xlim(0, total_dose)
			ylim_min = MIN_LOG_DW * 1.1
			ylim_max = 1.2
			ax.set_ylim(ylim_min, ylim_max)
			plt.tight_layout()
			plt.savefig(f'{CS_DW_FILE_IN[:-4]}_resbin{res_bin:04}.png', dpi=150, bbox_inches='tight')
			plt.close()
	
	print('Plotting some frames...')
	for frame in range(0, n_frames, PLOT_FREQUENCY):
		
		dose = frame_doses[frame]
		print(f'  {frame}: {dose:.2f} e/Å²')
		
		dw_original = np.log(dose_weights[frame, :])
		dw_merged = np.log(dose_weights_modified[frame, :])
		dw_nodownweightinitial = np.log(dose_weights_nodownweightinitial[frame, :])
		dw_grigorieff = np.log(dose_weights_grigorieff[frame, :])
		
		fig, ax = plt.subplots(figsize=(8, 5))
		ax.plot(res_bin_numbers, dw_original, label=OUTPUT_NAME_ORIGINAL, color=colors[0], marker='o', markersize=3)
		ax.plot(res_bin_numbers, dw_merged, label=OUTPUT_NAME_MODIFIED, color=colors[2], linestyle='--', linewidth=2)
		ax.plot(res_bin_numbers, dw_nodownweightinitial, label=OUTPUT_NAME_NOINITIALDOWNWEIGHTING, color=colors[3], linestyle='--', linewidth=2)
		ax.plot(res_bin_numbers, dw_grigorieff, label=OUTPUT_NAME_GRIGORIEFF, color=colors[4], linestyle='-.', linewidth=2)
		ax.set_title(f'Frame Number {frame} ({dose:.2f} e/Å²)', fontsize=13, fontweight='bold', pad=10)
		ax.set_xlabel('Resolution (Å)', fontsize=11, labelpad=6)
		ax.set_ylabel('Empirical Dose Weight (log)', fontsize=11, labelpad=6)
		ax.grid(True, alpha=0.3)
		plt.legend(loc='upper right')
		ax.set_xlim(0, n_res_bins - 1)
		ax.set_xticks(np.array(tick_pos))
		ax.set_xticklabels([f'{r:.2f}' for r in tick_res])
		ylim_min = MIN_LOG_DW * 1.1
		ylim_max = 1.2
		ax.set_ylim(ylim_min, ylim_max)
		plt.tight_layout()
		plt.savefig(f'{CS_DW_FILE_IN[:-4]}_frame{frame:04}.png', dpi=150, bbox_inches='tight')
		plt.close()
	
	print('Plotting heat map of all dose weights...')
	plot_dose_weights(dose_weights, frame_doses, ang_per_pixel, f'{CS_DW_FILE_IN[:-4]}_{OUTPUT_NAME_ORIGINAL}.png')
	plot_dose_weights(dose_weights_modified, frame_doses, ang_per_pixel, f'{CS_DW_FILE_IN[:-4]}_{OUTPUT_NAME_MODIFIED}.png')
	plot_dose_weights(dose_weights_nodownweightinitial, frame_doses, ang_per_pixel, f'{CS_DW_FILE_IN[:-4]}_{OUTPUT_NAME_NOINITIALDOWNWEIGHTING}.png')
	plot_dose_weights(dose_weights_grigorieff, frame_doses, ang_per_pixel, f'{CS_DW_FILE_IN[:-4]}_{OUTPUT_NAME_GRIGORIEFF}.png')
	plot_dose_weights(dose_weights_constant, frame_doses, ang_per_pixel, f'{CS_DW_FILE_IN[:-4]}_{OUTPUT_NAME_CONSTANT}.png')
	
	print(f'Saving dose weighting CryoSPARC files...')
	np.save(f'{CS_DW_FILE_IN[:-4]}_{OUTPUT_NAME_MODIFIED}.npy', dose_weights_modified)
	np.save(f'{CS_DW_FILE_IN[:-4]}_{OUTPUT_NAME_NOINITIALDOWNWEIGHTING}.npy', dose_weights_nodownweightinitial)
	np.save(f'{CS_DW_FILE_IN[:-4]}_{OUTPUT_NAME_GRIGORIEFF}.npy', dose_weights_grigorieff)
	np.save(f'{CS_DW_FILE_IN[:-4]}_{OUTPUT_NAME_CONSTANT}.npy', dose_weights_constant)
	shutil.copy(f'{CS_DW_FILE_IN[:-4]}_{output_name}.npy', CS_DW_FILE_IN)
	
	print('Done.')
