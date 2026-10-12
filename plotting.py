#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
from pprint import pformat
import sys, os, pickle, subprocess
from adaptiveTest import Sampling, Data

def main():
	if sys.argv[1] == "-":
		files = sys.stdin.read().rstrip().splitlines()
	else:
		files = sys.argv[1:]
	Plotting.plot(files)

class Plotting():
	def plot(files):
		rel_fig, rel_ax = plt.subplots()
		colours = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'grey', 'orange', 'purple', 'brown', 'lightgreen']
		
		# plot all the datasets
		for (i, file) in enumerate(files):
			Plotting.add_dataset(file, rel_ax, colours[i])
		
		rel_ax.set_xscale('log')
		
		# set labels and lims
		rel_ax.set_xticks([10**x for x in range(10)])
		rel_ax.set_xlabel("Rank", fontsize=10)
		rel_ax.xaxis.set_label_coords(.5, .04)
		rel_ax.yaxis.set_label_coords(.04, .5)
		rel_ax.set_ylabel("Relative Error", fontsize=10)
		rel_ax.legend(
			fontsize=3, 
			loc='upper left', 
			bbox_to_anchor=(.05, .95),
			labelspacing=1
		)
		rel_ax.tick_params(axis='both', which='major', labelsize=8)
		rel_ax.tick_params(axis='y', which='major', labelrotation=60)
		
		# Create offset transform by 5 points in x direction
		dx = 2/72.; dy = 9/72. 
		offset = matplotlib.transforms.ScaledTranslation(dx, dy, rel_fig.dpi_scale_trans)
		# apply offset transform to all x ticklabels.
		for label in rel_ax.yaxis.get_majorticklabels():
			label.set_transform(label.get_transform() + offset)
		
		# save the graph
		rel_fig.savefig("rel_log_temp.pdf")
		
		subprocess.run(["pdfcrop", "rel_log_temp.pdf", "rel_log.pdf"])
		os.remove("rel_log_temp.pdf")
		
	def add_dataset(file, rel_ax, colour):
		with open(file, mode='rb') as file:
			samples = pickle.load(file)
		points = samples.sample_points
		data = samples.data
		info = samples.info
		n = samples.n
		
		# plot relative
		x, y = Plotting.sparsify(points, [data.perc95[i]/max(points[i], 1) for i in range(len(points))])
		rel, = rel_ax.plot(x, y, c=colour, linestyle="-", lw=0.5)
		
		# set legend
		legend = pformat(info)
		rel.set_label(legend)
	
	def sparsify(x_in, y_in, s=50):
		x = []
		y = []
		i = 0
		while i <= len(x_in) - s:
			j, y_j = max(enumerate(y_in[ i : i+s ]), key=lambda t: t[1])
			x.append(x_in[i+j])
			y.append(y_j)
			i+=s
		return (x, y)

if __name__ == '__main__':
	main()
