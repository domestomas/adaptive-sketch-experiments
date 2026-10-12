#!/usr/bin/env python3
from relativeErrorSketch import RelativeErrorSketch
import argparse
from collections import namedtuple
import multiprocessing as mp
from multiprocessing import shared_memory
import pickle, os, random
import numpy as np
from pprint import pprint
from copy import deepcopy



Data = namedtuple('Data', ['perc68', 'perc95', 'perc99', 'avg', 'median'])
def run_the_sketch(n, order, k, merge, mem_name=''):
		sketch = RelativeErrorSketch(k=k)
		if mem_name != '':
			mem = shared_memory.SharedMemory(mem_name, track=False)
			stream = iter(np.ndarray((n,), dtype=np.int64, buffer=mem.buf))
		elif order == "sorted":
			stream = range(1, n+1)
		elif order == "reversed":
			stream = range(n, 0, -1)
		elif order == "random":
			stream = random.sample(range(1, n+1), k=n)
		
		# Streaming model
		if not merge:
			final = deepcopy(sketch)
			for item in stream:
				final.update(int(item)) # int(item) is important for PyPy3 to run fast
			
		# Balanced merging
		else:
			stack = list()
			used_items = 0
			while used_items < n:
				new_sketch = deepcopy(sketch)
				i = 0
				while used_items < n and i < 1000:
					i+=1
					used_items+=1
					new_sketch.update(stream.__next__())
				current = (1000, new_sketch)
				while len(stack) > 0 and stack[-1][0] == current[0]:
					last = stack.pop()
					current = (current[0]+last[0], RelativeErrorSketch.merge(last[1], current[1]))
				stack.append(current)
			final = stack.pop()[1]
			while len(stack) > 0:
				final = RelativeErrorSketch.merge(final, stack.pop()[1])
			
		if mem_name != '':
			mem.close()
		return final

def bisect(n, order, space, merge):
	cap = 0
	small = 3
	big = 50
	while big-small > 1:
		avg = (small + big) // 2
		s = run_the_sketch(n, 'sorted', 2*avg, merge)
		cap = sum(c.nomCapacity() for c in s.compactors)
		if cap < space:
			small = avg
		else:
			big = avg
	return 2*avg

def main():
	parser = argparse.ArgumentParser(description=
		'Program for testing Jagged Sketch.'
	)
	parser.add_argument(
		'-logn', type=int, default=6, 
		help='the approximate number of generated elements is 10^logn'
	)
	parser.add_argument(
		'-order', type=str, default='random',
		choices = ['sorted','reversed','random'], 
		help='the order of the streamed integers (in -random case stream must fit in memory)'
	)
	parser.add_argument(
		'-repeat', type=int, default=1,
		help='the number of times to repeat building the sketch (to calculate deviations)'
	)
	parser.add_argument(
		'-info', type=str, default='',
		help='additional info to include in the legend'
	)
	parser.add_argument(
		'-sectionsize', type=int, default=50, 
	)
	parser.add_argument(
		'-space', type=int, default=0, 
	)
	parser.add_argument(
		'-merge', type=bool, default=False,
	)
	args = parser.parse_args()
	
	# parse input
	n = 10**args.logn
	repeat = args.repeat
	order = args.order
	user_info = f"{args.info+'_' if args.info!='' else ''}"
	k = args.sectionsize
	space = args.space
	merge = args.merge
	
	if space != 0:
		k = bisect(n, order, space, merge)
	
	filename = (
		f"rs_{int(n/1000000)}mil_{order}"
		f"_eps{k}_M{merge}{'_'+user_info if user_info!= '' else ''}"
		)
	
	# check files and folders
	os.makedirs("samples", exist_ok=True)
	
	if repeat > 1:
		if os.path.isfile(f"samples/{filename}"):
			exit("file already exists")
		mem_name = ''
		if order == "random":
			rng = np.random.default_rng()
			a = np.array(rng.permutation(np.arange(1,n+1)))
			mem = shared_memory.SharedMemory(create=True, size=a.nbytes)
			random_stream = np.ndarray((n,), np.int64, buffer=mem.buf)
			random_stream[:] = a[:]
			mem_name = mem.name
		
		# run the sketch in paralel "repeat" times
		with mp.Pool() as pool:
			async_runs = [pool.apply_async(
					run_the_sketch, (n, order, k, merge, mem_name)
				) for _ in range(repeat)]
			runs = [x.get() for x in async_runs]
		
		if order == "random":
			mem.close()
			mem.unlink()
		
		s = runs[0]
		
	else:
		s = run_the_sketch(n, order, k, merge)
	
	# get skech info
	sketch_info = {
		"n":n, "repeat":repeat, "cap":sum(c.nomCapacity() for c in s.compactors), 
		"C":max(c.nomCapacity() for c in s.compactors), "H":s.H, 
		"k":k, "order":order, "sketch":"rs","M":merge, "user":user_info,
	}
	
	if repeat == 1:
		pprint(sketch_info)
		for c in s.compactors:
			print(f"sl: {c.sectionSize}, comp: {c.numCompactions} C: {c.nomCapacity()}")
	else:
		# dump the results to file
		ranks = [r.ranks() for r in runs]
		with open(f"samples/{filename}", mode='xb') as file:
			pickle.dump(Sampling(ranks, sketch_info, n), file)

class Sampling():

	def __init__(self, ranks, info, n):
		self.info = info
		self.repeat = len(ranks)
		self.sample_points = self.choose_sample_points(ranks)
		self.data = Data([], [], [], [], [])
		self.n = n
		self.prepare_data(ranks)
	
	def choose_sample_points(self, ranks):
		all_points = [x for run in ranks for (x, _) in run]
		all_points.sort()
		return [all_points[i*self.repeat+1] for i in range(len(all_points)//(self.repeat)-1)]

	def prepare_data(self, ranks):
		repeat = self.repeat
		current_index = [0 for _ in range(repeat)]
		sd_percents = (68*repeat//100, 95*repeat//100, 99*repeat//100)
		
		# find SD values for each point
		# go through all the lists in sync to achieve linear time
		for point in self.sample_points:
			errors = []
			error_sum = 0
			
			# for every run
			for i in range(repeat):
				run = ranks[i]
				curr = current_index[i]
				
				# move a little to find error in given point
				while curr < len(run) and run[curr][0] <= point: 
					curr += 1
				current_index[i] = curr
				
				# calculate error in given point
				err = run[curr-1][1] - point

				# save the error
				error_sum += err
				errors.append(err)
			
			# find an save the confidence intervals
			errors.sort()
			self.data.median.append(errors[len(errors)//2])
			self.data.avg.append(error_sum/repeat)
			errors = [abs(e) for e in errors]
			errors.sort()
			for i in range(3):
				self.data[i].append(errors[sd_percents[i]])

if __name__ == '__main__':
	main()
