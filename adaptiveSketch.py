from random import random
from math import log, ceil
from warnings import warn

GAP_BASE = 5
SECTION_SIZE = 32
SECTION_COUNT_BASE = 6
TIME_FACTOR = 2 # larger time factor should give faster time in larger space
NUM_PROTECTED = 2

def mylog(x):
	return 1 if x < 2 else ceil(log(x, 2))
	
class AdaptiveSketch(list):
	def __init__(self, gap=GAP_BASE, section_size=SECTION_SIZE,time_factor=TIME_FACTOR):
		if not (gap == int(gap) and gap > 0 and
			section_size == int(section_size) and section_size > 0):
			raise ValueError("gap and section size must be positive integers")
		if section_size % 2 != 0:
			raise ValueError("section size must be even")
		if time_factor * SECTION_COUNT_BASE < SECTION_COUNT_BASE + 1:
			raise ValueError("time_factor * section_count_base must be at least section_count_base + 1")
		if gap < 5 or section_size < 32:
			warn("The theoratical guarantees hold only for gap >= 5 and section_size >= 2^5")
		self.section_size = section_size # section size K from theory
		self.gap = gap # the gap parameter G from theory
		self.time_factor = time_factor
		self.max_section_count = 0
		self.max_space = 0 # for comparison in tests
		self.space = 0 # for comparison in tests
		self.append(AdaptiveCompactor(self, 0))
	
	# Merge two sketches into one, destroys the old sketches
	@staticmethod
	def merge(s, t):
		if s.H() < t.H():
			s, t = t, s
		s.gap = max(s.gap, t.gap)
		max_space = max(s.max_space, t.max_space)
		for i in range(t.H()):
			s[i] = AdaptiveCompactor.merge(s[i], t[i], s)
		s.update_parameters(t.H()-1)
		for (h, compactor) in enumerate(s):
			if compactor.is_full():
				if h == s.H() - 1:
					s.append(AdaptiveCompactor(s, h+1))
				s[h+1].extend(compactor.compact())
		s.update_parameters(s.H()-1)
		s.space = sum(len(c) for c in s)
		s.max_space = max(s.space, max_space)
		t.clear()
		return s
	
	def update_space(self, x):
		self.space += x
		self.max_space = max(self.max_space, self.space)
	
	def H(self):
		return len(self)
			
	# Adds new item to the skech
	def update(self, item):
		self[0].append(item)
		self.update_space(1)
		if self[0].is_full():
			first_unchanged = self.compress()
			self.update_parameters(first_unchanged)
			
	# Do the compaction on level zero and possibly on higher levels
	# Return the lowest level that was not compacted
	def compress(self):
		for (h, compactor) in enumerate(self):
			if compactor.is_full():
				if h == self.H() - 1:
					self.append(AdaptiveCompactor(self, h+1))
				self[h+1].extend(compactor.compact())
			else:
				return h
	
	# Update bounding levels and capacities for levels under self[top]
	def update_parameters(self, top):
		ctr_max = self[top].bounding_ctr
		for h in range(top-1, -1, -1):
			if h+self.gap <= self.H()-1 and ctr_max < self[h+self.gap].ctr:
				ctr_max = self[h+self.gap].ctr
			self[h].bounding_ctr = ctr_max
			self[h].update_section_count()
		new_gap = max(self.gap, ceil( 
			mylog(self.max_section_count)/2 + 3 * mylog(mylog(self.max_section_count))
		))
		if new_gap > self.gap:
			self.gap = new_gap
			self.update_parameters(self.H()-1)
	
	# Computes a list of items and their ranks
	def ranks(self):
		ranks_list = []
		items_and_weights = []
		for (h, items) in enumerate(self):
			items_and_weights.extend( (item, 2**h) for item in items )
		items_and_weights.sort()
		cum_weight = 0
		for (item, weight) in items_and_weights:
			cum_weight += weight
			ranks_list.append( (item, cum_weight) )
		return ranks_list

	# Computes cummulative distribution function (as a list of items 
	# and their ranks expressed as a number in [0,1])
	def cdf(self):
		cdf = []
		rank_list = self.ranks()
		_, total_weight = rank_list[-1]
		for (item, cum_weight) in rank_list:
			cdf.append( (item, cum_weight / total_weight) )
		return cdf

	# Returns an approximate rank of value
	def rank(self, value):
		return sum(c.rank(value)*2**h for (h, c) in enumerate(self))
		
	# Returns an input item which is approx. q-quantile 
 	# (i.e. has rank approx. q*self.N)
	def quantile(self, q):
		assert (q >= 0 and q <= 1), f"parameter q must be in [0, 1], but q = {q}"
		desired_rank = q*self.N
		ranks = self.ranks()
		i = 0
		j = len(ranks)
		while i < j:
			m = (i + j) // 2
			(item, rank) = ranks[m]
			if desired_rank > rank:
				i = m + 1
			else: j = m
		(item, rank) = ranks[i]
		return item

class AdaptiveCompactor(list):
	def __init__(self, sketch, h):
		self.num_compactions = 0
		self.sketch = sketch
		self.h = h
		self.ctr = 0
		self.bounding_ctr = 2
		self.markers = list()
		self.offset = 0
		self.section_count = 0
		self.update_section_count()
	
	# Merges two adaptive compactors and returns the result
	@staticmethod
	def merge(a, b, sketch):
		if len(a) < len(b):
			a, b = b, a
		a.extend(b)
		a.num_compactions += b.num_compactions
		a.sketch = sketch
		if len(b.markers) == 0 or (len(a.markers) > 0 and  b.markers[0] < a.markers[0]):
			a.ctr = b.ctr
		a.markers.extend(b.markers)
		a.offset = a.offset^b.offset
		b.clear()
		return a
	
	# Update section_count based on our ctr and bounding level's ctr
	def update_section_count(self):
		self.section_count = mylog(self.ctr**4) + SECTION_COUNT_BASE + 2*max(
			0, ceil(self.ctr - self.bounding_ctr)
		)
		if self.section_count > self.sketch.max_section_count:
			self.sketch.max_section_count = self.section_count
	
	def rank(self, value):
		return sum(1 for v in self if v <= value)

	def is_full(self):
		return len(self) >= (self.section_count + 1) * self.sketch.section_size 
		# ATTENTION: time factor not used to make a fair comparison with ReqSketch
		# return len(self) >= self.section_count * self.sketch.time_factor * self.sketch.section_size
	
	# Returns the canonical marking
	def compute_marking(self):
		m_i = len(self.markers) - 1
		sections_remaining = len(self) // self.sketch.section_size
		marking = []
		while m_i >= 0:
			if self[sections_remaining*self.sketch.section_size-1] < self.markers[m_i]:
				m_i -= 1
				marking.append(sections_remaining-1)
			sections_remaining -= 1
		return marking
	
	# Updates ctr based on given cannonical marking
	def update_ctr(self, marking):
		m_i = len(marking)-1
		marked_before_pivot = 0
		while m_i >= 0 and marking[m_i] == NUM_PROTECTED + marked_before_pivot:
			marked_before_pivot += 1
			m_i -= 1
		ctr = 0
		m_i = 0
		index_diff = 2 * marked_before_pivot + NUM_PROTECTED
		for section in range(0 if len(marking) == 0 else marking[0], marked_before_pivot+1, -1):
			# S[i] is marked
			if m_i < len(marking) and marking[m_i] == section:
				ctr += 2**( ( index_diff - section ) / 2 )
				m_i += 1
			# S[i] is unmarked
			else:
				ctr = max(0, ctr - 2**( (index_diff - section + 1) /2 ))
		ctr += (2**(marked_before_pivot/2) - 1) * (1 + 2**(1/2))
		self.ctr = ctr
	
	# Removes markers and returns the number of items to keep in the buffer
	def prepare_compaction(self, marking):
		keep_index = self.section_count - 1
		m = 0
		while m < len(marking) and marking[m] > keep_index:
			self.markers.pop()
			m+=1
		while m < len(marking) and marking[m] == keep_index:
			self.markers.pop()
			m+=1
			keep_index-=1
		return self.sketch.section_size * (keep_index+1)
	
	def compact(self):
		original_len = len(self)
		self.sort()
		self.markers.sort()
		marking = self.compute_marking()
		self.update_ctr(marking)
		self.update_section_count()
		items_to_keep = self.prepare_compaction(marking) # removes the markers
		
		# ensure we have even number of items in the buffer
		last = None
		if len(self) % 2 == 1:
			last = self.pop()
		
		if (self.h + self.sketch.gap >= self.sketch.H() - 1 or 
				self[items_to_keep] < self.sketch[self.h+self.sketch.gap].markers[0]):
			self.markers.append(self[items_to_keep])
		
		# standard trick for reducing randomness via anticorellation
		# random choice done every other time
		if self.num_compactions % 2 == 1:
			self.offset = 1 - self.offset
		else:
			self.offset = int(random() < 0.5)
			
		# yield half of non-protected and delete all of them from self
		for i in range(items_to_keep + self.offset, len(self), 2):
			yield self[i] # yield selected items
		del self[items_to_keep : ]
		if last is not None:
			self.append(last)
		self.num_compactions += 1
		self.update_ctr(self.compute_marking())
		self.sketch.update_space((len(self) - original_len)//2)
