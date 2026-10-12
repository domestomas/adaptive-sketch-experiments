# Adaptive sketch experiments
This is a proof-of-concept implementatnion of Adaptive sketch. The file `adaptiveSketch.py` contains all the source code of the algoritm; 
`adaptiveTest.py` and `plotting.py` together form a simple testing framework. 
The file `relativeErrorSketch.py` is an exact copy of the proof-of-concept implementation[^1] of ReqSketch writen by one of the authors Pavel Veselý used for comparison and
`relativeTest.py` is an analogy to `adpativeTest.py` for ReqSketch.

## Usage
First run `adpativeTest.py` to do the experiment. As command line parameters, you can control the logarithm of the size `-logn` of the input, 
the `-order` of the input items (sorted, reverse-sorted, random), the number of repetition of the experiment (`-repeat`),
and the `-gap` and `-section_size` parameters from theory. If you set the `space` parameter, 
the section size will be automatically set up for the sketch to take up approximately `space` items -- this is useful
for comparison with ReqSketch. If the option `-merge` is used, the sketch is build by repeatedly merging sketches summarizing 
similar number of items. If `-repeat` is large than 1, a sample file in subfolder *samples* is created; for sufficiently many items $y$ we 
save the average, median and the 68-th, 95-th and 99-th percentile of the error of the querry for rank($y$).

For example, this line

`./adaptiveTest.py -repeat 100 -logn 6 -order random -gap 7 -section_size 50`

creates a file `samples/as_10\^6_random_rep100_budget0_gap7_ss50_time2.0_M\=False`.

The file `relativeTest.py` can be used analogously to test ReqSketch.

Finally we use `./plotting.py` giving it paths to sapmple files to plot the 95-th percentiles of error from these files as a relative error to a log graph.

`./plotting.py samples/as_10\^6_random_rep100_budget0_gap7_ss50_time2.0_M\=False samples/rs_1mil_random_eps50_MFalse`

[^1]: <https://github.com/edoliberty/streaming-quantiles/blob/master/relativeErrorSketch.py>
