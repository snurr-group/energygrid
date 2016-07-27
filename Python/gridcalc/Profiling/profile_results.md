# Optimization notes
Ben Bucior, July 26, 2016

## General
line_profiler/kernprof is an annoying tool to use in code lacking functions.  Let's try a slightly newer package, pprofile.  It does not require the code annotations (which makes it slower), so it just needs the right bash invocation to prevent it from profiling all of the submodules (including numpy, etc.):

`pprofile --out test.out --exclude-syspath --exclude "aprdf*" scotty_map.py`

Note: exports to the changed directory instead of the root, which is also probably why the code doesn't export correctly to the profiling results


## Trials and Improvement
### Trial 1: Remove extra cell transformation
Baseline for these profile trials, based on the version sent to Arun and Scotty (16e9985 or c680c14).  The difference is that an erroneous fractional to Cartesian conversion has been removed.

### Trial 2: Replace six of the subtractions with simpler statements
As expected, this gives a speed boost by removing 4 of the 6 statements in the inner loop.  Since the calculation still takes time, this gives a 10-15% speed boost overall.  Removing loops or moving statements within them will be much more lucrative overall.  Metric.txt is also still consistent with the previous trial.

### Trial 3: Switch order of grid and atom loops
Move the "per-atom" calculations out of the inner-most loop to outside of the nx,ny,nz loop.  The rationale is that the atom properties do not need to be recalculated 20x20x20 times.  Since there are 4 statements moved, we would expect in the ballpark of 25-30% speed improvement.  As expected, ~24% improvement from Trial 2.  Including the inevitable performance penalty from using the profiler, the code takes about 16 minutes to run in its current iteration.

### Trial 4: See if the notation for assignment makes a difference
Test if there's a difference between something like `a=a+b` and `a+=b` for performance.  If there is, I'd only expect a few percent.  But apparently this is important for larger arrays to [avoid implicit copying](http://ipython-books.github.io/featured-01/).  Not really any change (results were slightly higher).  Overall, the code still shows ~99% of the execution time consumed by looping over the atoms (17,280,000 calls).

## Conclusion
With these minor changes to the code, I was able to reduce the profiled execution time from 1434 seconds to just under 1000.  Huge improvements would probably come by either writing the inner loop in C or turning towards GPU acceleration. (might be several libraries to do that from within Python)
