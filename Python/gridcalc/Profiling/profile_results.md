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

### Trial 5 (7/29): Encapsulating the loop within a function
Something I read suggested ~20% speed boosts by keeping critical code within functions (avoiding global variable lookups).  I doubt that, but I need to do this step anyway for the next Trial, so might as well see if it makes a difference.  The performance actually looks about 10% worse (overhead from translating the variables to the function?  low memory on my laptop?)

Using the `time` function in bash, this version took approximately 3m38s to run.

### Trial 6: Boost performance of the lj function with numba
I've heard great things about numba in its ease of use and orders-of-magnitude performance gain.  See [this blog post](http://quant-econ.net/py/need_for_speed.html) and the [official numba documentation](http://numba.pydata.org/) for usage instructions.  Basically, it looks like you just add an `from numba import jit` statement and `@jit` annotations wherever they're needed (or manually compile the functions without decorators).

Again using the `time` function, this version took approximately 2m44s to run (consistent between two runs).

### Trial 7: Try implementation with Cython (TODO)
See <http://nealhughes.net/cython1/> for more information.  It looks like there's some awesome visualization tools included for profiling the code.


## Conclusion 7/26
With these minor changes to the code, I was able to reduce the profiled execution time from 1434 seconds to just under 1000.  Huge improvements would probably come by either writing the inner loop in C or turning towards GPU acceleration. (might be several libraries to do that from within Python)
