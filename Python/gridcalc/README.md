# gridcalc
Generate a LJ potential energy surface 

## Dependencies
numpy, PyCIFRW or Pymatgen, pyevtk, Python 2.7

## Pymatgen notes
Pymatgen has compilation issues on Windows, so the code will automatically use a modification of the PyCIFRW package if Pymatgen does not work.  Pymatgen's cifio submodule also does not install cleanly with `pip install pymatgen` on a clean Ubuntu 16.04 for some reason.

As a temporary debugging workaround, one can use the [official docker test container](http://pymatgen.org/index.html#sample-docker-version).  This container skips pyevtk, so relevant lines must be removed from the source before usage.  To install docker, use the `docker.io` package from the Canonical repo in ubuntu.  [Add your user](https://github.com/docker/docker/issues/17645) to the docker group for the correct permissions (then log out).  After installing with `docker pull materialsvirtuallab/pymatgen`, run the container from your source code directory by `docker run -t -i -v $PWD:/opt/research materialsvirtuallab/pymatgen`.

You can then use pymatgen from within the docker environment by using the following commands in the iPython console:

```{python}
import os
os.chdir("/opt/research")
run no_evtk.py
```


