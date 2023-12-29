# SHGOFAST
A monkey-patched Cython/C++ backend for bottleneck functions and classes
of the `scipy.optimize.shgo` code using the simplicial sampling method.

## Install
The module can be installed using your preferred `pip` command from the
source archives root directory (where this README resides).

The installation requires a reasonably modern C++ compiler (tested with g++
13.1.1) and the boost software library (tested with version 1.81.0). If boost
is not available, defining the C++ macro `STD_SET_ONLY` will use only standard
library components.

## Use
Simply import the `patch` and `unpatch` functions from `shgofast`:
```python
from scipy.optimize import shgo
from shgofast import patch
patch()

# This should run faster now if cost_fun is not the bottleneck:
shgo(cost_fun, sampling_method='simplicial', **kwargs)
```

## Improvement
From a one-shot measurement, execution time of identical `shgo` calls has
been reduced from 13 seconds to 2.5 seconds for a four-dimensional problem
in which `shgo` was not limited by the function but foremost by the
`VertexCacheField.proc_minimisers` method.


## Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### [0.1.0] - 2023-12-29
- Initial version.
