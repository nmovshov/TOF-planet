# Set workspace (paths and/or common variables).

import sys, os
import numpy as np
if sys.platform == 'win32':
    import matplotlib.pyplot as plt
    plt.style.use('ndefault')
    from IPython import get_ipython
    get_ipython().run_line_magic('matplotlib','')

sys.path.append(os.getcwd())

import TOFPlanet
