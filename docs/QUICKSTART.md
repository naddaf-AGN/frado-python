# FRADO Python in Google Colab (5-minute start)

Nothing needs installing on your own laptop; Colab runs Python remotely. The installation cell below installs packages in Colab only.

### 1. Upload and extract

Open [Google Colab](https://colab.research.google.com/), choose **New notebook**, and execute:

```python
from google.colab import files
uploaded = files.upload()  # Select the downloaded FRADO Python ZIP file.
archive = next(iter(uploaded))
```

```python
!unzip -qo {archive}
%cd /content/frado-python
!pip -q install -r requirements-optimized.txt
```

### 2. Verify the installation

```python
!python -m pytest -q
!python run.py --config examples/optimized_sweep.json --dry-run
```

The first command runs the automated tests. The second prints the requested cloud radii without running the physics.

### 3. Run a small set of clouds

```python
!python run.py --config examples/optimized_sweep.json --output-dir OutPut/demo
```

This example covers 8,370, 8,380, and 8,390 Rg. It uses the accelerated Numba backend. The **first** run may include a compilation delay.

### 4. Inspect and download the results

```python
import pandas as pd
summary = pd.read_csv('OutPut/demo/summary.csv')
display(summary[['launch_radius_rg', 'termination_reason', 'numerical_status']])
```

```python
!python examples/plot_sweep_summary.py OutPut/demo
from IPython.display import Image, display
display(Image(filename='OutPut/demo/sweep_summary.png'))
```

```python
!zip -qr /content/frado_results.zip OutPut/demo
from google.colab import files
files.download('/content/frado_results.zip')
```

**Change physical inputs:** edit a copy of `examples/optimized_sweep.json` or pass CLI overrides (`!python run.py --help`). Use a **different output folder** for each run, because output files in an existing folder can be overwritten. Colab sessions are temporary; download or save results to Drive.

### Local equivalent

From the repository root, use the same commands *without* `!` or `%cd`, in a Python environment with the required packages.
