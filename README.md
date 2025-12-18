# Varlociraptor Visualization

Three Python functions to visualize Varlociraptor VCF records.

## Features

- **Event Probabilities**: Visualize PROB_* fields as bar charts
- **Allele Frequency Distribution**: Scatter plot with ML frequency highlighted
- **Observations**: Visualization of read-level evidence

## Installation

Using [pixi](https://pixi.sh):
```bash
# Install pixi
curl -fsSL https://pixi.sh/install.sh | bash

# Clone repository
git clone git@github.com:YOUR_USERNAME/varlociraptor-viz.git
cd varlociraptor-viz

# Install dependencies
pixi install
```

## Usage
```python
import pysam
from varlociraptor_viz import (
    visualize_event_probabilities,
    visualize_allele_frequency_distribution,
    visualize_observations
)

# Open VCF file
vcf = pysam.VariantFile("variants.vcf")
record = next(vcf)

# Create visualizations
chart1 = visualize_event_probabilities(record)
chart1.save("events.html")

chart2 = visualize_allele_frequency_distribution(record, "sample_name")
chart2.save("afd.html")

chart3 = visualize_observations(record, "sample_name")
chart3.save("obs.html")
```

## Running Tests
```bash
# Run all tests
pixi run test

# Run with coverage
pixi run test-cov
```

## Functions

### 1. visualize_event_probabilities(record)
Visualizes posterior probabilities for all events (PROB_* INFO fields).

### 2. visualize_allele_frequency_distribution(record, sample_name=None)
Visualizes allele frequency distribution (AFD field).

### 3. visualize_observations(record, sample_name=None)
Visualizes read observations (OBS field).

## License

MIT

## Author

Christy (christykanto1998@gmail.com)
