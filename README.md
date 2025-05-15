# SCAN: Sparse reCovery trAnsmitter detectioN

**SCAN** is an unsupervised framework for detecting and characterizing transmitters from power spectral density (PSD) traces using sparse dictionary coding. It accurately separates multiple overlapping transmitters and recovers their time-frequency activity, even under low signal-to-noise ratio (SNR) conditions. This repository contains the official MATLAB implementation of SCAN, as described in our [IEEE INFOCOM 2025 paper](https://doi.org/10.5281/zenodo.14618294).

---

## ✨ Features

- Unsupervised transmitter detection (no labeled data required)
- Robust to noise, overlap, and low SNR
- Supports narrowband and short-lived signals
- Validated on synthetic, controlled, and over-the-air traces

---

## 🧠 Method Overview

SCAN models PSD matrices as a sum of sparse outer products between temporal and frequency dictionary atoms. It operates iteratively, performing:
1. Greedy frequency band selection via dictionary matching,
2. Sparse temporal modeling via convex optimization,
3. Residual masking to iteratively isolate and extract transmitters.

See **Section V** of the [paper](https://doi.org/10.5281/zenodo.14618294) for technical details.

---

## 🚀 Getting Started

### 📦 Requirements

- MATLAB R2022b or later
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

### 📥 Input Format

The input to SCAN is a **2D PSD matrix** `X` of size **(time × frequency)**. Rows represent time steps, and columns represent frequency bins.

---

## ▶️ Running the Demo

1. Clone the repository:
   ```bash
   git clone https://github.com/okoroandrew/SCAN-code.git
 2. Open MATLAB and Navigate to SCAN-code/ directory
 3. run demo.m

### 📤 Output

Scan produces 2 outputs:
1. **Omega**: A binary matrix the same size as the input PSD. Zeros indicate the presence of a detected transmitter; ones indicate background.
2. **pred**: A cell array where each element is a binary matrix representing one detected transmitter (same dimensions as input). The number of elements corresponds to the number of transmitters detected.

### Notes and Tips

1.	**Input length**: SCAN expects the number of time steps (rows) to be a power of 2. If it is not, it automatically uses the first 2^N rows.
2.	**Input frequency size**: For best performance on machines without a GPU, limit the number of frequency bins (columns) to 512 or 256. Above 1024 will be slow.
3.	**Precomputed dictionaries**: The precomputedH folder includes precomputed frequency dictionaries (H) for common sizes. If a required size is missing, SCAN will compute it automatically. You can optionally save new ones for reuse.
4.	**Speed vs. accuracy**: In scan_wrapper.m, the beta parameter controls the optimization speed. For large datasets, set beta between 0.2 and 0.5 to reduce runtime. Note that lower beta values may reduce detection accuracy.


### 📚 Citation
If you use this code, please cite:
```
@inproceedings{okoro2025scan,
  title={SCAN: Sparse reCovery trAnsmitter detectioN},
  author={Okoro, Blessing and McNeil, Maxwell and Meka, Kavya and Doke, Karyn and Bogdanov, Petko and Zheleva, Mariya},
  booktitle={IEEE INFOCOM},
  year={2025}
}
```

### 📬 Contact

For questions, please contact Blessing Andrew Okoro at aaokoro@albany.edu


