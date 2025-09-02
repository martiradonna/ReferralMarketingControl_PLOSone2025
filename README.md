# ReferralMarketingControl\_PLOSone2025

This repository contains the MATLAB code associated with the manuscript:

> **Lacitignola D., Martiradonna A. (2025, forthcoming)**
> *Can we enhance trust in the circular economy through referral marketing control? PLOS One*

It includes the **model implementation**, **simulation scripts**, and **figure generation** for full reproducibility of the results presented in the manuscript.

---

## 🧱 Overview

The manuscript proposes a dynamic model of **trust** and **behavior diffusion** in a circular economy system, enhanced via **Z-control mechanisms** on two subpopulations:

* **Broadcasters** \$b(t)\$: agents that spread trust
* **Inerts** \$i(t)\$: non-active participants

The control strategy is evaluated under deterministic and stochastic settings, both in **forward** and **backward** policy scenarios.

This repository provides the full set of MATLAB codes used to simulate these scenarios and produce the manuscript figures.

---

## 📂 Code List and Figures

| Script Filename                       | Description                                                               | Figure    |
| ------------------------------------- | ------------------------------------------------------------------------- | --------- |
| `simulate_uncontrolled_model.m`       | Simulates uncontrolled model (forward and backward scenarios)             | Figs. 2–3 |
| `simulate_zcontrolled_model_on_b.m`   | Z-control applied to broadcasters \$b(t)\$                                | Fig. 4    |
| `simulate_zcontrolled_model_on_i.m`   | Z-control applied to inerts \$i(t)\$                                      | Fig. 5    |
| `simulate_sensitivity_GT_vs_alpha.m`  | Sensitivity of final gain \$G(T)\$ vs. \$\alpha\$ and % variation         | Fig. 6    |
| `simulate_uncontrolled_stochastics.m` | Stochastic evolution of state variables without control                   | Fig. 7    |
| `simulate_G_zcontrol_stochastic.m`    | Stochastic evolution of gain \$G(t)\$ under Z-control (backward scenario) | Fig. 8    |

All scripts generate `.eps` and `.tiff` versions of the figures in high resolution.

---

## ▶️ How to Run

1. **Clone the repository**:

   ```bash
   git clone https://github.com/<your-username>/ReferralMarketingControl_PLOSone2025.git
   cd ReferralMarketingControl_PLOSone2025
   ```

2. **Open MATLAB**, navigate to the cloned folder, and run the desired script.
   For example, to simulate the Z-control on broadcasters:

   ```matlab
   run('simulate_zcontrolled_model_on_b.m');
   ```

3. **Output**: The script will automatically generate and save the corresponding figure(s) in the current MATLAB working directory.

---

## 🛠 Requirements

* **MATLAB R2021a or later**
* No specific toolboxes required
* Scripts may run with minor modifications on **GNU Octave**
  *(Note: some figure export options may differ)*

---

## 🚪 License

All scripts are released under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.

You are free to:

* **Share** — copy and redistribute the material in any medium or format
* **Adapt** — remix, transform, and build upon the material for any purpose, even commercially

**Under the following terms**:

* **Attribution** — You must give appropriate credit by citing the original manuscript (see below)

📄 [Full license text](https://creativecommons.org/licenses/by/4.0/)

---

## 📚 Citation

If you use any script or result from this repository, please cite the following manuscript:

> **Lacitignola D., Martiradonna A. (2025)**
> *Can we enhance trust in the circular economy through referral marketing control?*
> *Manuscript submitted to PLOS ONE*, Ms. No. PONE-D-25-26382

**BibTeX**:

```bibtex
@unpublished{LacitignolaMartiradonna2025,
  author    = {Lacitignola, Deborah and Martiradonna, Angela},
  title     = {Can we enhance trust in the circular economy through referral marketing control?},
  note      = {Manuscript submitted to PLOS ONE, Ms. No. PONE-D-25-26382},
  year      = {2025}
}
```

---

## 👥 Authors

* **Deborah Lacitignola** – Università di Cassino e del Lazio Meridionale, Cassino (FR), Italy
* **Angela Martiradonna** – Università di Foggia, Foggia, Italy
  *All MATLAB code and simulations authored by Angela Martiradonna.*

📬 For academic correspondence: **[angela.martiradonna@unifg.it](mailto:angela.martiradonna@unifg.it)**

---
