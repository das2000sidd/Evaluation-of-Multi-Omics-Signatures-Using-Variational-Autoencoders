# Evaluation of Multi-Omics Signatures Using Variational Autoencoders

## Overview

This project investigates whether multi-omics data can be represented using modality-specific variational autoencoders (VAEs) and subsequently integrated for molecular classification.

Using TCGA-BRC* data, the project evaluates RNA-seq, copy-number variation (CNV), and DNA methylation data for estrogen receptor (ER) status classification.

Each omics modality is first compressed into a low-dimensional latent representation using a VAE. These latent representations are then combined through late latent fusion and evaluated using logistic regression.

The project also compares PCA, autoencoders (AEs), and VAEs as alternative representation-learning approaches.



## Objectives

The main objectives are to:

* Integrate multiple TCGA-BRCA molecular data types.
* Develop modality-specific VAEs for high-dimensional omics data.
* Compare PCA, AE, and VAE representations.
* Evaluate the predictive information contained in each omics modality.
* Investigate whether combining latent representations improves ER-status classification.
* Quantify uncertainty in model performance using bootstrap confidence intervals.
* Compare multimodal models against RNA-only models using paired bootstrap AUC differences.
* Examine individual samples where multimodal integration changes the classification relative to RNA alone.


## Dataset

The analysis uses data from The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA) cohort.

Three molecular modalities were used:

RNA-seq = Gene expression                  |
CNV = Gene-level copy-number variation |
DNA methylation = CpG methylation beta values      |

Clinical ER status was obtained from TCGA clinical data and used as the classification outcome.

ER status was encoded as:

* `0` = ER-positive
* `1` = ER-negative

After intersecting the three omics datasets and clinical information, **285 samples** with valid ER status were available for the VAE analysis.

---

## Data Processing

Samples were aligned across all molecular modalities using TCGA sample identifiers.

A fixed stratified train/validation/test split was used:

| Set              | Samples |
| ---------------- | ------: |
| Training         |     159 |
| Validation       |      40 |
| Independent test |      86 |
| Total            |     285 |

The split was fixed using `random_state=42`.

All preprocessing steps were fitted using the training data only and subsequently applied to validation and test data to avoid information leakage.

### RNA-seq

RNA-seq features were filtered using a training-set variance threshold and standardized before VAE training.

### CNV

CNV data were processed using:

* Removal of features with insufficient non-missing values
* Training-set median imputation
* Variance filtering
* Standardization

### DNA methylation

DNA methylation beta values were used as the third molecular modality and processed using the same train-only preprocessing principle.

---

## Model Architecture

A separate VAE was trained for each omics modality.

The architecture consists of:

```text
Input features
      │
      ▼
Linear layer
      │
     ReLU
      │
      ├───────────────┐
      ▼               ▼
   μ layer        log(σ²) layer
      │               │
      └───────┬───────┘
              ▼
       Reparameterization
              │
              ▼
        16-dimensional
        latent space
              │
              ▼
          Decoder
              │
              ▼
       Reconstructed input
```

Each VAE uses a **16-dimensional latent representation**.

The VAE objective combines reconstruction loss with KL divergence:

```text
VAE loss = reconstruction loss + β × KL divergence
```

The latent representation used for downstream classification is the deterministic **latent mean (`μ`)**, rather than a sampled latent vector.

---

## PCA, AE and VAE Comparison

Before evaluating multimodal fusion, RNA-seq representation learning was compared using:

* Principal Component Analysis (PCA)
* Autoencoder (AE)
* Variational Autoencoder (VAE)

Each method generated a 16-dimensional representation, which was subsequently evaluated using logistic regression.

This provides a comparison between a conventional linear dimensionality-reduction method and two neural-network-based representation-learning approaches.

---

## Multimodal Latent Fusion

After training separate VAEs, the latent representations were combined using **late fusion**.

```text
RNA-seq ─────────► RNA VAE ─────────► 16D latent
                                      │
CNV ─────────────► CNV VAE ──────────► 16D latent ──┐
                                      │             │
Methylation ─────► Methylation VAE ──► 16D latent ─┤
                                                    ▼
                                             Concatenation
                                                    │
                                                    ▼
                                          Logistic Regression
                                                    │
                                                    ▼
                                             ER classification
```

Seven configurations were evaluated:

1. RNA
2. CNV
3. DNA methylation
4. RNA + CNV
5. RNA + methylation
6. CNV + methylation
7. RNA + CNV + methylation

This is **late latent fusion** rather than a joint multimodal VAE.

---

## Results

### Independent Test Set

All models were evaluated on the same independent test cohort of 86 samples.

| Model       | Latent dimensions | Validation AUC | Validation accuracy |   Test AUC | Test accuracy |
| ----------- | ----------------: | -------------: | ------------------: | ---------: | ------------: |
| RNA         |                16 |          1.000 |               97.5% | **0.9937** |         93.0% |
| CNV         |                16 |         0.9749 |               92.5% |     0.9694 |         89.5% |
| Methylation |                16 |          1.000 |               95.0% |     0.9772 |         90.7% |
| RNA + CNV   |                32 |         0.9964 |                97.5 |            |               |

Bootstrap Confidence Intervals

Bootstrap 95% confidence intervals were estimated using 2,000 resamples.

```
Model	Test AUC	95% CI	Accuracy	95% CI
RNA	0.9937	0.9804–1.0000	93.0%	87.2–97.7%
CNV	0.9694	0.9275–0.9964	89.5%	82.6–95.3%
Methylation	0.9772	0.9446–0.9979	90.7%	83.7–96.5%
RNA + CNV	0.9945	0.9812–1.0000	93.0%	87.2–97.7%
RNA + Methylation	0.9898	0.9712–1.0000	94.2%	88.4–98.8%
CNV + Methylation	0.9804	0.9538–0.9986	91.9%	86.1–97.7%
RNA + CNV + Methylation	0.9890	0.9696–1.0000	94.2%	88.4–98.8%
```


Comparison With RNA Alone

RNA provided a highly informative latent representation for ER-status classification.

Paired bootstrap comparisons of test-set AUC were performed relative to RNA alone:

```
Comparison	ΔAUC	95% CI
CNV − RNA	−0.0244	−0.0641 to 0.0026
Methylation − RNA	−0.0165	−0.0421 to 0.0000
RNA + CNV − RNA	+0.0008	−0.0036 to 0.0063
RNA + Methylation − RNA	−0.0039	−0.0146 to 0.0018
CNV + Methylation − RNA	−0.0134	−0.0347 to 0.0020
RNA + CNV + Methylation − RNA	−0.0047	−0.0172 to 0.0030
```

All confidence intervals include zero.

Therefore, within this 86-sample test cohort, the observed differences in AUC do not provide clear evidence of an incremental discrimination benefit from multimodal latent fusion over RNA alone.

The results nevertheless demonstrate that CNV and methylation contain substantial ER-status-related information, while the additional modalities do not consistently improve classification when combined with the already highly informative RNA representation.

Interpretation

Several observations emerge from the analysis:

RNA-seq generated a highly predictive latent representation.
CNV and methylation independently contained useful predictive information.
Combining latent representations produced high-performing classifiers.
Multimodal fusion did not consistently improve AUC over RNA alone.
Differences between the models were small relative to the uncertainty estimated from the test cohort.
Multimodal integration changed predictions for individual samples, providing an opportunity to investigate cases where additional molecular modalities either supported or contradicted the RNA-based prediction.

Importantly, the objective of this project is to demonstrate representation learning and multimodal integration, rather than to develop a clinically validated ER-status classifier.

Per-Sample Analysis

In addition to aggregate performance metrics, individual test samples were examined to determine how multimodal fusion changed predictions relative to RNA alone.

For each sample, the analysis records:

True ER status
RNA probability
Multimodal probability
RNA prediction
Multimodal prediction
Whether each prediction was correct
Change in predicted probability after multimodal fusion

This allows individual rescues and errors introduced by multimodal integration to be examined rather than relying solely on aggregate AUC.

Reproducibility

The analysis uses a fixed random seed and sample-level split to maintain consistent train, validation and test cohorts.

The final test cohort is kept fixed across all seven VAE configurations to enable paired comparisons.

Technologies
Python
PyTorch
scikit-learn
NumPy
pandas
Matplotlib
TCGA/GDC data
Variational Autoencoders
Autoencoders
PCA
Logistic Regression
Bootstrap resampling
Limitations

**This analysis has several important limitations:**

The final multi-omics cohort contains only 285 samples.
The independent test set contains 86 samples.
TCGA-BRCA is a heterogeneous cohort containing multiple breast cancer subtypes.
ER status is used as a classification endpoint rather than a clinical outcome.
The three molecular modalities have substantially different feature structures and dimensions.
The test set has been used for exploratory comparison of multiple model configurations; therefore, the reported test comparisons should not be interpreted as a confirmatory model-selection procedure.
External validation on an independent cohort has not been performed.


**Conclusion**

This project demonstrates an end-to-end workflow for multi-omics representation learning using variational autoencoders.

Modality-specific VAEs generated informative low-dimensional representations from RNA-seq, CNV and DNA methylation data. RNA-seq produced a highly predictive representation for ER-status classification, while CNV and methylation also contained substantial predictive information.

Latent fusion successfully integrated the molecular representations, but did not consistently provide an incremental AUC benefit over RNA alone in this cohort. This highlights an important aspect of multimodal modelling: adding additional data modalities does not necessarily improve predictive performance when one modality already contains strong discriminatory information.

The project provides a practical framework for exploring representation learning, dimensionality reduction and multimodal integration in cancer genomics.
