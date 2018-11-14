# ECL2
ECL2 is an advanced version of [ECL](https://github.com/fcyu/ECL). It has a linear computational complexity. It also supports multi-thread computation and multiple variable modifications.

## Executable file:
Download the zipped file from https://github.com/fcyu/ECL2/releases/latest.

## How to use it?
Requirements: Java 1.8.

Usage:
```
java -Xmx128g -jar /path/to/ECL2.jar <parameter_file> <data_file>
```
1. <parameter_file>: parameter file. Can be downloaded along with ECL2.
2. <data_file>: spectra data file (mzXML or MGF).

Example: java -Xmx128g -jar ECL2.jar parameter.def data.mzxml

## An example of the result file
| scan_num | spectrum_id | spectrum_mz | spectrum_mass | peptide_mass | rt   | C13_correction | charge | score    | delta_C  | ppm      | peptide                                                            | protein                 | protein_annotation_1                                                          | protein_annotation_2                                                                       | e_value  | q_value | mgf_title |
|----------|-------------|-------------|---------------|--------------|------|----------------|--------|----------|----------|----------|--------------------------------------------------------------------|-------------------------|-------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|----------|---------|-----------|
| 47330    | 47330       | 1215.974    | 3642.894      | 3642.889     | 7401 | -2             | 3      | 2.136501 | 0.112479 | 1.352205 | n[34.063]QLTEMKGHETK[34.063]c-6-n[34.063]EYKLTYYTPEYETK[34.063]c-3 | AT5G37830.1-ATCG00490.1 | Symbols: OXP1 oxoprolinase 1 chr5:15056635-15060525 REVERSE   LENGTH=1266 | Symbols: RBCL ribulose-bisphosphate carboxylases chrC:54958-56397   FORWARD LENGTH=479 | 2.79E-06 | 0       |           |

### Explainations:
1. `scan_num`: Scan num (starts from 1) of the spectrum.
2. `spectrum_id`: Spectrum ID.
3. `spectrum_mz`: Observed spectrum precursor MZ value.
4. `spectrum_mass`: Observed neutral mass of the spectrum.
5. `peptide_mass`: Theoretical mass of the identified peptide.
6. `rt`: Retention time (in seconds) of the spectrum. Only applicable to mzXML format.
7. `C13_correction`: The number of C13 isotope correction.
8. `charge`: Precursor charge.
9. `score`: Identification score. It's similar to XCorr.
10. `delta_C`: The relative difference between the top score and the second best score.
11. `ppm`: PPM between `spectrum_mass` and `peptide_mass`.
12. `peptide`: identified peptide sequence.
13. `protein`: Protein accessions.
14. `protein_annotation_1`: The protein annotation of the first chain.
15. `protein_annotation_2`: The protein annotation of the second chain.
16. `e_value`: e-value of the peptide-spectrum match.
17. `q_value`: q-value (used as FDR) of the peptide-spectrum match.
18. `mgf_title`: MGF title if the spectral file is in MGF format.

## Cite
Yu, Fengchao, Ning Li, and Weichuan Yu. "Exhaustively Identifying Cross-Linked Peptides with a Linear Computational Complexity." Journal of Proteomics Research 16.10 (2017): 3942-3952.

[bibtex]
```bibtex
@article{yu2017exhaustively,
  title={Exhaustively Identifying Cross-Linked Peptides with a Linear Computational Complexity},
  author={Yu, Fengchao and Li, Ning and Yu, Weichuan},
  journal={Journal of Proteome Research},
  volume={16},
  number={10},
  pages={3942--3952},
  year={2017},
  publisher={ACS Publications},
  doi={10.1021/acs.jproteome.7b00338}
}
```
