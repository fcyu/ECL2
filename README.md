# ECL2
ECL2 is an advanced version of ECL. It has a linear computational complexity, supports multi-thread computation, and multiple variable modifications.

## How to use it?
Requirements: Java 1.8 or later

Usage:
```
java -Xmx128g -jar /path/to/ECL2.jar <parameter_file> <data_file>
```
1. <parameter_file>: parameter file. Can be downloaded along with ECL2.
2. <data_file>: spectra data file (mzXML or MGF).

example: java -Xmx128g -jar ECL2.jar parameter.def data.mzxml

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
  publisher={ACS Publications}
}
```
