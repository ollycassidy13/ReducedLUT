# ReducedLUT: Table Decomposition with "Don't Care" Conditions 
[![DOI](https://img.shields.io/badge/DOI-10.1145/3706628.3708823-orange)](https://doi.org/10.1145/3706628.3708823) [![arXiv](https://img.shields.io/badge/arXiv-2412.18579-b31b1b.svg?style=flat)](https://arxiv.org/abs/2412.18579) <a href="https://doi.org/10.5281/zenodo.14499541"><img src="https://zenodo.org/badge/874439825.svg" alt="DOI"></a> 

![ReducedLUT](img/logo.png)

Lookup tables (LUTs) are frequently used to efficiently store arrays of precomputed values for complex mathematical computations. When used in the context of neural networks, these functions exhibit a lack of recognizable patterns which presents an unusual challenge for conventional logic synthesis techniques. ReducedLUT is a tool for the compression of lookup tables and generation of their hardware files in Verilog for RTL designs, as we demonstrated across multiple machine learning applications where don't care conditions were leveraged for greater compression. This project is a derivative work of [CompressedLUT](https://github.com/kiabuzz/CompressedLUT).

This code is part of a publication in the ACM/SIGDA International Symposium on Field-Programmable Gate Arrays 2025, which is available in the [ACM digital library](https://dl.acm.org/doi/10.1145/3706628.3708823) and on [arXiv](https://arxiv.org/abs/2412.18579). We also present a [NeuraLUT-based toolflow](https://github.com/MartaAndronic/NeuraLUT/tree/reducedlut).
> Oliver Cassidy, Marta Andronic, Samuel Coward, and George A. Constantinides. 2025. ReducedLUT: Table Decomposition with "Don't Care" Conditions. In Proceedings of the 2025 ACM/SIGDA International Symposium on Field Programmable Gate Arrays (FPGA '25). Association for Computing Machinery, New York, NY, USA, 36–42. https://doi.org/10.1145/3706628.3708823

## ⚠️ Installation
```bash
git clone https://github.com/ollycassidy13/ReducedLUT
cd ReducedLUT
make
```
> Note: To run the makefile, you will need to have GNU Make and a g++ compiler supporting C++11 on your system
    
## 🌱 Getting Started
#### Lookup Table as a Text File
A text (.txt) file, containing the values of your lookup table, should be prepared as an input. The file must contain a power of 2 lines, each of which is a single hexadecimal value in ascending input order. An example of such a text file can be found in `table.txt`. The following command generates hardware files corresponding to the lookup table described in that text file.

```bash
./reducedlut -table table.txt
```

#### Introducing *Don't Cares*
A text (.txt) file, containing the values of your lookup table, should be prepared as before. Another text file containing the input data should be prepared too, containing the data to be used to determine *don't cares*. This should be arranged as binary values each on a separate line with a bitwidth equal to the table's input bitwidth. A rarity and exiguity threshold should also be specified (for more information on these parameters see the `help.txt` file).

```bash
./reducedlut -table table.txt -input input.txt -exiguity 4 -rarity 1
```

See the `help.txt` file for command line arguments in more detail.

## 🔢 ReducedLUT Flow

This section describes how to evaluate LUT-based neural networks using the ReducedLUT methodology as desribed in our publication. The process described in the publication is based on the [NeuraLUT publication](https://github.com/MartaAndronic/NeuraLUT/tree/reducedlut), where trained models are mapped to a series of L-LUTs. To get started:

1. **Train and Map**  
   First, train a NeuraLUT model following the instructions provided in the [NeuraLUT repository](https://github.com/MartaAndronic/NeuraLUT/tree/reducedlut). The trained model should then be mapped into L-LUTs, one for each neuron.

2. **Inference and Logging**  
   After the L-LUTs are obtained, run an inference pass on your training dataset. During this step, log each neuron's input in text files. 

3. **Run ReducedLUT**  
   Use the logged neuron inputs (text files) as well as each neuron's corresponding L-LUT to run the ReducedLUT script. Each neuron should be processed independently. ReducedLUT will evaluate the L-LUT and input data and provide the compressed L-LUT as a Verilog file.

4. **Results**  
    The model of compressed L-LUTs should be tested before a synthesis and place and route to obtain the final test accuracy, P-LUT and Fmax results.

Below is a diagram illustrating the flow:

![ReducedLUT Flow Diagram](img/flow.jpg)

> Note: The toolflow used in the publication is available on the [ReducedLUT branch of the NeuraLUT repository](https://github.com/MartaAndronic/NeuraLUT/tree/reducedlut).

## 🧪 Summary of Major Modifications
- The novelty of our work is through the integration of *don't cares* within the heart of the LUT decomposition to enhance its capabilities
- ReducedLUT is specifically optimised to implement *don't cares* within machine learning applications
- ReducedLUT takes 3 new inputs - an inputs file, rarity threshold and exiguity parameter 
- ReducedLUT introduces two new vectors storing unique sub-table indices in ascending and descending order of element value
- ReducedLUT integrates a new approach to manipulate the similarity matrix and vector to enhance compression when generating Verilog files
- ReducedLUT uses a new optimised calculation for the use of *don't cares* to determine the optimal compression
- ReducedLUT can store bit totals for a model within a text file for retrieval

## 📖 Citation
Should you find this work valuable, we kindly request that you consider referencing our paper as below:
```
@inproceedings{reducedlut,
author = {Cassidy, Oliver and Andronic, Marta and Coward, Samuel and Constantinides, George A.},
title = "{ReducedLUT: Table Decomposition with ``Don't Care'' Conditions}",
year = {2025},
isbn = {9798400713965},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
note = {doi: 10.1145/3706628.3708823},
booktitle = {Proceedings of the 2025 ACM/SIGDA International Symposium on Field Programmable Gate Arrays},
pages = {36–42},
location = {Monterey, CA, USA},
}
```

## 📜 Copyright & License Notice
It can be freely used for educational and research purposes by non-profit institutions and US government agencies only. Other organizations are allowed to use ReducedLUT only for evaluation purposes, and any further uses will require prior approval. The software may not be sold or redistributed without prior approval. One may make copies of the software for their use provided that the copies, are not sold or distributed, and are used under the same terms and conditions.


As unestablished research software, this code is provided on an "as is" basis without warranty of any kind, either expressed or implied. The downloading, or executing any part of this software constitutes an implicit agreement to these terms. These terms and conditions are subject to change at any time without prior notice.
