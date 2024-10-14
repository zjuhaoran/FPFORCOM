# FPFORCOM

Data and code to present the workflow and results of combinatorial mutation space exploration

## Package requirements

Simply install the following packages:
- Python >= 3.9.15
- numpy >= 1.23.5
- pandas >= 1.5.2
- scikit-learn >= 1.2.0
- scipy >= 1.10.0
- aaindex >= 1.0.5
- tqdm >= 4.64.1
- matplotlib >= 3.6.2

## Usage
`main.ipynb` -Jupyter notebook with our codes. The notebook has been used to demonstrate the process of how to perform combinatorial mutation activity prediction using the FFT_PLSR model in the article. Also, it introduces the execution process of the four related tasks described in the article and presents the model's prediction results.

## Directories
`src`- Directory for storing the required Python  files.

`input`- Directory for necessary data files.
* `data.xlsx`  : Contain the training sets for model construction.
* `seq_IFRS.txt` : Include the amino acid sequences of IFRS.
* `seq_Com1.txt` : Include the amino acid sequences of Com1-IFRS.
* `com1_11520_multi_m.csv`: The mutation space sequences generated during the model construction using data based on Com1-IFRS.
Since code generation is time-consuming, the complete mutation sequence file is provided here.

`output`- The data files saved  in the Jupyter notebook can be found in this folder.




##  System Requirements

The notebook is runs in Windows.

## Hardware requirements

Only a standard computer is required.

## Software requirements

Anaconda can be download by https://www.anaconda.com/products/individual

## Contact 
If you have any questions, please contact us by email: `yuhaoran@zju.edu.cn`

## License

This project is covered under the MIT licence. See `LICENSE` for more details.  

[Back to top](#TOP)
