# Paper code

Data and code to present the workflow and results of combinatorial mutation space exploration

## Software requirements

Anaconda can be download by https://www.anaconda.com/products/individual

## Package requirements

The package requirements can be found in the `fp4com.yaml` file.

To make a conda environment and  launch a Jupyter Notebook, execute the following commands in a bash shell:

```bash
conda env create -f fp4com.yaml #create the conda environment named fp4com
conda install -n fp4com ipykernel #install ipykernel
python -m ipykernel install --name=fp4com # put the conda environment named fp4com in the jupyter notebook
jupyter notebook #open the jupyter notebook
```
Then all you have to do is use the Jupyter notebook with the fp4com conda environment and execute the notebook `FP4COM.ipynb`.

## Usage
`FP4COM.ipynb` -Jupyter notebook with our codes. The notebook has been used to demonstrate the process of how to perform combinatorial mutation activity prediction using the FFT_PLSR model in the article. Also, it introduces the execution process of the four related tasks described in the article and presents the model's prediction results.

## Directories
`src`- Directory for storing the required Python  files.

`input`- Directory for necessary data files.
* `data.xlsx`  : Contain the training sets for model construction.
* `seq_IFRS.txt` : Include the amino acid sequences of IFRS.
* `seq_Com1.txt` : Include the amino acid sequences of Com1-IFRS.
* `com1_11520_multi_m.csv`: The mutation space sequences generated during the model construction using data based on Com1-IFRS.
Since code generation is time-consuming, the complete mutation sequence file is provided here.

`data`- The folder contains data files formatted in compliance with the EnzymeML standard, which include the relative activity data for all combinatorial mutations.

`output`- The data files saved  in the Jupyter notebook can be found in this folder.


##  System Requirements

The notebook is runs in Windows.

## Hardware requirements

Only a standard computer is required.


## Contact 
If you have any questions, please contact us by email: `yuhaoran@zju.edu.cn`

## License

This project is covered under the MIT licence. See `LICENSE` for more details.  

[Back to top](#TOP)