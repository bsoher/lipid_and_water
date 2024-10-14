# Requirements

**First, if a GPU is available:**

Install the specific versions of PyTorch and CUDA using conda. Replace `<cuda_version>` with your CUDA version (download and install CUDA from [here](https://developer.nvidia.com/cuda-11-7-0-download-archive)). Version 11.7 is the highest version supported by PyTorch 1.13.1.

`conda install pytorch==1.13.1 pytorch-cuda=<cuda_version> -c pytorch -c nvidia`

replacing <cuda_version> with your cuda (download and install from: https://developer.nvidia.com/cuda-11-7-0-download-archive), 11.7 is the highest version supported by torch 1.13.1.

**Second, install the rest of the requirements:**

Install the remaining requirements using pip: `pip install -r requirements.txt`

## Generating Datasets

To generate a dataset for testing:

`python generate_datasets.py -o <output filename.csv> --size <number of simulated spectra>`

The generated list of spectra will be used for benchmarking.

## Benchmarking time efficiency TensorFit and TDFDFit

The `run_tensorfit_times.py` script runs both TensorFit and TDFDFit, returning the time taken for each.

`python run_tensorfit_times.py -i <simulation filename> --size <number of spectra to fit>`

Optional flag: `--cpu`, to force the fitting to run on CPU when a GPU is available

## Benchmarking Accuracy

To run a benchmark on accuracy, use the `test_acc_simulated_data.py` script. It will compute the accuracy for both TensorFit and TDFDFit on the simulated dataset.

`python test_acc_simulated_data.py -i <simulation filename> --size <number of spectra to fit>`

## Citation:

If you find this repository useful for your research or work, please cite: 

**Turco F, Capiglioni M, Weng G, Slotboom J. TensorFit: A torch-based tool for ultrafast metabolite fitting of large MRSI data sets. Magn Reson Med. Published online March 12, 2024. doi:10.1002/mrm.30084**
