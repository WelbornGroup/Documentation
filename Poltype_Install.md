# Poltype Installation Guide

To get started you will need to download [Anaconda](https://www.anaconda.com/download#) (Linux Distribution) to your home directory on ARC

- To dowload it to your working directory use `wget "link"` using the link to the anaconda install 

- Use `chmod +x` to make the install script executable, run the install script

- Towards the end of the install it may prompt you to initialize anaconda for your bash shell: `conda init bash`
  This will modify your bashrc script to allow ARC to activate your base conda environment on start up of a new terminal. If you opt out of this step you will need to specify the absolute path to anaconda in the future

- For now activate your conda environment with `source activate` (or if that does not work `conda activate`) and you should see your base environment is activated


Once Conda is running change directories to: `/projects/welbornlab/Poltype2/master/Environments`

Now we will create our Poltype environment with the provided .yml file 

`conda env create -f environment-psi417.yml`

This will take a little while as the dependencies are downloading. Check the environment list to see if successful

`conda env list`

You can also check by activating the environment through `source activate poltype_psi417` (make sure to switch back to your base environment through `source deactivate` for the next step) 

Create a new environment and activate that environment to install the necessary software

`conda create --name xtbenv`

`conda activate xtbenv`

`conda install -c conda-forge xtb -yes`

This conda env will need its path 


