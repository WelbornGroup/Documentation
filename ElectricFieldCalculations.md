# Electric field calculations from Tinker (AMOEBA force field)

Open a Terminal window and log onto Cascades. You will need to install two pieces of codes to run electric field calculations. Create a new directory in your home directory to separate these elements from the rest of your files. This can be done doing:

```sh
ssh username@cascades1.arc.vt.edu
mkdir ElectricFields
```
In this case, you will keep the codes in the directory named "ElectricFields", located at `/home/username/ElectricFields/`

### Installation of MDI-enabled Tinker on Cascades
```sh
module load Anaconda/5.2.0
module load cmake
conda create -n tinker_dependencies -c conda-forge numpy pandas texinfo matplotlib
```
You will be prompted to confirm you want to proceed. Press "Y" or type enter. Finally, activate the numpy, panda and texinfo packages with:

```sh
source activate tinker_dependencies
```
The source code is available on Github, on an [MDI-enabled fork](https://github.com/taylor-a-barnes/Tinker) of Tinker. This can be acquired using:

```sh
cd ElectricFields 
git clone --branch mdi-ef https://github.com/taylor-a-barnes/Tinker.git
```


 
 ```sh
cd Tinker/dev
./full_build.sh
 ```
 Your compiled files (fftw, mdi and the Tinker source) are now in `/home/username/ElectricFields/Tinker/build/tinker/`. Only the "dynamic.x" executable is required by the driver. Check that the file is there by doing:
 
 ```sh
 ls /home/username/ElectricFields/Tinker/build/tinker/source/dynamic.x
 ```
 
 Now that you have downloaded and compiled MDI-enabled Tinker, you will need to download and compile the code to do electric field analysis. 
 
 ### Installation of the electric field analysis code on Cascades
 Download the [electric field analysis code](https://github.com/taylor-a-barnes/MDI_EF_Analysis) from Github by doing:
 
 ```sh
 cd /home/username/ElectricFields
 git clone https://github.com/taylor-a-barnes/MDI_EF_Analysis.git
 ```
 Build it using CMake:
 
 ```sh
 cd /home/username/ElectricFields/MDI_EF_Analysis
 cmake .
 make
 ```
 You will next need to tell the driver the location of the files you compiled. Edit:
 - the file `/home/username/ElectricFields/MDI_EF_Analysis/test/locations/MDI_EF_Analysis` with the full path to the "MDI_EF_Analysis.py" Python script (`/home/username/ElectricFields/MDI_EF_Analysis/MDI_EF_Analysis/MDI_EF_Analysis.py`)
 - the file `/home/username/ElectricFields/MDI_EF_Analysis/test/locations/Tinker` with the full path to the "dynamic.x" executable (`/home/username/ElectricFields/Tinker/build/tinker/source/dynamic.x`)
 
### Testing
You can now run a quick test of the driver by changing directory to the `/honme/username/ElectricFields/MDI_EF_Analysis/test/bench5` directory and running the `tcp.sh` script:
```sh
    ./tcp.sh
```
The driver's output should match the reference output file (`proj_totfield.csv`) in the `/home/username/MDI_EF_Analysis/sample_analysis` directory.

 
