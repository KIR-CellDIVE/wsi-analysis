# Whole-slide image analysis of Cell DIVE multiplex microscopy images

This work aims to facilitate and simplify the analysis of segmented and quantified whole-slide images for researchers using the Cell DIVE multiplex imaging platform. This analysis workflow closely follows the approach established by Korsunsky et al. ([paper](https://doi.org/10.1016%2Fj.medj.2022.05.002)/[github](https://github.com/immunogenomics/fibroblastatlas2022)) and aims to make it more generally applicable to new data. This analysis workflow is part of [DIVEMAP](https://github.com/KIR-CellDIVE/DIVE-MAP) and a STAR protocol publication (doi: TO BE ADDED).

## Installation

### Setup system for building the whole-slide image analysis container
Please follow the steps outlined here [here](https://github.com/KIR-CellDIVE/wsi-analysis) to prepare your computer for building and running the whole-slide image analysis container. You only have to setup you system once to build the various containers that are part of [DIVEMAP](https://github.com/KIR-CellDIVE/DIVE-MAP)

### Build whole-slide image analysis container

If you on Windows enter your previously created WSL virtual environment by typing `wsl -d Ubuntu -u ubuntu` (if you have not already done so) or if you on Linux open your favourite terminal emulator. To build the WSI segmentation container we start by creating a `builds` folder in the HOME `~` directory and cloning/downloading this repository from GitHub: 


```bash
mkdir -p ~/builds \
&& cd ~/builds \
&& git clone https://github.com/KIR-CellDIVE/wsi-analysis.git
```
Next, we build a Apptainer container called `wsi_analysis.sif` based on definition file `Apptainer`:

```bash
cd wsi-analysis/apptainer \
&& sudo apptainer build wsi_analysis.sif Apptainer
```

In order to make it easier to run the container in the future we create a bash scripts `wsi-analysis` in `~/.local/bin` that can simply be called from anywhere inside the console. Adapt these commands if you decide to download and build the container in a different directory. (Skip this step if you'd rather start the containers directly yourself). 

We make sure that `~/.local/bin` exists.
```bash
mkdir -p ~/.local/bin
```
Then, we create two bash scripts in `~/.local/bin` to make starting the container to run the analysis more straightforward.


```bash
echo "#! /bin/bash
## run wsi-analysis without GPU acceleration
[ -d "/mnt" ] && apptainer run --writable-tmpfs \"\$@\" --bind /mnt:/opt/analysis/drives --bind /:/opt/analysis/host $HOME/builds/wsi-analysis/apptainer/wsi_analysis.sif || apptainer run --writable-tmpfs \"\$@\" --bind /:/opt/analysis/host $HOME/builds/wsi-analysis/apptainer/wsi_analysis.sif" > ~/.local/bin/wsi-analysis
```
Lastly, we make these two bash scripts executable

```bash
chmod +x ~/.local/bin/wsi-analysis
```
and reload the `~/.profile` file to add `~/.local/bin` to `$PATH`.
```bash
source ~/.profile
```



## Run whole-slide image segmentation

If you have followed the installation step you should be able to run the whole-slide image segmentation Jupyter Notebook server now. If you are on `Windows` and you use `WSL`, first open `PowerShell` and enter the previously created WSL environment `Ubuntu` as the user `ubuntu` if you haven't already done so:

```bash
wsl -d Ubuntu -u ubuntu
```

Once you are in the `WSL` environment you can start the analysis container by typing
```bash
wsi-analysis
```

> You can pass additional Apptainer arguments if you want. For example, to bind a results folder to a directory `/data` to make it more easily accessible inside the notebook. In `WSL` the `C:` drive, `D:` drive, etc are mounted and located at `/mnt/c`, `/mnt/d`, etc, respectively. To mount your data folder to `/data` start the notebooks as follows:
>```bash 
> wsi-analysis --bind /path/to/result:/data
>```
>

You should now see a link similar to `http://127.0.0.1:9999/lab/workspaces/lab?reset?token=...`, copy it and open it in your preferred browser. Then, in the left sidebar navigate to the `notebooks` folder and open the `2_WSI_Analysis.ipnyb` notebook. Follow the instructions at the top of the notebook to save and open a copy of the notebook. Once done, you can start the analysis of your segmented images by using the template notebook that is part of this repository.


## macOS installation
`Apptainer` can also be installed under MacOS making use of virtualisation using `Vagrant`. However, we can not give any guarantees and support for running this container and segmentation notebook under macOS. Thus, please refer to the official [Apptainer Documentation](https://apptainer.org/docs/admin/main/installation.html#mac) for detailed installation instructions of the container environment. These installation instruction should provide you with a Linux environment, which you can use to build the whole-slide image segmentation container. However, at this moment in time this method does not support GPU-accelerated segmentation which will make it very slow for large Cell DIVE slides.


## References and Acknowledgments

The work in this repository and protocol paper was based on the work by Korsunsky et al. ([paper](https://doi.org/10.1016%2Fj.medj.2022.05.002)/[github](https://github.com/immunogenomics/fibroblastatlas2022)).


## How to cite

If you use this work as part of your analysis please cite this `wsi-analysis` repo directly (https://github.com/KIR-CellDIVE/wsi-analysis) as well as the accompanying publication: (**to be added**). Please also refer to the repositories acknowledged here and ensure compliance with all licensing requirements.

* Authors, Title, Journal, Year, DOI
