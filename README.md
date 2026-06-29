# 3D-PDR
Code for modelling the chemistry of photodissociation regions in 1D and 3D

Click [here](http://itamos.readthedocs.io/) to visit the manual of the ITAMOS project and the 3D-PDR code.

## Quick installation 

### SUNDIALS

You need to have `cmake` already installed in your system.

Install SUNDIALS using the following commands in a directory where you will generally have the 3D-PDR code.

```console
$ git clone https://github.com/LLNL/sundials.git sundials
```

Create and enter a build directory:

```console
$ mkdir sundials/build
$ cd sundials/build
```

Configure SUNDIALS using `cmake`:

```console
cmake -DCMAKE_INSTALL_PREFIX=/YOUR-HOMEPATH/3D-PDR/sundials ../
```

Build and install

```console
$ make
$ make install
```

Next, edit your shell configuration file (e.g. `~/.bashrc`) and add the following lines:

```console
export LD_LIBRARY_PATH=/home/USERNAME/3D-PDR/sundials/lib
export SUNDIALS_DIR=/home/USERNAME/3D-PDR/sundials
```

Then, reload your shell configuration:

```console
$ source ~/.bashrc
```

### PDR-studio

Go to the `PDR-studio/` directory you will find in the `3D-PDR/` folder you have cloned the code. You will see the `run.sh` file:

```console
$ chmod 755 run.sh
$ ./run.sh
```

This will install all requirements needed for the PDR-studio to run. This will be done only once. After the installation finishes, you will be prompted with IP addresses. If your browser does not start automatically, CTRL+Click on one of the two IP addresses.

Now run and analyse your models with **PDR-studio**.

Note that in the Chemical Analysis section, you will need to install the Ollama LLM model. Follow the instructions you will see on that section if you wish to generate AI summaries of the model you run. 
