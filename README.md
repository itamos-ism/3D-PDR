# 3D-PDR
Code for modelling the chemistry of photodissociation regions in 1D and 3D

Click [here](http://itamos.readthedocs.io/) to visit the manual of the ITAMOS project and the 3D-PDR code.

## Quick installation 

### SUNDIALS

You will need to have `cmake` already installed in your system.

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

Navigate to the `PDR-studio/` directory inside your cloned `3D-PDR/` repository and run the following:

```bash
$ chmod 755 run.sh
$ ./run.sh
```

This will install all required dependencies — this step only runs once. Once the installation completes, you will be prompted with IP addresses. If your browser does not open automatically, Ctrl+Click on one of the displayed addresses.

You are now ready to run and analyse your models with PDR-studio.

> **Note:** The Chemical Analysis section requires an Ollama LLM model to generate AI summaries of your results. Instructions for installing and configuring Ollama are provided within that section.
````

