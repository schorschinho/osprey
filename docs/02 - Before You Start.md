# Before You Start

This chapter intends to provide an overview of `LCGannet` requirements,
installation guidelines, and helpful suggestions for organizing MRS data in a
meaningful folder structure. The efficient use of `LCGannet` requires some very
basic familiarity with MATLAB. If you have never used MATLAB before, someone
else in your group probably has, and can assist with setting everything up for
you to get started.

## Installation

You will need the following software and code:

- **MATLAB** (version 2017a or newer) with the *Optimization Toolbox* and the
  *Statistics Toolbox*. You may need to ask your institutional system
  administrator for a separate license for the toolboxes.

- Download the latest **LCGannet** code from our [Github
  repository](https://github.com/schorschinho/LCGannet), then extract and add
  the entire folder (with subfolders) to your MATLAB path. Make sure to
  regularly check for updates, as we frequently commit new features, bug fixes,
  and improved functions.

- *(Optional)* To perform voxel co-registration and tissue segmentation,
  download **SPM12** [from their
  website](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/), then extract and
  add to your MATLAB path.

- *(COMING SOON)* If you want to use the `LCGannet` Graphical User Interface (GUI),
  please download the following toolboxes from the MATLAB File Exchange:

    - [GUI Layout
      Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox)
      (David Sampson)

    - [Widget
      Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/66235-widgets-toolbox)
      (Robin Jackey)

  Download both toolboxes in the MATLAB toolbox format (.mltbx). You can
  double-click to install. MATLAB will automatically add the toolboxes to its
  path.

## Structure and workflow

`LCGannet` was designed to have easy workflow with as little user input as
possible. It has many built-in routines to recognize data formats and sequence
types, and will be able to perform most required processing steps automatically.
However, some input will need to be provided by the user in the form of a *job*.

A job specifies the locations of the files containing MRS and structural imaging
data, the type of MRS sequence, and some basic control options over the data
modeling procedure.

While a job is progressing through the `LCGannet` pipeline, all raw and processed
data associated with this job are stored in a MATLAB structure, the LCGannet
*data container*. By default, this container is called `MRSCont`, but you are free to give it a more meaningful variable name.

`LCGannet`'s data handling is based in large parts on the free MATLAB toolbox FID-A. The `LCGannet` folder contains a library of, sometimes modified, FID-A functions and additions. Please make sure that you do not include an installation of the original FID-A suite in your MATLAB path.

The `MRSCont` container functions as a super-structure, containing
FID-A structures for each dataset and processing step, along with additional
information (e.g. a few basic QA metrics, quantification results, etc).

## How to organize your raw data

Raw MRS data come in an overwhelming variety of formats, each producing
different numbers of files. `LCGannet` does not make a lot of assumptions with regard to your folder structure, since you specify the exact location for each file in a job file prior to each analysis. It is **highly** recommended, however, that you store different acquisitions in separate folders. Organizing your raw data in a consistent and meaningful way will save you a lot of time and nerves.

To ensure optimal functioning of `LCGannet`, we suggest adapting the folder
hierarchy [proposed by the BIDS (Brain Imaging Data Structure)
initiative](https://github.com/bids-standard/bids-starter-kit/wiki/The-BIDS-folder-hierarchy):

> There are four main levels of the folder hierarchy, these are:
```
> project/
> └── subject
>     └── session
>         └── acquisition
```
> With the exception of the top-level `project` folder, all sub-folders have a specific structure to their name (described below). Here's an example of how this hierarchy looks:
```
> myProject/
> └── sub-01
>     └── ses-01
>         └── anat
```

An example adaptation of this hierarchy for Philips MRS data could look like this:

```
myProject/
├── sub-01
│   └── ses-01
│       ├── anat
│       │   └── sub-01_T1w.nii
│       └── mrs
│           ├── sub-01_mega-press_act
│           │   ├── sub-01_mega-press_act.sdat
│           │   └── sub-01_mega-press_act.spar
│           ├── sub-01_mega-press_ref
│           │   ├── sub-01_mega-press_ref.sdat
│           │   └── sub-01_mega-press_ref.spar
│           └── sub-01_press-water
│               ├── sub-01_press-water_act.sdat
│               └── sub-01_press-water_act.spar
├── sub-02
│   └── ses-01
│       ├── anat
│       │   └── sub-02_T1w.nii
│       └── mrs
│           ├── sub-02_mega-press_act
│           │   ├── sub-02_mega-press_act.sdat
│           │   └── sub-02_mega-press_act.spar
│           ├── sub-02_mega-press_ref
│           │   ├── sub-02_mega-press_ref.sdat
│           │   └── sub-02_mega-press_ref.spar
│           └── sub-02_press-water
│               ├── sub-02_press-water_act.sdat
│               └── sub-02_press-water_act.spar
...
...
```

While adhering to the BIDS standard is a recommendation, it is by no means
binding for most file formats. You can, for example, keep Siemens TWIX files
(`.dat`) all in the same folder, and `LCGannet` will be able to handle them, if
you define their paths in the job file accordingly.

**There is one exception: For data stored in single-average DICOM (`*.IMA`, `*.DCM`) or single-average Siemens RDA format, it is absolutely necessary to keep every scan in a separate folder.**
