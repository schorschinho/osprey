# How to organize your raw data

Raw MRS data come in an overwhelming variety of formats, each producing
different numbers of files. **Osprey** does not make a lot of assumptions with regard to your folder structure, since you specify the exact location for each file in a job file prior to each analysis. It is **highly** recommended, however, that you store different acquisitions in separate folders. Organizing your raw data in a consistent and meaningful way will save you a lot of time and nerves.

To ensure optimal functioning of **Osprey**, we suggest adapting the folder
hierarchy [proposed by the BIDS (Brain Imaging Data Structure)
initiative](https://github.com/bids-standard/bids-starter-kit/wiki/The-BIDS-folder-hierarchy):

> There are four main levels of the folder hierarchy, these are:
```
project/
└── subject
    └── session
        └── acquisition
```
> With the exception of the top-level `project` folder, all sub-folders have a specific structure to their name (described below). Here's an example of how this hierarchy looks:
```
myProject/
└── sub-01
    └── ses-01
        └── anat
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
(`.dat`) all in the same folder, and **Osprey** will be able to handle them, if
you define their paths in the job file accordingly.

!!! warning
    There is one exception: For data stored in single-average DICOM (`*.IMA`, `*.DCM`) or single-average Siemens RDA format, it is **absolutely necessary** to keep every scan in a separate folder.
