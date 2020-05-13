# Installing Osprey

Get the latest **Osprey** code from our [GitHub repository](https://github.com/schorschinho/osprey):

  - Clone the repository, or download and extract the contents of the .ZIP file to a directory on your drive.

    <img src="/img/02-installing-osprey-download-clone.png" alt="Downloading Osprey" width="300"/>

  - Add the entire `osprey` directory (with subfolders) to your MATLAB path.
  - Make sure to regularly check the repository for updates, as we frequently commit new features, bug fixes, and improved functions.

To perform voxel co-registration and tissue segmentation, download **SPM12** [from the University College London website](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/). You will need to provide the SPM team with some information before you can access the download link.

  - Extract the downloaded archive so that the `spm12` folder is on the same directory level as the `osprey` folder:

      ![Installing SPM12](/img/02-installing-osprey-spm-folder.png "Installing SPM12")

  - Add the `spm12` folder to your MATLAB path, but **without** subfolders - during testing, we have found that adding SPM subfolders to the MATLAB path can cause functions to fail.

If you want to use the **Osprey** Graphical User Interface (GUI), please download the following toolboxes from the MATLAB File Exchange:

  - [GUI Layout Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox)
    (David Sampson)

  - [Widget Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/66235-widgets-toolbox)
    (Robin Jackey)

    Download both toolboxes in the MATLAB toolbox format (.mltbx). You can
    double-click to install. MATLAB will automatically add the toolboxes to its path.
