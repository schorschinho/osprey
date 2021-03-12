function MRSCont = osp_senseRecon(MRSCont)
%% MRS_struct = osp_senseRecon(MRS_struct)
%   Loads complex coil sensitivity reference scans, converts them into
%   NIfTI format, co-registers with MRS voxels, and determines complex coil
%   sensitivities following the SENSE formalism.
%
%   This is usually only used to reconstruct multi-voxel excitation data,
%   e.g. PRIAM.
%
%   Input:
%       MRSCont     = Osprey MRS data container.
%
%   Output:
%       MRSCont     = Osprey MRS data container filled with information 
%                   about the SENSE reconstruction (in MRSCont.SENSE)
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-03-15)
%       goeltzs1@jhmi.edu
%
%   Credits:
%       This code is based on an initial PRIAM reconstruction routine.
%       Dr. Vincent O. Boer (vincentob@drcmr.dk)
%       Danish Research Centre for Magnetic Resonance (Hvidovre Hospital)
%
%   History:
%       2018-03-15: First version of the code.
%

% Generate sub-folder for the reconstruction in the output directory
outputFolder    = MRSCont.outputFolder;
if ~exist([outputFolder filesep 'SenseReconstruction'],'dir')
    mkdir([outputFolder filesep 'SenseReconstruction']);
end

% If existing coil sensitivity images are saved in the SenseReconstruction
% folder, load them here. Otherwise, start loading them.
for kk = 1:MRSCont.nDatasets    
    % Check whether the provided CPX file exists
    cpxFile = MRSCont.files_sense{kk};
    [cpxPath, cpxName, cpxExt] = fileparts(cpxFile);
    if exist(cpxFile, 'file')
        % If it exists, check whether it has previously been loaded and
        % converted to NIfTI
        sinFile = [cpxFile(1:end-4) '.sin'];
        niiSOSFile = [outputFolder filesep 'SenseReconstruction' filesep cpxName '_SOS.nii'];
        niiSENSEFile = [outputFolder filesep 'SenseReconstruction' filesep cpxName '_SENSE.nii'];
%         if exist(niiSOSFile, 'file') && exist(niiSENSEFile, 'file')
%             niiSOS      = nii_tool('load', niiSOSFile);
%             niiSENSE    = nii_tool('load', niiSENSEFile);
%         else
            [niiSOS, niiSENSE, noise_array] = osp_cpx2nii(cpxFile, sinFile, niiSOSFile, niiSENSEFile, 0);
            MRSCont.SENSE{kk}.noise_array = noise_array;
%         end
    else
        error('Philips complex coil sensitivity file %s not found!\n', cpxFile);
    end
    

    % Calculate the SENSE unfolding matrix based on MRS voxel locations
    % (This is where the bulk of the work happens)
    % Set the last argument for calcUnfoldingMatrix to 1 to view the
    % overlays (for debugging purposes)
    MRSCont = calcUnfoldingMatrix(MRSCont, niiSOSFile, niiSENSEFile, kk, 0);
end


% Apply the SENSE unfolding matrix
for kk = 1:MRSCont.nDatasets
    MRSCont.raw{kk} = osp_SenseUnfolding(MRSCont.raw_uncomb{kk},MRSCont.SENSE{kk});
    % Now do the same for the water reference signal
    if MRSCont.flags.hasRef
        MRSCont.raw_ref{kk} = osp_SenseUnfolding(MRSCont.raw_ref_uncomb{kk},MRSCont.SENSE{kk});
    end
    % Now do the same for the (short-TE) water signal
    if MRSCont.flags.hasWater
        MRSCont.raw_w{kk} = osp_SenseUnfolding(MRSCont.raw_w_uncomb{kk},MRSCont.SENSE{kk});
    end
end

% Set flags
MRSCont.flags.coilsCombined     = 1;

% Delete coil-un-combined data to free up memory
raw_fields = {'raw_uncomb','raw_ref_uncomb','raw_w_uncomb'};
for kk = 1:length(raw_fields)
    if isfield(MRSCont, raw_fields{kk})
        MRSCont = rmfield(MRSCont, raw_fields{kk});
    end
end

end
