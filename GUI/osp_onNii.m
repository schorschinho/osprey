function osp_onNii( ~, ~,gui)
%% osp_onSave
%   Callback function on nii button click. Opens the external Nifti Viewer.
%
%
%   USAGE:
%       osp_onNii( ~, ~ ,gui);
%
%   INPUT:      gui      = gui class containing all handles and the MRSCont 
%
%   OUTPUT:     Changes in gui parameters and MRSCont are written into the
%               gui class
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-01-16)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2021-02-06: First version of the code.
%

    fprintf('Opening external nii viever...\n');
    MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class 
    if ~(isfield(MRSCont.flags,'addImages') && (MRSCont.flags.addImages == 1) && MRSCont.flags.moved)
        if  ~((isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI))
            if exist(MRSCont.coreg.vol_mask{gui.controls.Selected}.fname, 'file')
                nii_viewer(MRSCont.files_nii{gui.controls.Selected}, MRSCont.coreg.vol_mask{gui.controls.Selected}.fname);
            else if exist([MRSCont.coreg.vol_mask{gui.controls.Selected}.fname, '.gz'], 'file')
                nii_viewer(MRSCont.files_nii{gui.controls.Selected}, [MRSCont.coreg.vol_mask{gui.controls.Selected}.fname, '.gz']);
                end
            end
        else
            if exist(MRSCont.coreg.vol_mask{gui.controls.Selected}{1}.fname, 'file')
                 nii_viewer(MRSCont.files_nii{gui.controls.Selected}, {MRSCont.coreg.vol_mask{gui.controls.Selected}{1}.fname,MRSCont.coreg.vol_mask{gui.controls.Selected}{2}.fname})
            else if exist([MRSCont.coreg.vol_mask{gui.controls.Selected}{1}.fname, '.gz'], 'file')
                 nii_viewer(MRSCont.files_nii{gui.controls.Selected}, {[MRSCont.coreg.vol_mask{gui.controls.Selected}{1}.fname, '.gz'],[MRSCont.coreg.vol_mask{gui.controls.Selected}{2}.fname, '.gz']})
                end
            end

        end
        fprintf('... done.\n');
    else
        fprintf('The MRS container has been moved and you have no access to the nifti files anymore.\n');
    end
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class

end % osp_onNii