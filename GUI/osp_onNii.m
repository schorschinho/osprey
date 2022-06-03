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
            file_code = [num2str(exist(MRSCont.coreg.vol_image{gui.controls.Selected}.fname,'file'))...
                        num2str(exist([MRSCont.coreg.vol_image{gui.controls.Selected}.fname, '.gz'],'file'))...
                        num2str(exist(MRSCont.coreg.vol_mask{gui.controls.Selected}.fname, 'file')) ...
                        num2str(exist([MRSCont.coreg.vol_mask{gui.controls.Selected}.fname, '.gz'], 'file'))];                    
            switch file_code
                case '2020'
                    nii_viewer(MRSCont.coreg.vol_image{gui.controls.Selected}.fname, MRSCont.coreg.vol_mask{gui.controls.Selected}.fname);
                case '2002'
                    nii_viewer(MRSCont.coreg.vol_image{gui.controls.Selected}.fname, [MRSCont.coreg.vol_mask{gui.controls.Selected}.fname, '.gz']);
                case '0220'
                    nii_viewer([MRSCont.coreg.vol_image{gui.controls.Selected}.fname, '.gz'], MRSCont.coreg.vol_mask{gui.controls.Selected}.fname);
                case '0202'
                    nii_viewer([MRSCont.coreg.vol_image{gui.controls.Selected}.fname, '.gz'], [MRSCont.coreg.vol_mask{gui.controls.Selected}.fname, '.gz']);
                case '2000'
                    nii_viewer(MRSCont.coreg.vol_image{gui.controls.Selected}.fname);
                case '0200'
                    nii_viewer([MRSCont.coreg.vol_image{gui.controls.Selected}.fname, '.gz']);
                otherwise
                    fprintf('Nothing to open here.\n');
            end
        else
            file_code = [num2str(exist(MRSCont.coreg.vol_image{gui.controls.Selected}{1}.fname,'file'))...
                        num2str(exist([MRSCont.coreg.vol_image{gui.controls.Selected}{1}.fname, '.gz'],'file'))...
                        num2str(exist(MRSCont.coreg.vol_mask{gui.controls.Selected}{1}.fname, 'file')) ...
                        num2str(exist([MRSCont.coreg.vol_mask{gui.controls.Selected}{1}.fname, '.gz'], 'file'))];                    
            switch file_code
                case '2020'
                    nii_viewer(MRSCont.coreg.vol_image{gui.controls.Selected}{1}.fname, MRSCont.coreg.vol_mask{gui.controls.Selected}{1}.fname, MRSCont.coreg.vol_mask{gui.controls.Selected}{2}.fname);
                case '2002'
                    nii_viewer(MRSCont.coreg.vol_image{gui.controls.Selected}{1}.fname, [MRSCont.coreg.vol_mask{gui.controls.Selected}{1}.fname, '.gz'], [MRSCont.coreg.vol_mask{gui.controls.Selected}{2}.fname, '.gz']);
                case '0220'
                    nii_viewer([MRSCont.coreg.vol_image{gui.controls.Selected}{1}.fname, '.gz'], MRSCont.coreg.vol_mask{gui.controls.Selected}{1}.fname, MRSCont.coreg.vol_mask{gui.controls.Selected}{2}.fname);
                case '0202'
                    nii_viewer([MRSCont.coreg.vol_image{gui.controls.Selected}{1}.fname, '.gz'], [MRSCont.coreg.vol_mask{gui.controls.Selected}{1}.fname, '.gz'], [MRSCont.coreg.vol_mask{gui.controls.Selected}{2}.fname, '.gz']);
                case '2000'
                    nii_viewer(MRSCont.coreg.vol_image{gui.controls.Selected}{1}.fname);
                case '0200'
                    nii_viewer([MRSCont.coreg.vol_image{gui.controls.Selected}{1}.fname, '.gz']);
                otherwise
                    fprintf('Nothing to open here.\n');
            end
        end
        fprintf('... done.\n');
    else
        fprintf('The MRS container has been moved and you have no access to the nifti files anymore.\n');
    end
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class

end % osp_onNii