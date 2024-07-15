function out = osp_add_nii_mrs_field(in,OspreyVersion)
%% nii_hdr_ext = osp_add_nii_mrs_field(in,OspreyVersion)
%   This function adds nii_mes field to input data matlab structure format
%   following the standard described in "Clarke WT, Bell TK, Emir UE,
%   Mikkelsen M, Oeltzschner G, Shamaei A, Soher BJ, Wilson M. 
%   NIfTI-MRS: A standard data format for magnetic resonance spectroscopy. 
%   Magn Reson Med. 2022. doi: 10.1002/mrm.29418.
%
%   USAGE:
%       out             = osp_add_nii_mrs_field(in,OspreyVersion)
%
%   INPUTS:
%       in              = input data in matlab structure format.
%       OspreyVersion   = Osprey version string from MRSCont.ver.Osp.
%
%   OUTPUTS:
%       out             = output data in matlab structure format.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-06)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%       This code is also based on numerous functions from the spec2nii
%       conversion tool by Dr. William Clarke (University of Oxford)
%       https://github.com/wtclarke/spec2nii
%       Clarke et al., Magn Reson Med 88:2358-2370 (2022)
%
%   HISTORY:
%       2023-03-06: First version of the code.

    out = in;

    % Initialize the NIfTI struct using the dicm2nii toolbox
    % (https://github.com/xiangruili/dicm2nii)
    try
        nii = nii_tool('init', out.fids);
    catch ME
        switch ME.identifier
            case 'MATLAB:UndefinedFunction'
                error(['Cannot find the function ''nii_tool.m''.' ...
                    ' Please ensure that you have downloaded the required', ...
                    ' dcm2nii toolbox (https://github.com/xiangruili/dicm2nii)', ...
                    ' and added it to your MATLAB path.']);
            otherwise
                rethrow(ME);
        end
    end

    % Add nii_mrs field to the FID-A struct
    if ~isfield(out,'nii_mrs')
        out.nii_mrs.hdr_ext = osp_generate_nii_hdr_ext(out,OspreyVersion);
        out.nii_mrs.hdr     = osp_generate_nii_hdr(out, nii.hdr,out.OriginalFile);  
    end

    % Add Processing software to header extension
    ver = split(OspreyVersion,' ');
    ver = ver{2};
    out.nii_mrs.hdr_ext.ProcessingSoftwareVersion.Value           = ver;
    out.nii_mrs.hdr_ext.ProcessingSoftwareVersion.Description     = 'The Osprey version used to process the data stored in this file.';
end