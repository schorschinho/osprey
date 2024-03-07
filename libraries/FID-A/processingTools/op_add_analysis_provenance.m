function out = op_add_analysis_provenance(in,fields)
%%  out = op_add_analysis_provenance(in,fields)
%   This function adds analysis provenance to the the NIfTI header 
%   extension for NIfTI-MRS following the standard described in "Clarke WT,
%   Bell TK, Emir UE, Mikkelsen M, Oeltzschner G, Shamaei A, Soher BJ, Wilson M. 
%   NIfTI-MRS: A standard data format for magnetic resonance spectroscopy. 
%   Magn Reson Med. 2022. doi: 10.1002/mrm.29418.
%
%   USAGE:
%       out = op_add_analysis_provenance(in,fields)
%
%   INPUTS:
%       in              = input data in matlab structure format.
%       nii_hdr         = NIfTI header.
%
%   OUTPUTS:
%       out             = input data in matlab structure format.
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

    if isfield(in,'nii_mrs') % Has NIfTI information
        out = in;
        
        % Get current NIFTI-MRS header extension
        hdr_ext = out.nii_mrs.hdr_ext;
    
        % Add shared fields
        fields.Time     = datestr(now,30);
        fields.Program  = 'Osprey';
        fields.Version  = hdr_ext.ProcessingSoftwareVersion;
    
        fields = orderfields(fields,{'Time','Program','Version','Method', 'Details'});;
    
        % Add to Processing Applied entries
        if ~isfield(hdr_ext,'ProcessingApplied')
            hdr_ext.ProcessingApplied(1) = fields;
        else
            hdr_ext.ProcessingApplied(end+1,1) = fields;
        end

        % Add updated field back into struct
        out.nii_mrs.hdr_ext = hdr_ext;

    else
        out = in;
    end
end