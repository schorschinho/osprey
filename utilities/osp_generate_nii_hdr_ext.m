function nii_hdr_ext = osp_generate_nii_hdr_ext(in,OspreyVersion)
%% nii_hdr_ext = osp_generate_nii_hdr_ext(in,OspreyVersion)
%   This function creates the NIfTI header extension for NIfTI-MRS
%   following the standard described in "Clarke WT, Bell TK, Emir UE,
%   Mikkelsen M, Oeltzschner G, Shamaei A, Soher BJ, Wilson M. 
%   NIfTI-MRS: A standard data format for magnetic resonance spectroscopy. 
%   Magn Reson Med. 2022. doi: 10.1002/mrm.29418.
%
%   USAGE:
%       nii_hdr_ext = osp_generate_nii_hdr_ext(in,OspreyVersion)
%
%   INPUTS:
%       in              = input data in matlab structure format.
%       OspreyVersion   = Osprey version string from MRSCont.ver.Osp.
%
%   OUTPUTS:
%       nii_hdr_ext     = NIfTI header extension.
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


    nii_hdr_ext.SpectrometerFrequency = in.txfrq/1e6;
    try
        nii_hdr_ext.ResonantNucleus = {in.nucleus};
    catch
        nii_hdr_ext.ResonantNucleus = {'1H'};
    end
    nii_hdr_ext.EchoTime            = double(in.te/1e3);
    nii_hdr_ext.RepetitionTime      = double(in.tr/1e3);
    nii_hdr_ext.TxOffset            = double(0);
    nii_hdr_ext.Manufacturer        = in.Manufacturer;
    nii_hdr_ext.SoftwareVersions    = in.software;
    nii_hdr_ext.SequenceName        = in.seq;
    nii_hdr_ext.PatientPosition     = in.PatientPosition;
    nii_hdr_ext.ConversionMethod    = OspreyVersion;
    nii_hdr_ext.ConversionTime      = datestr(now,30);
    nii_hdr_ext.kSpace              = logical(zeros(3,1));
    nii_hdr_ext.OriginalFile        = {in.OriginalFile};
    
end