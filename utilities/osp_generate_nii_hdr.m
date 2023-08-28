function nii_hdr = osp_generate_nii_hdr(in,nii_hdr,outfile)
%%  nii_hdr = osp_generate_nii_hdr(in,nii_hdr,outfile)
%   This function creates the NIfTI header extension for NIfTI-MRS
%   following the standard described in "Clarke WT, Bell TK, Emir UE,
%   Mikkelsen M, Oeltzschner G, Shamaei A, Soher BJ, Wilson M. 
%   NIfTI-MRS: A standard data format for magnetic resonance spectroscopy. 
%   Magn Reson Med. 2022. doi: 10.1002/mrm.29418.
%
%   USAGE:
%       nii_hdr = osp_generate_nii_hdr(in,nii_hdr,outfile)
%
%   INPUTS:
%       in              = input data in matlab structure format.
%       nii_hdr         = NIfTI header.
%       outfile         = outfile name (string)
%
%   OUTPUTS:
%       nii_hdr         = NIfTI header.
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
    
    % Set easy pixdim vector entries
    nii_hdr.pixdim(5) = in.dwelltime;
    nii_hdr.pixdim(6) = 1;
    nii_hdr.pixdim(7) = 1;
    nii_hdr.pixdim(8) = 1;
    
    % Get affine matrix from vendor native data
    [m44] = get_voxel_pos_for_nifti(in);

    % Fill s-form
    nii_hdr.srow_x = m44(1,:);
    nii_hdr.srow_y = m44(2,:);
    nii_hdr.srow_z = m44(3,:);

    % Calculate q-form
    [nii_hdr.quatern_b, nii_hdr.quatern_c, nii_hdr.quatern_d,...
    nii_hdr.qoffset_x, nii_hdr.qoffset_y, nii_hdr.qoffset_z,...
    nii_hdr.pixdim(2), nii_hdr.pixdim(3), nii_hdr.pixdim(4),...
    nii_hdr.pixdim(1)] = nifti_mat44_to_quatern(m44);

    % Remove unwanted fields from original nii_hdr
    nii_hdr = rmfield(nii_hdr,'db_name');
    nii_hdr = rmfield(nii_hdr,'extents');
    nii_hdr = rmfield(nii_hdr,'glmax');
    nii_hdr = rmfield(nii_hdr,'glmin');
    nii_hdr = rmfield(nii_hdr,'regular');
    nii_hdr = rmfield(nii_hdr,'session_error');
    nii_hdr = rmfield(nii_hdr,'data_type');

    %Other fields (have to check this)
    nii_hdr.sizeof_hdr = 540; % NIfTI-2
    nii_hdr.vox_offset = 992;
    nii_hdr.magic = 'n+2 \n';
    nii_hdr.datatype = 32;  % complex data
    nii_hdr.bitpix  = 64;
    nii_hdr.qform_code = 2;
    nii_hdr.sform_code = 2;
    nii_hdr.intent_name = 'mrs_v0_7';
    nii_hdr.unused_str = '';
    nii_hdr.extension = [1,0,0,0];
    nii_hdr.version = 2;
    nii_hdr.swap_endian = false;
    nii_hdr.file_name = outfile;
end

function m44 = get_voxel_pos_for_nifti(in)
%% m44 = get_voxel_pos_for_nifti(in)
% This function generates the affine transformation matrix for the
% NIfTI-MRS data for GE, Philips, and Siemens.
%   USAGE:
%       m44 = get_voxel_pos_for_nifti(in)
%
%   INPUTS:
%       in              = input data in matlab structure format.
%
%   OUTPUTS:
%       m44             = affine matrix.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-06)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is also based on numerous functions from the spec2nii
%       conversion tool by Dr. William Clarke (University of Oxford)
%       https://github.com/wtclarke/spec2nii
%       Clarke et al., Magn Reson Med 88:2358-2370 (2022)

    m44 =[];
    switch in.Manufacturer
        case 'Philips'
            % Get voxel information
            dimensions(1) = in.geometry.size.lr;
            dimensions(2) = in.geometry.size.ap;
            dimensions(3) = in.geometry.size.cc;
            euler_angles(1) = -in.geometry.rot.lr;
            euler_angles(2) = -in.geometry.rot.ap;
            euler_angles(3) = in.geometry.rot.cc;
            shift(1) = -in.geometry.pos.lr;
            shift(2) = -in.geometry.pos.ap;
            shift(3) = in.geometry.pos.cc;
             % Calcualte m44 matrix
            m44 = calcualte_affine_matrix(euler_angles, dimensions, shift);
        case 'Siemens'
%             res = [NaN vars.x_dim / vars.x_pts vars.y_dim / vars.y_pts
%             vars.z_dim / vars.z_pts 1 NaN 1/vars.fs];  %Keep for MRSI?
            % Get voxel information
            ima_norm = [in.geometry.rot.NormSag in.geometry.rot.NormCor in.geometry.rot.NormTra];
            ima_norm = ima_norm / norm(ima_norm);
            if isfield(in.geometry.pos,'TablePosTra')
                ima_pos = [in.geometry.pos.PosSag + in.geometry.pos.TablePosSag ...
                           in.geometry.pos.PosCor + in.geometry.pos.TablePosCor ...
                           in.geometry.pos.PosTra + in.geometry.pos.TablePosTra];
            else
                ima_pos = [in.geometry.pos.PosSag in.geometry.pos.PosCor in.geometry.pos.PosTra];
            end            
            rotation = in.geometry.rot.VoI_InPlaneRot;
            
            % Get voxel orientation from notrmal vectors
            [~, maxdir] = max([abs(in.geometry.rot.NormSag) ...
                               abs(in.geometry.rot.NormCor) ...
                               abs(in.geometry.rot.NormTra)]);

            [dGp, dGr] = fGSLCalcPRS(ima_norm, rotation,maxdir);

            col_vec = dGp;

            % Keep for MRSI case
%             if x_pts * y_pts * z_pts == 1
                row_vec = dGr;
%             else
%                 row_vec = -dGr;
%             end
            
            sli_vec = cross(row_vec, col_vec);
            ima_pos = ima_pos - row_vec * (1 / 2 - 0.5) * in.geometry.size.VoI_RoFOV / 1 - col_vec * (1 / 2 - 0.5) * in.geometry.size.VoI_PeFOV / 1;

            % Keep for MRSI case
%             ima_pos = ima_pos - row_vec * (vars.x_pts / 2 - 0.5) * vars.x_dim / vars.x_pts - col_vec * (vars.y_pts / 2 - 0.5) * vars.y_dim / vars.y_pts;
            
            m44 = [row_vec * in.geometry.size.VoI_RoFOV, 0; col_vec * in.geometry.size.VoI_PeFOV, 0; sli_vec * in.geometry.size.VoIThickness, 0; ima_pos, 1]';
            m44(1:2,:) = -m44(1:2,:);
        case 'GE'

            % Calculate the 4x4 affine matrix
            geom = in.geometry;
            
            % Convert from RAS to LPS
            VoxOffs  = geom.pos .* [-1 -1 1];
            tlhc_LPS = geom.rot.tlhc .* [-1 -1 1];
            trhc_LPS = geom.rot.trhc .* [-1 -1 1];
            brhc_LPS = geom.rot.brhc .* [-1 -1 1];
            
            e1_SVS_n = trhc_LPS - tlhc_LPS;
            e1_SVS_n = e1_SVS_n ./ norm(e1_SVS_n);
            e2_SVS_n = brhc_LPS - trhc_LPS;
            e2_SVS_n = e2_SVS_n ./ norm(e2_SVS_n);
            e3_SVS_n = -cross(e1_SVS_n, e2_SVS_n);
            
            [~,orientation_SVS] = max(abs(e3_SVS_n));
            
            if orientation_SVS == 3     % axial
                e1_SVS_n2 = e1_SVS_n;
                e2_SVS_n2 = e2_SVS_n;
                e3_SVS_n2 = e3_SVS_n;
            elseif orientation_SVS == 2 % coronal
                e1_SVS_n2 = e1_SVS_n;
                e2_SVS_n2 = e3_SVS_n;
                e3_SVS_n2 = e2_SVS_n;
            elseif orientation_SVS == 1 % sagittal
                e1_SVS_n2 = e3_SVS_n;
                e2_SVS_n2 = e1_SVS_n;
                e3_SVS_n2 = e2_SVS_n;
            end
            e1_SVS = geom.size.dim1 * e1_SVS_n2;
            e2_SVS = geom.size.dim2 * e2_SVS_n2;
            e3_SVS = geom.size.dim3 * e3_SVS_n2;

            m44 = zeros(4, 4);
            m44(1,1:3) = e1_SVS;
            m44(2,1:3) = e2_SVS;
            m44(3,1:3) = e3_SVS;
            m44(1:3, 4) = VoxOffs;
            m44(4, 4) = 1.0; 
            m44(1:2,:) = -m44(1:2,:);
        case 'LCModel RAW'
            m44 = zeros(4, 4);
            m44(1,1)=10000;
            m44(2,2)=10000;
            m44(3,3)=10000;
            m44(4,4)=1;

    end
end


function [dGp, dGr] = fGSLCalcPRS(dGs, dPhi,lOrientation)
%% [dGp, dGr] = fGSLCalcPRS(dGs, dPhi,lOrientation)
% This function generates the affine transformation matrix for the
% NIfTI-MRS data for GE, Philips, and Siemens.
%   USAGE:
%       [dGp, dGr] = fGSLCalcPRS(dGs, dPhi,lOrientation)
%
%   INPUTS:
%       dGs      = dGS vector (= slice normal vector)
%       dPhi     = the rotational angle around dGs
%
%   OUTPUTS:
%       dGp      = dGp vector (= phase direction vector).
%       dGr      = dGr vector (= read direction vector).
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-06)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is also based on numerous functions from the spec2nii
%       conversion tool by Dr. William Clarke (University of Oxford)
%       This is a MATLAB conversion fo the gslfunction.py to generate the
%       affine matrix for Siemens data
%       https://github.com/wexeee/spec2nii/blob/master/spec2nii/GSL/gslfunctions.py
%       https://github.com/wtclarke/spec2nii
%       Clarke et al., Magn Reson Med 88:2358-2370 (2022)

    dGp = zeros(1, 3);
    
    if lOrientation == 3        %TRANSVERSE
        dGp(1) = 0.0;
        dGp(2) = dGs(3) * sqrt(1. / (dGs(2)^2 + dGs(3)^2));
        dGp(3) = -dGs(2) * sqrt(1. / (dGs(2)^2 + dGs(3)^2));
    elseif lOrientation == 2    %CORONAL
        dGp(1) = dGs(2) * sqrt(1. / (dGs(1)^2 + dGs(2)^2));
        dGp(2) = -dGs(1) * sqrt(1. / (dGs(1)^2 + dGs(2)^2));
        dGp(3) = 0.0;
    elseif lOrientation == 1    %SAGITTAL
        dGp(1) = -dGs(2) * sqrt(1. / (dGs(1)^2 + dGs(2)^2));
        dGp(2) = dGs(1) * sqrt(1. / (dGs(1)^2 + dGs(2)^2));
        dGp(3) = 0.0;
    else
        error('Invalid slice orientation returned from fGSLClassOri');
    end
    
    % Calculate GR = GS x GP
    dGr = zeros(1, 3);
    dGr(1) = dGs(2) * dGp(3) - dGs(3) * dGp(2);
    dGr(2) = dGs(3) * dGp(1) - dGs(1) * dGp(3);
    dGr(3) = dGs(1) * dGp(2) - dGs(2) * dGp(1);
    
    
    dGp(1) = cos(dPhi) * dGp(1) - sin(dPhi) * dGr(1);
    dGp(2) = cos(dPhi) * dGp(2) - sin(dPhi) * dGr(2);
    dGp(3) = cos(dPhi) * dGp(3) - sin(dPhi) * dGr(3);
    
    % Calculate new GR = GS x GP
    dGr(1) = dGs(2) * dGp(3) - dGs(3) * dGp(2);
    dGr(2) = dGs(3) * dGp(1) - dGs(1) * dGp(3);
    dGr(3) = dGs(1) * dGp(2) - dGs(2) * dGp(1);
end

function R = eul2rotm_xyz(euler_angles)
% Converts Euler angles to a 3x3 rotation matrix using 'XYZ' convention
% euler_angles: 1x3 vector of [roll, pitch, yaw] angles in degrees

    phi = deg2rad(euler_angles(1));
    theta = deg2rad(euler_angles(2));
    psi = deg2rad(euler_angles(3));
    
    R_x = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    R_y = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    R_z = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    
    R = R_z * R_y * R_x;
end

function [roll, pitch, yaw] = rotm2eul_xyz(R)
% rotm2eul_xyz converts a rotation matrix to Euler angles using XYZ order
% The input must be a 3x3 rotation matrix.
% The output angles are in degree.

    % Extract the first two columns of the rotation matrix.
    x_axis = R(:, 1);
    y_axis = R(:, 2);
    
    pitch = rad2deg(atan2(-R(2,3),R(3,3)));
    
    roll = rad2deg(asin(R(1,3)));
    
    yaw = rad2deg(atan2(-R(1,2),R(1,1)));
end

function m44 = calcualte_affine_matrix(euler_angles, dimensions, shift)
% calcualte_affine_matrix converts euler angles, voxel dimesions, and voxel
% shifts to a rotation matrix
% The input angles are in degree.

    scalingMat = diag(dimensions);
    rot = eul2rotm_xyz(euler_angles);
    m33 = rot * scalingMat;
    m44 = eye(4);
    m44(1:3, 1:3) = m33;
    m44(1:3, 4) = shift;
end