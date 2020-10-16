function MRSCont = calcUnfoldingMatrix(MRSCont, ref_file, Ref, img_array, noise_array,spec_file,kk)
%% function MRSCont = calcUnfoldingMatrix(MRSCont, ref_file, Ref, img_array, noise_array,kk)
%   Calculates the unfolding matrix for a given coil reference image and
%   voxel locations specified in MRS_struct.
%
%   Input:
%       MRSCont     MRSCont filled with information from GannetLoad and
%                   senseRecon.
%       ref_file    cpx file used
%       img_array   4D stack of receiver coil images
%                   [number of coils, x, y, z]
%       noise_array 2D stack of receiver coil noise levels
%                   [number of coils, number of pixels estimated for noise]
%       spec_path   Folder containing the MRS data. calcUnfoldingMatrix assumes a
%                   subdirectory '/ref' containing the *.cpx data.
%
%   Output:
%       MRS_struct  MRS_struct filled with information about the SENSE
%                   reconstruction (in MRS_struct.p.SENSE)
%
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-03-18)
%       goeltzs1@jhmi.edu
%
%   Credits:
%       This code is based on an initial PRIAM reconstruction routine.
%       Dr. Vincent O. Boer (vincentob@drcmr.dk)
%       Danish Research Centre for Magnetic Resonance (Hvidovre Hospital)
%
%   History:
%       2018-03-18: First version of the code.
%


%% Extract voxel size and orientation data from *.sin files.

% Get dimensions, orientation and global offset of reference image scan
MRSCont.SENSE{kk}.refsc_sin_info = get_sin_info([ref_file(1:end-4) '.sin']);
refsc_voxel_size = MRSCont.SENSE{kk}.refsc_sin_info.voxel_size; % row - column - slice
refsc_offset = MRSCont.SENSE{kk}.refsc_sin_info.voxel_offsets; % ap - rl - fh
refsc_orientation = MRSCont.SENSE{kk}.refsc_sin_info.slice_orientation; % 0 = axial slices, 1 = sagittal slices, 2 = coronal slices

% Get dimensions, orientation and global offset of MRS scan
sin_file = strrep(spec_file,'.data','.sin');
if ~isfile(sin_file)
    msg = 'No sin file found. It needs to have the same name as the .data file.';
    error(msg);
end

MRSCont.SENSE{kk}.sin_info = get_sin_info(sin_file);


mrs_voxelsize = MRSCont.SENSE{kk}.sin_info.voxel_size; % row - column - slice
mrs_offset = MRSCont.SENSE{kk}.sin_info.voxel_offsets; % ap - rl - fh
mrs_orientation = MRSCont.SENSE{kk}.sin_info.slice_orientation; % not needed yet, 0 = axial slices, 1 = sagittal slices, 2 = coronal slices
mrs_angulation = [MRSCont.raw_uncomb{kk}.geometry.rot.ap MRSCont.raw_uncomb{kk}.geometry.rot.lr MRSCont.raw_uncomb{kk}.geometry.rot.cc]; % ap - rl - fh

% Calculate voxel offsets (this will need more work later)
mrs_offset1 = [ mrs_offset(1); mrs_offset(2); mrs_offset(3)];

% Calculate voxel offset in dependence of the orientation of the mrs voxel
switch refsc_orientation
    case 1 % saggital orientation ap - rl - fh
        mrs_offset2 = [ mrs_offset(1); ...
            mrs_offset(2) - MRSCont.SENSE{kk}.vox_sep * cos(mrs_angulation(1)*pi/180); ...
            mrs_offset(3) + MRSCont.SENSE{kk}.vox_sep * sin(mrs_angulation(1)*pi/180)]; % will depend on orientation and angulations extracted from SDAT

    case 2 % coronal orientation ap - rl - fh
        mrs_offset2 = [ mrs_offset(1)+ MRSCont.SENSE{kk}.vox_sep * sin(mrs_angulation(1)*pi/180); ...
            mrs_offset(2); ...
            mrs_offset(3) - MRSCont.SENSE{kk}.vox_sep * cos(mrs_angulation(1)*pi/180)]; % will depend on orientation and angulations extracted from SDAT
end
% Translate ref scan into 'real world' dimensions depending on orientation of
% the reference orientation. Arrange such that the order of dimensions is
% always AP-RL-FH.
switch refsc_orientation
    case 1 % sagittal slices
        % This arrives in the order
        % [ncoils fh ap rl]
        img_array = permute(img_array,[1 3 4 2]);
        dim_ap = size(img_array,2);
        dim_rl = size(img_array,3);
        dim_fh = size(img_array,4);

        img_array = flip(img_array, 3);
        % refscan voxel size arrives in order ap fh rl
        % swap to make ap rl fh
        refsc_voxel_size([1 2 3]) = refsc_voxel_size([1 3 2]);
%         mrs_offset([1 2 3]) = mrs_offset([1 3 2]);

    case 2 % coronal slices
        % This arrives in the order
        % [ncoils fh rl ap]
        img_array = permute(img_array,[1 4 3 2]);
        dim_ap = size(img_array,2);
        dim_rl = size(img_array,3);
        dim_fh = size(img_array,4);

        % refscan voxel size arrives in order fh rl ap % presumably!
        % swap to make ap rl fh
        refsc_voxel_size([1 2 3]) = refsc_voxel_size([3 2 1]);
%         mrs_offset([1 2 3]) = mrs_offset([3 2 1]);
%         refsc_offset([1 2 3]) = refsc_offset([3 2 1]);

    case 0 % axial slices
        dim_ap = size(img_array,2);
        dim_rl = size(img_array,3); % presumably!
        dim_fh = size(img_array,4); % presumably!

        % refscan voxel size arrives in order ap rl fh (presumably!)
        % swap to make ap rl fh
        refsc_voxel_size([1 2 3]) = refsc_voxel_size([1 2 3]);

end

% Correct by offset of reference scan
mrs_offset1 = mrs_offset1 - refsc_offset';
mrs_offset2 = mrs_offset2 - refsc_offset';

% This calculation below will only work for unrotated voxels
% Work on this!
lb1 = (mrs_offset1 - mrs_voxelsize'/2); % calculate coordinates of edges
ub1 = (mrs_offset1 + mrs_voxelsize'/2); % calculate coordinates of edges
lb1_vox = lb1 ./refsc_voxel_size'; % and transform into voxel numbers in the ref scan
ub1_vox = ub1 ./refsc_voxel_size'; % and transform into voxel numbers in the ref scan

ap_range = round([dim_ap/2 + 0.5 + lb1_vox(1):dim_ap/2 + 0.5 + ub1_vox(1)]);
rl_range = round([dim_rl/2 + 0.5 + lb1_vox(2):dim_rl/2 + 0.5 + ub1_vox(2)]);
fh_range = round([dim_fh/2 + 0.5 + lb1_vox(3):dim_fh/2 + 0.5 + ub1_vox(3)]);

% second voxel: fill in correct values!
lb2 = (mrs_offset2 - mrs_voxelsize'/2); % calculate coordinates of edges
ub2 = (mrs_offset2 + mrs_voxelsize'/2); % calculate coordinates of edges
lb2_vox = lb2 ./refsc_voxel_size'; % and transform into voxel numbers in the ref scan
ub2_vox = ub2 ./refsc_voxel_size'; % and transform into voxel numbers in the ref scan
ap_range2 = round(dim_ap/2 + 0.5 + lb2_vox(1):dim_ap/2 + 0.5 + ub2_vox(1));
rl_range2 = round(dim_rl/2 + 0.5 + lb2_vox(2):dim_rl/2 + 0.5 + ub2_vox(2));
fh_range2 = round(dim_fh/2 + 0.5 + lb2_vox(3):dim_fh/2 + 0.5 + ub2_vox(3));

%% Plot an overlay of the voxel on the body coil image
% Calculate the sensitivity maps with the sum of squares method
SOS = squeeze(sqrt(sum(img_array.*conj(img_array),1)));
%% Figure
% Define color map
gray2 = gray;
gray2(64,:) = [1 0 0];
gray2(1,:) = [0 1 0];

% Make new figure with subplots
figure(2)
clf
SOSdisplay = SOS;
SOSdisplay(ap_range,rl_range,fh_range) = 3; % set color of voxel 1 to red

% Plot coronal slice overlaid with voxel 1
subplot(2,3,1)
imagesc(rot90(squeeze(squeeze(SOSdisplay(round(mean(ap_range)),:,:)))),[-0.2 1.6])
title('vox1         COR          A/L')
% axis off;
axis square;

% Plot sagittal slice overlaid with voxel 1
subplot(2,3,2)
imagesc(rot90(squeeze(SOSdisplay(:,round(mean(rl_range)),:))),[-0.2 1.6])
title('vox1         SAG          H/L')
% axis off;
axis square;

% Plot axial slice overlaid with voxel 1
subplot(2,3,3)
imagesc(squeeze(SOSdisplay(:,:,round(mean(fh_range)))),[-0.2 1.6])
title('vox1         TRA          H/P')
% axis off;
axis square;

% Plot coronal slice overlaid with voxel 2
SOSdisplay = SOS;
SOSdisplay(ap_range2,rl_range2,fh_range2) =-3; % set color of voxel 2 to green
subplot(2,3,4)
imagesc(rot90(squeeze(SOSdisplay(round(mean(ap_range2)),:,:))),[-0.2 1.6])
title('vox2         COR          A/L')
% axis off;
axis square;

% Plot sagittal slice overlaid with voxel 2
subplot(2,3,5)
imagesc(rot90(squeeze(SOSdisplay(:,round(mean(rl_range2)),:))),[-0.2 1.6])
title('vox2         SAG          H/L')
% axis off;
axis square;

% Plot axial slice overlaid with voxel 2
subplot(2,3,6)
imagesc(squeeze(SOSdisplay(:,:,round(mean(fh_range2)))),[-0.2 1.6])
title('vox2         TRA          H/P')
% axis off;
axis square;

colormap(gray2)
%%
% Make mask with threshold
mask = (SOS/max(SOS(:)))>0.04;
% Save mask and sum-of-squares image
%save([spec_path filesep 'GannetRecon_output' filesep 'mask.mat'],'mask');

save(strrep( MRSCont.SENSE{kk}.senseCoilFile,'_ref_scan_sense_coil_img.mat', '_SOS_sense.mat'),'SOS');

%% generate sensitivity maps
% Now create the individual complex coil sensitivity maps by normalizing
% every individual coil image to the sum-of-squares image
sensitivities = zeros(size(img_array));
for c=1:size(sensitivities,1)
    sensitivities(c,:,:,:) = mask.*squeeze(img_array(c,:,:,:))./SOS;
    %sensitivities(c,:,:,:) = squeeze(img_array(c,:,:,:))./SOS;
end

%% generate SENSE matrix

sens1 = squeeze(mean(mean(mean(sensitivities(:,ap_range,rl_range,fh_range),2),3),4));
sens2 = squeeze(mean(mean(mean(sensitivities(:,ap_range2,rl_range2,fh_range2),2),3),4));

% Now pick only the coils that have been used in the actual MRS scan.
MRSCont.SENSE{kk}.nRecCoils = Ref.ncoils; % save number of receiver coils
for ii = 1:MRSCont.SENSE{kk}.nRecCoils
    for ll = 1:length(MRSCont.SENSE{kk}.sin_info.channel_names)
        MRSCont.SENSE{kk}.CoilsUsedInMRS(ii) = strcmp(MRSCont.SENSE{kk}.refsc_sin_info.channel_names{ii},MRSCont.SENSE{kk}.sin_info.channel_names{ll});
        if MRSCont.SENSE{kk}.CoilsUsedInMRS(ii) == 1
            break
        end
    end
end

toDelete = find(MRSCont.SENSE{kk}.CoilsUsedInMRS==0);
sens1(toDelete,:) = [];
sens2(toDelete,:) = [];
sensitivities(toDelete,:,:,:) = [];
noise_array(toDelete,:) = [];
%% Figure
%get sensitivities
if ~MRSCont.SENSE{kk}.sens_from_ref_scan

    ix = round(mean(ap_range));
    iy = round(mean(rl_range));
    iz = round(mean(fh_range));
    ix2 = round(mean(ap_range2));
    iy2 = round(mean(rl_range2));
    iz2 = round(mean(fh_range2));

    signal = sens1;
    signal = signal./sqrt(sum(abs(signal).^2));
    tmpimage = squeeze(sum(bsxfun(@times, sensitivities, conj(signal)),1));

    disp( ['orig. voxel positions; x1=' num2str(ix) '; y1=' num2str(iy) '; z1=' num2str(iz) ';'])
    [max_val, position] = max(abs(tmpimage(:)));
    [tx,ty,tz] = ind2sub(size(tmpimage),position);
    disp( ['opt. voxel positions; x1=' num2str(tx) '; y1=' num2str(ty) '; z1=' num2str(tz) ';'])
    disp(['similarity: ' num2str(max_val) ' vs ' num2str(abs(tmpimage(ix,iy,iz)))])

    figure(17)
    clf
    jet2 = jet;
    jet2(1,:) = [0;0;0];
    colormap(jet2)

    minpl = 0;
    maxpl = 1;
    
    tmpimage(ix,iy,1:end) = 0.3;
    tmpimage(ix,1:end,iz) = 0.3;
    tmpimage(1:end,iy,iz) = 0.3;

    subplot(2,3,1)

    imagesc(abs(rot90(squeeze(tmpimage(ix,:,:)))),[minpl maxpl])
    ylabel('y');ylabel('FH')
    xlabel('z');xlabel('LR')
    axis square;

    title('MPR')
    subplot(2,3,2)
    imagesc(abs(rot90(squeeze(tmpimage(:,iy,:)),3)),[minpl maxpl])
    ylabel('x');ylabel('AP')
    xlabel('z');xlabel('LR')
    axis square;
    set(gca,'YDir','normal');
    title('MPR')
    
    subplot(2,3,3)
    imagesc(abs(squeeze(tmpimage(:,:,iz))),[minpl maxpl])
    ylabel('x');ylabel('AP')
    xlabel('y');xlabel('FH')
    axis square;
    set(gca,'YDir','normal');
    title('MPR')

    signal = sens2;
    signal = signal./sqrt(sum(abs(signal).^2));
    tmpimage = squeeze(sum(bsxfun(@times, sensitivities, conj(signal)),1));
    tmpimage(ix2,iy2,1:end) = 0.3;
    tmpimage(ix2,1:end,iz2) = 0.3;
    tmpimage(1:end,iy2,iz2) = 0.3;

    disp( ['orig. voxel positions; x1=' num2str(ix2) '; y1=' num2str(iy2) '; z1=' num2str(iz2) ';'])
    [max_val, position] = max(abs(tmpimage(:)));
    [tx,ty,tz] = ind2sub(size(tmpimage),position);
    disp( ['opt. voxel positions; x2=' num2str(tx) '; y2=' num2str(ty) '; z2=' num2str(tz) ';'])
    disp(['similarity: ' num2str(max_val) ' vs ' num2str(abs(tmpimage(ix2,iy2,iz2)))])
    disp(' ')

    subplot(2,3,4)
    imagesc(abs(rot90(squeeze(tmpimage(ix2,:,:)))),[minpl maxpl])
    axis square;
    ylabel('y')
    xlabel('z')

    title('MPR')
    subplot(2,3,5)
    imagesc(abs(rot90(squeeze(tmpimage(:,iy2,:)),3)),[minpl maxpl])
    axis square;
    ylabel('x')
    xlabel('z')
    set(gca,'YDir','normal');

    title('MPR')
    subplot(2,3,6)
    imagesc(abs(squeeze(tmpimage(:,:,iz2))),[minpl maxpl])
    axis square;
    ylabel('x')
    xlabel('y')
    set(gca,'YDir','normal');
    title('MPR')

end
%%
% Normalize
sens1 = sens1./squeeze(sqrt(sum(abs(sens1).^2,1)));
sens2 = sens2./squeeze(sqrt(sum(abs(sens2).^2,1)));

% Build sensitivity matrix S
S = [ sens1 ...
    sens2 ...
    ];

PSI = noise_array*noise_array';
U = inv(S'*inv(PSI)*S)*S'*inv(PSI);

ga = inv(S'*inv(PSI)*S);
gb = S'*inv(PSI)*S;
MRSCont.SENSE{kk}.gfact = sqrt(ga.*gb);
MRSCont.SENSE{kk}.U = U;
MRSCont.SENSE{kk}.S = S;
MRSCont.SENSE{kk}.PSI = PSI;

disp(['gfactors: ' num2str(diag(real(MRSCont.SENSE{kk}.gfact))')])

end
