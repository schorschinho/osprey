function [MRSCont] = create_quickMaps(MRSCont)
%% [MRSCont] = create_quickMaps(MRSCont)
%   This function generates amplitude integral maps for quick inspection. You can
%   define different regions and spectra to be used.
%
%   USAGE:
%       MRSCont = create_quickMaps(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-10-31)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-10-31: First version of the code.
%% Genreate integral maps

% Pick target spectra
switch MRSCont.opts.MRSI.quickMaps.target
    case 'raw'
        % Loop over spectra listed in options to create dummy maps
        for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
            for ll = 1 : length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss})  )
                MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss} ).(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll})= zeros(MRSCont.(MRSCont.opts.MRSI.quickMaps.specs{ss} ){1}.sz(2:end));
            end
        end
        
        % Loop over datasets
        for kk =1 : MRSCont.nDatasets(1)
            % Loop over spectra listed in options
            for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
                % Magnitude or real part?
                if MRSCont.opts.MRSI.quickMaps.abs
                    specs = abs(MRSCont.(MRSCont.opts.MRSI.quickMaps.specs{ss}){kk}.specs);
                else
                    specs = real(MRSCont.(MRSCont.opts.MRSI.quickMaps.specs{ss}){kk}.specs);
                end
                % Create sum spectrum
                if MRSCont.(MRSCont.opts.MRSI.quickMaps.specs{ss}){kk}.dims.subSpecs > 0
                    specs = squeeze(sum(specs,MRSCont.(MRSCont.opts.MRSI.quickMaps.specs{ss}){kk}.dims.subSpecs));
                end
                % Do the integration across the different regions
                for ll = 1 : length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss})  )
                    MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll}) = ...
                    squeeze(sum(specs(MRSCont.(MRSCont.opts.MRSI.quickMaps.specs{ss}){kk}.ppm > MRSCont.opts.MRSI.quickMaps.limits.(MRSCont.opts.MRSI.quickMaps.specs{ss} )(ll,1) & ...
                              MRSCont.(MRSCont.opts.MRSI.quickMaps.specs{ss}){kk}.ppm < MRSCont.opts.MRSI.quickMaps.limits.(MRSCont.opts.MRSI.quickMaps.specs{ss} )(ll,2),:,:,:),1));
                end  
            end
        end
    case 'processed'
        % Loop over spectra listed in options
        for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
            for ll = 1 : length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss})  )
                MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss} ).(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll})= zeros(MRSCont.processed.(MRSCont.opts.MRSI.quickMaps.specs{ss} ){1}.sz(2:end));
            end
        end
        
        % Loop over datasets
        for kk =1 : MRSCont.nDatasets(1)
            % Loop over spectra listed in options
            for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
                % Magnitude or real part?
                if MRSCont.opts.MRSI.quickMaps.abs
                    specs = abs(MRSCont.processed.(MRSCont.opts.MRSI.quickMaps.specs{ss}){kk}.specs);
                else
                    specs = real(MRSCont.processed.(MRSCont.opts.MRSI.quickMaps.specs{ss}){kk}.specs);
                end
                % Do the integration across the different regions
                for ll = 1 : length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss})  )
                    MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll}) = ...
                    squeeze(sum(specs(MRSCont.processed.(MRSCont.opts.MRSI.quickMaps.specs{ss}){kk}.ppm > MRSCont.opts.MRSI.quickMaps.limits.(MRSCont.opts.MRSI.quickMaps.specs{ss} )(ll,1) & ...
                              MRSCont.processed.(MRSCont.opts.MRSI.quickMaps.specs{ss}){kk}.ppm < MRSCont.opts.MRSI.quickMaps.limits.(MRSCont.opts.MRSI.quickMaps.specs{ss} )(ll,2),:,:,:),1));
                end  
            end
        end
end

%% Export the results as NIfTI files
switch MRSCont.opts.MRSI.quickMaps.target
    case 'raw'
        files = [];
        outputFolder = MRSCont.outputFolder;
        data = MRSCont.raw{1};
        if (MRSCont.raw{1}.nZvoxels > 1) && ~MRSCont.opts.MRSI.pseudo3D 
                reorder = flip(1:data.nZvoxels);
                for ll = 1 : data.nZvoxels
                        ToExport = data;
                        if isfield(ToExport.geometry,'slice_distance')
                            VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                            ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice 
                        end
                        mkdir(fullfile(outputFolder,['slice_' num2str(reorder(ll))],'quickMaps'))
                        out.hdr = ToExport.nii_mrs.hdr;
                        out.hdr.dim(1) = 3;
                        out.hdr.dim(2) = 1;
                        out.hdr.pixdim(5) = 1;
                        out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes 
         
                        for mm = 1 : length(metab_names)
                                out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,ll))./MRSI_model_water.amplitudes(:,:,ll) *55500);         
                            
                                out.img(isnan(out.img)) =0;
                                out.img(isinf(out.img)) =0;
                                nii_tool('save', out, fullfile(outputFolder,['slice_' num2str(reorder(ll))],'concs','rawWaterScaled',[metab_names{mm}  '.nii.gz']));
                            
                        end
                           
                    shift = shift - 1;
                end
            else
                ToExport = data;
                mkdir(fullfile(outputFolder,'quickMaps'))
                out.hdr = ToExport.nii_mrs.hdr;
                out.hdr.dim(1) = 3;
                out.hdr.dim(2) = 1;
                out.hdr.pixdim(5) = 1;
                out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes 
        
                    for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
                        for ll = 1 : length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss})  )
                            out.img = MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll}); 
                            out.img(isnan(out.img)) =0;
                            out.img(isinf(out.img)) =0;
                            out.img = flip(out.img,3);
                            out.img = flip(out.img,1);
                            MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll})=out.img;
                            nii_tool('save', out, fullfile(outputFolder,'quickMaps',[MRSCont.opts.MRSI.quickMaps.specs{ss} '_' MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll}  '.nii.gz']));
                            files{end+1} = fullfile(outputFolder,'quickMaps',[MRSCont.opts.MRSI.quickMaps.specs{ss} '_' MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll}  '.nii.gz']);
                        end
                    end
        end
        
        % Do interpolation if needed
        if MRSCont.opts.MRSI.quickMaps.interpolation > 1
            for ff = 1 : length(files)
                gunzip(files{ff});
                resize_img(files{ff}(1:end-3), files{1}(1:end-3), 0, MRSCont.opts.MRSI.quickMaps.interpolation); 
                [outFolder,outFile,~]= fileparts(files{ff});
                struct_name_parts = strsplit(outFile, '_');
                if length(struct_name_parts) == 3
                    struct_name_parts{1} = [struct_name_parts{1} '_' struct_name_parts{2}];
                    struct_name_parts{2} = struct_name_parts{3};
                    struct_name_parts = struct_name_parts(1:2);
                end
                metab_img = spm_vol(fullfile(outFolder,['r' outFile]));
                [metab_img,~]    = spm_read_vols(metab_img);
                MRSCont.quickMapsInt.(struct_name_parts{1}).(struct_name_parts{2}(1:end-4)) =double(metab_img);
                gzip(fullfile(outFolder,['r' outFile]));
                delete(fullfile(outFolder,['r' outFile]));
                if ff > 1
                    delete(files{ff}(1:end-3));
                end
            end
            delete(files{1}(1:end-3));
        end
    case 'processed'
        files = [];
        outputFolder = MRSCont.outputFolder;
        data = MRSCont.processed.A{1};
        if (MRSCont.raw{1}.nZvoxels > 1) && ~MRSCont.opts.MRSI.pseudo3D 
                reorder = flip(1:data.nZvoxels);
                for ll = 1 : data.nZvoxels
                        ToExport = data;
                        if isfield(ToExport.geometry,'slice_distance')
                            VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                            ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice 
                        end
                        mkdir(fullfile(outputFolder,['slice_' num2str(reorder(ll))],'quickMaps'))
                        out.hdr = ToExport.nii_mrs.hdr;
                        out.hdr.dim(1) = 3;
                        out.hdr.dim(2) = 1;
                        out.hdr.pixdim(5) = 1;
                        out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes 
         
                        for mm = 1 : length(metab_names)
                                out.img = squeeze(squeeze(MRSI_model.amplitudes(mm,:,:,ll))./MRSI_model_water.amplitudes(:,:,ll) *55500);         
                            
                                out.img(isnan(out.img)) =0;
                                out.img(isinf(out.img)) =0;
                                MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll})=out.img;
                                nii_tool('save', out, fullfile(outputFolder,'quickMaps',[MRSCont.opts.MRSI.quickMaps.specs{ss} '_' MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll}  '.nii.gz']));
                                files{end+1} = fullfile(outputFolder,'quickMaps',[MRSCont.opts.MRSI.quickMaps.specs{ss} '_' MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll}  '.nii.gz']);
                            
                        end
                           
                    shift = shift - 1;
                end
            else
                ToExport = data;
                mkdir(fullfile(outputFolder,'quickMaps'))
                out.hdr = ToExport.nii_mrs.hdr;
                out.hdr.dim(1) = 3;
                out.hdr.dim(2) = 1;
                out.hdr.pixdim(5) = 1;
                out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes 
        
                for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
                    for ll = 1 : length(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss})  )
                        out.img = MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll}); 
                        out.img(isnan(out.img)) =0;
                        out.img(isinf(out.img)) =0;
                        out.img = flip(out.img,3);
                        out.img = flip(out.img,1);
                        MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).(MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll})=out.img;
                        nii_tool('save', out, fullfile(outputFolder,'quickMaps',[MRSCont.opts.MRSI.quickMaps.specs{ss} '_' MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll}  '.nii.gz']));
                        files{end+1} = fullfile(outputFolder,'quickMaps',[MRSCont.opts.MRSI.quickMaps.specs{ss} '_' MRSCont.opts.MRSI.quickMaps.names.(MRSCont.opts.MRSI.quickMaps.specs{ss}){ll}  '.nii.gz']);
                    end
                end
                
                % Export the SNR/FWHM maps
                for z = 1 : MRSCont.processed.A{1}.nZvoxels
                    for x = 1 : MRSCont.processed.A{1}.nXvoxels
                        for y = 1 : MRSCont.processed.A{1}.nYvoxels
                            for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
                                MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).SNR(x,y,z) = MRSCont.QM{x,y,z}.SNR.(MRSCont.opts.MRSI.quickMaps.specs{ss});
                                MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).FWHM(x,y,z) = MRSCont.QM{x,y,z}.FWHM.(MRSCont.opts.MRSI.quickMaps.specs{ss});
                            end
                        end
                    end
                end             
                
                MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).SNR = flip(MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).SNR,3);
                MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).SNR = flip(MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).SNR,1);
                MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).FWHM = flip(MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).FWHM,3);
                MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).FWHM = flip(MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).FWHM,1);

                for ss = 1 : length(MRSCont.opts.MRSI.quickMaps.specs)
                    out.img = MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).SNR; 
                    out.img(isnan(out.img)) =0;
                    out.img(isinf(out.img)) =0;
                    out.img = flip(out.img,3);
                    out.img = flip(out.img,1);
                    nii_tool('save', out, fullfile(outputFolder,'quickMaps',[MRSCont.opts.MRSI.quickMaps.specs{ss} '_SNR.nii.gz']));
                    files{end+1} = fullfile(outputFolder,'quickMaps',[MRSCont.opts.MRSI.quickMaps.specs{ss} '_SNR.nii.gz']);
                    out.img = MRSCont.quickMaps.(MRSCont.opts.MRSI.quickMaps.specs{ss}).FWHM; 
                    out.img(isnan(out.img)) =0;
                    out.img(isinf(out.img)) =0;
                    out.img = flip(out.img,3);
                    out.img = flip(out.img,1);
                    nii_tool('save', out, fullfile(outputFolder,'quickMaps',[MRSCont.opts.MRSI.quickMaps.specs{ss} '_FWHM.nii.gz']));
                    files{end+1} = fullfile(outputFolder,'quickMaps',[MRSCont.opts.MRSI.quickMaps.specs{ss} '_FWHM.nii.gz']);
                end

                % Export the abs SNR
                specs = abs(MRSCont.processed.A{kk}.specs);

                MRSCont.quickMaps.A.abstNAASNR = ...
                squeeze(max(specs(MRSCont.processed.A{kk}.ppm > 1.95 & ...
                          MRSCont.processed.A{kk}.ppm < 2.1,:,:,:),[],1)) ./ ...
                squeeze(std(specs(MRSCont.processed.A{kk}.ppm > -2 & ...
                          MRSCont.processed.A{kk}.ppm < 0 ,:,:,:),1));
                
                MRSCont.quickMaps.A.abstNAASNR = flip(MRSCont.quickMaps.A.abstNAASNR,3);
                MRSCont.quickMaps.A.abstNAASNR = flip(MRSCont.quickMaps.A.abstNAASNR,1);

                if isfield(MRSCont.processed,'AFID') % Also do for FID only data for maximum echo MRSI
                    specs = abs(MRSCont.processed.AFID{kk}.specs);
    
                    MRSCont.quickMaps.AFID.abstNAASNR = ...
                    squeeze(max(specs(MRSCont.processed.AFID{kk}.ppm > 1.95 & ...
                              MRSCont.processed.AFID{kk}.ppm < 2.1,:,:,:),[],1)) ./ ...
                    squeeze(std(specs(MRSCont.processed.A{kk}.ppm > -2 & ...
                              MRSCont.processed.AFID{kk}.ppm < 0 ,:,:,:),1));

                    MRSCont.quickMaps.AFID.abstNAASNR = flip(MRSCont.quickMaps.AFID.abstNAASNR,3);
                    MRSCont.quickMaps.AFID.abstNAASNR = flip(MRSCont.quickMaps.AFID.abstNAASNR,1);
                end
             
                
                out.img = MRSCont.quickMaps.A.abstNAASNR; 
                out.img(isnan(out.img)) =0;
                out.img(isinf(out.img)) =0;
                out.img = flip(out.img,3);
                out.img = flip(out.img,1);
                nii_tool('save', out, fullfile(outputFolder,'quickMaps','A_abstNAA_SNR.nii.gz'));
                files{end+1} = fullfile(outputFolder,'quickMaps','A_abstNAA_SNR.nii.gz');

                if isfield(MRSCont.processed,'AFID') % Also do for FID only data for maximum echo MRSI
                    out.img = MRSCont.quickMaps.AFID.abstNAASNR; 
                    out.img(isnan(out.img)) =0;
                    out.img(isinf(out.img)) =0;
                    out.img = flip(out.img,3);
                    out.img = flip(out.img,1);
                    nii_tool('save', out, fullfile(outputFolder,'quickMaps','AFID_abstNAA_SNR.nii.gz'));
                    files{end+1} = fullfile(outputFolder,'quickMaps','AFID_abstNAA_SNR.nii.gz');
                end


        end
        
         % Do interpolation if needed
        if MRSCont.opts.MRSI.quickMaps.interpolation > 1
            for ff = 1 : length(files)
                gunzip(files{ff});
                resize_img(files{ff}(1:end-3), files{1}(1:end-3), 0, MRSCont.opts.MRSI.quickMaps.interpolation); 
                [outFolder,outFile,~]= fileparts(files{ff});
                struct_name_parts = strsplit(outFile, '_');
                if length(struct_name_parts) == 3
                    struct_name_parts{1} = [struct_name_parts{1} '_' struct_name_parts{2}];
                    struct_name_parts{2} = struct_name_parts{3};
                    struct_name_parts = struct_name_parts(1:2);
                end
                metab_img = spm_vol(fullfile(outFolder,['r' outFile]));
                [metab_img,~]    = spm_read_vols(metab_img);
                MRSCont.quickMapsInt.(struct_name_parts{1}).(struct_name_parts{2}(1:end-4)) =double(metab_img);
                gzip(fullfile(outFolder,['r' outFile]));
                delete(fullfile(outFolder,['r' outFile]));
                if ff > 1
                    delete(files{ff}(1:end-3));
                end
            end
            delete(files{1}(1:end-3));
        end
    end
end
