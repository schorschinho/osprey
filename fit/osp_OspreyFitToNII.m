function out = osp_OspreyFitToNII(in, fitParams,resBasisSet,scale,fitOpts,isMRSI,whichSubspec,mask);
    if nargin < 7    
        whichSubspec = 'off';
        if nargin < 6
                isMRSI = 0;
        end
    end

    out = in;

    % Prepare input for parameter to model conversion   
    inputSettings.scale                 = scale;
    if ~strcmp(whichSubspec,'w')
        inputSettings.fitRangePPM           = fitOpts.range;
    else
        inputSettings.fitRangePPM           = fitOpts.rangeWater;
    end
    inputSettings.minKnotSpacingPPM     = fitOpts.bLineKnotSpace;
    try
        inputSettings.GAP = fitOpts.GAP;
    catch
    end
    inputData.dataToFit = out;
    inputData.basisSet = resBasisSet;

    if ~isMRSI
        % Get model results
        [ModelOutput] = fit_OspreyParamsToModel(inputData, inputSettings, fitParams);
    
        % Phase data according to fit
        out = op_zeropad(out,2);
        out=op_addphase(out,-fitParams.ph0);
        pivotPoint = 4.68;
        multiplier = ModelOutput.ppm - pivotPoint;
        ph1 = fitParams.ph1 * pi/180;
        start_index = find(out.ppm==ModelOutput.ppm(1));
        end_index = find(out.ppm==ModelOutput.ppm(end));
        out.specs(start_index:end_index,1)   = out.specs(start_index:end_index,1) .* (exp(-1i*ph1*multiplier)).';
        out.fids = ifft(fftshift(out.specs,1),[],1);
    
        % Prepare model
        out = addModelResults(out, ModelOutput,scale);
    else
        inputData.dataToFit =op_takeVoxel(in,[1 1 1]);
        zero_fill_sz = out.sz;        
        if ~strcmp(whichSubspec,'w') 
            zero_fill_sz(1) = zero_fill_sz(1)*2;
            [ModelOutput] = fit_OspreyParamsToModel(inputData, inputSettings, fitParams{1,1,1}.(whichSubspec).fitParams{1});      
            pivotPoint = 4.68;
            multiplier = ModelOutput.ppm - pivotPoint;
        else
            [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams{1,1,1}.(whichSubspec).fitParams{1});
        end
        
        if ~strcmp(whichSubspec,'w')
            zero_fill_sz = [zero_fill_sz(1) 1 1 size(ModelOutput.indivMets,2)+4 zero_fill_sz(2:end)];
        else
            zero_fill_sz = [zero_fill_sz(1) 1 1 4 zero_fill_sz(2:end)];
        end
        out.fids = zeros(zero_fill_sz);
        out.specs = zeros(zero_fill_sz);
        out.sz = size(out.fids);
        out.dims.averages = 2;
        out.dims.subSpecs = 3;
        out.dims.extras = 4;
        out.dims.Xvoxels = 5;
        out.dims.Yvoxels = 6;
        out.dims.Zvoxels = 7;
        for z = 1 : out.nZvoxels
            for y = 1 : out.nYvoxels
                for x = 1 : out.nXvoxels                    
                    temp = op_takeVoxel(in,[x y z]);
                    if ~strcmp(whichSubspec,'w') 
                        temp = op_zeropad(temp,2);
                    end
                    fprintf('\nTransforming model to NIfTI-MRS (%d, %d, %d) of the MRSI dataset...\n', x, y, z); 
                    if mask(x,y,z)>0
                        if ~strcmp(whichSubspec,'w') 
                            inputData.dataToFit =op_takeVoxel(in,[x y z]);  
                            temp   = op_ampScale(temp, 1/scale);
                            temp = op_freqshift(temp, -fitParams{x,y,z}.(whichSubspec).fitParams{1}.refShift);
                            temp=op_addphase(temp,-fitParams{x,y,z}.(whichSubspec).fitParams{1}.ph0);
                            ph1 = fitParams{x,y,z}.(whichSubspec).fitParams{1}.ph1 * pi/180;
                            % Get model results
                            [ModelOutput] = fit_OspreyParamsToModel(inputData, inputSettings, fitParams{x,y,z}.(whichSubspec).fitParams{1});
                            start_index = find(temp.ppm==ModelOutput.ppm(1));
                            end_index = find(temp.ppm==ModelOutput.ppm(end));
                            temp.specs(start_index:end_index,1)   = temp.specs(start_index:end_index,1) .* (exp(-1i*ph1*multiplier)).';
                            % Prepare model
                            temp = addModelResults(temp, ModelOutput,1);
                        else
                            inputData.dataToFit =op_takeVoxel(in,[x y z]);  
                            temp   = op_ampScale(temp, 1/scale);
                            % Get model results
                            [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams{x,y,z}.(whichSubspec).fitParams{1});
                            % Prepare model
                            temp = addModelResults(temp, ModelOutput,1);
                        end
                        out = op_addVoxel(out,temp,[x y z]);
                    else
                        temp   = op_ampScale(temp, 1/scale);
                        out.fids(:,1,1,1,x,y,z) = temp.fids(:);
                        out.specs(:,1,1,1,x,y,z) = temp.specs(:);
                    end
                end
            end
        end
        out.fids = ifft(fftshift(out.specs,1),[],1);
        out.t = temp.t;
        out.ppm = temp.ppm;
    end

        % figure, plot(ModelOutput.ppm,ModelOutput.data), hold on, plot(out.ppm, squeeze(real(out.specs(:,1,1,1))))


end

function out = addModelResults(in, Model,scale)
    out = in;
    start_index = find(out.ppm==Model.ppm(1));
    end_index = find(out.ppm==Model.ppm(end));

    TransModel.baseline = zeros(in.sz(1),1);
    TransModel.completeFit = zeros(in.sz(1),1);
    TransModel.residual = zeros(in.sz(1),1);
    if isfield(Model,'indivMets')
        TransModel.indivMets = zeros(in.sz(1),size(Model.indivMets,2));
    end

    if isfield(Model,'baseline')
       TransModel.baseline(start_index:end_index) = Model.baseline;
    end    
    TransModel.completeFit(start_index:end_index) = Model.completeFit;
    TransModel.residual(start_index:end_index) = Model.residual;
    if isfield(Model,'indivMets')
    TransModel.indivMets(start_index:end_index,:) = Model.indivMets;
    end

    TransModel.baseline = hilbert(TransModel.baseline);
    TransModel.completeFit = hilbert(TransModel.completeFit);
    TransModel.residual = hilbert(TransModel.residual);
    if isfield(Model,'indivMets')
        for rr = 1:size(TransModel.indivMets,2)
            TransModel.indivMets(:,rr) = hilbert(TransModel.indivMets(:,rr));
        end
    end

    out.specs(:,1,1,2)= TransModel.completeFit;
    out.specs(:,1,1,3)= TransModel.residual;
    out.specs(:,1,1,4)= TransModel.baseline;
    if isfield(Model,'indivMets')
        out.specs(:,1,1,5:4+size(TransModel.indivMets,2))= TransModel.indivMets + TransModel.baseline;
    end
    out.fids = ifft(fftshift(out.specs,1),[],1);
    out.sz = size(out.fids);
    out   = op_ampScale(out, 1/scale);
end