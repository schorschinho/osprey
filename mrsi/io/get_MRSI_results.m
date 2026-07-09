function MRSI_model = get_MRSI_results(FitMatrix)
%% MRSI_model = get_MRSI_results(FitMatrix)
%   Gets specta  from fit object into a struct of
%   matricies
%
%   USAGE:
%       MRSI_model = get_MRSI_results(FitMatrix)
%
%   INPUTS:
%       FitMatrix = Matrix of Osprey 3.0 fit objects
%
%   OUTPUTS:
%       out     = MRSI model struct with spectra
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.
%% Export data

    % Get the matrix dimensions
    dims = size(FitMatrix);
    if length(dims) == 2
        dims(end+1) = 1;
    end
    
    % Get fit object to setup loops
    non_zero = find(~cellfun('isempty', FitMatrix));
    tempModel = FitMatrix{non_zero(1)};

    
    if size(tempModel.Data.fids,2) == 1
        % Initialize empty matrices
        MRSI_model.data = zeros(size(tempModel.Data.fids,1),size(FitMatrix,1),size(FitMatrix,2),size(FitMatrix,3));
        MRSI_model.fit = MRSI_model.data;
        MRSI_model.baseline = MRSI_model.data;
        MRSI_model.residual = MRSI_model.data;
        MRSI_model.metabs = zeros(size(tempModel.Data.fids,1),size(tempModel.Model{tempModel.step}.fit.metabs,2),size(FitMatrix,1),size(FitMatrix,2),size(FitMatrix,3));
        
        % Loop over results to fill the struct
        for z = 1: dims(3)
            for x = 1: dims(1)
                for y = 1: dims(2)
                    if ~isempty(FitMatrix{x,y,z})
                        temp = FitMatrix{x,y,z}.returnModel(tempModel.step);
                        MRSI_model.data(:,x,y,z) = temp.data;
                        MRSI_model.fit(:,x,y,z) = temp.fit;
                        MRSI_model.baseline(:,x,y,z) = temp.baseline;
                        MRSI_model.residual(:,x,y,z) = temp.residual;
                        MRSI_model.metabs(:,:,x,y,z) = temp.metabs;
                    end
                end
            end
        end
    else
        % Initialize empty matrices
        MRSI_model.data = zeros(size(tempModel.Data.fids,1),size(FitMatrix,1),size(FitMatrix,2),size(FitMatrix,3),size(tempModel.Data.fids,2));
        MRSI_model.fit = MRSI_model.data;
        MRSI_model.baseline = MRSI_model.data;
        MRSI_model.residual = MRSI_model.data;
        MRSI_model.metabs = zeros(size(tempModel.Data.fids,1),size(tempModel.Model{tempModel.step}.fit.metabs,2),size(FitMatrix,1),size(FitMatrix,2),size(FitMatrix,3),size(tempModel.Data.fids,2));
        
        % Loop over results to fill the struct
        for z = 1: dims(3)
            for x = 1: dims(1)
                for y = 1: dims(2)
                    if ~isempty(FitMatrix{x,y,z})
                        temp = FitMatrix{x,y,z}.returnModel(FitMatrix{x,y,z}.step,NaN);
                        MRSI_model.data(:,x,y,z,:) = temp.data;
                        MRSI_model.fit(:,x,y,z,:) = temp.fit;
                        MRSI_model.baseline(:,x,y,z,:) = temp.baseline;
                        MRSI_model.residual(:,x,y,z,:) = temp.residual;
                        MRSI_model.metabs(:,:,x,y,z,:) = temp.metabs;
                    end
                end
            end
        end
    end

end