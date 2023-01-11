function excludeBasisFunctionFromFit(obj, input)
    % Sets the 'includeInFit' flag in the 'BasisSets' property to
    % 0, meaning that the corresponding basis function will be
    % excluded from in the fit.
    
    if ischar(input)
        input = {input};
    end
    
    if strcmpi(input, 'all')
        input = obj.BasisSets.names;
    end
    
    for ii = length(input)
        % Check which metabolites are available in the basis set
        % and match the input
        [metsToInclude, ~, ~] = intersect(obj.BasisSets.names, input, 'stable');
        for rr = 1:length(metsToInclude)
            idxToInclude = find(strcmp(metsToInclude{rr}, obj.BasisSets.names));
            obj.BasisSets.includeInFit(obj.step+1,idxToInclude) = 0;
        end
        
    end
end
        