function tableOut = showBasisFunctions(obj)
    % Returns a table showing the basis functions included in the
    % basis set, and whether they will be used in the fit.
    tableOut = table(obj.BasisSets.names', obj.BasisSets.includeInFit', ...
        'VariableNames', {'Basis function', 'Include in fit?'});
end