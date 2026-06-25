function callLCModel(controlFile, pathLCModelBinary)
% Wrapper function for LCModel binary

callLCMCommand = ['"' pathLCModelBinary '" < "' controlFile '"'];
[status, result] = system(callLCMCommand);

if status % If LCModel throws an errpr, capture it and report in the Matlab command window
    error('LCModel error while executing the following control file:\n%s/n/n%s',controlFile,result);
end

end