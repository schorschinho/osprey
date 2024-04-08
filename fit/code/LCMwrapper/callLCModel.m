function callLCModel(controlFile, pathLCModelBinary)
% Wrapper function for LCModel binary

callLCMCommand = ['"' pathLCModelBinary '" < "' controlFile '"'];
system(callLCMCommand);

end