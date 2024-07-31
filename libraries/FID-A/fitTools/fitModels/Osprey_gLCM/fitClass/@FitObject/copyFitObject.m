function obj = copyFitObject(obj)
% Method to copy a OspreyFitObj
    try
        % R2010b or newer - directly in memory (faster)
        objByteArray = getByteStreamFromArray(obj);
        newObj = getArrayFromByteStream(objByteArray);
    catch
        % R2010a or earlier - serialize via temp file (slower)
        fname = [tempname '.mat'];
        save(fname, 'obj');
        newObj = load(fname);
        newObj = newObj.obj;
        delete(fname);
    end
end

