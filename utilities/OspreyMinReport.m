function [MRSCont] = OspreyMinReport(MRSCont)
%% [MRSCont] = OspreyMinReport(MRSCont)
%   This function generates a markdown report with the minimal reporting standards defined in
%   MRSinMRS (doi.org/10.1002/nbm.4448) by Lin et al. (2021).
%
%   USAGE:
%       MRSCont = OspreyMinReport(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2021-04-16)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2021-04-16: First version of the code.

% Fall back to defaults if not provided
if ~MRSCont.flags.didQuantify
    msg = 'Trying to create a MRSinMRS markdown report, but data have not been fully processed yet. Run all modules up to OspreyQuantify first.';
    fprintf(msg);
    error(msg);
end
outputFolder    = MRSCont.outputFolder;
fid=fopen(fullfile(outputFolder,'SummaryMRSinMRS.md'),'w+');
fprintf(fid,'\n # Summary following minimum reporting standards in MRS generated in Osprey');
fprintf(fid,'\n See Lin et al. ''Minimum Reporting Standards for in vivo Magnetic Resonance Spectroscopy (MRSinMRS): Experts'' consensus recommendations. NMR in Biomedicine. 2021;e4484.');
fprintf(fid,'\n [doi.org/10.1002/nbm.4448](https://doi.org/10.1002/nbm.4448)');
fprintf(fid,'\n \n');
fprintf(fid,'\n You may need to add information that was not avaialble in the raw data.');
fprintf(fid,'\n \n');
%% 1.Hardware part
fprintf(fid,'|1. Hardware|  | \n');
fprintf(fid,'|--|--| \n');
fprintf(fid,'|a. Field strength [T]| %d T| \n', round(MRSCont.raw{1}.Bo,2));
fprintf(fid,'|b. Manufacturer| %s| \n', MRSCont.vendor  );
if isfield(MRSCont.raw{1},'software')
    fprintf(fid,'|c. Model (software version if available)| %s| \n', MRSCont.raw{1}.software);
else
    fprintf(fid,'|c. Model (software version if available)| %s| \n', '-');
end
if isfield(MRSCont.raw{1},'nucleus')
    if contains(MRSCont.raw{1}.nucleus,'H')
        nucleus = '1H';
    else
        nucleus = MRSCont.raw{1}.nucleus;
    end
else %Let us assume it is aproton scan
    nucleus = '1H';
end
if MRSCont.raw{1}.dims.coils ~= 0
    fprintf(fid,'|d. RF coils: nuclei (transmit/receive), number of channels, type, body part|%s %d channel| \n',nucleus, MRSCont.raw{1}.sz(MRSCont.raw{1}.dims.coils));
else
    fprintf(fid,'|d. RF coils: nuclei (transmit/receive), number of channels, type, body part| %s | \n',nucleus);
end
fprintf(fid,'|e. Additional hardware| -| \n');
fprintf(fid,'\n \n');

%% 2. Acquisition part
fprintf(fid,'|2. Acquisition|  | \n');
fprintf(fid,'|--|--| \n');
fprintf(fid,'|a. Pulse sequence | %s| \n', MRSCont.fit.basisSet.seq{1}); % Here we need something for cases without fitting was done
fprintf(fid,'|b. Volume of interest (VOI) locations | %s| \n', '-');
if isfield(MRSCont.raw{1}, 'geometry')
    dim_names = fieldnames(MRSCont.raw{1}.geometry.size);
    dim1 = MRSCont.raw{1}.geometry.size.(dim_names{1});
    dim2 = MRSCont.raw{1}.geometry.size.(dim_names{2});
    dim3 = MRSCont.raw{1}.geometry.size.(dim_names{3});
    fprintf(fid,'|c. Nominal VOI size [mm<sup>3</sup>]| %d x %d x %d mm<sup>3</sup>| \n', dim1,  dim2,  dim3);
else
    fprintf(fid,'|c. Nominal VOI size [mm<sup>3</sup>]| - x - x - mm<sup>3</sup>| \n');
end
fprintf(fid,'|d. Repetition time (TR), echo time (TE) [ms]| TR %d ms, TE %d ms| \n', MRSCont.raw{1}.tr, MRSCont.raw{1}.te);
fprintf(fid,'|e. Total number of averages per spectrum <br> i. Number of averaged specra per subspectrum | %d total averages with %d averages per subspectrum| \n', MRSCont.raw{1}.rawAverages,MRSCont.raw{1}.rawAverages/MRSCont.raw{1}.rawSubspecs);
if MRSCont.flags.isUnEdited
    fprintf(fid,'|f. Additional sequence parameters | F1: %d Hz, %d points| \n', MRSCont.raw{1}.spectralwidth,MRSCont.raw{1}.sz(1));
else if MRSCont.flags.isMEGA
        switch MRSCont.opts.editTarget{1}
            case 'GABA'
               ON = 1.9;
               OFF = 7.5;
            case 'GSH'
               ON = 4.56;
               OFF = 7.5;
             case 'Lac'
               ON = 4.1;
               OFF = 7.5;
             case 'PE322'
               ON = 3.22;
               OFF = 7.5;
             case 'PE398'
               ON = 3.98;
               OFF = 7.5;               
        end
        fprintf(fid,'|f. Additional sequence parameters <br> i. editing pulse frequencies | F1: %d Hz, %d points <br> ppm<sub>ON</sub> = %.2f ppm<sub>OFF</sub> = %.2f | \n', ...
                MRSCont.raw{1}.spectralwidth,MRSCont.raw{1}.sz(1),ON,OFF);
    else
        fprintf(fid,'|f. Additional sequence parameters | F1: %d Hz, %d points| \n', MRSCont.raw{1}.spectralwidth,MRSCont.raw{1}.sz(1));
    end
end
fprintf(fid,'|g. Water suppression method | %s| \n', '-');
fprintf(fid,'|h. Shimming method, reference peak, and threshold of acceptance of shim chosen | %s| \n', '-');
fprintf(fid,'|i. Trigger or motion correction| %s| \n', '-');
fprintf(fid,'\n \n');
%% 3. Data analysis part
fprintf(fid,'|3. Data analysis methods and outputs|  | \n');
fprintf(fid,'|--|--| \n');
fprintf(fid,'|a. Analysis software | %s| \n', MRSCont.ver.Osp); % Here we need something for cases without fitting was done
fprintf(fid,'|b. Processing steps deviating from Osprey | %s| \n', 'None');
if MRSCont.flags.isUnEdited
    try
        strings = fieldnames(MRSCont.quantify.tables.A);
    catch
        strings = fieldnames(MRSCont.quantify.tables.off);
    end
end
if MRSCont.flags.isMEGA
    strings = fieldnames(MRSCont.quantify.tables.off);
end
if  MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
    strings = fieldnames(MRSCont.quantify.tables.sum);
end
outs = '';
for ss = 1 : length(strings)
    outs = [outs strings{ss}];
    if ss ~= length(strings)
        outs = [outs ', '];
    end
end
fprintf(fid,'|c. Output measure | %s \n', outs);
outs = '';
br = 1;
if MRSCont.flags.isUnEdited
    includeMetabs = MRSCont.fit.resBasisSet.off.(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name;
end
if MRSCont.flags.isMEGA
    includeMetabs = MRSCont.fit.resBasisSet.diff1.(MRSCont.info.diff1.unique_ndatapoint_spectralwidth{1}).name;
end
if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
    includeMetabs = MRSCont.fit.resBasisSet.diff1.(MRSCont.info.diff1.unique_ndatapoint_spectralwidth{1}).name;
end
for ss = 1 : length(includeMetabs)
    if br < 10
        outs = [outs includeMetabs{ss}];        
    else
        outs = [outs includeMetabs{ss}  ',<br>']; 
        br = 1;
    end
    if ss ~= length(includeMetabs) && br ~= 1
        outs = [outs ','];
    end
    br = br + 1;
end
fprintf(fid,'|d. Quantification references and assumptions, fitting model assumptions| Basis set list:<br> %s <br>Fitting method: %s basline knot spacing %.2f ppm\n', outs,MRSCont.opts.fit.method,MRSCont.opts.fit.bLineKnotSpace );
fprintf(fid,'\n \n');

%% 4. Data quality part
fprintf(fid,'|4. Data quality|  | \n');
fprintf(fid,'|--|--| \n');
if MRSCont.flags.isUnEdited || MRSCont.flags.isMEGA
    fprintf(fid,'|a. SNR (NAA), linewidth (NAA) [Hz] | SNR: %.0f +- %.0f, linewidth %.2f +- %.2f Hz| \n', round(mean(MRSCont.QM.SNR.A)),round(std(MRSCont.QM.SNR.A),2),round(mean(MRSCont.QM.FWHM.A),2),round(std(MRSCont.QM.FWHM.A),2)); % Here we need something for cases without fitting was done
end
if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
    fprintf(fid,'|a. SNR (NAA), linewidth (NAA) [Hz] | SNR: %.0f +- %.0f, linewidth %.2f +- %.2f Hz| \n', round(mean(MRSCont.QM.SNR.sum)),round(std(MRSCont.QM.SNR.sum),2),round(mean(MRSCont.QM.FWHM.sum),2),round(std(MRSCont.QM.FWHM.sum),2)); % Here we need something for cases without fitting was done
end
fprintf(fid,'|b. Data exclusion criteria | %s| \n', 'None');
if MRSCont.flags.isUnEdited
    fprintf(fid,'|c. Quality measures of postporcessing model fitting (Mean Relative Amplitude Residual) | %.2f %% \n', round(mean(MRSCont.QM.relAmpl.A),2));
end
if MRSCont.flags.isMEGA
    fprintf(fid,'|c. Quality measures of postporcessing model fitting (Mean Relative Amplitude Residual (Residual/Noise)) <br> off <br> diff1 |<br> %.2f +- %.2f  <br> %.2f +- %.2f \n', round(mean(MRSCont.QM.relAmpl.A),2), round(std(MRSCont.QM.relAmpl.A),2),round(mean(MRSCont.QM.relAmpl.diff1),2),round(std(MRSCont.QM.relAmpl.diff1),2));
end
if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
    fprintf(fid,'|c. Quality measures of postporcessing model fitting (Mean Relative Amplitude Residual (Residual/Noise)) <br> sum <br> diff1 <br> diff2 |<br> %.2f +- %.2f  <br> %.2f +- %.2f  <br> %.2f +- %.2f \n', round(mean(MRSCont.QM.relAmpl.sum),2),round(std(MRSCont.QM.relAmpl.sum),2),round(mean(MRSCont.QM.relAmpl.diff1),2),round(std(MRSCont.QM.relAmpl.diff1),2),round(mean(MRSCont.QM.relAmpl.diff2),2),round(std(MRSCont.QM.relAmpl.diff2),2));
end
fprintf(fid,'|d. Mean spectrum created with OspreyOverview| %s \n', 'Figure 1' );

fclose(fid);
end