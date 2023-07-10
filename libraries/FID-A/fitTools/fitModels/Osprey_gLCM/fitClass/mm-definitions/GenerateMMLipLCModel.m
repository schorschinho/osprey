% Generator for default LCModel MM/Lipid parametrizations
% http://s-provencher.com/pub/LCModel/manual/manual.pdf
% Chapter 11.7
% Georg Oeltzschner, Johns Hopkins University 2023

% Create empty struct
MMLipLCModel = struct;

% Create Lip13a
MMLipLCModel.Lip13a.ppm = 1.28;
MMLipLCModel.Lip13a.fwmin = 0.15;
MMLipLCModel.Lip13a.amp = 2;
MMLipLCModel.Lip13a.sdppm = 0.01;
MMLipLCModel.Lip13a.fwex = 0.2;
MMLipLCModel.Lip13a.sdfw = 0.035;

% Create Lip13b
MMLipLCModel.Lip13b.ppm = 1.28;
MMLipLCModel.Lip13b.fwmin = 0.089;
MMLipLCModel.Lip13b.amp = 2;
MMLipLCModel.Lip13b.sdppm = 0.01;
MMLipLCModel.Lip13b.fwex = 0.09;
MMLipLCModel.Lip13b.sdfw = 0.035;

% Create Lip13c
MMLipLCModel.Lip13c.ppm = 1.30;
MMLipLCModel.Lip13c.fwmin = 0.089;
MMLipLCModel.Lip13c.amp = 2;
MMLipLCModel.Lip13c.sdppm = 0.01;
MMLipLCModel.Lip13c.fwex = 0.09;
MMLipLCModel.Lip13c.sdfw = 0.035;

% Create Lip13d
MMLipLCModel.Lip13d.ppm = 1.26;
MMLipLCModel.Lip13d.fwmin = 0.089;
MMLipLCModel.Lip13d.amp = 2;
MMLipLCModel.Lip13d.sdppm = 0.01;
MMLipLCModel.Lip13d.fwex = 0.09;
MMLipLCModel.Lip13d.sdfw = 0.035;

% Create Lip09
MMLipLCModel.Lip09.ppm = 0.89;
MMLipLCModel.Lip09.fwmin = 0.14;
MMLipLCModel.Lip09.amp = 3;
MMLipLCModel.Lip09.sdppm = 0.02;
MMLipLCModel.Lip09.fwex = 0.19;
MMLipLCModel.Lip09.sdfw = 0.035;

% Create MM09
MMLipLCModel.MM09.ppm = 0.91;
MMLipLCModel.MM09.fwmin = 0.14;
MMLipLCModel.MM09.amp = 3;
MMLipLCModel.MM09.sdppm = 0.02;
MMLipLCModel.MM09.fwex = 0.17;
MMLipLCModel.MM09.sdfw = 0.015;

% Create Lip20
MMLipLCModel.Lip20.ppm = [2.04; 2.25; 2.8];
MMLipLCModel.Lip20.fwmin = [0.15; 0.15; 0.2];
MMLipLCModel.Lip20.amp = [1.33; 0.67; 0.87];
MMLipLCModel.Lip20.sdppm = 0.005;
MMLipLCModel.Lip20.fwex = 0.2;
MMLipLCModel.Lip20.sdfw = 0.025;

% Create MM20
MMLipLCModel.MM20.ppm = [2.08; 2.25; 1.95; 3];
MMLipLCModel.MM20.fwmin = [0.15; 0.2; 0.15; 0.2];
MMLipLCModel.MM20.amp = [1.33; 0.33; 0.33; 0.4];
MMLipLCModel.MM20.sdppm = 0.005;
MMLipLCModel.MM20.fwex = 0.18;
MMLipLCModel.MM20.sdfw = 0.01;

% Create MM12
MMLipLCModel.MM12.ppm = 1.21;
MMLipLCModel.MM12.fwmin = 0.15;
MMLipLCModel.MM12.amp = 2;
MMLipLCModel.MM12.sdppm = 0.01;
MMLipLCModel.MM12.fwex = 0.2;
MMLipLCModel.MM12.sdfw = 0.02;

% Create MM14
MMLipLCModel.MM14.ppm = 1.43;
MMLipLCModel.MM14.fwmin = 0.17;
MMLipLCModel.MM14.amp = 2;
MMLipLCModel.MM14.sdppm = 0.02;
MMLipLCModel.MM14.fwex = 0.2;
MMLipLCModel.MM14.sdfw = 0.02;

% Create MM17
MMLipLCModel.MM17.ppm = 1.67;
MMLipLCModel.MM17.fwmin = 0.15;
MMLipLCModel.MM17.amp = 2;
MMLipLCModel.MM17.sdppm = 0.03;
MMLipLCModel.MM17.fwex = 0.17;
MMLipLCModel.MM17.sdfw = 0.02;

% Write to JSON
jsonOut = jsonencode(MMLipLCModel, PrettyPrint=true);
fid = fopen('MMLipLCModel.json','w');
fprintf(fid,'%s',jsonOut);
fclose(fid);

% Test if the written version is identical
MMLipLCModelTest = jsonToStruct('MMLipLCModel.json');
if isequal(MMLipLCModel, MMLipLCModelTest)
    disp('Test passed.');
else
    error('Something went wrong when the JSON file was written.');
end
