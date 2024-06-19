% Run LCM similar to old Osprey model
% Change basis set files in
% 2Step_Spline_invivo_GroupingPars_soft_constraint.json create one file for
% each TI

% To change polarity of basis sets modify basis set in Osprey_gLCM.m after line 135 

for kk = 1 : MRSCont.nDatasets(1)
    MRSCont.processed.metab{kk}.nucleus={'1H'};
end
[OldOsprey] = Osprey_gLCM(MRSCont.processed.metab{2},'which(fullfile(''Osprey_gLCM'',''fitClass'',''model-procedures'',''2Step_Spline_invivo_GroupingPars_soft_constraint.json''))',0,0,1);
MRSCont.fit.metab = OldOsprey;
MRSCont = OspreyQuantify(MRSCont);