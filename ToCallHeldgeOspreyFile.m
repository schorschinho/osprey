MRSCont = OspreyJob('jobPhilipsMPPRIAM.m');
MRSCont = OspreyLoad(MRSCont); 
MRSCont = OspreyProcess(MRSCont); 
MRSCont = OspreyFit(MRSCont);


MRSCont = OspreyJob('jobSDAT.m');


out = osp_plotLoad(MRSCont, 1, 'mets',1); %(data_number, voxel_1)
out = osp_plotProcess(MRSCont, 1, 'diff1',1);
out = osp_plotFit(MRSCont, 1, 'diff1',1); %(data_number, voxel_1)


out = osp_plotLoad(MRSCont, 1, 'mets',2); %(data_number, voxel_2)
out = osp_plotProcess(MRSCont, 1, 'diff1',2);
out = osp_plotFit(MRSCont, 1, 'diff1',2); % (data_number, voxel_2)