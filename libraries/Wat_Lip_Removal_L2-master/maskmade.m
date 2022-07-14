% csi: the csi data for generation of lipi mask
% num: determine range of the mask
function [metamask lipidmask] = maskmade(csi,num)
NN = size(csi); % size of csi data
lipstart = round(0.700 * NN(3)); % start point of lipid resonance
lipend = round(0.800 * NN(3)); % end point of lipid resonance
lipid = sum(abs(csi(:,:,lipstart:lipend)),3); % intergration of the lipid resonance
lipidmin = 0; % for finding voxel with lipid content greater than lipidmin
lipid_mask = lipid>lipidmin; % a matrix for labeling voxels with lipid content greater than lipidmin
lipidfind = find(lipid_mask>0); % find labeled voxels
NL = length(lipidfind); % the number of labeled voxels
while NL > num
    lipidmin=lipidmin+2;
    lipid_mask=lipid>lipidmin;
    lipidfind=find(lipid_mask>0);
    NL=length(lipidfind);
end
meta_mask=zeros(NN(1),NN(2));
add_mask=lipid_mask+meta_mask;
maskfind=find(add_mask>1);
N2=length(maskfind);
metamin=0;
while N2<8
    metamin=metamin+5;
    for i=1:NN(1)
        for j=1:NN(2)
            if 1*(abs(i-NN(1)/2))^2+1.4*(abs(j-NN(2)/2))^2<metamin || 1*(abs(i-NN(1)/2))^2+1.4*(abs(j-NN(2)/2-1))^2<metamin...
                    || 1*(abs(i-NN(1)/2-1))^2+1.4*(abs(j-NN(2)/2))^2<metamin || 1*(abs(i-NN(1)/2-1))^2+1.4*(abs(j-NN(2)/2-1))^2<metamin
            meta_mask(i,j)=1;
            end
        end
    end
    add_mask=lipid_mask+meta_mask;
    maskfind=find(add_mask==2);
    N2=length(maskfind);
end
metamask=meta_mask;
lipidmask=lipid_mask;
end

