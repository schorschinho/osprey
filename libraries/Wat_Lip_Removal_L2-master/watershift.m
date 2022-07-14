function [csi11, csi21, csi31]= watershift(csi1, csi2, csi3)
%B0 correction with non-water suppression csi scan as reference
csi00 = abs(csi1);  %% csi1 the water reference scan
N=size(csi1);
csi11=zeros(N); %% frequency shifted reference data
csi21=zeros(N); %% frequency shifted csi data
csi31=zeros(N);
for k=1:N(1)
    for l=1:N(2)
        aa=csi00(k,l,:);
        [a b]=max(aa);
%         b = b - 10;
        if abs(b-N(3)/2) < (N(3)/10) %% restrict the shifting range
            csi11(k,l,:) = circshift(csi1(k,l,:),[0 0 -b+N(3)/2]);
            if nargin == 2
                csi21(k,l,:) = circshift(csi2(k,l,:),[0 0 -b+N(3)/2]);
            end
            if nargin == 3
                csi21(k,l,:) = circshift(csi2(k,l,:),[0 0 -b+N(3)/2]);
                csi31(k,l,:) = circshift(csi3(k,l,:),[0 0 -b+N(3)/2]);
            end
        else
            csi11(k,l,:) = csi1(k,l,:);
            if nargin == 2
                csi21(k,l,:) = csi2(k,l,:);
            end
            if nargin == 3
                csi21(k,l,:) = csi2(k,l,:);
                csi31(k,l,:) = csi3(k,l,:);
            end
        end
    end
end
end

