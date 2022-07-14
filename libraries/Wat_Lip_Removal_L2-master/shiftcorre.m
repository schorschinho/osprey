function [ csi_sc ] = shiftcorre( csi )
csi_abs = abs(csi); %% absolute value spectra
csi_abs = mrs_baselinecorre(csi_abs,400); %% baseline correction
N=size(csi_abs); %% csi data size
csi_sc=zeros(N); %% csi with shift corrected
start1=round(0.578*N(3)); %% start searching point
end1=round(0.6055*N(3)); %% end searching point
x=start1:end1;
distance1=round(0.0115*N(3)); %% first to second point
distance2=round(0.0764*N(3)); %% first to third point
distance3=round(0.02*N(3));
LL = length(N); %% data dimensions 3/4
if LL == 3
    for k=1:N(1)
        for l=1:N(2)
            aa=csi_abs(k,l,x)+csi_abs(k,l,x+distance1)+0.3*csi_abs(k,l,x+distance2);
            [a b]=max(aa);
            csi_sc(k,l,:)=circshift(csi(k,l,:),[0 0 -b+distance3]);
        end
    end
end
if LL == 4
    for k=1:N(1)
        for l=1:N(2)
            for m=1:N(4)
                aa=csi_abs(k,l,x,m)+csi_abs(k,l,x+distance1,m)+0.3*csi_abs(k,l,x+distance2,m);
                [a b]=max(aa);
                csi_sc(k,l,:,m)=circshift(csi(k,l,:,m),[0 0 -b+distance3 0]);
            end
        end
    end
end
end

