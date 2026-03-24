function csi_bs= mrs_baselinecorre(csi, ws)
% for mrsi baseline correction
% csi: mrsi data before baseline correction
% ws: window size for the baseline correction

csi = abs(csi);
N = size(csi);
csi_bs = zeros(N);

if length(N) == 3
    csi = permute(csi,[3 2 1]);
    x=1:N(3); x=x';
    csi = reshape(csi,N(3),N(2)*N(1));
    csi = msbackadj(x,csi,'WINDOWSIZE',ws);
    csi = reshape(csi,N(3),N(2),N(1));
    csi_bs = permute(csi,[3 2 1]);
    csi_bs = abs(csi_bs);
end
if length(N) == 4
    csi = permute(csi,[3 2 1 4]);
    x=1:N(3); x=x';
    for ll = 1:N(4)
    csi1 = reshape(csi(:,:,:,ll),N(3),N(2)*N(1));
    csi1 = msbackadj(x,csi1);
    csi(:,:,:,ll) = reshape(csi1,N(3),N(2),N(1));
    end
    csi_bs = permute(csi,[3 2 1 4]);
    csi_bs = abs(csi_bs);
end
end
