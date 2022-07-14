function csiws = watersup_sim( csi, wat, beta)
% redisual water/lipid removal by L2 regularization
% csi: mrsi data with redisual water
% wat: water/lipid matrix
% delta: regularization parameter
beta = 10^(-beta);
NN = size(csi);
LL = length(NN);

%% single-voxel mrs
if LL == 1 || LL == 2
    N = size(csi);
    Water = wat;
    Water_inv = inv(eye(N(1))+beta*(Water*Water'));
    csiws = Water_inv*csi;
end

%% one-slice MRSI
if LL == 3
    csi = permute(csi,[3 2 1]);
    N = size(csi);
    Water = wat;
    Water_inv = inv(eye(N(1))+beta*(Water*Water'));
    csi = reshape(csi,N(1),N(2)*N(3));
    csi = Water_inv*csi;
    csi = reshape(csi,N(1),N(2),N(3));
    csi = permute(csi,[3 2 1]);
    csiws = csi;
end

%% multiple-slice MRSI
if LL == 4
    Water = wat;
    Water_inv = inv( eye(NN(3)) + beta * (Water * Water'));
    for ii = 1:NN(4)
        csi1 = permute(csi(:,:,:,ii),[3 2 1 4]);
        csi1 = csi(:,:,:,ii);
        N = size(csi1);
        csi1 = reshape(csi1,N(1),N(2)*N(3));
        csi1 = Water_inv * csi1;
        csi1 = reshape(csi1,N(1),N(2),N(3));
        csi1 = permute(csi1,[3 2 1]);
        csi(:,:,:,ii) = csi1;
    end
    csiws = csi;
end
end

