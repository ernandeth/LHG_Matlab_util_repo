function data_grid = spiralGrid(ks,data,kx,ky,kz)
    data_grid = griddata(ks(:,2), ks(:,1), ks(:,3), data, kx, ky, kz);
end

