function out =  reorder_slice_coord(vox)
    switch vox
        case 1
            out = 3;
        case 2
            out = 2;
        case 3
            out = 1;
    end
end