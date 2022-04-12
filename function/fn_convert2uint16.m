function mat = fn_convert2uint16(mat)
    matMax = max(mat(:));
    if matMax <= 1
        mat = uint16(matMax * 65535);
    elseif matMax <= 255
        mat = uint16(matMax / 255 * 65535);
    else
        mat = uint16(mat);
    end
end