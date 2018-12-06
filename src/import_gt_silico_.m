%--------------------------------------------------------------------------
function import_gt_silico_(vcFile_mat)
    if matchFileExt_(vcFile_mat, '.mat')
        % catalin's raster function
        S = load(vcFile_mat);
        vnSpk = cellfun(@numel, S.a);
        viClu = int32(jrclust.utils.neCell2mat(arrayfun(@(n)n*ones(1, vnSpk(n)), 1:numel(vnSpk), 'UniformOutput', 0)));
        viTime = int32(jrclust.utils.neCell2mat(S.a) * 20); % Convert to sample # (saved in ms unit & sampling rate =20KHZ)
    elseif matchFileExt_(vcFile_mat, '.mda')
        % import Jeremy Magland format
        mr = readmda_(vcFile_mat)';
        viClu = int32(mr(:,3));
        viTime = int32(mr(:,2));
    end
    S_gt = makeStruct_(viClu, viTime);
    vcFile_gt = subsFileExt_(vcFile_mat, '_gt.mat');
    write_struct_(vcFile_gt, S_gt);
end %func
