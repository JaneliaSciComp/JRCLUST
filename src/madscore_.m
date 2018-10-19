%--------------------------------------------------------------------------
function mr = madscore_(mr)
    % maximum absolute difference transformation

    mr = bsxfun(@minus, mr, median(mr));
    vr = median(abs(mr));
    mr = bsxfun(@rdivide, mr, vr);
end %func
