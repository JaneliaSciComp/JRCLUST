function vals = madScore(vals)
    %MADSCORE Transform vals into MAD units
    vals = bsxfun(@minus, vals, median(vals));
    madsc = median(abs(vals));

    vals = bsxfun(@rdivide, vals, madsc);
end
