%--------------------------------------------------------------------------
function [mr, vi_shuffle] = shuffle_static_(mr, dimm)
    % dimm = 1 or 2 (dimension to shuffle
    if nargin<2, dimm=1; end
    fStatic = 1;
    if fStatic
        s = RandStream('mt19937ar','Seed',0); %always same shuffle order
        switch dimm
            case 1
            vi_shuffle = randperm(s, size(mr,1));
            mr = mr(vi_shuffle, :);
            case 2
            vi_shuffle = randperm(s, size(mr,2));
            mr = mr(:, vi_shuffle);
        end
    else
        switch dimm
            case 1
            mr = mr(randperm(size(mr,1)), :);
            case 2
            mr = mr(:, randperm(size(mr,2)));
        end
    end
end %func
