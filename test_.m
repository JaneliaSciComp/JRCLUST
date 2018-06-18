%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and added test ouput
function S_out = test_(vcFunc, cell_Input, nOutput, fVerbose, fDeleteEmpty)
    % S_out = test_(vcFunc, {input1, input2, ...}, nOutput)

    if nargin<2, cell_Input = {}; end
    if nargin<3, nOutput = []; end
    if nargin<4, fVerbose = ''; end
    if nargin<5, fDeleteEmpty = ''; end

    if isempty(fDeleteEmpty), fDeleteEmpty = 1; end
    if isempty(nOutput), nOutput = 1; end
    if ~iscell(cell_Input), cell_Input = {cell_Input}; end
    if isempty(fVerbose), fVerbose = 1; end
    if fDeleteEmpty, delete_empty_files_(); end

    try
        switch nOutput
            case 0
            feval(vcFunc, cell_Input{:});
            S_out = [];
            case 1
            [out1] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1);
            case 2
            [out1, out2] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1, out2);
            case 3
            [out1, out2, out3] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1, out2, out3);
            case 4
            [out1, out2, out3, out4] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1, out2, out3, out4);
        end %switch
        if fVerbose
            if nOutput>=1, fprintf('[%s: out1]\n', vcFunc); disp(S_out.out1); end
            if nOutput>=2, fprintf('[%s: out2]\n', vcFunc); disp(S_out.out2); end
            if nOutput>=3, fprintf('[%s: out3]\n', vcFunc); disp(S_out.out3); end
            if nOutput>=4, fprintf('[%s: out4]\n', vcFunc); disp(S_out.out4); end
        end
    catch
        disperr_();
        S_out = [];
    end
end %func
