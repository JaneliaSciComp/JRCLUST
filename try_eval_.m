%--------------------------------------------------------------------------
function try_eval_(vcEval1, fVerbose)
    if nargin<2, fVerbose = 1; end
    try
        eval(vcEval1);
        fprintf('\t%s\n', vcEval1);
    catch
        if fVerbose
            fprintf(2, '\tError evaluating ''%s''\n', vcEval1);
        end
    end
end %func
