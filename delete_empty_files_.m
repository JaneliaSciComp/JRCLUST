%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and test
function delete_empty_files_(vcDir)
    if nargin<1, vcDir=[]; end
    delete_files_(find_empty_files_(vcDir));
end % function
