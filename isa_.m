%--------------------------------------------------------------------------
% 8/6/17 JJJ: Tested and documented
function flag = isa_(vr, vcClass)
    % subtract mean for mr or tr
    try
        if isa(vr, 'gpuArray')
            flag = strcmpi(vcClass, classUnderlying(vr));
        else
            flag = isa(vr, vcClass);
        end
    catch
        disperr_();
        flag = 0;
    end
end % function
