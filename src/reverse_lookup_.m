%--------------------------------------------------------------------------
function [viA, vl] = reverse_lookup_(viB, viA2B)
    % viB must belong to viA2B
    % viB=[3 1 1 5 3 3 2]; viA2B=[1 3 5];

    viB = int32(viB);
    viA2B = int32(viA2B);
    [vl, viA] = ismember(int32(viB), int32(viA2B));

    % vl = ismember(viB, viA2B);
    % assert_(all(vl), 'reverse_lookup_: all viB must belong to viA2B');
    % viA = arrayfun(@(i)find(viA2B==i), viB(vl), 'UniformOutput', 1);
end %func
