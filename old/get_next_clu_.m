%--------------------------------------------------------------------------
function iClu_next = get_next_clu_(S_clu, iClu1)
    %  get a next cluster number reading from Sclu.csNote_clu "=Clu#"
    iClu_next = [];
    vcNote_clu1 = S_clu.csNote_clu{iClu1};
    if isempty(vcNote_clu1), return; end
    iStart = find(vcNote_clu1 == '=', 1, 'first');
    if isempty(iStart), return; end
    vcNote_clu1 = vcNote_clu1(iStart+1:end);
    iEnd = find(vcNote_clu1 == ',' | vcNote_clu1 == ';' | vcNote_clu1 == ' ', 1, 'first');
    if ~isempty(iEnd), vcNote_clu1 = vcNote_clu1(1:iEnd-1); end
    iClu_next = str2double(vcNote_clu1);
    if isnan(iClu_next), iClu_next = []; return; end
end %func
