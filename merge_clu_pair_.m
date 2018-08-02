%--------------------------------------------------------------------------
function Sclu = merge_clu_pair_(Sclu, iClu1, iClu2)
    % if iClu1>iClu2, [iClu1, iClu2] = swap(iClu1, iClu2); end

    % update vnSpk_clu, viClu, clusterSites. move iClu2 to iClu1
    n1 = Sclu.vnSpk_clu(iClu1);
    n2 = Sclu.vnSpk_clu(iClu2);
    Sclu.vnSpk_clu(iClu1) = n1 + n2;
    Sclu.vnSpk_clu(iClu2) = 0;
    Sclu.viClu(Sclu.viClu == iClu2) = iClu1;
    Sclu.cviSpk_clu{iClu1} = find(Sclu.viClu == iClu1);
    Sclu.cviSpk_clu{iClu2} = [];
    try
        Sclu.csNote_clu{iClu1} = '';
        Sclu.csNote_clu{iClu2} = '';
    catch
    end
end %func
