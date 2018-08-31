function [mrMiss, mrFp, vnCluGt, miCluMatch, S_score_clu] = clusterVerify(viCluGt, viTimeGt, viClu, viTime, jitter)
% very rapid

% fNorm_mode = 3; %normalize by the ground truth unit

if nargin<5, jitter = 12; end

viTime = double(viTime);
viTimeGt = double(viTimeGt);
viClu = int32(viClu);
viCluGt = int32(viCluGt);
viCluGt_unique = unique(viCluGt);
% viClu_unique = unique(viClu);
nCluGt = numel(viCluGt_unique);
% nClu = numel(viClu_unique);
nClu = max(viClu);
[mrMiss, mrFp, mrAccuracy] = deal(nan(nClu, nCluGt, 'single'));
[vnDetected, vnCluGt] = deal(zeros(nCluGt,1,'int32'));
cviTime = arrayfun(@(iClu)unique(int32(viTime(viClu==iClu)/jitter)), 1:nClu, 'UniformOutput', 0);
% cviTimeGt = arrayfun(@(iClu)gpuArray(viTimeGt(viCluGt==iClu)/jitter), viCluGt_unique, 'UniformOutput', 0);
[cviTimeGt, cviHit_gt, cviMiss_gt, cviHit_clu, cviMiss_clu, ...
    cviSpk_gt_hit, cviSpk_gt_miss, spikesByCluster_hit, spikesByCluster_miss] = ...
    deal(cell(1, nCluGt));
viTime0 = (int32(viTime/jitter));
fprintf('Validating cluster\n\t');
t1 = tic;
parfor iCluGt1=1:nCluGt
    viTimeGt1 = viTimeGt(viCluGt == viCluGt_unique(iCluGt1));
    rGT1 = int32(single(viTimeGt1) / jitter);
    rGT1 = unique(rGT1);
    cviTimeGt{iCluGt1} = rGT1;
    if isempty(rGT1), continue; end
    vnCluGt(iCluGt1) = numel(rGT1);

    % detection check
    vnDetected(iCluGt1) = count_overlap_(rGT1, viTime0);
    [vrMiss_, vrFp_, vrAccuracy_] = deal(zeros(nClu, 1, 'single'));
    for iClu=1:nClu
        rComp1 = cviTime{iClu};
        n3 = count_overlap_(rGT1, rComp1);
        n2 = numel(rComp1);
        n1 = numel(rGT1);
        vrMiss_(iClu) = (n1-n3)/n1;
        vrFp_(iClu) = (n2-n3)/n2;   
        vrAccuracy_(iClu) = n3 / (n1+n2-n3);
    end
    mrMiss(:,iCluGt1) = vrMiss_;
    mrFp(:,iCluGt1) = vrFp_;    
    mrAccuracy(:,iCluGt1) = vrAccuracy_;
    fprintf('.');
end
fprintf('\n\ttook %0.1fs.\n', toc(t1));
[mrMiss, mrFp, vnDetected, vnCluGt] = multifun_(@gather, mrMiss, mrFp, vnDetected, vnCluGt);
vrDetected = double(vnDetected) ./ double(vnCluGt);
mrScore = 1-mrFp-mrMiss;
[vrScore, viCluMatch] = max(mrScore, [], 1);
[~, miCluMatch] = sort(mrScore, 'descend');
viFp = sub2ind(size(mrMiss), viCluMatch, 1:numel(viCluMatch));
vrMiss = mrMiss(viFp);
vrFp = mrFp(viFp);
vrScore = 1-vrMiss-vrFp;
vrScore(vrScore<0)=0;
for iCluGt=1:nCluGt
    viSpk_gt1 = find(viCluGt == viCluGt_unique(iCluGt));
    viTime_gt1 = viTimeGt(viSpk_gt1);    
    iClu1 = viCluMatch(iCluGt);
    viSpk_clu1 = find(viClu==iClu1);
    viTime_clu1 = viTime(viSpk_clu1);
    [vlSpk_gt1, vlSpk_clu1, viiSpk_gt1, viiSpk_clu1] = time_match2_(viTime_gt1, viTime_clu1, jitter);    
    cviHit_gt{iCluGt} = viTime_gt1(viiSpk_gt1);
    cviMiss_gt{iCluGt} = viTime_gt1(~vlSpk_gt1);
    cviHit_clu{iCluGt} = viTime_clu1(viiSpk_clu1);
    cviMiss_clu{iCluGt} = viTime_clu1(~vlSpk_clu1);
    cviSpk_gt_hit{iCluGt} = viSpk_gt1(vlSpk_gt1);
    cviSpk_gt_miss{iCluGt} = viSpk_gt1(~vlSpk_gt1);
    spikesByCluster_hit{iCluGt} = viSpk_clu1(vlSpk_clu1);
    spikesByCluster_miss{iCluGt} = viSpk_clu1(~vlSpk_clu1);  
end
[vrAccuracy, viCluMatch_accuracy] = max(mrAccuracy);
S_score_clu = makeStruct(vrScore, vrMiss, vrFp, viCluMatch, cviHit_gt, ...
    cviMiss_gt, cviHit_clu, cviMiss_clu, cviSpk_gt_hit, cviSpk_gt_miss, ...
    spikesByCluster_hit, spikesByCluster_miss, vrAccuracy, viCluMatch_accuracy);
% viCluMatch1 = viClu_unique(viCluMatch);
% vrScore = 1-vrMiss-vrFp;
func1=@(x)quantile(x, [.25,.5,.75])*100;
fprintf('Validation summary\n');
fprintf('\tEvents detected (<%0.1f%%> %0.1f %0.1f %0.1f): %s\n', ...
    nanmean(vrDetected)*100, func1(vrDetected), sprintf('%0.1f ', vrDetected*100));
fprintf('\tFalse-positives (<%0.1f%%> %0.1f %0.1f %0.1f): %s\n', ...
    nanmean(vrFp)*100, func1(vrFp), sprintf('%0.1f ', vrFp*100));
fprintf('\tFalse-negatives (<%0.1f%%> %0.1f %0.1f %0.1f): %s\n', ...
    nanmean(vrMiss)*100, func1(vrMiss), sprintf('%0.1f ', vrMiss*100));
fprintf('\tAccuracy (<%0.1f%%> %0.1f %0.1f %0.1f): %s\n', ...
    nanmean(vrAccuracy)*100, func1(vrAccuracy), sprintf('%0.1f ', vrAccuracy*100));
fprintf('\tScore (<%0.1f%%> %0.1f %0.1f %0.1f): %s\n', ...
    nanmean(vrScore)*100, func1(vrScore), sprintf('%0.1f ', vrScore*100));
fprintf('\tCluster-size: %s\n', sprintf('%d, ', vnCluGt));
fprintf('\tMatching clu: %s\n', sprintf('%d, ', viCluMatch));
end % function


%--------------------------------------------------------------------------
function nOverlap = count_overlap_(viGt, viTest)
nGt = numel(viGt);
nOverlap = nGt - numel(setdiff(setdiff(setdiff(viGt, viTest), viTest+1), viTest-1));
end % function


%--------------------------------------------------------------------------
function [vlA, vlB, viA1, viB1] = time_match2_(viA, viB, jitter)
% A: ground truth, B: matching unit
if nargin<3, jitter=25; end

if jitter>0
    viA = (double(viA)/jitter);
    viB = (double(viB)/jitter);
end
viA = int32(viA);
viB = int32(viB);
vlA = false(size(viA));
vlB = false(size(viB));
for i1=-1:1
    for i2=-1:1
        vlA = vlA | ismember(viA+i1, viB+i2);
    end
end
if nargout==1, return; end

%viA_match = find(vlA);
viA1 = find(vlA);
viB1 = zeros(size(viA1));
viA11 = viA(viA1);
for iA=1:numel(viA1)
    [~, viB1(iA)] = min(abs(viB - viA11(iA)));
end
vlB(viB1) = 1;
end % function


%--------------------------------------------------------------------------
function [ S ] = makeStruct( varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name

S=[];
for i=1:nargin
    S = setfield(S, inputname(i), varargin{i});
end
end % function


%--------------------------------------------------------------------------
function varargout = multifun_(hFun, varargin)
% apply same function to the input, unary function only

if nargout ~= numel(varargin), error('n arg mismatch'); end
for i=1:nargout
    try
        varargout{i} = hFun(varargin{i});
    catch
        varargout{i} = varargin{i};
    end
end
end % function