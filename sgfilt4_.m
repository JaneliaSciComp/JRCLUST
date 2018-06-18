%--------------------------------------------------------------------------
% 17/11/20 JJJ: created
% function vr1 = fft_align_mean_(mr)
% % take the median phase
% if ndims(mr)==3
%     mr1 = zeros(size(mr,1), size(mr,2), 'like', mr);
%     for i1=1:size(mr,2)
%         mr1(:,i1) = fft_align_mean_(squeeze(mr(:,i1,:)));
%     end
%     vr1 = mr1;
%     return;
% end
% mrF = fft(mr);
% vrAbs1 = mean(abs(mrF),2);
% vrAng1 = median(unwrap(angle(mrF)),2);
% vr1 = real(ifft(complex(vrAbs1.*cos(vrAng1), vrAbs1.*sin(vrAng1))));
% end %func


%--------------------------------------------------------------------------
function [via1, via2, via3, via4, vib1, vib2, vib3, vib4] = sgfilt4_(n1, fGpu)
    persistent n1_prev_ via1_ via2_ via3_ via4_ vib1_ vib2_ vib3_ vib4_
    if nargin<2, fGpu=0; end

    % Build filter coeff
    if isempty(n1_prev_), n1_prev_ = 0; end
    try a = size(via1_); catch, n1_prev_ =0; end
    if n1_prev_ ~= n1 %rebuild cache
        vi0 = int32(1:n1);
        vi0 = gpuArray_(vi0, fGpu);
        via4_ = vi0+4; via4_(end-3:end)=n1;   vib4_ = vi0-4; vib4_(1:4)=1;
        via3_ = vi0+3; via3_(end-2:end)=n1;   vib3_ = vi0-3; vib3_(1:3)=1;
        via2_ = vi0+2; via2_(end-1:end)=n1;   vib2_ = vi0-2; vib2_(1:2)=1;
        via1_ = vi0+1; via1_(end)=n1;         vib1_ = vi0-1; vib1_(1)=1;
        n1_prev_ = n1;
    end

    % Copy from cache
    via4 = via4_;   vib4 = vib4_;
    via3 = via3_;   vib3 = vib3_;
    via2 = via2_;   vib2 = vib2_;
    via1 = via1_;   vib1 = vib1_;
end %func
