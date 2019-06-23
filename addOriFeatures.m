% 
%   Copyright (C) 2016  Starsky Wong <sununs11@gmail.com>
% 
%   Note: The SIFT algorithm is patented in the United States and cannot be
%   used in commercial products without a license from the University of
%   British Columbia.  For more information, refer to the file LICENSE
%   that accompanied this distribution.

function [feat_index] = addOriFeatures(extreme_index,feat_index,ddata,hist)
% Function: Add good orientation for keypoints
global features;
global base_sigma;
global layer_para;
global ori_hist_bins
%omax = dominantOri(hist,ori_hist_bins);?????????????????
omax=max(hist);
% orientation magnitude relative to max that results in new feature
global ori_peak_ratio;

for i = 1:ori_hist_bins
    k=mod(i-2,ori_hist_bins)+1;
    r=mod(i,ori_hist_bins)+1;
    if ( hist(i) > hist(k) && hist(i) > hist(r) && hist(i) >= ori_peak_ratio*omax )
        bin = i + interp_hist_peak(hist(k),hist(i),hist(r))-1;%需要表达为[0-36)浮点数
        bin=mod(bin,36);
        accu_layer = ddata.layer + ddata.x_hat(3);
        features(feat_index).extreme_index = extreme_index;
        % first octave is double size
        features(feat_index).x = (ddata.x + ddata.x_hat(1))*2^(ddata.octv-2);
        features(feat_index).y = (ddata.y + ddata.x_hat(2))*2^(ddata.octv-2);
        features(feat_index).scl = base_sigma * power(2,ddata.octv-2 + (accu_layer-1)/layer_para);        
        features(feat_index).ori = bin/ori_hist_bins*2*pi - pi;
        features(feat_index).weight=hist(i);
        feat_index = feat_index + 1;
    end
end
end

function [position] = interp_hist_peak(p,c,r)
    position = 0.5*(p-r)/(p-2*c+r);
end