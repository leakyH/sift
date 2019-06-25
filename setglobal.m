global base_sigma;
base_sigma=1.6;
global layer_para;
layer_para=3;%larger than 3
global max_interp_steps;
max_interp_steps=5;
% number of bins in histogram
global ori_hist_bins
ori_hist_bins = 36;

% orientation magnitude relative to max that results in new feature
global ori_peak_ratio
ori_peak_ratio = 0.8;
% 特征点对比度阈值
global contr_thr;
contr_thr = 0.03;
% 边缘筛选阈值
global curv_thr;
curv_thr = 5;d