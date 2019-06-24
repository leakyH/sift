function [descrs,locs] = sift_mark(input)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
%input=im2double(imread('lena.bmp'));
%input=rgb2gray(im2double(imread('Iris.tif')));
%only gray images
%input=imrotate(imresize(input,0.8),60,'bilinear','crop');
setglobal();
global base_sigma;
global layer_para;
scale_para=5;



%% 构建Gaussian_pyr DOG_pyr
sigma_rate_k=2^(1/layer_para);
sigma_matrix=zeros(scale_para,layer_para+3);
for octave=1:scale_para
    for layer=1:layer_para+3
        sigma_matrix(octave,layer)=base_sigma.*2^(octave-1)*sigma_rate_k^(layer-1);
    end
end
%sigma矩阵构建完成
image_bundle=cell(1,scale_para);
image_bundle{1}=imresize(input,2);
Gaussian_pyr=cell(1,scale_para);
DOG_pyr=cell(1,scale_para);

for octave=1:scale_para
    image_scale=image_bundle{octave};
    gauss_bundle=zeros([size(image_scale),layer_para+3]);
    gauss_bundle(:,:,1)=image_scale;
    for i=2:layer_para+3
        gauss_bundle(:,:,i)=imgaussfilt(image_scale,sigma_matrix(octave,i));%高斯滤波
    end
    DOG_pyr{octave}=gauss_bundle(:,:,1:layer_para+2)-gauss_bundle(:,:,2:layer_para+3);%差分
    Gaussian_pyr{octave}=gauss_bundle;
    image_bundle{octave+1}=imresize(gauss_bundle(:,:,layer_para+1),0.5);%根据上一个octave的倒数第三张图下采样
end
%DOG构建完成


%% 关键点的精确定位
% 忽略图像边缘的点
img_border = 5;
% 特征点对比度阈值
contr_thr = 0.02;
% high threshold on feature ratio of principal curvatures
curv_thr = 5;
prelim_contr_thr = 0.5*contr_thr;
global extreme_array
extreme_array = struct('x',0,'y',0,'octv',0,'layer',0,'x_hat',[0,0,0],'scl_octv',0);
extreme_index = 1;
for i = 1:scale_para
    dog_imgs = DOG_pyr{i};
    [height, width] = size(dog_imgs(:,:,1));
    % find extrema in middle intvls
    for j = 2:layer_para+1
        dog_img = dog_imgs(:,:,j);
        for x = img_border+1:height-img_border
            for y = img_border+1:width-img_border
                % preliminary check on contrast
                if(abs(dog_img(x,y)) > prelim_contr_thr)%如果最大/最小值小于阈值一半，那么对比度一定不到阈值
                    % check 26 neighboring pixels
                        value = dog_imgs(x,y,j);
                        block = dog_imgs(x-1:x+1,y-1:y+1,j-1:j+1);
                        if ( value > 0 && value == max(block(:)) )
                            flag_MAXorMIN = 1;
                        elseif ( value == min(block(:)) )
                            flag_MAXorMIN = 1;
                        else
                            flag_MAXorMIN = 0;
                        end
                    if(flag_MAXorMIN)
                        extreme = interpLocation(dog_imgs,height,width,i,j,x,y,img_border,contr_thr);
                        if(~isempty(extreme))
                            if(~isEdgeLike(dog_img,extreme.x,extreme.y,curv_thr))
                                 extreme_array(extreme_index) = extreme;
                                 extreme_index = extreme_index + 1;
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Orientation Assignment
% number of detected points
n = size(extreme_array,2);
% determines gaussian sigma for orientation assignment
ori_sig_factr = 1.5;
% number of bins in histogram
ori_hist_bins = 36;
% array of feature
global features
features = struct('extreme_index',0,'x',0,'y',0,'scl',0,'ori',0,'descr',[],'weight',0);
feat_index = 1;
for i = 1:n
    extreme = extreme_array(i);
    ori_sigma = ori_sig_factr * extreme.scl_octv;
    % generate a histogram for the gradient distribution around a keypoint
    hist = oriHist(Gaussian_pyr{extreme.octv}(:,:,extreme.layer),extreme.x,extreme.y,ori_hist_bins,round(3*ori_sigma),ori_sigma);
    for j = 1:2
        smoothOriHist(hist,ori_hist_bins);
    end
    % generate feature from extreme and orientation hist peak
    % add orientations greater than or equal to 80% of the largest orientation magnitude
    feat_index = addOriFeatures(i,feat_index,extreme,hist);
end
%% 描述子
% number of features
n = size(features,2);
% width of 2d array of orientation histograms
descr_hist_d = 4;
% bins per orientation histogram
descr_hist_obins = 8;
% threshold on magnitude of elements of descriptor vector
descr_mag_thr = 0.2;
descr_length = descr_hist_d*descr_hist_d*descr_hist_obins;
local_features = features;
local_ddata_array = extreme_array;
for feat_index = 1:n
    feat = local_features(feat_index);
    extreme = local_ddata_array(feat.extreme_index);
    gauss_img = Gaussian_pyr{extreme.octv}(:,:,extreme.layer);
% computes the 2D array of orientation histograms that form the feature descriptor
    hist_width = 3*extreme.scl_octv;
    radius = round( hist_width * (descr_hist_d + 1) * sqrt(2) / 2 );
    feat_ori = feat.ori;
    ddata_x = extreme.x;
    ddata_y = extreme.y;
    hist = zeros(1,descr_length);
    for i = -radius:radius
        for j = -radius:radius
            j_rot = j*cos(feat_ori) - i*sin(feat_ori);
            i_rot = j*sin(feat_ori) + i*cos(feat_ori);
            r_bin = i_rot/hist_width + descr_hist_d/2 - 0.5;
            c_bin = j_rot/hist_width + descr_hist_d/2 - 0.5;
            if (r_bin > -1 && r_bin < descr_hist_d && c_bin > -1 && c_bin < descr_hist_d)
                mag_ori = calcGrad(gauss_img,ddata_x+i,ddata_y+j);
                if (mag_ori(1) ~= -1)
                    ori = mag_ori(2);
                    ori = ori - feat_ori;
                    ori=mod(ori,2*pi);
                    o_bin = ori * descr_hist_obins / (2*pi);
                    w = exp( -(j_rot*j_rot+i_rot*i_rot) / (2*(0.5*descr_hist_d*hist_width)^2) );
                    hist = interpHistEntry(hist,r_bin,c_bin,o_bin,mag_ori(1)*w,descr_hist_d,descr_hist_obins);
                end
            end
        end
    end
    local_features(feat_index) = hist2Descr(feat,hist,descr_mag_thr);
end
% sort the descriptors by descending scale order
features_scl = [local_features.scl];
[~,features_order] = sort(features_scl,'descend');
% return descriptors and locations
descrs = zeros(n,descr_length);
weights=zeros(n,1);
locs = zeros(n,2);
for i = 1:n
    weights(i,:) = local_features(features_order(i)).weight;
    descrs(i,:) = local_features(features_order(i)).descr;
    locs(i,2) = local_features(features_order(i)).x;
    locs(i,1) = local_features(features_order(i)).y;
end


%% 在图上标出
figure;
imshow(input);
viscircles(locs,weights*10,'LineWidth',1);
hold on
%annotation('arrow',round(locs(1,:)),round(locs(1,:)+[weights(1).*cos(features(1).ori),weights(1).*sin(features(1).ori)]));

quiver(locs(:,1),locs(:,2),weights.*cos([features.ori]'),weights.*sin([features.ori]'))
end

