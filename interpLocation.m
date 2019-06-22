% 
%   Copyright (C) 2016  Starsky Wong <sununs11@gmail.com>
% 
%   Note: The SIFT algorithm is patented in the United States and cannot be
%   used in commercial products without a license from the University of
%   British Columbia.  For more information, refer to the file LICENSE
%   that accompanied this distribution.

function [ ddata ] = interpLocation( dog_imgs, height, width, octv, layer, x, y, img_border, contr_thr )
% Function: Interpolates a scale-space extremum's location and scale
global base_sigma;
global layer_para;
global max_interp_steps;
sigma_rate_k=2^(1/layer_para);
i = 1;
while (i <= max_interp_steps)
    dD = deriv3D(layer,x,y);
    H = hessian3D(layer,x,y);
    [U,S,V] = svd(H);
    T=S;
    T(S~=0) = 1./S(S~=0);
    svd_inv_H = V * T' * U';
    x_hat = - svd_inv_H*dD;
    if( abs(x_hat(1)) < 0.5 && abs(x_hat(2)) < 0.5 && abs(x_hat(3)) < 0.5)
        break;
    end
    x = x + round(x_hat(1));
    y = y + round(x_hat(2));
    layer = layer + round(x_hat(3));
    if (layer < 2 || layer > layer_para+1 || x <= img_border || y <= img_border || x > height-img_border || y > width-img_border)
        ddata = [];
        return;
    end
    i = i+1;
end
if (i > max_interp_steps)
    ddata = [];
    return;
end
contr = dog_imgs(x,y,layer) + 0.5*dD'*x_hat;
if ( abs(contr) < contr_thr )
    ddata = [];
    return;
end
ddata.x = x;
ddata.y = y;
ddata.octv = octv;
ddata.layer = layer;
ddata.x_hat = x_hat;
ddata.scl_octv = base_sigma * power(sigma_rate_k,(layer+x_hat(3)-1)/layer_para);

function [ result ] = deriv3D(layer, x, y)
    dx = (dog_imgs(x+1,y,layer) - dog_imgs(x-1,y,layer))/2;
    dy = (dog_imgs(x,y+1,layer) - dog_imgs(x,y-1,layer))/2;
    ds = (dog_imgs(x,y,layer+1) - dog_imgs(x,y,layer-1))/2;
    result = [dx,dy,ds]';
end

function [ result ] = hessian3D(layer, x, y)
    center = dog_imgs(x,y,layer);
    dxx = dog_imgs(x+1,y,layer) + dog_imgs(x-1,y,layer) - 2*center;
    dyy = dog_imgs(x,y+1,layer) + dog_imgs(x,y-1,layer) - 2*center;
    dss = dog_imgs(x,y,layer+1) + dog_imgs(x,y,layer-1) - 2*center;

    dxy = (dog_imgs(x+1,y+1,layer)+dog_imgs(x-1,y-1,layer)-dog_imgs(x+1,y-1,layer)-dog_imgs(x-1,y+1,layer))/4;
    dxs = (dog_imgs(x+1,y,layer+1)+dog_imgs(x-1,y,layer-1)-dog_imgs(x+1,y,layer-1)-dog_imgs(x-1,y,layer+1))/4;
    dys = (dog_imgs(x,y+1,layer+1)+dog_imgs(x,y-1,layer-1)-dog_imgs(x,y-1,layer+1)-dog_imgs(x,y+1,layer-1))/4;

    result = [dxx,dxy,dxs;dxy,dyy,dys;dxs,dys,dss];
end

end


