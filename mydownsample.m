function [output] = mydownsample(data_0,rate)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
[W,L,D]=size(data_0);
data=zeros(W+1,L+1,D);%补零操作
data(1:W,1:L,:)=data_0;
output=zeros(ceil(W/rate(1)),ceil(W/rate(2)),D);
step_W=(W-1)/(ceil(W/rate(1))-1);
step_L=(L-1)/(ceil(L/rate(2))-1);
for i=0:ceil(W/rate(1))-1
    m=floor(i*step_W);%每一行 的m和m+1的位置是固定的
    for j=0:ceil(L/rate(1))-1
        n=floor(j*step_L);
    output(i+1,j+1,:)=((m+1-i*step_W)*(n+1-j*step_L)*data(m+1,n+1,:)+...
                                (m+1-i*step_W)*(j*step_L-n)*data(m+1,n+2,:)+...
                                (i*step_W-m)*(n+1-j*step_L)*data(m+2,n+1,:)+...
                                (i*step_W-m)*(j*step_L-n)*data(m+2,n+2,:));
    end
end
end

