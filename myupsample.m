function [output] = myupsample(data_0,rate)
%UNTITLED5 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[W,L,D]=size(data_0);
data=zeros(W+1,L+1,D);
data(1:W,1:L,:)=data_0;%�������
output=zeros(W*rate(1),L*rate(2),D);
step_W=(W-1)/(W*rate(1)-1);%<1
step_L=(L-1)/(L*rate(2)-1);
for row=0:W*rate(1)-1
    m=floor(row*step_W);%ÿһ���й̶���m��m+1
    for col=0:L*rate(2)-1
                n=floor(col*step_L);
                output(row+1,col+1,:)=(m+1-row*step_W)*(n+1-col*step_L)*data(m+1,n+1,:)+...
                            (m+1-row*step_W)*(col*step_L-n)*data(m+1,n+2,:)+...
                            (row*step_W-m)*(n+1-col*step_L)*data(m+2,n+1,:)+...
                            (row*step_W-m)*(col*step_L-n)*data(m+2,n+2,:);
    end
end
end


