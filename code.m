main;%sift_mark(input);
input=im2double(imread('lena.bmp'));
a=round(rand(1,4)*256);%�ü�
input2=input(a(1):a(2)+256,a(3):a(4)+256);
match_mark(input,input2)
b=rand();%��С
input3=imresize(input,b);
match_mark(input,input3);
c=1./rand();%�Ŵ�
input4=imresize(input,c);
match_mark(input,input4);
d=rand().*360;%��ת
input5=imrotate(input,d,'bilinear','crop');
match_mark(input,input5);
input6=input;%��������
input6(rand(512,512)>0.98)=1;
input6(rand(512,512)<0.02)=0;
match_mark(input,input6);
input7=input+randn(512,512)*0.1;%��˹������
match_mark(input,input7);