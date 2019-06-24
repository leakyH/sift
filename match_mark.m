function [] = match_mark(input1,input2)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
[des1,loc1]=sift_mark(input1);
[des2,loc2]=sift_mark(input2);
match_thr = 0.7;
% for each descriptor in the first image, select its match to second image.
des2t = des2';
n = size(des1,1);
matched = zeros(1,n);
for i = 1 : n
   dotprods = des1(i,:) * des2t;
    [values,index] = sort(acos(dotprods));
%    [values,index] = sort(dotprods,'descend');
   disp(values(1)./values(2));
   if (values(1) < match_thr * values(2))
      matched(i) = index(1);
   else
      matched(i) = 0;
   end
%       if (values(2) < match_thr * values(1))
%       matched(i) = index(1);
%    else
%       matched(i) = 0;
%    end
end
drawMatched(matched, input1, input2, loc1, loc2);

end

