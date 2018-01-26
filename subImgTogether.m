function [Img] = subImgTogether(subAperImg, ML, SH)

Img = zeros(SH.resolution);
for idx=1:length(ML.coor)
    Img(ML.coor(idx,1):ML.coor(idx,2), ML.coor(idx,3):ML.coor(idx,4)) =  subAperImg{idx};
end