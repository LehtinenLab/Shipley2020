%put in perivascular boundary, the stromal mask, the vessel mask the image to measure, and the width of
%dilation

function [fluor, mask_all] = MeasurePeriProfile(peri, vess, stroma, im, w)
    for j = 1:w
        mask_zero = peri;
        mask_vess(:,:,j) = (imdilate(peri,strel('disk',j)) & ~imdilate(peri,strel('disk',j-1))) & vess;
        mask_stroma(:,:,j) = (imdilate(peri,strel('disk',j)) & ~imdilate(peri,strel('disk',j-1))) & stroma;
    end  
    mask_all = cat(3,flip(mask_stroma,3),mask_zero,mask_vess);

    mask_i = reshape(mask_all,[],size(mask_all,3));
    im_i = im(:);
    fluor = (im_i' * mask_i) ./ sum(mask_i,1);
end