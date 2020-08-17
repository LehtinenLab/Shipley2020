function [unwarped_chunk] = ApplyOptotuneWarp(chunk, tforms_optotune)

for j = 1:size(chunk,4)
    for i = 1:size(chunk,3)
        slice = squeeze(chunk(:,:,i,j));
        unwarped_chunk(:,:,i,j) = imwarp(slice,tforms_optotune(i),'OutputView',imref2d(size(slice)));   
    end
end