fdir = 'C:\Users\LehtinenLab\Dropbox\Fred\Shipley2019\Figures\LPS supplement\LPS+';
fname_raw = '571 LVCP-big 2 GFP PECAM-568 ZO1-647 20x-2.tif';
fname_manual = 'ROI12_manual.tiff';

%% load manual mask
mask_manual = load_tiff(strcat(fdir,filesep,fname_manual));
mask_manual = logical(mask_manual);
%% load image and 
raw_im = load_tiff(strcat(fdir,filesep,fname_raw));
r_raw = mat2gray(raw_im(:,:,2));
r_lo = prctile(r_raw(:),20);
r_hi = prctile(r_raw(:),90);
r_adj = mat2gray(r_raw, [r_lo,r_hi]);

%get vesselness
r_v = vesselness2D_fbs(r_adj,[8:15],[1;1],1,true);
%binarize to find vessels
r_v_b = imbinarize(r_v);
%take boundary and dilate to get perivascular space
r_v_b_e = boundarymask(r_v_b);
%     peri = imdilate(r_v_b_e, strel('disk',w));
peri = r_v_b_e;
% find vessel space
mask_auto = r_v_b & ~peri;

%%
similarity = dice(mask_manual,mask_auto);
figure

subplot(2,2,1);
imshow(r_adj);
title('image of vasculature');

subplot(2,2,2);
imshowpair(r_adj,mask_manual);
title('Manual mask overlay');

subplot(2,2,3);
imshowpair(r_adj,mask_auto);
title('Automatic segmentation overlay');

subplot(2,2,4);
imshowpair(mask_manual, mask_auto);
title(['Dice Index = ' num2str(similarity)]);



