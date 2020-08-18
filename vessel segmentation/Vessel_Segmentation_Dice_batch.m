fdir = 'C:\Users\LehtinenLab\Dropbox\Fred\Shipley2019\Figures\LPS supplement\dice overlap';
fbase = dir(strcat(fdir,filesep,'*.png'));

%%
for i = 1:numel(fbase)
    [~, tempbase{i}, ~] = fileparts(fbase(i).name);
    fname_raw = strcat(tempbase{i},'_raw.tif');
    fname_manual = strcat(tempbase{i},'_manual.tif');
    
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
    dice_similarity(i) = dice(mask_manual,mask_auto);
    bf_similarity(i) = bfscore(mask_manual,mask_auto);
end

T = table(tempbase',dice_similarity',bf_similarity','VariableNames',{'filename','dice_index','bf_score'});
writetable(T,strcat(fdir,filesep,'similarity_score.csv'));
fclose all;