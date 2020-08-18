clear all;

fdir_list = {'C:\Users\LehtinenLab\Dropbox\Fred\Shipley2019\Figures\LPS supplement\LPS+';
    'C:\Users\LehtinenLab\Dropbox\Fred\Shipley2019\Figures\LPS supplement\LPS-'};

% fdir = 'C:\Users\LehtinenLab\Dropbox\Fred\Shipley2019\Figures\LPS supplement\LPS+';
allnames = struct;
for L = 1:numel(fdir_list)
    allnames_temp = dir(strcat(fdir_list{L},filesep,'*.tif'));
    if L == 1
        allnames = allnames_temp;
    else
        allnames = [allnames; allnames_temp];
    end
%     allnames = cat(1,allnames,dir(strcat(fdir_list{L},filesep,'*.tif')));
end
%     [~,condition,~] = fileparts(fdir);
%%
w = 8;
for i = 1:numel(allnames)
    path = strcat(allnames(i).folder, filesep, allnames(i).name);
    [~,condition,~] = fileparts(allnames(i).folder);
    im = load_tiff(path);
    r_raw = mat2gray(im(:,:,2));
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
    % find stromal, non-perivascular space
    stroma = ~(r_v_b | peri);

    % find vessel space
    vess = r_v_b & ~peri;

    %
%     for j = 1:w
%         mask_zero = peri;
%         mask_vess(:,:,j) = (imdilate(peri,strel('disk',j)) & ~imdilate(peri,strel('disk',j-1))) & vess;
%         mask_stroma(:,:,j) = (imdilate(peri,strel('disk',j)) & ~imdilate(peri,strel('disk',j-1))) & stroma;
%     end  
%     mask_all = cat(3,flip(mask_stroma,3),mask_zero,mask_vess);
% 
%     mask_i = reshape(mask_all,[],size(mask_all,3));
%     r_adj_i = r_adj(:);
%     fluor = (r_adj_i' * mask_i) ./ sum(mask_i,1);
%     plot(fluor);

    [fluor, mask_all] = MeasurePeriProfile(peri, vess, stroma, r_adj, w);
    [fluor_rot, ~] = MeasurePeriProfile(imrotate(peri,90), imrotate(vess,90), imrotate(stroma,90), r_adj,w);
    
    data(i).fname = allnames(i).name;
    data(i).condition = condition;
    data(i).r_adj = r_adj;
    data(i).peri = peri;
    data(i).stroma = stroma;
    data(i).vess = vess;
    data(i).masks = mask_all;
    data(i).fluor = fluor';
    data(i).fluor_rot = fluor_rot';
end

%%
X = [-w:w];

fluor_all = [data(:).fluor];
fluor_mean = mean(fluor_all,2);

fluor_rot_all = [data(:).fluor_rot];
fluor_rot_mean = mean(fluor_rot_all,2);

figure,
hold on;
for i = 1:size(fluor_all,2)
    plot(X,fluor_rot_all(:,i),'Color',[1,.5,.5],'LineWidth',0.5);
end
p1 = plot(X,fluor_rot_mean,'r','LineWidth',2);

for i = 1:size(fluor_all,2)
    plot(X,fluor_all(:,i),'Color',[0.5,0.5,0.5],'LineWidth',0.5);
end
p2 = plot(X,fluor_mean,'k','LineWidth',2);

hold off;
ylim([0,1]);
axis square;
legend([p2,p1],{'periluminal','rotated periluminal'});



