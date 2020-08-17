function tform_cum_ordered = MultiStackReg_Fiji_affine(mean_vol,fdir,Nz)

% fdir = 'H:\VF284\20180518\';
% fname = 'VF284_180518_001';
% cd([fdir,'chunks\']);
% chunk = load_tiff('chunk1.tif');
% Nz = size(mean_vol,3);
% chunk = reshape(chunk,size(chunk,1),size(chunk,2),2,Nz,[]);
% tproj = mean(chunk,5);
% mean_im = squeeze(mean_vol(:,:,1,:));

B = sum(reshape(mean_vol,[],Nz),1);

[~,I] = max(B);

Miji;
MIJ.createImage(mean_vol);
MIJ.setSlice(I);
MIJ.run('MultiStackReg', ['stack_1=[Import from Matlab] action_1=Align file_1=',fdir,'TransformationMatricesAffine.txt stack_2=None action_2=Ignore file_2=[] transformation=Affine save']);
MIJ.run('Close');
MIJ.exit;

%%
transforms = LoadTransforms([fdir,'TransformationMatricesAffine.txt']);

M = repmat([true; true; true; false; false; false],Nz-1,2);

movingPoints = transforms(M);
movingPoints = reshape(movingPoints,[],2);

%% get step-wise transformation (NOT cumulative)
idx = 1:3:size(movingPoints,1);
fixedPoints = transforms(4:6,:);
tform(1) = fitgeotrans(fixedPoints,fixedPoints,'affine');
for i = 1:length(idx)
    A = movingPoints(idx(i):idx(i)+2,:);

    tform(i+1) = fitgeotrans(A,fixedPoints,'affine');
end

tform_ordered = [fliplr(tform(1:I)),tform(I+1:end)];

%HAVE to invert the transformation for those coming BEFORE reference point
for i = 1:I
    tform_ordered(i) = invert(tform_ordered(i));
end

%% get cumulative transforms
M_cum = eye(3);
for i = 1:length(tform_ordered)
    M_cum = M_cum * tform_ordered(i).T;
    tform_cum_ordered(i) = affine2d(M_cum);
end

end
















