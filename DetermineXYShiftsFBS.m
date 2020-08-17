function[RowShifts,ColumnShifts] = DetermineXYShiftsFBS(full_vol,...
    BlurFactor,KeepingFactor,ReferenceVolume)

%   DETERMINEXYSHITS: determine XY shifts per plane with a parfor loop
%   Credits: Alex Fratzl, edited by Fred Shipley
%
%   Inputs:
%     full_vol -- 4D matrix of uint16 or other, dim (x,y,z,t)
%     BlurFactor -- width of the gaussian filter (ex: 1)
%     KeepingFactor -- 0 < KeepingFactor < 1, % of FOV taken into account,
%       removes edges to determine shifts (ex: 0.95)
%     ReferenceVolume -- 4D matrix of uint, moving volume reference, 
%       dim (x,y,z,t/n) with n = nb avg frames per reference volume
%   Outputs:
%     RowShifts -- 2D matrix of doubles, dim (z,t)
%     ColumnShifts -- 2D matrix of doubles, dim (z,t)

    tStartDXYS = tic;

    Size = size(full_vol);
    S3 = size(full_vol, 3); S4 = size(full_vol, 4);
    Keep = KeepingFactor;
    red_vol = full_vol(ceil(Size(1)*(1-Keep)/2):...
        ceil(Size(1)*(1-(1-Keep)/2)),...
        ceil(Size(2)*(1-Keep)/2):...
        ceil(Size(2)*(1-(1-Keep)/2)),:,:);
    chunck = floor(S4/size(ReferenceVolume, 4));
    ReferenceVolume = ReferenceVolume(...
        ceil(Size(1)*(1-Keep)/2):...
        ceil(Size(1)*(1-(1-Keep)/2)),...
        ceil(Size(2)*(1-Keep)/2):...
        ceil(Size(2)*(1-(1-Keep)/2)),:,:);
    RowShifts = zeros(S3, S4);
    ColumnShifts = zeros(S3, S4);
    
    for t = 1:S4
        reft = ReferenceVolume(:,:,:,ceil(t/chunck));
        for i = 1:S3
            ref = reft(:,:,i);
            output = dftregistrationAlex(fft2(imgaussfilt(...
                ref,BlurFactor)),...
                fft2(imgaussfilt(red_vol(:,:,i,t),BlurFactor)),100);
            RowShifts(i,t) = output(1);
            ColumnShifts(i,t) = output(2);
        end
    end

    tEndDXYS = toc(tStartDXYS);
end