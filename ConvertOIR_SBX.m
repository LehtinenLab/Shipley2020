function shift = ConvertOIR_SBX(mouse,date,run,fdir,fbase,Nx,Ny,Nz,Nt,Nchan,varargin)
    p = inputParser;
    addOptional(p,'lineshift',false);
    parse(p,varargin{:});
    p=p.Results;
%     EXAMPLE PARAMETERS:
%     mouse = '512-725-1';
%     date = '190919';
%     run = 1;
%     ftype = 'sbx';
%     server = 'Mummy';
%     fbase = '512-725-1_20190919_WAY';
% 
%     Nx = 512;
%     Ny = 512;
%     Nz = 81;
%     Nt = 1200;
%     Nchan = 2;
    
    p = inputParser;
    addOptional(p,'lineshift',false);
    parse(p,varargin{:});
    p=p.Results;

    % get different path names
    tifdir = strcat(fdir,filesep,fbase,'.tif.frames');
    savename = strcat(mouse,'_',date,'_',sprintf('%03d',run));
    path = strcat(fdir, filesep, savename, '.sbx');
    
    %check to see if tifs have already been converted to sbx. If so, end
    %code and move on
    if ~exist(tifdir) && exist(path)
        disp('Tiff directory does not exist. May have been deleted already');
        return
    end
    
    
    info = SpoofSBXinfo3D(Ny, Nx, Nz, Nt, Nchan); %generate info file
    save(strcat(fdir, filesep, savename,'.mat'),'info', '-v7.3'); %save info mat file
    
    %look for info txt file and copy it to main directory
    try
        path_txt = strcat(tifdir,filesep,fbase,'.txt');
        copyfile(path_txt,fdir);
    end
    
    %% determine line shift atomatically
    shift = [];
    if p.lineshift
        A = zeros(Nchan, Nx,Ny,Nz,10);
        for t = 1:10
            for z = 1:Nz
                for c = 1:Nchan
                    tempname = strcat(fbase,'_C',sprintf('%03d',c),'Z',sprintf('%03d',z),'T',sprintf('%03d',t),'.tif');
                    A(c,:,:,z,t) = load_tiff_nobar(strcat(tifdir,filesep,tempname));
                end
            end
        end
        A = reshape(A,Nchan,Nx,Ny,[]);
        
%         A = pipe.imread(path,1,10*Nz,-1,[]);
        R = -5:5;
        for i = 1:numel(R)
            temp = A;
            temp(:,1:2:end,:,:) = circshift(temp(:,1:2:end,:,:),R(i),3);
            temp = mean(temp,4);

            I(1,:,:) = rescale(temp(1,:,:));
            I(2,:,:) = rescale(temp(2,:,:));
            Im = squeeze(mean(I,1));
            Im = rescale(Im);
            B(:,:,i) = abs(diff(Im,1,1));
        end

        trace = squeeze(sum(sum(B,1),2));

%         figure,plot(R,trace);

        [~,pkidx] = min(trace);
        shift = R(pkidx);
    end
    %% 
    % convert all of the tif frames (exported from FluoViewer) into .sbx format
    % ONLY NEED TO DO THIS ONCE; delete the ".tif.frames" folder afterwards
    % sa        ve .sbx and .oir (should be redundant anyway)
    rw = pipe.io.RegWriter(path,info,'.sbx',true);

    H = parfor_progressbar(Nt, 'converting to sbx...');
%     for t = 1:800
    for t = 1:Nt
        vol = zeros(Nchan, Nx, Ny, Nz);
        parfor z = 1:Nz %loads a volume in parallel into memory
            for c = 1:Nchan
                tempname = strcat(fbase,'_C',sprintf('%03d',c),'Z',sprintf('%03d',z),'T',sprintf('%03d',t),'.tif');
%                 disp(tempname);
                vol(c,:,:,z) = load_tiff_nobar(strcat(tifdir,filesep,tempname));
            end
        end
        rw.write(uint16(vol));
        H.iterate(1);
    end
    rw.delete;
    close(H);

    %% deletes tif file folder
    rmdir(tifdir, 's');
end