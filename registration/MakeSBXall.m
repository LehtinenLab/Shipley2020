function zproj_mean = MakeSBXall(path,shiftpath, varargin)
    
    p = inputParser;
    addOptional(p,'edges',[0,0,0,0]);
    addOptional(p,'Nt',[]); %how many volumes to register (if empty, do all)
    addOptional(p,'pmt',[]); %leave empty to run all channels.
    addOptional(p,'lineshift',0);
    addOptional(p,'optotune','true');
    addOptional(p,'refchannel',2);
    
    info_sbx = pipe.io.sbxInfo(path);
    Nz = size(info_sbx.otwave,2);    
    addOptional(p,'proj_range',round(0.25*Nz):round(0.75*Nz));
    
    parse(p,varargin{:});
    p = p.Results;
    

    
    fdir = fileparts(path);
    
    Nx = info_sbx.sz(1) - p.edges(3) - p.edges(4);
    Ny = info_sbx.sz(2) - p.edges(1) - p.edges(2);

    
    load(shiftpath,'-mat');
    
    ZS_total = ZS+ZS_chunk;
    ZS_total = ZS_total - median(ZS_total);

    RS_total = RS+RS_chunk;
    RS_total = RS_total - median(RS_total);

    CS_total = CS+CS_chunk;
    CS_total = CS_total - median(CS_total);
    
    if isempty(p.Nt)
        p.Nt = size(CS,2);
    elseif p.Nt > size(CS,2)
        fprintf('Trying to save more volumes than have been registered');
        return
    end
    
    if isempty(p.pmt)
        Nc = info_sbx.nchan;
        if Nc == 2
            p.pmt = -1;
        else
            p.pmt = 1;
        end
    elseif p.pmt == -1
        Nc = 2;
    elseif (p.pmt == 1 || p.pmt == 2)
        Nc = 1;
    else
        fprintf('please enter valid pmt');
        return
    end
    
    
    if ~exist('tforms_optotune_full')
        tforms_optotune_full = affine2d(eye(3));
    end
    
    w = waitbar(0,'doing zproj');

    tic
    
    %work one volume at a time, do the zproj registration
    zproj_raw = zeros(Nc,Nx,Ny,p.Nt);
    
    for i = 1:p.Nt
        %read the volume
        raw_vol = pipe.imread(path,Nz*(i-1)+1, Nz,p.pmt,[]);
        raw_vol = reshape(raw_vol,Nc,Nx,Ny,Nz);
        %apply line shift
        raw_vol(:,1:2:end,:,:) = circshift(raw_vol(:,1:2:end,:,:),p.lineshift,3);

        %crop it based on edges
        raw_vol = raw_vol(:,p.edges(3)+1:end-p.edges(4),p.edges(1)+1:end-p.edges(2),:);
        %if using optitune correction, perform correction
        if p.optotune
            parfor j = 1:Nz
                for c = 1:Nc
                    slice = squeeze(raw_vol(c,:,:,j));
                    warpvol(c,:,:,j) = imwarp(slice,tforms_optotune_full(j),'OutputView',imref2d(size(slice)));
                end
            end
        else %if not using optitune correction, just skip this step and go straight to applying DFT shifts
            warpvol = raw_vol;
        end
        
        reg_vol = zeros(size(warpvol));
        for c = 1:Nc
            %first do XY shifts 
            A = cell(Nz,1);
            parfor j = 1:Nz
                R = RS_total(j,i);
                C = CS_total(j,i);
                slice = squeeze(warpvol(c,:,:,j));
                A{j} = imtranslate(slice,[C,R]);
            end
            
            %reassemble parallel shifts
            for j = 1:Nz
                reg_vol(c,:,:,j) = A{j};
            end
            
        end
        %then do z shift
        Z = round(ZS_total(i));
        
        reg_vol = circshift(reg_vol,Z,4);
        if Z>0
            reg_vol(:,:,:,1:Z) = 0;
        elseif Z <0
            reg_vol(:,:,:,end+Z:end) = 0;
        end
        
        zproj_raw(:,:,:,i) = mean(reg_vol(:,:,:,p.proj_range),4);
        
    waitbar(i/p.Nt);
    end
    
    [zproj_mean,R_zproj,C_zproj] = zproj_reg(1, p.Nt, p.pmt, p.proj_range, 'refchan',p.refchannel,'zproj_raw',zproj_raw);
    
    R_zproj = transpose(repmat(R_zproj,1,Nz));
    C_zproj = transpose(repmat(C_zproj,1,Nz));
    
    RS_total = RS_total + R_zproj;
    CS_total = CS_total + C_zproj;
    
    % generate the "info" file again fro the .sbxz file
    info = SpoofSBXinfo3D(Nx,Ny,Nz,p.Nt,Nc);
    
    %open regwriter object for writing z-interpolated figure
    rw = pipe.io.RegWriter(path,info,'.sbxall',true);
    w = waitbar(0,'writing sbxall');
    for i = 1:p.Nt
        %read the volume
        raw_vol = pipe.imread(path,Nz*(i-1)+1, Nz,p.pmt,[]);
        raw_vol = reshape(raw_vol,Nc,Nx,Ny,Nz);
        %apply line shift
        raw_vol(:,1:2:end,:,:) = circshift(raw_vol(:,1:2:end,:,:),p.lineshift,3);

        %crop it based on edges
        raw_vol = raw_vol(:,p.edges(3)+1:end-p.edges(4),p.edges(1)+1:end-p.edges(2),:);
        %if using optitune correction, perform correction
        if p.optotune
            parfor j = 1:Nz
                for c = 1:Nc
                    slice = squeeze(raw_vol(c,:,:,j));
                    warpvol(c,:,:,j) = imwarp(slice,tforms_optotune_full(j),'OutputView',imref2d(size(slice)));
                end
            end
        else %if not using optitune correction, just skip this step and go straight to applying DFT shifts
            warpvol = raw_vol;
        end
        
        reg_vol = zeros(size(warpvol));
        for c = 1:Nc
            %first do XY shifts 
            A = cell(Nz,1);
            parfor j = 1:Nz
                R = RS_total(j,i);
                C = CS_total(j,i);
                slice = squeeze(warpvol(c,:,:,j));
                A{j} = imtranslate(slice,[C,R]);
            end
            
            %reassemble parallel shifts
            for j = 1:Nz
                reg_vol(c,:,:,j) = A{j};
            end
            
        end
        %then do z shift
        Z = round(ZS_total(i));
        
        reg_vol = circshift(reg_vol,Z,4);
        if Z>0
            reg_vol(:,:,:,1:Z) = 0;
        elseif Z <0
            reg_vol(:,:,:,end+Z:end) = 0;
        end
        
        %crop the registered image based on edge-sizes
        reg_vol_masked = zeros(size(reg_vol,1),size(reg_vol,2)+p.edges(3)+p.edges(4),size(reg_vol,3)+p.edges(1)+p.edges(2),size(reg_vol,4));
        reg_vol_masked(:,p.edges(3)+1:p.edges(3)+size(reg_vol,2),p.edges(1)+1:p.edges(1)+size(reg_vol,3),:) = reg_vol;
        reg_vol_masked = uint16(reg_vol_masked);

        %write volume to .sbxz file
        rw.write(reg_vol_masked);
        waitbar(i/p.Nt);

        %display an update every 10 volumes written
        if mod(i,10) == 0
            t = toc;
            fprintf('written %d volumes to .sbxall in %4.2f seconds\n',i,t);
            tic
        end
    end
    close(w);
end