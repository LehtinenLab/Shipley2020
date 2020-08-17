function out = zproj(path, k, Nt, pmt, otrange, varargin)

info = pipe.io.sbxInfo(path);
%if not specified, average first thirty volumes
if nargin < 2, Nt = floor(info.nframes / info.otlevels) - 1; end
%if not specified, load both PMTs
if nargin < 3, pmt = -1; end

p = inputParser;
addOptional(p, 'proj_type', 'mean');
addOptional(p, 'write_tiff', false); %check whether to save the tproj as a tiff, or just return it
addOptional(p, 'mtype', '.sbx'); %movie type suffix (other options are .sbxz or .sbxreg)
addOptional(p, 'tiff_suffix', '_zproj'); %name appended to tproj tiff file
addOptional(p, 'reg_suffix', []); %add a registration path
addOptional(p, 'registration', true); %do the registration
addOptional(p, 'line_shift', 0);
parse(p, varargin{:});
p = p.Results;
    
if (pmt == -1) & (info.nchan == 2)
    Nc = 2;
else
    Nc = 1;
end

Nx = info.width;
Ny = info.height;
Nz = info.otlevels;

[base,name,~] = fileparts(path);
movpath = fullfile(base,strcat(name,p.mtype));

if p.registration
    regpath = fullfile(base,strcat(name,p.reg_suffix));
else
    regpath = [];
end
    
A = cell(k+Nt-1,1);

%load each volume (parallelized) and take projection
h = parfor_progressbar(Nt,'doing z projection...');

parfor i = k:k+Nt-1
    
    if p.registration
        vol = pipe.imread(movpath,(i-1)*Nz+1,Nz,pmt,[], ...
        'register',true,'registration_path',regpath,'mtype',p.mtype);
    else
        vol = pipe.imread(movpath,(i-1)*Nz+1,Nz,pmt,[],'mtype',p.mtype);
    end
    
    if Nc == 1
        vol = vol(:,:,otrange);
    else
        vol = vol(:,:,:,otrange);
    end
    
    vol(vol==0) = NaN;
    
    if strcmp(p.proj_type,'mean')
        if pmt == -1
            slice = squeeze(mean(vol,4,'omitnan'));
        else
            slice = squeeze(mean(vol,3,'omitnan'));
        end
    elseif strcmp(p.proj_type,'max')
        if pmt == -1            
            slice = squeeze(max(vol,[],4,'omitnan'));
        else         
            slice = squeeze(max(vol,[],3,'omitnan'));
        end
    elseif strcmp(p.proj_type,'median')
        if pmt == -1
            slice = squeeze(median(vol,4,'omitnan'));
        else
            slice = squeeze(median(vol,3,'omitnan'));
        end
    else
        disp('invalid projection type');
    end
    slice(isnan(slice)) = 0;
    A{i,:} = slice;
    h.iterate(1);
end
out = zeros(Nc,Nx,Ny,Nt);
out = squeeze(out);


%concatenate each projected volume
for j = k:k+Nt-1
    if ndims(A{end}) == 3
        out(:,:,:,j-k+1) = A{j,1};
    else
        out(:,:,j-k+1) = A{j,1};
    end
end


clear A;

%write the file
if p.write_tiff
    if ndims(out) == 4
        out = permute(out, [2,3,1,4]);
    end
    out = reshape(out,size(out,1),size(out,2),[]);
    [base,name,~] = fileparts(path);
    fn = strcat(base,'\',name,p.tiff_suffix);
    pipe.io.writeTiff(uint16(out),fn);
end

close(h);
end