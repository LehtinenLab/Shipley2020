% path = [];
% fdir = 'L:\1120-0115-962\200227\';
% fbase = '20200227_1120-0115-962_002_aCSF 2ul.min_1h';
% 
% fdir = 'Z:\Fred\adult_in_vivo_cx3cr1\512-725-3\20190920';
% fbase = '512-725-3_20190920_LPS';

function [Nchan, Nx, Ny, Nz, Nt] = GetDimensions(path,fdir,fbase)

    try %look for an sbx info file if it's already been converted to sbx
        info = pipe.io.sbxInfo(path);
        Nchan = info.nchan;
        Nx = info.sz(1);
        Ny = info.sz(2);
        Nz = info.otlevels; 
        Nt = floor(info.nframes / info.otlevels);

    catch % fill this in if this is your first time analyzing/converting data. 
        tifdir = strcat(fdir,filesep,fbase,'.tif.frames');
        tiflist = dir(strcat(tifdir,'\*.tif'));
        first_im = load_tiff(strcat(tiflist(1).folder,filesep,tiflist(1).name));

        tiflist_sort = sort_nat({tiflist.name});
        last_name = tiflist_sort{end};
        
        
        last_name_split = split(last_name,["_","."]);
        A = last_name_split{end-1};
        limits = regexp(A,'[\d.]+','match');

        Nchan = str2num(limits{1});
        Nx = size(first_im,1);
        Ny = size(first_im,2);
        Nz = str2num(limits{2});
        Nt = str2num(limits{3});

    end

end