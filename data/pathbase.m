function base = pathbase(server)
%BASEPATH Hard codes the base data directories for the Andermann lab based
%   on servers and host names.

    % Cache the hostname to accurately get the server
    persistent cached_hostname
    if ~isempty(cached_hostname)
        hn = cached_hostname;
    else
        [success, syshostname] = system('hostname');

        % some error checks
        assert(success == 0, 'Error running hostname');
        assert(~any(syshostname == '.'), 'Dots found in hostname: is it a fqdn?');

        hn = deblank(syshostname);
        cached_hostname = hn;
    end

    % Set base path depending on server
    if nargin < 1 || isempty(server) || strcmpi(hn, server)
        if strcmp(hn, 'Megatron')
            base = 'D:\twophoton_data\2photon\scan\';
        elseif strcmp(hn, 'Atlas')
            base = 'E:\twophoton_data\2photon\raw\';
        elseif strcmp(hn, 'BeastMode')
            base = 'S:\twophoton_data\2photon\scan\';
        elseif strcmp(hn, 'Sweetness')
            base = 'D:\2p_data\scan\';
        elseif strcmpi(hn, 'santiago')
            base = 'D:\2p_data\scan\';
        %following for IRONWOLF drives (fbs)
        elseif strcmpi(hn, 'Werewolf')
            base = 'E:\';
        elseif strcmpi(hn, 'Vampire')
            base = 'F:\';
        elseif strcmpi(hn, 'Frankenstein')
            base = 'H:\';
        elseif strcmpi(hn, 'Banshee')
            base = 'G:\';
        elseif strcmpi(hn, 'Mummy')
            base = 'J:\';
        elseif strcmpi(hn, 'Zombie')
            base = 'L:\';
        elseif strcmpi(hn,'DeepFreezeCX3CR1')
            base = 'Z:\Fred\adult_in_vivo_cx3cr1\';
        elseif strcmpi(hn,'Ahram')
            base = 'Z:\Ahram\2p_imaging_videos\';
        elseif strcmpi(hn,'Euclid')
            base = 'F:\';
        elseif strcmpi(hn, 'Pythagoras')
            base = 'G:\';
        else strcmpi(hn, 'Hippocrates');
            base = 'H:\';
        end
        
    else
        if strcmpi(server, 'santiago');
            base = '\\santiago\2p_data\scan\';
        elseif strcmpi(server, 'sweetness');
            base = '\\sweetness\2p_data\scan\';
        elseif strcmpi(server, 'megatron');
            base = '\\megatron\2photon\scan\';
        elseif strcmpi(server, 'storage') && strcmp(hn, 'Megatron');
            base = 'E:\scan\';
        elseif strcmpi(server, 'storage');
            base = '\\megatron\E\scan\';
        elseif strcmpi(server, 'anastasia');
            base = '\\anastasia\data\2p\';
        %following for IRONWOLF drives (fbs)
        elseif strcmpi(server, 'Werewolf')
            base = 'E:\';
        elseif strcmpi(server, 'Vampire')
            base = 'F:\';
        elseif strcmpi(server, 'Frankenstein')
            base = 'H:\';
        elseif strcmpi(server, 'Mummy')
            base = 'J:\';
        elseif strcmpi(server, 'Banshee')
            base = 'G:\';
        elseif strcmpi(server, 'Zombie')
            base = 'L:\';
        elseif strcmpi(server,'DeepFreezeCX3CR1')
            base = 'Z:\Fred\adult_in_vivo_cx3cr1\';
        elseif strcmpi(server,'Ahram')
            base = 'Z:\Ahram\2p_imaging_videos\';
        elseif strcmpi(server,'Euclid')
            base = 'F:\';
        elseif strcmpi(server, 'Pythagoras')
            base = 'G:\';
        else strcmpi(server, 'Hippocrates');
            base = 'H:\';
        end
        
    end
end

