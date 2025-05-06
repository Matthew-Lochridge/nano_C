% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
% 
% Function called by select_param() to parse nanostructure input for carbon allotrope and lattice parameters (nanoribbons and nanotubes).
% Inputs:
%   param = container for nanostructure parameters
%   nanostructure = input from main()
% Outputs:
%   updated param with lattice parameters (nanoribbon or nanotube)
%   allotrope = char vector defining the type of structure as 'graphene, 'agnr', 'zgnr', 'gnr', or 'cnt'

function [param, allotrope] = parse_input(param)
    nanostructure = param.nanostructure;
    % default output
    allotrope = nanostructure;

    if ~isempty(nanostructure)
        % define error message
        err = "Carbon nanostructure not properly specified: accepted formats are 'n-AGNR', 'n-ZGNR', '(n,m)-GNR' and '(n1,m1)@(n2,m2)@...-CNT' for integer n's and m's.";

        if strcmpi(nanostructure,'graphene')
            return;

        elseif strcmpi(nanostructure,'ppp')
            return;

        else
            decomp = split(nanostructure,'-');
            input_param = decomp{1};
            allotrope = decomp{2};
            if ~strcmpi(allotrope,'agnr') && ~strcmpi(allotrope,'zgnr') && ~strcmpi(allotrope,'gnr') && ~strcmpi(allotrope,'cnt')
                disp(err)
                return;
            end
        
            if strcmpi(allotrope,'agnr') || strcmpi(allotrope,'zgnr') % N_w-AGNR or N_w-ZGNR
                if ~isempty(str2num(input_param))
                    param.N_dw = str2num(input_param); % width (transverse in-plane) of a ribbon within the supercell in dimer lines
                    param.N_dl = 2; % infinite ribbons
                else
                    disp(err)
                    return;
                end

            elseif strcmpi(allotrope,'gnr') && strcmp(input_param(1),'(') && strcmp(input_param(end),')') % (N_dl,N_dw)-GNR
                input_param = input_param(2:end-1);
                param_decomp = split(input_param,',');
                if ~isempty(str2num(param_decomp{1}))
                    param.N_dl = str2num(param_decomp{1}); % length (axial) of a ribbon within the supercell in dimer lines
                else
                    disp(err)
                    return;
                end
                if ~isempty(str2num(param_decomp{2}))
                    param.N_dw = str2num(param_decomp{2}); % width (transverse in-plane) of a ribbon within the supercell in dimer lines
                else
                    disp(err)
                    return;
                end

            elseif strcmpi(allotrope,'cnt') && strcmp(input_param(1),'(') && strcmp(input_param(end),')') % (N_1(1),N_2(1))@(N_1(2),N_2(2))@...-CNT
                input_param = split(input_param,'@');
                param.N_1 = [];
                param.N_2 = [];
                for t = 1:length(input_param)
                    input_param_t = input_param{t};
                    if strcmp(input_param_t(1),'(') && strcmp(input_param_t(end),')')
                        input_param_t = input_param_t(2:end-1);
                        input_param_t_decomp = split(input_param_t,',');
                        if ~isempty(str2num(input_param_t_decomp{1}))
                            param.N_1(t) = str2num(input_param_t_decomp{1}); % number of one primitive graphene translation for a nanotube wrapping vector
                        else
                            disp(err)
                            return;
                        end
                        if ~isempty(str2num(input_param_t_decomp{2}))
                            param.N_2(t) = str2num(input_param_t_decomp{2}); % number of the other primitive graphene translation for a nanotube wrapping vector
                        else
                            disp(err)
                            return;
                        end
                    end
                end

            else
                disp(err)
                return;
            end
        end
    end
end