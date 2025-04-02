% Author: Matthew Lochridge
% Term Project for MSEN 5377 (Spring 2025)
% 
% Function called by select_param to parse structure for selected nanoribbon or nanotube parameters
% parameters.
% Input:
%   structure = name of selected nanostructure as a char vector,
%               e.g. '9-AGNR', '8-ZGNR', '(10,10)@(15,15)-CNT', etc.
%               If unspecified, 'graphene' is called by default.
% Outputs:
%   allotrope = char vector defining the type of structure as 'graphene, 'agnr', 'zgnr', 'gnr', or 'cnt'
%   N_dl, N_dw = length (axial) and width (transverse in-plane) of a ribbon within the supercell in dimer lines
%   N_1, N_2 = numbers of each primitive graphene translation for a nanotube wrapping vector

function [allotrope, N_dl, N_dw, N_1, N_2] = parse_structure(structure)
    % default outputs
    allotrope = structure;
    N_dl = 0;
    N_dw = 0;
    N_1 = 0;
    N_2 = 0;

    if ~isempty(structure)
        % define error message
        err = "Carbon nanostructure not properly specified: accepted formats are 'n-AGNR', 'n-ZGNR', '(n,m)-GNR' and '(n1,m1)@(n2,m2)@...-CNT' for integer n's and m's.";

        if strcmpi(allotrope,'graphene')
            return;

        else
            decomp = split(structure,'-');
            param = decomp{1};
            allotrope = decomp{2};
            if ~strcmpi(allotrope,'agnr') && ~strcmpi(allotrope,'zgnr') && ~strcmpi(allotrope,'gnr') && ~strcmpi(allotrope,'cnt')
                disp(err)
                return;
            end
        
            if strcmpi(allotrope,'agnr') || strcmpi(allotrope,'zgnr') % N_w-AGNR or N_w-ZGNR
                if ~isempty(str2num(param))
                    N_dw = str2num(param);
                    N_dl = 2; % infinite ribbons
                else
                    disp(err)
                    return;
                end

            elseif strcmpi(allotrope,'gnr') && strcmp(param(1),'(') && strcmp(param(end),')') % (N_dl,N_dw)-GNR
                param = param(2:end-1);
                param_decomp = split(param,',');
                if ~isempty(str2num(param_decomp{1}))
                    N_dl = str2num(param_decomp{1});
                else
                    disp(err)
                    return;
                end
                if ~isempty(str2num(param_decomp{2}))
                    N_dw = str2num(param_decomp{2});
                else
                    disp(err)
                    return;
                end

            elseif strcmpi(allotrope,'cnt') && strcmp(param(1),'(') && strcmp(param(end),')') % (N_1(1),N_2(1))@(N_1(2),N_2(2))@...-CNT
                param = split(param,'@');
                for t = 1:length(param)
                    param_t = param{t};
                    if strcmp(param_t(1),'(') && strcmp(param_t(end),')')
                        param_t = param_t(2:end-1);
                        param_t_decomp = split(param_t,',');
                        if ~isempty(str2num(param_t_decomp{1}))
                            N_1(t) = str2num(param_t_decomp{1});
                        else
                            disp(err)
                            return;
                        end
                        if ~isempty(str2num(param_t_decomp{2}))
                            N_2(t) = str2num(param_t_decomp{2});
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