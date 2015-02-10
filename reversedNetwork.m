% Only used for directed networks
classdef reversedNetwork
    methods(Static = true)
        % Read the network frome file. Return the network object and the
        % total weights of the network
        function [invnet] = getReversedNetwork(networkFile, isUnweighted, isUndirected)
            %import java.util.*;
            invnet = java.util.HashMap;
            fid=fopen(networkFile);
            while 1
                tline = fgetl(fid);
                % Exit when empty
                if ~ischar(tline), break, end
                weight = 1;
                % unweighted network
                if isUnweighted
                    [src dst] = strread(tline, '%d%d', 1, 'delimiter', ' |\t');
                else
                    [src dst weight] = strread(tline, '%d%d%f', 1, 'delimiter', ' |\t');
                end
                
                if src == dst
                    continue;
                end
                
                if ~invnet.containsKey(dst)
                    invnet.put(dst, java.util.HashMap);
                end
                neighbors = invnet.get(dst);
                %disp(neighbors);
                if ~neighbors.containsKey(src)
                    neighbors.put(src, weight);
                end
                
                if ~invnet.containsKey(src)
                    invnet.put(src, java.util.HashMap);
                end
                % undirected network
                if isUndirected
                    neighbors = invnet.get(src);
                    if ~neighbors.containsKey(dst)
                        neighbors.put(dst, weight);
                    end
                end
            end
            fclose(fid);
        end
    end
end