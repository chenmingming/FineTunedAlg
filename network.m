classdef network
    methods(Static = true)
        % Read the network frome file. Return the network object and the
        % total weights of the network
        function [net totalEdge totalWeight] = getNetwork(networkFile, isUnweighted, isUndirected)
            %import java.util.*;
            net = java.util.HashMap;
            totalEdge = 0;
            totalWeight = 0;
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
                
                if ~net.containsKey(src)
                    net.put(src, java.util.HashMap);
                end
                neighbors = net.get(src);
                %disp(neighbors);
                if ~neighbors.containsKey(dst)
                    neighbors.put(dst, weight);
                    totalEdge = totalEdge + 1;
                    totalWeight = totalWeight + weight;
                end
                
                if ~net.containsKey(dst)
                    net.put(dst, java.util.HashMap);
                end
                % undirected network
                if isUndirected
                    neighbors = net.get(dst);
                    if ~neighbors.containsKey(src)
                        neighbors.put(src, weight);
                        totalEdge = totalEdge + 1;
                        totalWeight = totalWeight + weight;
                    end
                end
            end
            fclose(fid);
        end
    end
end