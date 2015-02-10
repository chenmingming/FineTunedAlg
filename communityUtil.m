classdef communityUtil
    methods(Static = true)
        function [communities]=getCommunities(communityFile)
            import java.util.*;
            communities=ArrayList;
            fid=fopen(communityFile);
            while 1
                tline = fgetl(fid);
                % Exit when empty
                if ~ischar(tline), break, end
                c=textscan(tline,'%d');
                c=c{1};
                len=size(c,1);
                community=HashSet;
                for i=1:len
                    community.add(c(i,1));
                end
                communities.add(community);
            end
            fclose(fid);
        end
    end
end