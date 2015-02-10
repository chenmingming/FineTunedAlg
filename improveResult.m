import java.util.*;

% Please change the following two variables according to the properties the
% networks
isUnweighted = true;
isUndirected = true;

% The path and filename of the networks
networkFile = 'example_data/network1.txt';
% The community file of other community detection algorithms
communityFile = 'example_data/result.txt';
% The community file improved with Fine-tuned Qds
communityResultFile = 'example_data/improved_result.txt';

[net totalEdge totalWeight] = network.getNetwork(networkFile, isUnweighted, isUndirected);
        
if ~isUndirected
    revnet=reversedNetwork.getReversedNetwork(networkFile, isUnweighted, isUndirected);
else
    revnet=HashMap;
end

% Read the existed communities of other algorithms
communities=communityUtil.getCommunities(communityFile);
disp(communities.size);

for i=0:communities.size-1
   disp(communities.get(i));
end

disp('After improved with Qds.');

communities = fineTuneQds_sort(net,revnet,totalEdge,totalWeight,isUndirected,communities);
%communities = fineTuneQ_sort(net,revnet,totalWeight,isUndirected,communities);

% Output community detection result
fid = fopen(communityResultFile,'w');
communitySize=communities.size;
for i=0:communitySize-1
   community=communities.get(i);
   disp(community);
   comNodeIter=community.iterator;
   while comNodeIter.hasNext
       nodeId=comNodeIter.next;
       fprintf(fid,'%d ',nodeId);
   end
   fprintf(fid,'\n');
end
fclose(fid);