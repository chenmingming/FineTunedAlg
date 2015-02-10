import java.util.*;

% Please change the following two variables according to the properties the
% networks
isUnweighted = true;
isUndirected = true;

% The path and filename of the networks
networkFile = 'example_data/network1.txt';
% The file to store the detected communities
communityResultFile = 'example_data/result.txt';

[net totalEdge totalWeight] = network.getNetwork(networkFile, isUnweighted, isUndirected);

disp(['Number of nodes = ' num2str(net.size) ', number of edges = ' num2str(totalEdge)]);

if ~isUndirected
    revnet=reversedNetwork.getReversedNetwork(networkFile, isUnweighted, isUndirected);
else
    revnet=0;
end

% The communities
% Get communities from file or treat all nodes as a whole community
communities = ArrayList;
community = HashSet;
nodeIter = net.entrySet.iterator;
while nodeIter.hasNext
    item = nodeIter.next;
    community.add(item.getKey);
end
communities.add(community);

% Fine-tuned Q
%communities = fineTuneQ_sort(net,revnet,totalWeight,isUndirected,communities);
% Fine-tuned Qds
communities = fineTuneQds_sort(net,revnet,totalEdge,totalWeight,isUndirected,communities);

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




