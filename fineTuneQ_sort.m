function [communities] = fineTuneQ_sort(net,revnet,totalWeight,isUndirected,communities)
    splitSize = 0;
    mergeSize = 0;
    count=1;

    communitySize = communities.size;
    % When no splitting and no merging, exit.
    while communitySize~=splitSize||communitySize~=mergeSize
         disp(['Iteration ' num2str(count)]);
         count=count+1;
         communitySize = communities.size;
         
         % Split
         communities = splitQSort(net,revnet,totalWeight,isUndirected,communities);
         splitSize = communities.size;
         
         % Merge
         communities = mergeQ(net,totalWeight,communities);
         mergeSize = communities.size;
     end
     
     disp(['iterations ' num2str(count)]);
end

%---------------------------------------------------------------------------------------------------------------------------%
function [nodeCommunities communityWeights wouts]=initial(net,totalWeight,communities)
    import java.util.*;
    communitySize = communities.size;
    communityWeights=zeros(communitySize,communitySize);
	wouts=zeros(communitySize,1);
    nodeCommunities = HashMap;
    % Java ArrayList's index starts from 0
    for i=0:communitySize-1
        community = communities.get(i);
        comIter = community.iterator;
        while comIter.hasNext
            nodeCommunities.put(comIter.next, i);
        end
    end
        
    if communitySize==1
        communityWeights(1,1) = totalWeight;
    else 
        % Traverse the network to get info between communities
        nodeIter = net.entrySet.iterator;
        while nodeIter.hasNext
            nodeItem = nodeIter.next;
            nodeId = nodeItem.getKey;
            % MATLAB matrix index starts from 1 while Java ArrayList index
            % starts from 0
            communityId = nodeCommunities.get(nodeId)+1;
            neighbors = nodeItem.getValue;
            neighborIter = neighbors.entrySet.iterator;
            while neighborIter.hasNext
                neighborItem = neighborIter.next;
                neighborNodeId = neighborItem.getKey;
                neighborComId = nodeCommunities.get(neighborNodeId)+1;
                weight = neighborItem.getValue;
                communityWeights(communityId,neighborComId)=communityWeights(communityId,neighborComId)+weight;
				if communityId ~= neighborComId
				    wouts(communityId)=wouts(communityId)+weight;
				end
            end % neighbor while
        end % node while
    end % if
end

%-----------------------------------------------------------------------------------------------------------------------------------------------%
function [Qes] = getCommunityQ(communitySize,totalWeight,communityWeights,wouts)
    Qes=zeros(communitySize,1);
    for i=0:communitySize-1
        Q=communityWeights(i+1,i+1)/totalWeight-((communityWeights(i+1,i+1)+wouts(i+1))/totalWeight)^2;
        Qes(i+1)=Q;
    end
end

%----------------------------------------------------------------------%
function [indices]=laplacianMatrixSortCut(net,isUndirected,nodeCommunities,community)
    import java.util.*;
    % key: network node id; value: new subnetwork node id.
	nodeIdsMap = HashMap;
    % key: new subnetwork node id; value: network node id.
	newNodeIdsMap = HashMap;
    % Get subnetwork nodes' id
    comIter=community.iterator;
    count=1;
    while comIter.hasNext
        nodeId=comIter.next;
        nodeIdsMap.put(nodeId,count);
        newNodeIdsMap.put(count,nodeId);
        count=count+1;
    end
        
    comIter=community.iterator;
    laplacianMatrix=zeros(community.size,community.size);
    while comIter.hasNext
        nodeId=comIter.next;
        communityId=nodeCommunities.get(nodeId)+1;
        newNodeId=nodeIdsMap.get(nodeId);
        neighbors=net.get(nodeId);
        neighborIter=neighbors.entrySet.iterator;
        sumWeight=0;
        while neighborIter.hasNext
            neighborItem=neighborIter.next;
            nbId=neighborItem.getKey;
            nbComId=nodeCommunities.get(nbId)+1;
            if communityId==nbComId
                weight=neighborItem.getValue;
                newNbId=nodeIdsMap.get(nbId);
                laplacianMatrix(newNodeId,newNbId)=-weight;
                sumWeight=sumWeight+weight;
            end
        end % neighbor while
            
        %The diagonal degree matrix
        laplacianMatrix(newNodeId,newNodeId)=sumWeight;
    end % com node while

    % Get Fiedler Vector by using Lanczos Iterative Method
    laplacianMatrix=sparse(laplacianMatrix);
    opts.tol = 1e-6;
    oneVector=ones(community.size,1);
    %opts.v0=ones(community.size,1)/community.size;
    opts.v0=oneVector/norm(oneVector,2);
    if isUndirected
        opts.issym=1;
        [fiedlerVector d]=eigs(laplacianMatrix,2,'SA',opts);
    else
        opts.issym=0;
        [fiedlerVector d]=eigs(laplacianMatrix,2,'SR',opts);
    end
    fiedlerVector=fiedlerVector(:,2);
	%disp(fiedlerVector');
    % Sort fiedler vector in ascending order
	[fiedlerVector indices]=sort(fiedlerVector,'descend');
	%disp(fiedlerVector');
	for i=1:community.size
	    newNodeId=indices(i);
		nodeId=newNodeIdsMap.get(newNodeId);
		indices(i)=nodeId;
    end
    %disp(indices');
end

%----------------------------------------------------------------------------------------------------------------------%
% Split communities with Fiedler Vector and only when split increases the Modularity Density
function [obtainedCommunities] = splitQSort(net,revnet,totalWeight,isUndirected,communities)
    import java.util.*;
    % O(m): compute the weight matrix k*k
    [nodeCommunities communityWeights wouts]=initial(net,totalWeight,communities);
    communitySize=communities.size;
    % O(k): Compute Q of each community
    Qes=getCommunityQ(communitySize,totalWeight,communityWeights,wouts);
    % The new communities after splitting
    obtainedCommunities=ArrayList;
    % O(m+n): To see whether it is necessary to split a community
    for i=0:communitySize-1
        community=communities.get(i);
        if community.size <= 1
            obtainedCommunities.add(community);
            continue;
        end
        
		maxQ=-Inf;
		maxPosition=-1;
		% O(cnlog(cn))
        indices=laplacianMatrixSortCut(net,isUndirected,nodeCommunities,community);
		splittedComOne=HashSet;
		splittedComTwo=HashSet;
		for j=1:community.size
		    splittedComTwo.add(indices(j));
		end
       
        oneWin=0;
		oneWout=0;
        twoWin=communityWeights(i+1,i+1);
		twoWout=wouts(i+1);
        
		% Move 1 node to the first split community step by step, n-2 steps
		for j=1:community.size-1
		    nodeId=indices(j);
			% Move nodeId from the second split community to the first one
			splittedComOne.add(nodeId);
			splittedComTwo.remove(nodeId);
			% Outgoing
            neighbors=net.get(nodeId);
            neighborIter=neighbors.entrySet.iterator;
            while neighborIter.hasNext
                neighborItem=neighborIter.next;
                nbId=neighborItem.getKey;
                weight=neighborItem.getValue;
				% O(1): HashSet.contains
				if splittedComOne.contains(nbId)
				    oneWin=oneWin+weight;
                    twoWout=twoWout-weight;
                    if isUndirected
                        oneWin=oneWin+weight;
                        oneWout=oneWout-weight;
                    end
                elseif splittedComTwo.contains(nbId)
				    twoWin=twoWin-weight;
                    oneWout=oneWout+weight;
                    if isUndirected
                        twoWin=twoWin-weight;
                        twoWout=twoWout+weight;
                    end
                else
                    oneWout=oneWout+weight;
                    twoWout=twoWout-weight;
                end
            end
            
            % Incoming
            if ~isUndirected
                neighbors=revnet.get(nodeId);
                neighborIter=neighbors.entrySet.iterator;
                while neighborIter.hasNext
                    neighborItem=neighborIter.next;
                    nbId=neighborItem.getKey;
                    weight=neighborItem.getValue;
				    % O(1): HashSet.contains
				    if splittedComOne.contains(nbId)
				        oneWin=oneWin+weight;
                        oneWout=oneWout-weight;
                    elseif splittedComTwo.contains(nbId)
				        twoWin=twoWin-weight;
                        twoWout=twoWout+weight;
                    end
                end
            end
            
            %disp([oneWin oneWout twoWin twoWout]);
			oneQ=oneWin/totalWeight-((oneWin+oneWout)/totalWeight)^2;
			twoQ=twoWin/totalWeight-((twoWin+twoWout)/totalWeight)^2;
			sumQ=oneQ+twoQ;
            %disp([oneQ twoQ sumQ Qes(i+1)]);
		    if sumQ>maxQ
		        maxQ=sumQ;
			    maxPosition=j;
		    end
		end
        
        splitFlag=false;
        
        if maxQ>Qes(i+1)
            splitFlag=true;
            if (maxPosition==1)||(maxPosition==community.size-1) 
                splitFlag=false;
                comId=nodeCommunities.get(nodeId)+1;
                inComEdges=0;
                
                % Outgoing
                outNeighbors=net.get(nodeId);
                outNeighborIter=outNeighbors.entrySet.iterator;
                while outNeighborIter.hasNext
                    neighborItem=outNeighborIter.next;
                    nbId=neighborItem.getKey;
                    nbComId=nodeCommunities.get(nbId)+1;
                    if comId==nbComId
                        inComEdges=inComEdges+1;
                    end
                end
                
                % Incoming
                if ~isUndirected
                    inNeighbors=revnet.get(nodeId);
                    inNeighborIter=inNeighbors.entrySet.iterator;
                    while inNeighborIter.hasNext
                        neighborItem=inNeighborIter.next;
                        nbId=neighborItem.getKey;
                        nbComId=nodeCommunities.get(nbId)+1;
                        if comId==nbComId
                            inComEdges=inComEdges+1;
                        end
                    end
                end
                
                if inComEdges==0
                    splitFlag=true;
                end
            end
        end
        
		if splitFlag
            splittedComOne.clear;
			splittedComTwo.clear;
			for j=1:community.size
			    if j<=maxPosition
				    splittedComOne.add(indices(j));
                else
				    splittedComTwo.add(indices(j));
                end
            end
            
            if splittedComOne.size>splittedComTwo.size
                obtainedCommunities.add(splittedComOne);
                obtainedCommunities.add(splittedComTwo);
            else
                obtainedCommunities.add(splittedComTwo);
                obtainedCommunities.add(splittedComOne);
            end
        else
            obtainedCommunities.add(communities.get(i));
        end % end if max
    end % community for
    
    % Clear the collections to save the space
	communities.clear();
	nodeCommunities.clear();
end

%-------------------------------------------------------------------------------------------------------------------------%
function [obtainedCommunities] = mergeQ(net,totalWeight,communities)
    import java.util.*;
    % O(m): compute the weight matrix k*k
    [nodeCommunities communityWeights wouts]=initial(net,totalWeight,communities);
    communitySize=communities.size;
    % O(k): Compute Q of each community
    Qes=getCommunityQ(communitySize,totalWeight,communityWeights,wouts);
    
    % O(k^2+k^2logk^2): Compute Q of two combined communities
    combinedCommunities=TreeMap;
    for i=0:communitySize-1
        for j=i+1:communitySize-1
            % Do not need to consider to merge 2 communities without any edge between them
            if communityWeights(i+1,j+1)==0&&communityWeights(j+1,i+1)==0
                continue;
            end
            win=communityWeights(i+1,i+1)+communityWeights(j+1,j+1)+communityWeights(i+1,j+1)+communityWeights(j+1,i+1);
			wout=wouts(i+1)+wouts(j+1)-communityWeights(i+1,j+1)-communityWeights(j+1,i+1);
            Q=win/totalWeight-((win + wout)/totalWeight)^2;
            sumQ=Qes(i+1)+Qes(j+1);
            if Q>sumQ
                increase=Q-sumQ;
                %disp([i j Q sumQ]);
                % If key exists, add eps to generate a new key
                while combinedCommunities.containsKey(-increase)
                    increase=increase+eps;
                end
                combinedCommunities.put(-increase, [i j]);
            end
        end %j for 
    end %i for
    
    totalIncrease=0;
    combinedLists=HashMap;
    combinedComIter=combinedCommunities.entrySet.iterator;
    while combinedComIter.hasNext
        combinedComItem=combinedComIter.next;
        increase=combinedComItem.getKey;
        twoComs=combinedComItem.getValue;
        firstCom=twoComs(1);
        secondCom=twoComs(2);
        if ~combinedLists.containsKey(firstCom)&&~combinedLists.containsKey(secondCom)
            combinedLists.put(firstCom,secondCom);
            combinedLists.put(secondCom,firstCom);
            totalIncrease=totalIncrease-increase;
        end
    end
    
    % O(k): obtain the combined communities
    obtainedCommunities=ArrayList;
    for i=0:communitySize-1
        community=communities.get(i);
        if combinedLists.containsKey(i)
            combinedComId=combinedLists.get(i);
            if i<combinedComId
                community.addAll(communities.get(combinedComId));
                obtainedCommunities.add(community);
            end
        else
            obtainedCommunities.add(community);
        end
    end
    
    communities.clear();
    nodeCommunities.clear();
	combinedCommunities.clear();
	combinedLists.clear();
end