function [communities] = fineTuneQds_sort(net,revnet,totalEdge,totalWeight,isUndirected,communities)
    splitSize = 0;
    mergeSize = 0;
    count=1;
    communitySize = communities.size;
    % When no splitting and not merging, exit.
    while communitySize~=splitSize||communitySize~=mergeSize
        disp(['Iteration ' num2str(count)]);
        count=count+1;
        communitySize = communities.size;
        %disp('split:');
        communities = splitQdsSort(net,revnet,totalEdge,totalWeight,isUndirected,communities);
        splitSize = communities.size;
        
        %disp('merge:');
        communities = mergeQds(net,totalEdge,totalWeight,communities);
        mergeSize = communities.size;
    end
    
    disp(['iterations ' num2str(count)]);
end

%---------------------------------------------------------------------------------------------------------------------------%
function [nodeCommunities communityWeights communityEdges communityDensities]=initial(net,totalEdge,totalWeight,communities)
    import java.util.*;
    communitySize = communities.size;
    communityWeights=zeros(communitySize,communitySize);
    communityEdges=zeros(communitySize,communitySize);
    communityDensities=zeros(communitySize,communitySize);
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
        numNodes = net.size;
        communityWeights(1,1) = totalWeight;
        communityEdges(1,1) = totalEdge;
        communityDensities(1,1)=totalEdge/(numNodes*(numNodes-1));
    else 
        % Traverse the network to get info between communities
        nodeIter = net.entrySet.iterator;
        while nodeIter.hasNext
            nodeItem = nodeIter.next;
            nodeId = nodeItem.getKey;
            % MATLAB matrix index starts from 1 while Java ArrayList index
            % starts from 0
            communityId = nodeCommunities.get(nodeId)+1;
            comSize = communities.get(communityId-1).size;
            neighbors = nodeItem.getValue;
            neighborIter = neighbors.entrySet.iterator;
            while neighborIter.hasNext
                neighborItem = neighborIter.next;
                neighborNodeId = neighborItem.getKey;
                neighborComId = nodeCommunities.get(neighborNodeId)+1;
                neighborComSize = communities.get(neighborComId-1).size;
                weight = neighborItem.getValue;
                communityWeights(communityId,neighborComId)=communityWeights(communityId,neighborComId)+weight;
                communityEdges(communityId,neighborComId)=communityEdges(communityId,neighborComId)+1;
                if communityId == neighborComId
                    communityDensities(communityId,communityId)=communityDensities(communityId,communityId)+1/(comSize*(comSize-1));
                else
                    communityDensities(communityId,neighborComId)=communityDensities(communityId,neighborComId)+1/(comSize*neighborComSize);
                end
            end % neighbor while
        end % node while
    end % if
end

%-----------------------------------------------------------------------------------------------------------------------------------------------%
function [Qdses] = getCommunityQds(communitySize,totalWeight,communityWeights,communityDensities)
    Qdses=zeros(communitySize,1);
    for i=0:communitySize-1
        wout=0;
        splitPenalty=0;
        for j=0:communitySize-1
            if j~=i
                wout=wout+communityWeights(i+1,j+1);
                splitPenalty=splitPenalty+(communityWeights(i+1,j+1)/totalWeight)*communityDensities(i+1,j+1);
            end
        end
        Qds=(communityWeights(i+1,i+1)/totalWeight)*communityDensities(i+1,i+1)-...
            (((communityWeights(i+1,i+1)+wout)/totalWeight)*communityDensities(i+1,i+1))^2-splitPenalty;
        Qdses(i+1)=Qds;
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
    %laplacianMatrix=[];
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
	
	%laplacianMatrix
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
	%d
	%fiedlerVector
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

%---------------------------------------------------------------------------------------------------------------%
% Get the weights and densities information of the split community to other communities
function [weights edges densities]=getSplitCommunityInfo(net,communities,nodeCommunities,communitySize,splitCommunity,twoComSize)
    weights=zeros(communitySize+1,1);
	edges=zeros(communitySize+1,1);
	densities=zeros(communitySize+1,1);
	comSize=splitCommunity.size;
	comIter=splitCommunity.iterator;
	while comIter.hasNext
        nodeId=comIter.next;
        communityId=nodeCommunities.get(nodeId)+1;
        neighbors=net.get(nodeId);
        neighborIter=neighbors.entrySet.iterator;
        while neighborIter.hasNext
            neighborItem=neighborIter.next;
            nbId=neighborItem.getKey;
            nbComId=nodeCommunities.get(nbId)+1;
            weight=neighborItem.getValue;
            if communityId==nbComId
                if splitCommunity.contains(nbId)
                    weights(communityId)=weights(communityId)+weight;
					edges(communityId)=edges(communityId)+1;
                    densities(communityId)=densities(communityId)+1/(comSize*(comSize-1));
                else
                    weights(communitySize+1)=weights(communitySize+1)+weight;
					edges(communitySize+1)=edges(communitySize+1)+1;
                    densities(communitySize+1)=densities(communitySize+1)+1/(comSize*twoComSize);
                end
            else
                nbComSize=communities.get(nbComId-1).size;
                weights(nbComId)=weights(nbComId)+weight;
				edges(nbComId)=edges(nbComId)+1;
                densities(nbComId)=densities(nbComId)+1/(comSize*nbComSize);
            end
        end % neighbor while
    end % com one while
end

%----------------------------------------------------------------------------------------------------------------------%
% Split communities with Fiedler Vector and only when split increases the Modularity Density
function [obtainedCommunities] = splitQdsSort(net,revnet,totalEdge,totalWeight,isUndirected,communities)
    import java.util.*;
    % O(m): compute the weight, edge, and density matrix k*k
    [nodeCommunities communityWeights communityEdges communityDensities]=initial(net,totalEdge,totalWeight,communities);
    communitySize=communities.size;
    % O(k^2): Compute Modularity Density of each community
    Qdses=getCommunityQds(communitySize,totalWeight,communityWeights,communityDensities);
    % The new communities after splitting
    obtainedCommunities=ArrayList;
    % O(m+kn+nlgn): To see whether it is necessary to split a community
    for i=0:communitySize-1
        community=communities.get(i);
		comSize=community.size;
        if comSize <= 1
            obtainedCommunities.add(community);
            continue;
        end
        
        maxQds=-Inf;
		maxPosition=-1;
		% O((cn)log(cn)): cn is the size of the subcommunity
        indices=laplacianMatrixSortCut(net,isUndirected,nodeCommunities,community);
        %disp('indices');
        %disp(indices');
        splittedComOne=HashSet;
		splittedComTwo=HashSet;
        % O(cn)
	    for j=1:comSize
		    splittedComTwo.add(indices(j));
        end
        oneWeights=zeros(communitySize+1,1);
	    oneEdges=zeros(communitySize+1,1);
        % Get weight and density matrix for the second splitted community
        [twoWeights twoEdges twoDensities]=getSplitCommunityInfo(net,communities,nodeCommunities,communitySize,splittedComTwo,0);
        % The change of Split Penalty of other communities for directed networks
		if ~isUndirected
            revOneWeights=zeros(communitySize+1,1);
            revOneEdges=zeros(communitySize+1,1);
			[revTwoWeights revTwoEdges revTwoDensities]=getSplitCommunityInfo(revnet,communities,nodeCommunities,communitySize,splittedComTwo,0);
            %disp('initial');
            %disp(revTwoWeights');
        end
        
        % Move 1 node to the first split community step by step, n-2 steps
		for j=1:comSize-1
            sumQds=0;
		    nodeId=indices(j);
            communityId=nodeCommunities.get(nodeId)+1;
			% Move nodeId from the second split community to the first one
			splittedComOne.add(nodeId);
			splittedComTwo.remove(nodeId);                        
            oneComSize=splittedComOne.size;
            twoComSize=splittedComTwo.size;
			% Outgoing
            neighbors=net.get(nodeId);
            neighborIter=neighbors.entrySet.iterator;
            while neighborIter.hasNext
                neighborItem=neighborIter.next;
                nbId=neighborItem.getKey;
                nbComId=nodeCommunities.get(nbId)+1;
                weight=neighborItem.getValue;
                % O(1): HashSet.contains
				if splittedComOne.contains(nbId)
				    oneWeights(communityId)=oneWeights(communityId)+weight;
                    oneEdges(communityId)=oneEdges(communityId)+1;
                    twoWeights(communitySize+1)=twoWeights(communitySize+1)-weight;
                    twoEdges(communitySize+1)=twoEdges(communitySize+1)-1;
                    if isUndirected
                        oneWeights(communityId)=oneWeights(communityId)+weight;
                        oneEdges(communityId)=oneEdges(communityId)+1;
                        oneWeights(communitySize+1)=oneWeights(communitySize+1)-weight;
                        oneEdges(communitySize+1)=oneEdges(communitySize+1)-1;
                    end
                elseif splittedComTwo.contains(nbId)
				    twoWeights(communityId)=twoWeights(communityId)-weight;
                    twoEdges(communityId)=twoEdges(communityId)-1;
                    oneWeights(communitySize+1)=oneWeights(communitySize+1)+weight;
                    oneEdges(communitySize+1)=oneEdges(communitySize+1)+1;
                    if isUndirected
                        twoWeights(communityId)=twoWeights(communityId)-weight;
                        twoEdges(communityId)=twoEdges(communityId)-1;
                        twoWeights(communitySize+1)=twoWeights(communitySize+1)+weight;
                        twoEdges(communitySize+1)=twoEdges(communitySize+1)+1;
                    end
                else
                    oneWeights(nbComId)=oneWeights(nbComId)+weight;
                    oneEdges(nbComId)=oneEdges(nbComId)+1;
                    twoWeights(nbComId)=twoWeights(nbComId)-weight;
                    twoEdges(nbComId)=twoEdges(nbComId)-1;
                end
            end % while
            
            % Incoming
            if ~isUndirected
                neighbors=revnet.get(nodeId);
                neighborIter=neighbors.entrySet.iterator;
                while neighborIter.hasNext
                    neighborItem=neighborIter.next;
                    nbId=neighborItem.getKey;
                    nbComId=nodeCommunities.get(nbId)+1;
                    weight=neighborItem.getValue;
				    % O(1): HashSet.contains
				    if splittedComOne.contains(nbId)
				        oneWeights(communityId)=oneWeights(communityId)+weight;
                        oneEdges(communityId)=oneEdges(communityId)+1;
                        oneWeights(communitySize+1)=oneWeights(communitySize+1)-weight;
                        oneEdges(communitySize+1)=oneEdges(communitySize+1)-1;
                    elseif splittedComTwo.contains(nbId)
				        twoWeights(communityId)=twoWeights(communityId)-weight;
                        twoEdges(communityId)=twoEdges(communityId)-1;
                        twoWeights(communitySize+1)=twoWeights(communitySize+1)+weight;
                        twoEdges(communitySize+1)=twoEdges(communitySize+1)+1;
                    else
                        revOneWeights(nbComId)=revOneWeights(nbComId)+weight;
                        revOneEdges(nbComId)=revOneEdges(nbComId)+1;
                        revTwoWeights(nbComId)=revTwoWeights(nbComId)-weight;
                        revTwoEdges(nbComId)=revTwoEdges(nbComId)-1;
                    end
                end
            end
            
            % Compute the Modularity Density of first splitted community
            oneWin=oneWeights(i+1);
            oneEin=oneEdges(i+1);
            oneWout=0;
            oneSplitPenalty=0;
            for k=1:communitySize+1
                if (k-1)~=i
                    if k~=communitySize+1
                        nbCommunity=communities.get(k-1);
						nbComSize=nbCommunity.size;
                        interDensity=oneEdges(k)/(oneComSize*nbComSize);
                        sp=(oneWeights(k)/totalWeight)*interDensity;
			            % The change of Split Penalty of other communities
			            if isUndirected      
				            sumQds=sumQds+(communityWeights(k,i+1)/totalWeight)*(communityEdges(k,i+1)/(comSize*nbComSize))-sp;
                            %disp([((communityWeights(k,i+1)/totalWeight)*(communityEdges(k,i+1)/(comSize*nbComSize))) sp]);
                        end
                    else
                        interDensity=oneEdges(k)/(oneComSize*twoComSize);
                        sp=(oneWeights(k)/totalWeight)*interDensity;
                    end
                    oneWout=oneWout+oneWeights(k);
                    oneSplitPenalty=oneSplitPenalty+sp;
                end
            end
            if oneComSize<=1
                inDensity=0;
                %inDensity = 1;
            else
                inDensity=oneEin/(oneComSize*(oneComSize-1));
            end
            oneQds=(oneWin/totalWeight)*inDensity-(((oneWin+oneWout)/totalWeight)*inDensity)^2-oneSplitPenalty;
            sumQds=sumQds+oneQds;
            
            % Compute the Modularity Density of first splitted community
            twoWin=twoWeights(i+1);
            twoEin=twoEdges(i+1);
            twoWout=0;
            twoSplitPenalty=0;
            for k=1:communitySize+1
                if (k-1)~=i
                    if k~=communitySize+1
                        nbCommunity=communities.get(k-1);
						nbComSize=nbCommunity.size;
                        interDensity=twoEdges(k)/(twoComSize*nbComSize);
                        sp=(twoWeights(k)/totalWeight)*interDensity;
                        % The change of Split Penalty of other communities
			            if isUndirected
				            sumQds=sumQds-sp;
                            %disp(sp);
                        end
                    else
                        interDensity=twoEdges(k)/(twoComSize*oneComSize);
                        sp=(twoWeights(k)/totalWeight)*interDensity;
                    end
                    twoWout=twoWout+twoWeights(k);
                    twoSplitPenalty=twoSplitPenalty+sp;
                end
            end
            if twoComSize<=1
                inDensity=0;
                %inDensity = 1;
            else
                inDensity=twoEin/(twoComSize*(twoComSize-1));
            end
            twoQds=(twoWin/totalWeight)*inDensity-(((twoWin+twoWout)/totalWeight)*inDensity)^2-twoSplitPenalty;
            sumQds=sumQds+twoQds;
            
            % The change of Split Penalty of other communities for directed networks
		    if ~isUndirected
			    for k=1:communitySize
			        if (k-1)~=i 
				        nbCommunity=communities.get(k-1);
				        nbComSize=nbCommunity.size;
                        %tmp=sumQds;
				        sumQds=sumQds+(communityWeights(k,i+1)/totalWeight)*(communityEdges(k,i+1)/(comSize*nbComSize))...
					        -(revOneWeights(k)/totalWeight)*(revOneEdges(k)/(oneComSize*nbComSize))...
                            -(revTwoWeights(k)/totalWeight)*(revTwoEdges(k)/(twoComSize*nbComSize));
                        %disp(sumQds-tmp);
                    end
                end
            end
            
            if sumQds>maxQds
                maxQds=sumQds;
                maxPosition=j;
            end
        end  % for
        
        splitFlag=false;
        
        if maxQds>Qdses(i+1)
            splitFlag=true;
            if (maxPosition==1)||(maxPosition==comSize-1) 
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
			for j=1:comSize
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
function [obtainedCommunities] = mergeQds(net,totalEdge,totalWeight,communities)
    import java.util.*;
    % O(m): compute the weight, edge, and density matrix k*k
    [nodeCommunities communityWeights communityEdges communityDensities]=initial(net,totalEdge,totalWeight,communities);
    communitySize=communities.size;
    % O(k^2): Compute Modularity Density of each community
    Qdses=getCommunityQds(communitySize,totalWeight,communityWeights,communityDensities);
    
    % O(k^3+k^2logk^2): Compute Modularity Density of two combined communities
    combinedCommunities=TreeMap;
    for i=0:communitySize-1
        for j=i+1:communitySize-1
            % Do not need to consider to merge 2 communities without any edge between them
            if communityWeights(i+1,j+1)==0&&communityWeights(j+1,i+1)==0
                continue;
            end
			deltaQds=0;
            win=communityWeights(i+1,i+1)+communityWeights(j+1,j+1)+communityWeights(i+1,j+1)+communityWeights(j+1,i+1);
			ein=communityEdges(i+1,i+1)+communityEdges(j+1,j+1)+communityEdges(i+1,j+1)+communityEdges(j+1,i+1);
			icomSize=communities.get(i).size;
			jcomSize=communities.get(j).size;
			csize=icomSize+jcomSize;
			wout=0;
			splitPenalty=0;
            
            for k=0:communitySize-1
                if k~=i&&k~=j
				    kcomSize=communities.get(k).size;
                    if communityWeights(i+1,k+1)>0||communityWeights(j+1,k+1)>0
                        cout=communityWeights(i+1,k+1)+communityWeights(j+1,k+1);
                        wout=wout+cout;
                        outDensity=(communityEdges(i+1,k+1)+communityEdges(j+1,k+1))/(csize*kcomSize);
						splitPenalty=splitPenalty+(cout/totalWeight)*outDensity;
                    end
					
					% The change of Split Penalty of other communities
					if communityWeights(k+1,i+1)>0||communityWeights(k+1,j+1)>0
						deltaQds=deltaQds-((communityWeights(k+1,i+1)+communityWeights(k+1,j+1))/totalWeight)*((communityEdges(k+1,i+1)+communityEdges(k+1,j+1))/(kcomSize*csize))...
						    +(communityWeights(k+1,i+1)/totalWeight)*(communityEdges(k+1,i+1)/(kcomSize*icomSize))...
							+(communityWeights(k+1,j+1)/totalWeight)*(communityEdges(k+1,j+1)/(kcomSize*jcomSize));
				    end
                end
            end
            
            inDensity=ein/(csize*(csize-1));
            Qds=(win/totalWeight)*inDensity-(((win + wout)/totalWeight)*inDensity)^2-splitPenalty;
            %communities.get(i)
            %communities.get(j)
            %disp([deltaQds Qds Qdses(i+1) Qdses(j+1)]);
            deltaQds=deltaQds+Qds-Qdses(i+1)-Qdses(j+1);
            %disp(deltaQds);
            if deltaQds>0
                %disp('merge:')
                %disp([Qds Qdses(i+1) Qdses(j+1) (Qds-Qdses(i+1)-Qdses(j+1)) deltaQds]);
                increase=deltaQds;
                % If key exists, add eps to generate a new key
                while combinedCommunities.containsKey(-increase)
                    disp('exist');
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
                %disp(community);
                %disp(communities.get(combinedComId));
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