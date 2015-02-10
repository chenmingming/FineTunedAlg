--------- To run the algorithms ---------
1. Please look at "fineTuneAlgorithm.m" for how to run Fine-tuned Q and Fine-tuned Qds algorithms.
2. Please look at "imporveResult.m" for how to use Fine-tuned Qds to improve the community deteciton results that are found with other community detection algorithms.
3. You should change the varaibles "isUnweighted" and "isUndirected" to corresponds to the properties of the netowrks (undirected, unweighted, directed, weighted).
4. Fine-tuned Q is implemented in "fineTuneQ_sort.m", and Fine-tuned Qds is implemented in "fineTuneQds_sort.m".
5. The algorithms are implemented with MATLAB, incorporating some features from Java, so use MATLAB to run these scripts or functions.

--------- Input network and community format --------- 
1. Network format: each line in the file should have a format "srcId dstId" for unweighted networks or "srcId dstId weight" for weighted networks. The delimiter can be a single space or a tab. please look at "example_data/network1.txt" for an example.
2. Community format: each line in the file represents a community. The nodes in each line (or community) is seperated with a single space. Please look at file "example_data/result.txt" for an example.

---------- Output community format ----------
1. The detected communities are outputed into a file. Each line in the file represents a community. The nodes in each line (or community) is seperated with a single space.

---------- Note ----------
Note: The computational complexity of both algorithms is roughly O(|E|+|V|log|V|). However, the number of iterations that each algorithm takes may be a little large for large networks, so these two algorithms cannot scale to too large networks. We are proposing faster algorithms, please check in the near future if you are interested in.