# Grivan-Newman Algorithm

I will be implementing a "Divisive Method" of the Girvan-Newman Algorithm
In divisive methods, we start with the complete graph and take off the edges iteratively. The edge with the highest weight is removed first. At every step, the edge-weight calculation is repeated, since the weight of the remaining edges changes after an edge is removed. After a certain number of steps, we get clusters of densely connected nodes.

https://www.analyticsvidhya.com/blog/2020/04/community-detection-graphs-networks/ 

I also needed to implement a helper algorithm to find th number of shortest paths for each vertex from a single source:
https://www.baeldung.com/cs/graph-number-of-shortest-paths
