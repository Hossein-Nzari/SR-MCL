# SR-MCL
[MCL(Markov clustering algorithm) and R-MCL](https://sci-hub.se/https:/doi.org/10.1145/1854776.1854812) are two of the most famous graph clustering algorithms. 
However, they have their downsides. For instance, MCL does not scale well, even to moderate-sized graphs. 
And R-MCL cannot identify hierarchical modules,  despite being effective and efficient in hard clustering. 
Hence [SR-MCL](https://academic.oup.com/bioinformatics/article/28/18/i473/243788?login=true) was developed to improve the R-MCL performance in recognizing overlapping clusters in 
large graphs, like protein protein interaction networks. In this project, the implementation and performance of 
SR-MCL will be explored. 
 
1. Implementation of SR-MCL algorithm 
2. Deployment of the algorithm on [WI-PHI](https://www.yeastgenome.org/reference/S000120766) 
3. Defining an evaluation metric for graph clustering and studying and explaining the effects of SR-MCL’s 
parameters on the algorithm performance based on the evaluation metric.

# SR-MCL (Self-Regulative Markov Cluster Algorithm) for Network Clustering

SR-MCL is a bioinformatics algorithm designed for network clustering. It detects clusters within a network graph by iteratively transforming and analyzing the graph. This algorithm can be particularly useful in identifying functional modules or communities in biological networks. In this GitHub repository, you'll find a Python implementation of SR-MCL for network clustering.

## Getting Started

To use SR-MCL, follow these steps:

1. **Prerequisites**
   - Python 3.x
   - Required Python libraries (NumPy, SciPy, NetworkX, Pandas, and Matplotlib)

2. **Installation**
   - Clone the repository to your local machine:
     ```
     git clone <[repository_url](https://github.com/Hossein-Nzari/SR-MCL.git)>
     ```
   
3. **Data Preparation**
   - Prepare your network data in an Excel file (`.xls`) with three columns: 'Protein A,' 'Protein B,' and 'Score.' This data represents the edges and edge weights in your network.

4. **Run SR-MCL**
   - Open the provided Python script and make sure to update the file path for your Excel data.
   - Customize algorithm parameters like `r`, `p`, and `ω` to suit your specific dataset and clustering requirements.
   - Execute the script to run SR-MCL and obtain network clusters.

## Algorithm Parameters

- `r`: This parameter controls the inflation factor, affecting the granularity of clusters. Higher values make clusters more compact.
- `p`: The parameter for post-processing, determining the minimum overlap required between clusters to merge them. It filters out smaller, less significant clusters.
- `ω` (omega): The threshold for pruning matrix values. Values below this threshold are set to zero, effectively reducing noise in the network.

## Results and Clusters

After running SR-MCL, the algorithm will generate network clusters. The `c` variable contains the clusters as lists of nodes. You can visualize these clusters using NetworkX or export the results for further analysis.

## Post-processing

To further refine the clusters, a post-processing step is included. It filters out clusters with fewer than three nodes and merges clusters based on a similarity threshold `p`. The average quality function (`qf`) of clusters is also calculated.

## Example

For a detailed example of how to use SR-MCL, please refer to the provided Python script in this repository.

## Tuning Parameters

As mentioned in the code comments, adjusting the parameters `r`, `p`, and `ω` can significantly impact the results. Increasing these parameters will lead to a decrease in the number of clusters and an increase in the average quality function (`qf`), potentially giving you more consolidated and meaningful clusters.
