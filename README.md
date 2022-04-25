### Description
FilterReads is a quality-based filtering program for paired-end sequencing data. Given read files 1 and 2 and the total number of filtered reads (k), this program keeps k middle read pairs of the input data.

### Assumptions:
- Loading the entire input data is not allowed due to the limited memory space for such big data.
- Total number of input read pairs is not known apriori.
- Input fastq.gz files are locally available.

### Time complexity:
This program comprises of the following steps:
1. Collect the read pair quality scores: **O(n)**
2. Random shuffle for unbiased selection in case of ties: **O(n)**
3. Sort read indices based on their quality score: **O(nlog(n))**
4. Load read indices into a set in order to look up membership in average constant time: **O(n)**
5. Read and filter input for the selected read pairs: **O(n)**

Therefore, the overall time complexity is **O(nlog(n))**.

### Space complexity:
This program keeps an array of (index, quality) pairs of all input read pairs. So, the overall space complexity is **O(n)**.
