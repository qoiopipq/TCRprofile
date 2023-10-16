# TCRprofile
TCRprofile R package is designed for pairing alpha beta chains to label functional T cells from 10X immune profile outputs. 

## Aim 

The aim of TCRprofile package is to select functional T cells with paired TCR alpha and beta chains.

The selection criteria of paired TCR alpha and beta chains are either a cell with:

  -  1 alpha chain + 1 beta chain

  -  2 alpha chains + 1 beta chain

As a result of these, the summary of clonotype pairing table will provide a brief summary.

  -   A cell with functional alpha beta chains, will be classified as paired

  -   A cell with more than two chains (except for 2 alpha chains + 1 beta chain) will be classified as multi (cell with multiple chains).

  -   A cell with only 1 alpha chain or 1 beta chain, will be classified as tra or trb.

Details are in vignette. 
