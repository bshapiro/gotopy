# GOTOpy

GOTOpy provides a collection of functions for computing the Gene Ontology Term Overlap of collections of genes or methylation probes (hereafter referred to as entities). It includes functions that allow the comparison of pairs of entities, clusters of entities or clusterings of entities. Additionally, we include a script fetch_go_terms.py which generates maps from genes or probes to corresponding sets of Gene Ontology terms for use with GOTOpy. 

GOTOpy is also configurable for use with different methylation annotation sets, gene ontology categories (Biological Process, Molecular Function, Cellular Component, or all of these) and options for filtering the gene ontology terms to use in overlap calculations. 

The GO term fetcher script is documented further using optparse. Please run 'python fetch_go_terms.py -h' for more detailed information on the various possible configurations of your GOTOpy analysis. 

## Setup 
GOTOPY requires the following packages:

- numpy
- pandas
- GOATOOLS
- mygene

NB: Before running, please ensure that the global variable GOTOPY_HOME in gotopy.py is set to the root directory of your GOTOpy checkout. 

## More information 

GOTOpy was developed as part of the work presented in academic work on Bayesian Multi-View Clustering. If you have any questions regarding its use, please feel free to submit issues here, or reach out to the authors. 

GOTOpy is available for public use under a BSD-3 license. 
