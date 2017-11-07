# complenet17-graph-generator
Implementation of the method described in the paper: 

Miguel E. P. Silva, P. Paredes, and P. Ribeiro, “Network motifs detection using random networks with prescribed subgraph frequencies,” Proceedings of the 8th International Conference on Complex Networks (COMPLENET'17), 2017.

PDF file of the paper can be found at: http://www.dcc.fc.up.pt/~msilva/Pubs/complenet17.pdf

Usage:
1. make
2. make -f Makefile.1
3. ./genk -h

This Source Code uses the work developed in the following articles:

* ASONAM'2013 -http://asonam.cpsc.ucalgary.ca/ : Pedro Paredes and Pedro Ribeiro - Towards a Faster
Network-Centric Subgraph Census: http://dl.acm.org/citation.cfm?doid=2492517.2492535.
* (Submitted to) SNAM - https://www.springer.com/computer/database+management+%26+information+retrieval/journal/13278 : Pedro Paredes and Pedro Ribeiro - FaSE: Fast Exact and Approximate Subgraph Census

This software uses the nauty program version 2.4 by Brendan McKay. Therefore, nauty's
license restrictions apply to the usage of FaSE and GENK.
