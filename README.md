# An Exact and Polynomial Approach for a Bi-Objective Integer Programming Problem Regarding Network Flow Routing -- Source code

This is the implementation of an algorithm based on the ϵ-constraint approach for finding a minimal complete set of Pareto-optimal solutions for the biojective problem (P) in which the first (bottleneck) objective function minimizes the load in the most congested links and the second objective function minimizes the total cost of routing the flows. In each iteration of the algorithm, a single-objective problem, (P_ϵ), is solved. The integer and positive parameter ϵ is used to control the values of the bottleneck objective function and to ensure the optimality of the algorithm.


## BibTex

@Inbook{hansen:80,

    author={Hansen, Pierre},
    title={{Bicriterion Path Problems}},
    bookTitle={Multiple Criteria Decision Making Theory and Application: Proceedings of the Third Conference Hagen/K{\"o}nigswinter, West Germany, August 20--24, 1979},
    year={1980},
    publisher={Springer Berlin Heidelberg},
    address={Berlin, Heidelberg},
    pages={109--127},
    OPTdoi={10.1007/978-3-642-48782-8_9}
}


@Article{climaco:82,

    author = {J. C. N. Cl\'{i}maco and E. Q. V. Martins},
    title = {{A bicriterion shortest path algorithm}},
    journal = {European Journal of Operational Research},
    year = {1982},
    volume = {11},
    number = {4},
    pages = {399--404},
    OPTdoi = {10.1016/0377-2217(82)90205-3}
}

@Article{klein:82,

    author = {Dieter Klein and Edward Hannan},
    title = {An algorithm for the multiple objective integer linear programming problem},
    journal = {European Journal of Operational Research},
    volume = {9},
    number = {4},
    pages = {378--385},
    year = {1982},
    OPTdoi = {10.1016/0377-2217(82)90182-5}
}

@Article{chalmet:86,

    author = {L.G. Chalmet and L. Lemonidis and D.J. Elzinga},
    title = {An algorithm for the bi-criterion integer programming problem},
    journal = {European Journal of Operational Research},
    volume = {25},
    number = {2},
    pages = {292--300},
    year = {1986},
    OPTdoi = {10.1016/0377-2217(86)90093-7}
}

@Book{schrijver:86,

    author = {Schrijver, Alexander},
    title = {{Theory of Linear and Integer Programming}},
    year = {1986},
    publisher = {John Wiley \& Sons, Inc.},
    address = {New York, NY, USA},
    OPTisbn = {0-471-90854-1}
} 


@Article{berman:90,

    author = {Berman, O. and Einav, D. and Handler, G.},
    title = {{The Constrained Bottleneck Problem in Networks}},
    journal = {Operations Research},
    year = {1990},
    volume = {38},
    number = {1},
    pages = {178--181},
    OPTdoi = {10.1287/opre.38.1.178}
 } 

@Article{barabasi:99,

	author = {Barab{\'a}si, Albert-L{\'a}szl{\'o} and Albert, R{\'e}ka},
	title = {Emergence of Scaling in Random Networks},
	journal = {Science},
	year = {1999},
	volume = {286},
	number = {5439},
	pages = {509--512},
	OPTpublisher = {American Association for the Advancement of Science},
	OPTdoi = {10.1126/science.286.5439.509}
}

@Article{pham:04,

    author = {Peter P. Pham and Sylvie Perreau},
    title = {{Increasing the network performance using multi-path routing mechanism with load balance}},
    journal = {Ad Hoc Networks},
    year = {2004},
    volume = {2},
    number = {4},
    pages = {433--459},
    OPTdoi = {10.1016/j.adhoc.2003.09.003}
}

@book{ehrgott2005multicriteria,

    title={Multicriteria optimization},
    author={Ehrgott, Matthias},
    volume={491},
    year={2005},
    publisher={Springer Science \& Business Media}
}

@Article{ehrgott:06,

    author={Ehrgott, Matthias},
    title={A discussion of scalarization techniques for multiple objective integer programming},
    journal={Annals of Operations Research},
    year={2006},
    volume={147},
    number={1},
    pages={343--360},
    OPTdoi={10.1007/s10479-006-0074-z}
}

@book{chankong2008multiobjective,

    title={Multiobjective decision making: theory and methodology},
    author={Chankong, Vira and Haimes, Yacov Y},
    year={2008},
    publisher={Courier Dover Publications}
}

@Article{bornstein:12,

    author = {Cl\'{a}udio T. Bornstein and Nelson Maculan and Marta Pascoal and Leizer L. Pinto},
    title = {{Multiobjective combinatorial optimization problems with a cost and several bottleneck objective functions: An algorithm with reoptimization}} ,
    journal = {Computers \& Operations Research},
    year = {2012},
    volume = {39},
    number = {9},
    pages = {1969--1976},
    OPTdoi = {10.1016/j.cor.2011.09.006}
}

@Article{liu:12, 

    author={Q. Liu and J. Yin and V. C. M. Leung and Z. Cai}, 
    title={{ISAR: Improved Situation-Aware Routing Method for Wireless Mesh Backbones}},
    journal={IEEE Communications Letters}, 
    year={2012},
    month={Sep},
    volume={16},
    number={9},
    pages={1404--1407},
    OPTdoi={10.1109/LCOMM.2012.072012.120502}
}

@Article{galvez:13,

	author   = {Juan J. G\'{a}lvez and Pedro M. Ruiz},
	title    = {{Efficient Rate Allocation, Routing and Channel Assignment in Wireless Mesh Networks Supporting Dynamic Traffic Flows}},
	journal  = {Ad Hoc Networks},
	year     = {2013},
	volume   = {11},
	number   = {6},
	pages    = {1765--1781},
	OPTmonth = {August},
	OPTdoi   = {10.1016/j.adhoc.2013.04.002}
}

@Article{lokman:13,

     author={Lokman, Banu and K{\"o}ksalan, Murat},
     title={Finding all nondominated points of multi-objective integer programs},
     journal={Journal of Global Optimization},
     year={2013},
     volume={57},
     number={2},
     pages={347--365},
     OPTdoi={10.1007/s10898-012-9955-7}
}

@Article{bhushan:14, 

     author={N. Bhushan and J. Li and D. Malladi and R. Gilmore and D. Brenner and A. Damnjanovic and R. T. Sukhavasi and C. Patel and S. Geirhofer},
     title={{Network densification: the dominant theme for wireless evolution into 5G}},
     journal={IEEE Communications Magazine},
     year={2014},
     month={Feb},
     volume={52},
     number={2},
     pages={82--89},
     OPTdoi={10.1109/MCOM.2014.6736747}
}

@Article{ozlen:14,

     author={Ozlen, Melih and Burton, Benjamin A. and MacRae, Cameron A. G.},
     title={{Multi-Objective Integer Programming: An Improved Recursive Algorithm}},
     journal={Journal of Optimization Theory and Applications},
     year={2014},
     volume={160},
     number={2},
     pages={470--482},
     OPTdoi={10.1007/s10957-013-0364-y}
}

@Article{boland:16,

     author={Boland, Natashia and Charkhgard, Hadi and Savelsbergh, Martin},
     title={{The L-shape search method for triobjective integer programming}},
     journal={Mathematical Programming Computation},
     year={2016},
     volume={8},
     number={2},
     pages={217--251},
     OPTdoi={10.1007/s12532-015-0093-3}
}

@Article{mello:16,

    author = {Micael O.M.C. de Mello and Vinicius C.M. Borges and Leizer L. Pinto and Kleber V. Cardoso},
    title = {Improving load balancing, path length, and stability in low-cost wireless backhauls},
    journal = {Ad Hoc Networks},
    year = {2016},
    volume = {48},
    number = {},
    pages = {16--28},
    OPTdoi = {10.1016/j.adhoc.2016.05.002}
}

@article{gadegaard2016bi,

    title={A bi-objective approach to discrete cost-bottleneck location problems},
    author={Gadegaard, Sune Lauth and Klose, Andreas and Nielsen, Lars Relund},
    journal={Annals of Operations Research},
    pages={1--23},
    year={2016},
    publisher={Springer}
}

@book{ahuja2017network,

    title={Network flows: theory, algorithms, and applications},
    author={Ahuja, Ravindra K},
    year={2017},
    publisher={Pearson Education}
}
