Atomic models
=============

Atomic models have to be defined in  in order to carry out a
calculation. This section describes the model atom file in detail by
using the example ``helium.mod`` that is included in the present version
of Hazel. The atom model for He I is given by:
     
::

     2
     5
     1        0
                        0.00
     2        2
                        0.00
                       -0.987913
                       -1.064340
     3        0
                        0.00
     4        2
                        0.00
                       -0.270647
                       -0.292616
     5        4
                        0.00
                       -0.044187
                       -0.046722
     4
     1    1    2    1.022d7    10829.0911    1.0000000    1.0000000    0.0000000
     2    1    4    9.478d6    3888.6046    0.2000000    1.0000000    0.0000000
     3    2    3    2.780d7    7065.7085    1.0000000    1.0000000    0.0000000
     4    2    5    7.060d7    5875.9663    1.0000000    1.0000000    0.0000000

The first two numbers define the general properties of the atom. The
first line of the file is equal to :math:`2S`, where :math:`S` is the
value of the spin of the terms. In the example, :math:`S=1`. At present,
the code does not treat transitions between terms of different
multiplicity which are, otherwise, of reduced importance due to their
small transition probability. The second line contains the number of
terms included in the model atom. This example represents the triplet
system of He i with the lowest five terms, 2s\ :math:`^3`\ S,
3s\ :math:`^3`\ S, 2p\ :math:`^3`\ P, 3p\ :math:`^3`\ P and
3d\ :math:`^3`\ D


Then, the following lines define the term levels included in the model. The
information for each term consist of a line with an index (0,1,2,…) that
is used just to label each term and the value of :math:`2L`, where
:math:`L` is the value of the electronic orbital angular momentum. Then,
for each term, we must supply a list containing the energy separation in
cm\ :math:`^{-1}` between each :math:`J`-level and the level with the
smallest absolute value of :math:`J`. In case only one value of
:math:`J` is possible in the term, just put 0 in the energy difference.
 

Finally, the list of transitions has to be supplied. The first number
indicates the number of radiative transitions included in the model.
Then, the list contains the following numbers for each transition: index
number, index of lower level, index of upper level, Einstein coefficient
for spontaneous emission :math:`A_{ul}` of the transition, modification
factor :math:`f(\bar{n})`, modification factor :math:`f(w)` and value of
:math:`J^1_0/J^0_0`. The modification factors :math:`f(\bar{n})` and
:math:`f(w)` are multiplied by the mean number of photons per mode
:math:`\bar{n}` and the anisotropy factor :math:`w`, respectively. Since
 uses the value of :math:`\bar{n}` and :math:`w` calculated from the
tabulated solar CLV and taking into account geometrical effects, these
factors can be used to analyze the behavior of the emergent Stokes
profiles when, for some reason, the anisotropy or the intensity of the
radiation field is increased or decreased by an arbitrary factor.
Finally, if the radiation illuminating the atoms has non-zero net
circular polarization, it is possible to include its effect in the
statistical equilibrium equations by giving the value of
:math:`J^1_0/J^0_0`.