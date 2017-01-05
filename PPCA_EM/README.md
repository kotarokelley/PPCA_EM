# Project Name

PPCA_EM

## Synopsis
An implementation of PPCA as derived in Tipping, Bishop, 2006 Neural Computation, to unsupervised classification of unaligned 2-D electron microscopy images. Still a work in progress.

## Motivation


Many methods for classifying single particle 2-D electron microscopy images are broken up into two steps. To globally align all particles, a synthetic fiduciary or randomly selected particles from the data set are used as references. This aligned data set is then subject to a classification procedure. As argued in Zhao and Singer 2014 Journal of Structural Biology, this procedure may results in unsatisfactory results whenever there are many different views; owing to the so called "hairy ball theorem". Therefore, a method that avoids the alignment process and directly classifies the particles may be useful. 

## General Idea of PPCA_EM algorithm
Probabalistic principal component analysis (PPCA) as described in Tipping, Bishop, 2006 Neural Computation allows many local PPCA models, each with its own mean, to be fit to the data by an expectation maximization procedure. Each data point or particle in this case contribues to each local model through a prior probabiliby density, which can be used for classification, for example by its maximum. In this way, a soft classification criterion can be used. Instead of an initial alignment procedure, the original unaligned data set will be expanded with rotated and translated copies. A tradeoff between accuracy and speed is controlled by a user supplied number of components to keep and local models to assign. 

## Current State of Project
The PPCA expectation maximization algorithm as described in Tipping, Bishop, 2006 Neural Computation has been written and tested. The aligned zip code data set from US postal is being used as an initial test set. As expected, matrix inversions are a major bottleneck in terms of running time as it is O(n^3), where n is the number of particles. Before proceeding to an expanded data set, the use of approximate inverses that take less time, are being explored. Parallelizing parameter updates for each particle for each iteration of expectation maximization will need to be worked out as well. A rotationally and translationally preturbed zip code data set will then be used as a test set. To rotationally and translationally expand the input data, linear interpolation will be used. The makefiles works for OSX, but not for linux machines. To run tests on a faster 8core linux machine, a separate makefile is used. A better solution would be to cross-compile. 

## Future Plans
Once the algorithm is tested and sped up, a synthetic single particle EM data set will be used as a test set. This will allow investigation of the effects of elongated shape, flexibility, and noise on the stability of the algorithm. Through simulations, a heuristic for choosing a good number for the components to keep and local models to assign will be investigated. 

## Contributing

1. Please suggest any improvements, feature requests. 
2. Please notify of any bugs. 
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D