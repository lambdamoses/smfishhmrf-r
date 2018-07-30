\name{smfishHmrf-package}
\alias{smfishHmrf-package}
\docType{package}
\title{Single Molecule FISH Spatial Pattern Discovery}
\description{
  Use HMRF to dissect spatial patterns in smFISH data such as seqFISH.
}
\section{Introduction}{This package provides tools for spatial pattern discovery
  using multivariate Gaussian models and hidden Markov random field fitted by 
  expectation-maximization.

  smFISH, such as recently popularized seqFISH approach, is used to spatially
  profile cells at single-cell resolution. Using sequential multiplexed barcoding
  schemes, seqFISH can profile hundreds of genes' transcript count for 
  individual cells with high sensitivity. 

  Using gene expression from seqFISH, we can infer the spatial pattern that might
  exist in the spatial profile. This is analogous to image segmentation in computer
  vision. In the single-cell imaging, spatial pattern looks like a contiguous
  domain of cells, characterized by the expression of spatially dependent genes, and
  may suggest something about the cell's local environment. 
  This spatial state of the cell is unknown and need to be identified from seqFISH data.
  
  As the spatial states are assumed to be independent, we can use mixture of multivariate
  Gaussian distribution to model the observed expression with parameters estimated by the EM
  algorithm. Since nearby cells tend to be of the same state, a Markov
  random field model (we used Potts' model for our case) can be used to capture 
  the spatial similarity of cells by making homogeneous the relationship 
  between neighboring cells.

  As the typical package for doing HMRF on magnetic resonance imaging (MRI) data is 
  restricted to 1-channel (grayscale) and fixed 4-neighbor regular grid (of pixels), 
  the model in smfishHmrf is more general. First it extends simple 1D Gaussian to
  multivariate Gaussian with multidimensional mean and covariance matrices estimated by
  EM. Second, as the cells do not fall on a regular grid, we adapt to a local neighborhood
  graph instead of a grid. This allows for a more flexible specification of number of neighbors 
  and neighbor structure. To illustrate the method, we have tested it on 
  mouse brain visual cortex seqFISH imaged cells. 

}

\keyword{package}
