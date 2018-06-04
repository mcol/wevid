##=============================================================================
##
## Copyright (c) 2018 Paul McKeigue
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##=============================================================================

##
## wevid.R
##
## Package documentation.
##

#' Quantifying performance of a diagnostic test using the sampling distribution of the weight of evidence
#' favouring case over noncase status
#'
#' This package provides functions for quantifying the performance of a diagnostic test
#' (or any other binary classifier) by calculating and plotting the distributions in cases
#' and noncases of the weight of evidence favouring case over noncase status.
#' 
#' To use it, you should have computed on a test dataset (or on test folds used for
#' cross-validation:
#' 
#' 1. The prior probability of case status (this may be just the frequency of cases in the
#' training data.
#' 
#' 2. The posterior probability of case status (using the model learned on the training data
#' to predict on the test data)
#' 
#' 3. The observed case status (coded as 0=noncase, 1=case). 
#'
#' @author
#' Paul McKeigue \email{paul.mckeigue@@ed.ac.uk}
#'
#' Citation for the statistical methods used in this package:
#' McKeigue P. Quantifying performance of a diagnostic test as the expected information
#' for discrimination: relation to the C-statistic.
#' Statistical Methods for Medical Research 2018, in press.
#' 
#' @docType package
#' @import ggplot2
"_PACKAGE"
