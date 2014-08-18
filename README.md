bar1_project
============

Source code for the the publication [Yeast Mating and Image-Based Quantification of Spatial Pattern Formation](http://dx.doi.org/10.1371/journal.pcbi.1003690).
Since it is a pipeline it comprises various individual tasks which are melt together by several scripts. Here you will find the entire
source code as well as directions to additional programs in order to reproduce the published results.

Some of the code is already documented, and we are working on documenting it much betters as it currently is. However, in order
to fully understand all of the decision we made in the implementation we also recommend to read the Methods section of the publication and 
particularly the [Supporting Information](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003690#s5). 

The pipeline is summarized in the following parts.

Image analysis and cell tracking
================================

**found in the cellid subfolder**

For details on image acquisition and preprocessing please see our [Supplemental Text](http://dx.doi.org/10.1371/journal.pcbi.1003690.s002). Cells are tracked using
[CellID and VCellID](http://dx.doi.org/10.1002/0471142727.mb1418s100), where we modified the original version of CellID in order to
give additional output such as the location of membrane pixels and interior pixel-wise fluorescence values. For your convenience you will find our 
modified CellID code as well as the origianl VCellID source here. 

Installation
------------

**Requisites**
* libtiff


Automatic coupling of tracked cells to models
=============================================

**is mesh_tools.R**

Most of the actual *glueing* part is done using [R](http://www.r-project.org). We created a convenience library in R to analyze CellID output
and couple it to our modeling and analysis pipeline.

Solving and fitting model using DUNE and PDElab
===============================================

**found in the dune sub-folder**

Work tasks
==========

**found in the tasks subfolder**
