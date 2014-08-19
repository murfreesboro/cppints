#!/bin/sh
##
## create.sh
## 
## Made by liufenglai
## Login   <liufenglai@phoenix>
## 
## Started on  Fri May 29 07:46:28 2009 liufenglai
## Last update Fri May 29 07:49:36 2009 liufenglai
##
##
## used to create the reference.bib from the ../../../ref folder
##
## 

currentDir=`pwd`
refDir=/home/fenglai/personal/fenglai-cvs-archive/ref
cd $refDir
cat *.bib > $currentDir/ref.bib
