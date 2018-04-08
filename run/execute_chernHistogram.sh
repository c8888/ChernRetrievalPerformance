#!/bin/sh

#PBS -N nyx

#PBS -q normal

#PBS -l cput=12:00:00

#PBS -m abe

cd $HOME/ChernRetrievalPerformance/run/
math -script chernHistogram.m