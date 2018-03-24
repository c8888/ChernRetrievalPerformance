#!/bin/sh

#PBS -N nyx

#PBS -q long

#PBS -l cput=72:00:00

#PBS -m abe

cd $HOME/ChernRetrievalPerformance/run/
math -script chernHistogram.m