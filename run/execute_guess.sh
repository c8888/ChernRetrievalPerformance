#!/bin/sh

#PBS -N nesoi

#PBS -q mp64

#PBS -l cput=1200:00:00

#PBS -m abe

cd $HOME/ChernRetrievalPerformance/run/
math -script guess.m