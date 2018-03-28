#!/bin/sh

#PBS -N phanes

#PBS -q mp16

#PBS -l cput=900:00:00

#PBS -m abe

cd $HOME/ChernRetrievalPerformance/run/
math -script startHint.m