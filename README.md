Heat-Modeling-Project


A Worked in pairs to create a program simulating heat flow and aging in a multicore computer processor.

Uses the Runge-Kutta (rk) method to solve a system of differential equations

To Run


compile with command 'gcc thermalModel.c -o thermalModel.out'

run with command "./thermalModel.out paraFile.txt powerFile.txt outputFile.txt"

parafile gives thermal resistances and capacitances

powerfile gives the power used by each core

output is printed to the outputfile, includes heat of each core

optional ambientinputfile which gives initial temperature of each core can be written
