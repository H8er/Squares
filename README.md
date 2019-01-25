# Squares
circle.cpp - functions creates circle in (0,0)
easyf.cpp - ring in (0,0)
figure2.cpp - figure from pdf
parallel_figure.cu - crap
parallel_grid.cu - 1st version. w/o borders. computations on grid size
parallel_grid_v2.cu - 2nd version. w/o borders. computations on threadIds
parallel_grid_v3.cu - 3rd. with borders. computations by groups on each thread. keeps coords of each square. need to sync accuracy and gridsize
parallel_grid_v4_ind.cu - keeps indexes of squares instead of 3rd. in cout computation of coords.
