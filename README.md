2D CDFFAR



the input is a RDM Matric

the fist to for-loops implement the CUT iterator for the RDM matrix, the iterators consider the bound so that the side Training and Guard Cells are skipped.



the next to for-loops implement the moving window relative to the CUT, where it sum every taring data cells (as power), and counts the amount of Training cells.

after the moving window is ready it compute the average value, convert to db add the offset. IF  the result is bigger then the CUT it writes a 0 in the Result Matrix otherwise a 1.

After the filter is ready it zeros every cell that was not tested (was not a CUT) . 



