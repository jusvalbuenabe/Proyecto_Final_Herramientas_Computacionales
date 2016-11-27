g++ -std=c++11 -Wall 2soliton-Tex.cpp -o 2soliton-Tex.x 

./2soliton-Tex.x | gnuplot 

latex SolitonP2_3D.tex
dvipdf SolitonP2_3D.dvi SolitonP2_3D.pdf

mv SolitonP2_3D.pdf imgs/SolitonP2_3D.pdf
rm SolitonP2_3D.tex SolitonP2_3D.log SolitonP2_3D.aux SolitonP2_3D-inc.eps SolitonP2_3D.dvi

latex SolitonP2_Contour.tex
dvipdf SolitonP2_Contour.dvi SolitonP2_Contour.pdf

mv SolitonP2_Contour.pdf imgs/SolitonP2_Contour.pdf
rm SolitonP2_Contour.tex SolitonP2_Contour.aux SolitonP2_Contour.log SolitonP2_Contour-inc.eps SolitonP2_Contour.dvi
