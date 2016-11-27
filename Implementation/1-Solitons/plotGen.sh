g++ -std=c++11 -Wall -Werror solitonV1-Tex.cpp -o solitonV1-Tex.x 

./solitonV1-Tex.x | gnuplot 

latex SolitonP1_3D.tex
dvipdf SolitonP1_3D.dvi SolitonP1_3D.pdf

mv SolitonP1_3D.pdf imgs/SolitonP1_3D.pdf
rm SolitonP1_3D.tex SolitonP1_3D.log SolitonP1_3D.aux SolitonP1_3D-inc.eps SolitonP1_3D.dvi

latex SolitonP1_Contour.tex
dvipdf SolitonP1_Contour.dvi SolitonP1_Contour.pdf

mv SolitonP1_Contour.pdf imgs/SolitonP1_Contour.pdf
rm SolitonP1_Contour.tex SolitonP1_Contour.aux SolitonP1_Contour.log SolitonP1_Contour-inc.eps SolitonP1_Contour.dvi
