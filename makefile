all:
	pdflatex thesis
	biber thesis
	pdflatex thesis
	pdflatex thesis
	#dvips paper1
	#ps2pdf paper1
	#rm *-eps-converted-to.pdf
	rm *.{aux,bbl,blg,d,dvi,log,llt,toc,bcf,out,run.xml}
