file: hw1.tex
	@pdflatex hw1.tex > latex_temp.txt
	@makeindex hw1.nlo -s nomencl.ist -o hw1.nls > latex_temp.txt
	@pdflatex hw1.tex > latex_temp.txt
	@pdflatex hw1.tex > latex_temp.txt
	@cp -r *.log logs/.
	@rm *.aux *.ilg *.lof *.nlo *.nls *.out *.toc latex_temp.txt *.log *.fls *.fdb_latexmk

clean_after_vscode:
	@test ! $(wildcard *.log) || cp -r *.log logs/.
	@rm *.fdb_latexmk *.fls *.synctex.gz *.aux *.bbl *.log *.blg *.toc *.out *.lof *.nlo

test:
	@test ! $(wildcard *.pdf) || echo File does exist