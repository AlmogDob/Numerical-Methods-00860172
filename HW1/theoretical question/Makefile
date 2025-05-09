file: NM_hw1_TQ.tex
	@pdflatex NM_hw1_TQ.tex > latex_temp.txt
	@makeindex NM_hw1_TQ.nlo -s nomencl.ist -o NM_hw1_TQ.nls > latex_temp.txt
	@pdflatex NM_hw1_TQ.tex > latex_temp.txt
	@pdflatex NM_hw1_TQ.tex > latex_temp.txt
	@cp -r *.log logs/.
	@rm *.aux *.ilg *.lof *.nlo *.nls *.out *.toc latex_temp.txt *.log *.fls *.fdb_latexmk *.lol

clean_after_vscode:
	@test ! $(wildcard *.log) || cp -r *.log logs/.
	@rm *.fdb_latexmk *.fls *.synctex.gz *.aux *.bbl *.log *.blg *.toc *.out *.lof *.nlo *.lol

test:
	@test ! $(wildcard *.pdf) || echo File does exist