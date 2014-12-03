.PHONY : FORCE_MAKE
all : main.pdf
%.pdf : %.tex FORCE_MAKE
	latexmk -norc -lualatex -bibtex $<
clean :
	latexmk -norc -lualatex -bibtex -CA main

# Create BibTeX .bib file containing only cited references of a bigger file
# biber main.bcf --output_format=bibtex
