.PHONY : FORCE_MAKE
all : main.pdf
%.pdf : %.tex FORCE_MAKE
	latexmk -norc -pdf $<
clean :
	latexmk -norc -CA main

# Create BibTeX .bib file containing only cited references of a bigger file
# bibexport -o main_extracted.bib main.aux
