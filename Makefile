# Create minimal bib file from $HOME/text/papers/main.bib:
# biber main.bcf --output_format=bibtex

.PHONY : FORCE_MAKE
all : main.pdf
%.pdf : %.tex FORCE_MAKE
	latexmk -norc -lualatex $<
clean :
	latexmk -norc -CA main
