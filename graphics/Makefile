all: main_aux_grids.pdf variable_step_integration.pdf

main_aux_grids.pdf : main_aux_grids.svg
	inkscape -z -D --file=$< --export-latex --export-pdf=$@

variable_step_integration.pdf : variable_step_integration.svg
	inkscape -z -D --file=$< --export-latex --export-pdf=$@

.PHONY: clean
clean:
	rm main_aux_grids.pdf main_aux_grids.pdf_tex variable_step_integration.pdf variable_step_integration.pdf_tex
