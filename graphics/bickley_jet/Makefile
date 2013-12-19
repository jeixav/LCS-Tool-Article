.PHONY: all
all:	ftle_overview.pdf \
	lambda_lcs_full.pdf \
	strain_lcs \
	stretch_lcs

ftle_overview.pdf: ftle_overview.m
	matlab -nodisplay -r "ftle_overview, exit" >/dev/null

lambda_lcs_full.pdf: lambda_lcs_full.m
	matlab -nodisplay -r "lambda_lcs_full, exit" >/dev/null

strain_lcs.pdf: strain_lcs.m
	matlab -nodisplay -r "strain_lcs, exit" >/dev/null

stretch_lcs.pdf: stretch_lcs.m
	matlab -nodisplay -r "strain_lcs, exit" >/dev/null

.PHONY: clean
clean:
	rm -f ftle_overview.pdf \
	lambda_lcs_full.pdf \
	strain_lcs \
	stretch_lcs
