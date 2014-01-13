.PHONY: all
all:	ftle_overview.pdf \
	lambda_lcs_full.pdf \
	strain_lcs.pdf \
	stretch_lcs.pdf \
	lambda_lcs_convergence_400.pdf \
	lambda_lcs_convergence_500.pdf \
	lambda_lcs_convergence_600.pdf

ftle_overview.pdf: ftle_overview.m
	matlab -nodisplay -r "ftle_overview, exit" >/dev/null

lambda_lcs_full.pdf: lambda_lcs_full.m
	matlab -nodisplay -r "lambda_lcs_full, exit" >/dev/null

strain_lcs.pdf: strain_lcs.m
	matlab -nodisplay -r "strain_lcs, exit" >/dev/null

stretch_lcs.pdf: stretch_lcs.m
	matlab -nodisplay -r "strain_lcs, exit" >/dev/null

lambda_lcs_convergence_400.pdf: lambda_lcs_convergence_subdomain.m
	matlab -nodisplay -r "lambda_lcs_convergence_subdomain, exit" >/dev/null

lambda_lcs_convergence_500.pdf: lambda_lcs_convergence_subdomain.m
	matlab -nodisplay -r "lambda_lcs_convergence_subdomain, exit" >/dev/null

lambda_lcs_convergence_600.pdf: lambda_lcs_convergence_subdomain.m
	matlab -nodisplay -r "lambda_lcs_convergence_subdomain, exit" >/dev/null

.PHONY: clean
clean:
	rm -f ftle_overview.pdf \
	lambda_lcs_full.pdf \
	strain_lcs.pdf \
	stretch_lcs.pdf \
	lambda_lcs_convergence_400.pdf \
	lambda_lcs_convergence_500.pdf \
	lambda_lcs_convergence_600.pdf
