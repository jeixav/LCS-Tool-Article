.PHONY: all
all:	poincare_return_map.pdf \
	ftle.pdf \
	lambda_lcs_convergence_500.pdf \
	lambda_lcs_convergence_750.pdf \
	lambda_lcs_convergence_1000.pdf \
	lambda_strain_lcs.pdf \
	lambda_stretch_lcs.pdf

poincare_return_map.pdf: poincare_return_map.m
	matlab -nodisplay -r "poincare_return_map, exit" >/dev/null

ftle.pdf: ftle_.m
	matlab -nodisplay -r "ftle_, exit" >/dev/null

lambda_lcs_convergence_500.pdf: lambda_lcs_convergence.m
	matlab -nodisplay -r "lambda_lcs_convergence, exit" >/dev/null

lambda_lcs_convergence_750.pdf: lambda_lcs_convergence.m
	matlab -nodisplay -r "lambda_lcs_convergence, exit" >/dev/null

lambda_lcs_convergence_1000.pdf: lambda_lcs_convergence.m
	matlab -nodisplay -r "lambda_lcs_convergence, exit" >/dev/null

lambda_strain_lcs.pdf: lambda_strain_lcs.m
	matlab -nodisplay -r "lambda_strain_lcs, exit" >/dev/null

lambda_stretch_lcs.pdf: lambda_stretch_lcs.m
	matlab -nodisplay -r "lambda_stretch_lcs, exit" >/dev/null

.PHONY: clean
clean:
	rm -f poincare_return_map.pdf \
	ftle_pdf \
	lambda_lcs_convergence_500.pdf \
	lambda_lcs_convergence_750.pdf \
	lambda_lcs_convergence_1000.pdf \
	lambda_strain_lcs.pdf \
	lambda_stretch_lcs.pdf
