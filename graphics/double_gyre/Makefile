.PHONY: all
all:	hyperbolic_shear_lcs_forward.pdf \
	hyperbolic_shear_lcs_backward.pdf \
	hyperbolic_shear_lcs_details_forward.pdf \
	hyperbolic_shear_lcs_details_backward.pdf

hyperbolic_shear_lcs_forward.pdf: hyperbolic_shear_lcs_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_pdf, catch, end, exit" >/dev/null

hyperbolic_shear_lcs_backward.pdf: hyperbolic_shear_lcs_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_pdf, catch, end, exit" > /dev/null

hyperbolic_shear_lcs_details_forward.pdf: hyperbolic_shear_lcs_details_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_details_pdf, catch, end, exit" >/dev/null

hyperbolic_shear_lcs_details_backward.pdf: hyperbolic_shear_lcs_details_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_details_pdf, catch, end, exit" >/dev/null

.PHONY: clean
clean:
	rm -f hyperbolic_shear_lcs_forward.pdf \
	hyperbolic_shear_lcs_backward.pdf \
	hyperbolic_shear_lcs_details_forward.pdf \
	hyperbolic_shear_lcs_details_backward.pdf \
	hyperbolic_shear_lcs_forward.fig \
	hyperbolic_shear_lcs_backward.fig
