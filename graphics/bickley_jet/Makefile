.PHONY: all
all:	hyperbolic_shear_lcs_strainline.pdf \
	hyperbolic_shear_lcs_backward.pdf \
	hyperbolic_shear_lcs_details_strainline.pdf \
	hyperbolic_shear_lcs_details_stretchline.pdf

hyperbolic_shear_lcs_strainline.pdf: hyperbolic_shear_lcs_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_pdf, catch, end, exit" >/dev/null

hyperbolic_shear_lcs_stretchline.pdf: hyperbolic_shear_lcs_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_pdf, catch, end, exit" > /dev/null

hyperbolic_shear_lcs_details_strainline.pdf: hyperbolic_shear_lcs_details_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_details_pdf, catch, end, exit" >/dev/null

hyperbolic_shear_lcs_details_stretchline.pdf: hyperbolic_shear_lcs_details_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_details_pdf, catch, end, exit" >/dev/null

.PHONY: clean
clean:
	rm -f hyperbolic_shear_lcs_strainline.pdf \
	hyperbolic_shear_lcs_stretchline.pdf \
	hyperbolic_shear_lcs_details_strainline.pdf \
	hyperbolic_shear_lcs_details_stretchline.pdf \
	hyperbolic_shear_lcs_strainline.fig \
	hyperbolic_shear_lcs_stretchline.fig
