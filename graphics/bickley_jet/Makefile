.PHONY: all
all:	hyperbolic_shear_lcs_strainline.pdf \
	hyperbolic_shear_lcs_stretchline.pdf \
	hyperbolic_shear_lcs_details_strainline.pdf \
	hyperbolic_shear_lcs_details_stretchline.pdf \
	ftle_100_39.pdf \
	ftle_250_99.pdf \
	ftle_500_199.pdf \
	ftle_1000_398.pdf

hyperbolic_shear_lcs_strainline.pdf: hyperbolic_shear_lcs_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_pdf, catch, end, exit" >/dev/null

hyperbolic_shear_lcs_stretchline.pdf: hyperbolic_shear_lcs_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_pdf, catch, end, exit" > /dev/null

hyperbolic_shear_lcs_details_strainline.pdf: hyperbolic_shear_lcs_details_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_details_pdf, catch, end, exit" >/dev/null

hyperbolic_shear_lcs_details_stretchline.pdf: hyperbolic_shear_lcs_details_pdf.m
	matlab -nodisplay -r "try, hyperbolic_shear_lcs_details_pdf, catch, end, exit" >/dev/null

ftle_100_39.pdf: ftle_resolutions.m
	matlab -nodisplay -r "try, ftle_resolutions, catch, end, exit" >/dev/null

ftle_250_99.pdf: ftle_resolutions.m
	matlab -nodisplay -r "try, ftle_resolutions, catch, end, exit" >/dev/null

ftle_500_199.pdf: ftle_resolutions.m
	matlab -nodisplay -r "try, ftle_resolutions, catch, end, exit" >/dev/null

ftle_1000_398.pdf: ftle_resolutions.m
	matlab -nodisplay -r "try, ftle_resolutions, catch, end, exit" >/dev/null

.PHONY: clean
clean:
	rm -f hyperbolic_shear_lcs_strainline.pdf \
	hyperbolic_shear_lcs_stretchline.pdf \
	hyperbolic_shear_lcs_details_strainline.pdf \
	hyperbolic_shear_lcs_details_stretchline.pdf \
	hyperbolic_shear_lcs_strainline.fig \
	hyperbolic_shear_lcs_stretchline.fig \
	ftle_100_39.pdf \
	ftle_250_99.pdf \
	ftle_500_199.pdf \
	ftle_1000_398.pdf
