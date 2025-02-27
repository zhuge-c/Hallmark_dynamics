function df_area=distribution_cross_area(nor_xi,nor_pdf,cancer_xi,cancer_pdf)

[left_p,lindex]=min(cancer_xi);
[right_p,rindex]=max(nor_xi);

nor_xi_common=nor_xi(nor_xi>=left_p);
nor_yi_common=nor_pdf(nor_xi>=left_p);

cancer_xi_common=cancer_xi(cancer_xi<=right_p);
cancer_yi_common=cancer_pdf(cancer_xi<=right_p);

common_range = linspace(min([nor_xi_common, cancer_xi_common]), max([nor_xi_common, cancer_xi_common]), 100);
nor_interp = interp1(nor_xi_common, nor_yi_common, common_range, 'pchip');
cancer_interp = interp1(cancer_xi_common, cancer_yi_common, common_range, 'pchip');

cross_area_pdf = min(nor_interp, cancer_interp);
df_area= trapz(common_range, cross_area_pdf);

