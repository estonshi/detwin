/*
 * detwin.c
 *
 * Assemble and process FEL Bragg intensities
 * Process and deal with indexing ambiguity
 *
 * Copyright © 2014-2019 Beijing Computational Science Research Center (CSRC)
 * Copyright © 2012-2016 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2015      Keitaro Yamashita <k.yamashita@spring8.or.jp>
 *   2009-2016 Thomas White <taw@physics.org>
 *   2011      Andrew Martin <andrew.martin@desy.de>
 *   2012      Lorenzo Galli <lorenzo.galli@desy.de>
 *   2014      Chunhong Yoon <chun.hong.yoon@desy.de>
 *   2014      Haiguang Liu <hgliu@csrc.ac.cn>
 *   2019      Yingchen Shi <shiyc12@csrc.ac.cn> 
 *
 * This file is a 3rd-party patch for CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdarg.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>

#include "utils.h"
#include "statistics.h"
#include "reflist-utils.h"
#include "symmetry.h"
#include "stream.h"
#include "reflist.h"
#include "image.h"
#include "crystal.h"
#include "thread-pool.h"
#include "geometry.h"
#include "cell-utils.h"

#include "detwin.h"

#include <math.h>
#include <gsl/gsl_statistics.h>

#define epsilon 1.e-15

static int N_TWINS = 0;
static int space_group_num = -1;
static int ADD_N_TWINS = 0;
static char *added_operators[99][3];

struct image_sm
{
	int n_ref;
	unsigned int *serial;
	float *intensity;
	float *sigma;
	UnitCell* cell;
};

static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"read from stream file and convert the data to miller indexed FEL Bragg intensities.\n"
"output all crystals HKL, even numbered crystals HKL, odd numbered crystals HKL\n"
"\n"
"  -h, --help                Display this help message.\n"
"      --version             Print CrystFEL version number and exit.\n"
"  -i, --input=<filename>    Specify input stream or hkl filename (\"-\" for stdin).\n"
"  -o, --output=<filename>   Specify output hkl filename for merged intensities\n"
"                             Default: processed.hkl).\n"
"      --stat=<filename>     Specify output filename for merging statistics.\n"
"  -y, --symmetry=<sym>      Merge according to point group <sym>.\n"
"  -k, --spacegroupNum=<k>   Use space group number to determine twinning operator.\n"
"  -m, --max-niter=<t>       Number of iterations for de-twinning, default is 30\n"
"\n"
"  -s  --start-after=<n>     Skip <n> crystals at the start of the stream.\n"
"  -f  --stop-after=<n>      Stop after merging <n> crystals.\n"
"  -g, --histogram=<h,k,l>   Calculate the histogram of measurements for this\n"
"                             reflection.\n"
"  -z, --hist-parameters     Set the range for the histogram and the number of\n"
"          =<min,max,nbins>   bins. \n"
"\n"
"      --scale               Scale each pattern for best fit with the current\n"
"                             model. Operated before detwining.\n"
"      --winner-takes-all    If set, among all possible twinning modes, only the one \n"
"                             with higest cc to reference model will merge.\n"
"                             Default: all twinning modes merge in a weighted way\n"
"      --highres=<n>         Reject reflections with resolution (A) higher than n \n"
"                            while calculating CC in detwinning, default is Inf.\n"
"      --lowres=<n>          Reject reflections with resolution (A) lower than n \n"
"                            while calculating CC in detwinning, default is 0.\n"
"      --no-polarisation     Disable polarisation correction.\n"
"      --cc-only             Only calculate and display CC between -i and -o hkl file\n"
"      --min-measurements=<n> Require at least <n> measurements before a\n"
"                             reflection appears in the output.  Default: 2\n"
"      --min-snr=<n>         Require individual intensity measurements to\n"
"                             have I > n * sigma(I).  Default: -infinity.\n"
"      --min-cc=<n>          Reject frames with CC less than n. Default: infinity.\n"
"      --max-adu=<n>         Maximum peak value.  Default: infinity.\n"
"      --min-res=<n>         Merge only crystals which diffract above <n> A.\n"
"      --push-res=<n>        Integrate higher than apparent resolution cutoff.\n"
"      --write-assignments   Write reindexed results of the crystals in original stream\n"
"              =<filename>   file to filename.\n"
"      --add-operators=<op>  Add twinning operators manually, such as '--add-operators=\n"
"                            -h,-k,-l/k,h,l'. DO NOT give space between characters !\n"
"      --write-stream=<fn>   Write reindexed stream file to fn\n"
);
}


static int get_crystals_num(char *filename)
{
	int num = 0;
	char line[1024];
	char *rval = NULL;

	FILE *fh;
	fh = fopen(filename, "r");

	do{
		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) break;
		if (strcmp(line, CRYSTAL_START_MARKER"\n") == 0) num++;
	} while(1);

	fclose(fh);

	return num;
}


static void parse_operators(char *operator_str)
{
	int num = 0, snum = 0, i = 0;
	char *str_copy, *pch, *spch;
	char *op_tmp[99];
	str_copy = strdup(operator_str);

	pch = strtok(str_copy, "/");
	while(pch != NULL){
		op_tmp[num] = pch;
		num += 1;
		pch = strtok(NULL, "/");
	}

	ADD_N_TWINS = num;
	N_TWINS += num;

	for( i=0; i<ADD_N_TWINS; i++){
		pch = op_tmp[i];
		snum = 0;
		spch = strtok(pch, ",");
		while(spch != NULL){
			if(snum >= 3) ERROR("Value of option --add-operators is invalid");
			if(strcmp(spch,"h") && strcmp(spch,"-h") &&
			   strcmp(spch,"k") && strcmp(spch,"-k") &&
			   strcmp(spch,"l") && strcmp(spch,"-l"))
				ERROR("Value of option --add-operators is invalid");
			added_operators[i][snum] = spch;
			snum += 1;
			spch = strtok(NULL, ",");
		}
	}
}


static void plot_histogram(double *vals, int n, float hist_min, float hist_max,
                           int nbins)
{
	int i;
	double max = -INFINITY;
	double min = +INFINITY;
	double step;
	int histo[nbins];
	FILE *fh;

	fh = fopen("histogram.dat", "w");
	if ( fh == NULL ) {
		ERROR("Couldn't open 'histogram.dat'\n");
		return;
	}

	if ( hist_min == hist_max ) {
		for ( i=0; i<n; i++ ) {
			if ( vals[i] > max ) max = vals[i];
			if ( vals[i] < min ) min = vals[i];
		}
	} else {
		min = hist_min;
		max = hist_max;
	}
	STATUS("min max nbins: %f %f %i\n", min, max, nbins);
	min--;  max++;

	for ( i=0; i<nbins; i++ ) {
		histo[i] = 0;
	}

	step = (max-min)/nbins;

	for ( i=0; i<n; i++ ) {
		int bin;
		if ( (vals[i] > min) && (vals[i] < max) ) {
			bin = (vals[i]-min)/step;
			histo[bin]++;
		}
	}

	for ( i=0; i<nbins; i++ ) {
		fprintf(fh, "%f %i\n", min+step*i, histo[i]);
	}

	fclose(fh);
}


static double scale_intensities(RefList *reference, RefList *new,
                              const SymOpList *sym)
{
        double s;
        double top = 0.0;
        double bot = 0.0;
        Reflection *refl;
        RefListIterator *iter;

        for ( refl = first_refl(new, &iter);
              refl != NULL;
              refl = next_refl(refl, iter) )
        {

                double i1, i2;
                signed int hu, ku, lu;
                signed int h, k, l;
                Reflection *reference_version;

                get_indices(refl, &h, &k, &l);
                get_asymm(sym, h, k, l, &hu, &ku, &lu);

                reference_version = find_refl(reference, hu, ku, lu);
                if ( reference_version == NULL ) continue;

                i1 = get_intensity(reference_version);
                i2 = get_intensity(refl);

                /* Calculate LSQ estimate of scaling factor */
                top += i1 * i2;
                bot += i2 * i2;

        }

        s = top / bot;

        return s;
}


static double cc_intensities(RefList *reference, RefList *new,
                             const SymOpList *sym)
{
	/* "x" is "reference" */
	float s_xy = 0.0;
	float s_x = 0.0;
	float s_y = 0.0;
	float s_x2 = 0.0;
	float s_y2 = 0.0;
	int n = 0;
	float t1, t2;

	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(new, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{

		double i1, i2;
		signed int hu, ku, lu;
		signed int h, k, l;
		Reflection *reference_version;

		get_indices(refl, &h, &k, &l);
		get_asymm(sym, h, k, l, &hu, &ku, &lu);

		reference_version = find_refl(reference, hu, ku, lu);
		if ( reference_version == NULL ) continue;

		i1 = get_intensity(reference_version);
		i2 = get_intensity(refl);


		s_xy += i1 * i2;
		s_x += i1;
		s_y += i2;
		s_x2 += i1 * i1;
		s_y2 += i2 * i2;
		n++;

	}

	t1 = s_x2 - s_x*s_x / n;
	t2 = s_y2 - s_y*s_y / n;

	if ( (t1 <= 0.0) || (t2 <= 0.0) ) return 0.0;

	return (s_xy - s_x*s_y/n) / sqrt(t1*t2);
}


static double *check_hist_size(int n, double *hist_vals)
{
	int ns;
	double *tryMe;

	if ( n % 1000 ) return hist_vals;

	ns = n / 1000;
	ns = (ns+1)*1000;

	tryMe = realloc(hist_vals, ns*sizeof(double));
	if ( tryMe == NULL ) {
		ERROR("Failed to allocate space for histogram.\n");
	}
	return tryMe;
}


static void display_progress(int n_images, int n_crystals, int n_crystals_used)
{
        if ( !isatty(STDERR_FILENO) ) return;
        if ( tcgetpgrp(STDERR_FILENO) != getpgrp() ) return;

        pthread_mutex_lock(&stderr_lock);
        fprintf(stderr, "\r%i images processed, %i crystals, %i crystals used.",
                n_images, n_crystals, n_crystals_used);
        pthread_mutex_unlock(&stderr_lock);

        fflush(stdout);
}


static unsigned char *flags_from_list(RefList *list)
{
        Reflection *refl;
        RefListIterator *iter;
        unsigned char *out = new_arr_flag();

        for ( refl = first_refl(list, &iter);
              refl != NULL;
              refl = next_refl(refl, iter) ) {

                signed int h, k, l;

                get_indices(refl, &h, &k, &l);

                set_arr_flag(out, h, k, l, 1);

        }

        return out;

}


static double *intensities_from_list(RefList *list)
{
        Reflection *refl;
        RefListIterator *iter;
        double *out = new_arr_intensity();

        for ( refl = first_refl(list, &iter);
              refl != NULL;
              refl = next_refl(refl, iter) ) {

                signed int h, k, l;
                double intensity = get_intensity(refl);

                get_indices(refl, &h, &k, &l);

                set_arr_intensity(out, h, k, l, intensity);

        }

        return out;
}


static double sym_lookup_intensity(const double *intensities,
                                   const unsigned char *flags,
                                   const SymOpList *sym,
                                   signed int h, signed int k, signed int l)
{
        int i;
        double ret = 0.0;

        for ( i=0; i<num_equivs(sym, NULL); i++ ) {

                signed int he;
                signed int ke;
                signed int le;
                double f, val;

                get_equiv(sym, NULL, i, h, k, l, &he, &ke, &le);

                f = (double)lookup_arr_flag(flags, he, ke, le);
                val = lookup_arr_intensity(intensities, he, ke, le);

                ret += f*val;

        }

        return ret;
}


static void copy_to_image_sm( RefList *model, struct image_sm* model_sm )
{
	RefListIterator *iter;
	Reflection *refl;
	int n = 0;
	signed int h, k, l;

	for ( refl = first_refl(model, &iter);
		  refl != NULL;
		  refl = next_refl(refl, iter) )
	{
		double var, esd;
		int red;

		get_indices(refl, &h, &k, &l);
		model_sm->serial[n] = SERIAL(h,k,l);
		model_sm->intensity[n] = get_intensity(refl);

		red = get_redundancy(refl);
		if( red<1 ) { esd = 0; }
		else {
			var = get_temp2(refl) / get_temp1(refl);
			esd = sqrt(var)/sqrt(red);
		}
		model_sm->sigma[n] = esd;

		n++;
	}

	model_sm->n_ref = n;

}


static int add_crystal(RefList *model, struct image *image, Crystal *cr,
                       RefList *reference, const SymOpList *sym,
                       struct image_sm *this_image_sm,
                       double **hist_vals, signed int hist_h,
                       signed int hist_k, signed int hist_l, int *hist_n,
                       int config_nopolar, double min_snr, double max_adu,
                       double push_res, double min_cc, int do_scale,
                       FILE *stat)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *new_reflist;
	RefList *this_image;
	double scale, cc;
	int n_ref = 0;

	new_reflist = crystal_get_reflections(cr);
	this_image = reflist_new();

	/* First, correct for polarisation */
	if ( !config_nopolar ) {
		polarisation_correction(new_reflist, crystal_get_cell(cr), image);
	}

	if ( reference != NULL ) {
		if ( do_scale ){
			scale = scale_intensities(reference, new_reflist, sym);
		} else {
			scale = 1.0;
		}
		cc = cc_intensities(reference, new_reflist, sym);
		if ( cc < min_cc ) return 1;
		if ( isnan(scale) ) return 1;
		if ( scale <= 0.0 ) return 1;
		if ( stat != NULL ) {
			fprintf(stat, "%s %s %f %f\n", image->filename,
			        get_event_string(image->event), scale, cc);
		}
	} 
	else {
		scale = 1.0;
	}

	double refl_intensity, refl_sigma, refl_pk;
	signed int h, k, l;
	int model_redundancy;
	Reflection *model_version;
	double w;
	double temp, delta, R, mean, M2, sumweight;
	double res, max_res;

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{

		refl_intensity = scale * get_intensity(refl);
		refl_sigma = scale * get_esd_intensity(refl);
		refl_pk = get_peak(refl);
		w = 1.0;

		if ( (min_snr > -INFINITY) && isnan(refl_sigma) ) continue;
		if ( refl_intensity < min_snr * refl_sigma ) continue;
		if ( refl_pk > max_adu ) continue;

		get_indices(refl, &h, &k, &l);
		/* resolution limit */
		max_res = push_res + crystal_get_resolution_limit(cr);
		res = 2.0*resolution(crystal_get_cell(cr), h, k, l);
		if ( res > max_res ) continue;

		/* Put into the asymmetric unit for the target group */
		get_asymm(sym, h, k, l, &h, &k, &l);

		model_version = find_refl(model, h, k, l);
		if ( model_version == NULL ) {
			model_version = add_refl(model, h, k, l);
		}
		
		// 'add_refl()' initiate all data in 'model_version' to be 0
		mean = get_intensity(model_version);
		sumweight = get_temp1(model_version);
		M2 = get_temp2(model_version);

		temp = w + sumweight;
		delta = refl_intensity - mean;
		R = delta * w / temp;
		set_intensity(model_version, mean + R);
		set_temp2(model_version, M2 + sumweight * delta * R);
		set_temp1(model_version, temp);

		model_redundancy = get_redundancy(model_version);
		set_redundancy(model_version, ++model_redundancy);

		if ( *hist_vals != NULL ) {

			if ( (h==hist_h) && (k==hist_k) && (l==hist_l) ) {

				*hist_vals = check_hist_size(*hist_n,
							     *hist_vals);

				/* Check again because realloc might have
				 * failed */
				if ( *hist_vals != NULL ) {
					(*hist_vals)[*hist_n] = refl_intensity;
					*hist_n += 1;
				}

			}

		}

		if( this_image_sm == 0 ) continue;

	    // add this reflection to list ('this_image') 
        model_version = find_refl(this_image, h, k, l);
        if ( model_version == NULL ) {
			model_version = add_refl(this_image, h, k, l);
			n_ref++;
        }
        // 'add_refl()' initiate all data in 'model_version' to be 0
        mean = get_intensity(model_version);
		sumweight = get_temp1(model_version);
		M2 = get_temp2(model_version);
		temp = w + sumweight;
		delta = refl_intensity - mean;
		R = delta * w / temp;
		// set
		set_intensity( model_version, mean + R );
		set_temp2(model_version, M2 + sumweight * delta * R);
		set_temp1(model_version, temp);
		model_redundancy = get_redundancy(model_version);
		set_redundancy( model_version, ++model_redundancy);
		
	}

	if(this_image_sm != NULL && num_reflections( this_image ) == 0) return 1;

	if( this_image_sm != NULL ){
		/* allocate memory for this_image_sm */
		this_image_sm->n_ref = n_ref;
		this_image_sm->serial = malloc(n_ref*sizeof(unsigned int));
		this_image_sm->intensity = malloc(n_ref*sizeof(float));
		this_image_sm->sigma = malloc(n_ref*sizeof(float));
		this_image_sm->cell = cell_new_from_cell(crystal_get_cell(cr));

		/* set this_image_sm */
		copy_to_image_sm(this_image, this_image_sm);
	}

	reflist_free(this_image);
	return 0;
}


static int add_all(Stream *st, RefList *model, RefList *reference,
                     const SymOpList *sym, struct image_sm **image_array,
                     double **hist_vals, 
                     signed int hist_h, signed int hist_k, signed int hist_l,
                     int *hist_i, int config_nopolar, int min_measurements,
                     double min_snr, double max_adu, 
                     int start_after, int stop_after, double min_res, 
                     double push_res, double min_cc, int do_scale, 
                     char *stat_output, int *assignments)
{
	/*
		image_array and assignments should be NULL or allocated together
	*/

	int rval, i;
	int n_images = 0;
	int n_crystals_used = 0, n_crystals_seen = 0;

	Reflection *refl;
	RefListIterator *iter;

	FILE *stat = NULL;
	if ( stat_output != NULL ) {
		stat = fopen(stat_output, "w");
		if ( stat == NULL ) {
			ERROR("Failed to open statistics output file %s\n",
			      stat_output);
		}
	}

	do {

		struct image image;
		//int i;

		image.det = NULL;

		/* Get data from next chunk */
		rval = read_chunk_2(st, &image, STREAM_READ_UNITCELL | STREAM_READ_REFLECTIONS);
		if ( rval ) break;

		image_feature_list_free(image.features);

		n_images++;

		for ( i=0; i<image.n_crystals; i++ ) {

			int r;
			Crystal *cr = image.crystals[i];
			n_crystals_seen++;

			if(n_crystals_seen <= start_after) continue;
			if(crystal_get_resolution_limit(cr) < min_res) continue;

			if( image_array != NULL ){
				r = add_crystal(model, &image, cr, reference, sym, 
								image_array[n_crystals_seen-1], 
								hist_vals, hist_h, hist_k, hist_l, hist_i,
								config_nopolar, min_snr, max_adu, 
								push_res, min_cc, do_scale, stat);
				if ( r == 0 )
					assignments[n_crystals_seen-1] = 0;
			}
			else{
				r = add_crystal(model, &image, cr, reference, sym, 
								NULL, 
								hist_vals, hist_h, hist_k, hist_l, hist_i,
								config_nopolar, min_snr, max_adu, 
								push_res, min_cc, do_scale, stat);
			}

			if ( r == 0 ) {
				n_crystals_used++;
			}

			reflist_free(crystal_get_reflections(cr));
			cell_free(crystal_get_cell(cr));
			crystal_free(cr);

			if ( n_crystals_used == stop_after ) break;

		}

		free(image.filename);
		free(image.crystals);

        display_progress(n_images, n_crystals_seen, n_crystals_used);

		if ( (stop_after>0) && (n_crystals_used == stop_after) ) break;

	} while ( rval == 0 );

	for ( refl = first_refl(model, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double var;
		int red;

		red = get_redundancy(refl);
		if ( red < min_measurements ) {
			set_redundancy(refl, 0);
			continue;
		}

		var = get_temp2(refl) / get_temp1(refl);
		set_esd_intensity(refl, sqrt(var)/sqrt(red));
	}

	if ( stat != NULL ) {
		fclose(stat);
	}

	return n_crystals_seen;
}


void get_twin_num( int space_group_num )
{
	// http://www.ccp4.ac.uk/html/twinning.html
	// N_TWINS is a global var defined at top of this file, default is 2
	if((space_group_num>=168 && space_group_num<=173) || (space_group_num>=195 && space_group_num<=199) || space_group_num==146 || (space_group_num>=75 && space_group_num<=80)) N_TWINS += 2;
	else if(space_group_num>=143 && space_group_num<=145) N_TWINS += 4;
	else if(space_group_num==149 || space_group_num==151 || space_group_num==153) N_TWINS += 3;
	else if(space_group_num==150 || space_group_num==152 || space_group_num==154) N_TWINS += 3;
	else N_TWINS += 1;
}


void get_ith_twin( int h, int k, int l, int *he, int *ke, int *le, int ith )
{
	char * op_tmp = NULL;
	int space_group_op_num = 0, i = 0;

	// space group twinning op
	// http://www.ccp4.ac.uk/html/twinning.html
	if((space_group_num>=168 && space_group_num<=173) || (space_group_num>=195 && space_group_num<=199) || space_group_num==146 || (space_group_num>=75 && space_group_num<=80)){
		assert(ith < 2 + ADD_N_TWINS);
		space_group_op_num = 2;
		if(ith == 0) {*he= h; *ke= k; *le= l; return; }
		if(ith == 1) {*he= k; *ke= h; *le=-l; return; }
	}
	else if(space_group_num>=143 && space_group_num<=145){
		assert(ith < 4 + ADD_N_TWINS);
		space_group_op_num = 4;
		if(ith == 0) {*he= h; *ke= k; *le= l; return; }
		if(ith == 1) {*he=-h; *ke=-k; *le= l; return; }
		if(ith == 2) {*he= k; *ke= h; *le=-l; return; }
		if(ith == 3) {*he=-k; *ke=-h; *le=-l; return; }
	}
	else if(space_group_num==149 || space_group_num==151 || space_group_num==153){
		assert(ith < 3 + ADD_N_TWINS);
		space_group_op_num = 3;
		if(ith == 0) {*he= h; *ke= k; *le= l; return; }
		if(ith == 1) {*he=-h; *ke=-k; *le= l; return; }
		if(ith == 2) {*he= k; *ke= h; *le=-l; return; }
	}
	else if(space_group_num==150 || space_group_num==152 || space_group_num==154){
		assert(ith < 3 + ADD_N_TWINS);
		space_group_op_num = 3;
		if(ith == 0) {*he= h; *ke= k; *le= l; return; }
		if(ith == 1) {*he=-h; *ke=-k; *le= l; return; }
		if(ith == 2) {*he=-k; *ke=-h; *le=-l; return; }
	}
	else if(space_group_num > 0){
		assert(ith < 1 + ADD_N_TWINS);
		space_group_op_num = 1;
		if(ith == 0) {*he= h; *ke= k; *le= l; return; }
	}

	// added op
	for ( i=0; i<3; i++){
		op_tmp = added_operators[ith - space_group_op_num][i];
		if (strcmp(op_tmp, "h") == 0){
			if ( i == 0 ) {*he= h; continue;}
			if ( i == 1 ) {*ke= h; continue;}
			if ( i == 2 ) {*le= h; continue;}
		}
		if (strcmp(op_tmp, "-h") == 0){
			if ( i == 0 ) {*he=-h; continue;}
			if ( i == 1 ) {*ke=-h; continue;}
			if ( i == 2 ) {*le=-h; continue;}
		}
		if (strcmp(op_tmp, "k") == 0){
			if ( i == 0 ) {*he= k; continue;}
			if ( i == 1 ) {*ke= k; continue;}
			if ( i == 2 ) {*le= k; continue;}
		}
		if (strcmp(op_tmp, "-k") == 0){
			if ( i == 0 ) {*he=-k; continue;}
			if ( i == 1 ) {*ke=-k; continue;}
			if ( i == 2 ) {*le=-k; continue;}
		}
		if (strcmp(op_tmp, "l") == 0){
			if ( i == 0 ) {*he= l; continue;}
			if ( i == 1 ) {*ke= l; continue;}
			if ( i == 2 ) {*le= l; continue;}
		}
		if (strcmp(op_tmp, "-l") == 0){
			if ( i == 0 ) {*he=-l; continue;}
			if ( i == 1 ) {*ke=-l; continue;}
			if ( i == 2 ) {*le=-l; continue;}
		}
		ERROR("I can't recognize the input ambiguity operator (No.%d) !", i+1);
	}
}


SymOpList* get_ith_twin_op(int ith)
{
	char *op_str = NULL, op_str_2[99];
	SymOpList *amb = NULL;
	int space_group_op_num = 0, i = 0;

	assert(ith >= 0);
	if(ith == 0) return NULL;

	// space group
	if((space_group_num>=168 && space_group_num<=173) || (space_group_num>=195 && space_group_num<=199) || space_group_num==146 || (space_group_num>=75 && space_group_num<=80)){
		assert(ith < 2 + ADD_N_TWINS);
		space_group_op_num = 2;
		if(ith == 1) op_str = strdup("k,h,-l");
	}
	else if(space_group_num>=143 && space_group_num<=145){
		assert(ith < 4 + ADD_N_TWINS);
		space_group_op_num = 4;
		if(ith == 1) op_str = strdup("-h,-k,l");
		else if(ith == 2) op_str = strdup("k,h,-l");
		else if(ith == 3) op_str = strdup("-k,-h,-l");
	}
	else if(space_group_num==149 || space_group_num==151 || space_group_num==153){
		assert(ith < 3 + ADD_N_TWINS);
		space_group_op_num = 3;
		if(ith == 1) op_str = strdup("-h,-k,l");
		else if(ith == 2) op_str = strdup("k,h,-l");
	}
	else if(space_group_num==150 || space_group_num==152 || space_group_num==154){
		assert(ith < 3 + ADD_N_TWINS);
		space_group_op_num = 3;
		if(ith == 1) op_str = strdup("-h,-k,l");
		else if(ith == 2) op_str = strdup("-k,-h,-l");
	}
	else if(space_group_num > 0){
		assert(ith < 1 + ADD_N_TWINS);
		space_group_op_num = 1;
	}

	// added op
	i = ith - space_group_op_num;
	if(i >= 0){
		sprintf(op_str_2, "%s,%s,%s", added_operators[i][0], added_operators[i][1], added_operators[i][2]);
		op_str = strdup(op_str_2);
	}

	// parse amb
	if(op_str != NULL){
		amb = parse_symmetry_operations(op_str);
		return amb;
	}
	else
		return NULL;
}


void stat_pearson_i_sp(struct image_sm *image, RefList *full_list, double *val, 
						const SymOpList *sym, double rmin, double rmax)
{

	double **vec3, **vec4;
	int ni = image->n_ref;
	int *nacc_twin = calloc(N_TWINS, sizeof(int));
	int i, tw;
	// init vec
	vec3 = (double**)malloc(N_TWINS*sizeof(double*));
	vec4 = (double**)malloc(N_TWINS*sizeof(double*));
	for(i=0; i<N_TWINS; i++){
		vec3[i] = calloc(ni, sizeof(double));
		vec4[i] = calloc(ni, sizeof(double));
	}

	double i1, i2, res;
	signed int h, k, l;
	signed int hp, kp, lp;
	signed int he, ke, le;

	int refl1_idx;
	Reflection *refl2;
	UnitCell *cell;

	for ( refl1_idx=0; refl1_idx<ni; refl1_idx++)
	{

		h = GET_H(image->serial[refl1_idx]);
		k = GET_K(image->serial[refl1_idx]);
		l = GET_L(image->serial[refl1_idx]);
		cell = image->cell;

		// judge resolution
		if( cell != NULL ){
			res = 2.0*resolution(cell, h, k, l);
			if ( res < rmin ) continue;
			if ( res > rmax ) continue;
		}

		i1 = image->intensity[refl1_idx];
		if( i1<=0 ) continue;

		for(tw=0; tw<N_TWINS; tw++){

			get_ith_twin(h, k, l, &he, &ke, &le, tw);
			get_asymm(sym, he, ke, le, &hp, &kp, &lp);
			refl2 = find_refl(full_list, hp, kp, lp);

			if ( refl2 != NULL && get_redundancy(refl2) > 0 ) /* This is a common reflection */
			{
				i2 = get_intensity(refl2);
				if( i2<=0 ) continue;

				vec3[tw][nacc_twin[tw]] = i1;
				vec4[tw][nacc_twin[tw]] = i2;
				nacc_twin[tw]++;
			}
		}

	}

	for(tw=0; tw<N_TWINS; tw++){
		if (nacc_twin[tw] < 2 ) val[tw] = 0;
		else                    val[tw] = gsl_stats_correlation(vec3[tw], 1, vec4[tw], 1, nacc_twin[tw]);
	}

	// free
	for(i=0; i<N_TWINS; i++){
		free(vec3[i]);
		free(vec4[i]);
	}
	free(vec3);
	free(vec4);
	free(nacc_twin);

	return;
}


void stat_pearson_i_sp2(RefList *image, RefList *full_list, double *val, 
						const SymOpList *sym, double rmin, double rmax, UnitCell *cell)
{

	double **vec3, **vec4;
	int ni = num_reflections(image);
	int *nacc_twin = calloc(N_TWINS, sizeof(int));
	int i, tw;
	// init vec
	vec3 = (double**)malloc(N_TWINS*sizeof(double*));
	vec4 = (double**)malloc(N_TWINS*sizeof(double*));
	for(i=0; i<N_TWINS; i++){
		vec3[i] = calloc(ni, sizeof(double));
		vec4[i] = calloc(ni, sizeof(double));
	}

	double i1, i2, res;
	signed int h, k, l;
	signed int hp, kp, lp;
	signed int he, ke, le;

	Reflection *refl1;
	Reflection *refl2;
	RefListIterator *iter;

	for ( refl1 = first_refl(image, &iter);
			refl1 != NULL;
			refl1 = next_refl(refl1, iter) )
	{

		get_indices(refl1, &h, &k, &l);

		// judge resolution
		if( cell != NULL ){
			res = 2.0*resolution(cell, h, k, l);
			if ( res < rmin ) continue;
			if ( res > rmax ) continue;
		}

		i1 = get_intensity(refl1);
		if( i1<=0 ) continue;

		for(tw=0; tw<N_TWINS; tw++){

			get_ith_twin(h, k, l, &he, &ke, &le, tw);
			get_asymm(sym, he, ke, le, &hp, &kp, &lp);
			refl2 = find_refl(full_list, hp, kp, lp);

			if ( refl2 != NULL && get_redundancy(refl2) > 0 ) /* This is a common reflection */
			{
				i2 = get_intensity(refl2);
				if( i2<=0 ) continue;

				vec3[tw][nacc_twin[tw]] = i1;
				vec4[tw][nacc_twin[tw]] = i2;
				nacc_twin[tw]++;
			}
		}

	}

	for(tw=0; tw<N_TWINS; tw++){
		if (nacc_twin[tw] < 2 ) val[tw] = 0;
		else                    val[tw] = gsl_stats_correlation(vec3[tw], 1, vec4[tw], 1, nacc_twin[tw]);
	}

	// free
	for(i=0; i<N_TWINS; i++){
		free(vec3[i]);
		free(vec4[i]);
	}
	free(vec3);
	free(vec4);
	free(nacc_twin);

	return;
}


int index_of_max_value( double* cc, int n_twins )
{
	int i;
	int winner=0;
	double max_value = cc[0];
	for(i=1;i<n_twins;i++)
	{
		if(cc[i] >= max_value) { max_value = cc[i]; winner=i; }
	}
	return winner;
}


void merge_image_at_winner_orientation( RefList* model, struct image_sm *image, int winner,
										double weight, const SymOpList *sym, double min_snr)
{

	int refl_idx;
	Reflection *model_version;
	int h, k, l;
	int model_redundancy;
	double refl_intensity, refl_sigma;
	double w;
	double temp, delta, R, mean, M2, sumweight;
	w = weight;

	for ( refl_idx = 0; refl_idx < image->n_ref; refl_idx++ )
	{
		
		refl_sigma = image->sigma[refl_idx];
		refl_intensity = image->intensity[refl_idx];

		if ( (min_snr > -INFINITY) && isnan(refl_sigma) ) continue;
		if ( refl_intensity < min_snr * refl_sigma ) continue;

		h = GET_H(image->serial[refl_idx]);
		k = GET_K(image->serial[refl_idx]);
		l = GET_L(image->serial[refl_idx]);
		get_ith_twin(h, k, l, &h, &k, &l, winner);
		get_asymm(sym, h, k, l, &h, &k, &l);

		// merge to model
		model_version = find_refl(model, h, k, l);
		if ( model_version == NULL ) {
			model_version = add_refl(model, h, k, l);
		}

		mean = get_intensity(model_version);
		sumweight = get_temp1(model_version);
		M2 = get_temp2(model_version);

		temp = w + sumweight;
		delta = refl_intensity - mean;
		R = delta * w / temp;
		set_intensity(model_version, mean + R);
		set_temp2(model_version, M2 + sumweight * delta * R);
		set_temp1(model_version, temp);

		model_redundancy = get_redundancy(model_version);
		set_redundancy(model_version, ++model_redundancy);
	}
}


void compute_weights( double *cc, double *weights, int iter )
{
	double sum_w = 0;
	// N_TWINS is a global var
	for(int i=0; i<N_TWINS; i++)
	{
		cc[i] = exp( cc[i]*(iter+1) );
		sum_w += cc[i] ;		
	}
	for(int i=0; i<N_TWINS; i++)
	{
		weights[i] = cc[i] / sum_w;
	}
}


void merge_image_to_model( RefList* model, struct image_sm *image, double* cc, int iter, 
							const SymOpList *sym, double min_snr )
{
	int twin;
	double weights[ N_TWINS ];
	
	compute_weights( cc, weights, iter);

	for( twin=0; twin<N_TWINS; twin++)
		merge_image_at_winner_orientation( model, image, twin, weights[twin], sym, min_snr);	
}


RefList* make_reflections_for_uc_from_asymm( RefList* asymm, bool random_intensity, const SymOpList *sym, 
											int min_measurements)
{
	Reflection *refl;
	Reflection *new_refl;
	RefListIterator *iter;
	int h,k,l, index;
	int he,ke,le;
	int i_twin, this_N_TWINS;
	double refl_intensity; 
	RefList* this_ref_model;
	this_ref_model=reflist_new();
	/* the following lines are for intensity generating for UC from Asymmetric Unit */
	double *intensities;
	unsigned char *flags;
	intensities = intensities_from_list(asymm);
	flags = flags_from_list(asymm);

	if (random_intensity ) this_N_TWINS = N_TWINS; // generating twin pairs
	else                   this_N_TWINS = 1; // only itself, strict to original data

	index = 0;
	srand(time(NULL));

	for ( refl = first_refl(asymm, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		if ( get_redundancy(refl) < min_measurements ) continue;
		get_indices(refl, &h, &k, &l);

		for(i_twin=0;i_twin<this_N_TWINS;i_twin++)
		{
			get_ith_twin(h, k, l, &he, &ke, &le, i_twin);
			get_asymm(sym, he, ke, le, &he, &ke, &le);
		
			if( random_intensity ) refl_intensity = rand();
			else refl_intensity = sym_lookup_intensity(intensities, flags, sym, he, ke, le);

			new_refl = find_refl( this_ref_model, he, ke, le );
			if (new_refl == NULL ) new_refl = add_refl( this_ref_model, he,ke,le );

			set_intensity( new_refl, refl_intensity );
			set_redundancy( new_refl, 1 );
		}
	}

	return this_ref_model;
	
}


RefList* emc( struct image_sm **image_array, RefList* full_list, 
			RefList* even_list, RefList* odd_list, const SymOpList *sym, 
			int n_crystals_all, int max_n_iter, int WINNER_TAKES_ALL, double min_snr, 
			int min_measurements, double rmin, double rmax, int *assignments )
{

	int image_index;
	int winner, n_iter;
	struct image_sm *image_a;
	RefList *this_full_list;

	Reflection *refl;
	RefListIterator *iter;

	double cc[ N_TWINS ];
	double winner_cc;
	int counter[ N_TWINS ], ii, n_crystals_used;

	n_iter = 0;
	
	printf("iteration    winner_cc    crystals_for_each_mode\n");

	while(n_iter<max_n_iter)
	{

		this_full_list = reflist_new();

		for( ii=0; ii<N_TWINS; ii++ ) counter[ii] = 0;
		winner_cc = 0.0;
		n_crystals_used = 0;

		for( image_index=0; image_index<n_crystals_all; image_index++)
		{

			// good crystal ?
			if ( assignments[image_index] < 0 ) continue;

			image_a = image_array[ image_index ];
			stat_pearson_i_sp(image_a, full_list, cc, sym, rmin, rmax);	
			/* winner takes all *
			 * find the largest cc, and assign the orientation to it
			 * add this to the diffraction volume for this round
			 */
			winner = index_of_max_value( cc, N_TWINS);
			winner_cc += cc[winner];
			counter[ winner ] = counter[ winner ] + 1;

			if ( n_iter == max_n_iter-1 ){
				// full lsit
				merge_image_at_winner_orientation( this_full_list, image_a, winner, 1, sym, min_snr);
				// even / odd partial list
				if ( image_index % 2 == 1 )
					merge_image_at_winner_orientation( odd_list, image_a, winner, 1, sym, min_snr);
				else if ( image_index % 2 == 0 )
					merge_image_at_winner_orientation( even_list, image_a, winner, 1, sym, min_snr);
				// update assignments
				assignments[image_index] = winner;
			}
			else{
				if ( WINNER_TAKES_ALL )
					merge_image_at_winner_orientation( this_full_list, image_a, winner, 1, sym, min_snr);
				else
					merge_image_to_model( this_full_list, image_a, cc, n_iter, sym, min_snr);
			}

			n_crystals_used++;
		}

		winner_cc /= n_crystals_used;

		printf( "%-13d%-13.5f", n_iter, winner_cc);
		for(ii=0;ii<N_TWINS;ii++) printf("%d ", counter[ii]);
		printf("\n");
		
		/* set min_measurements and esd */
		for ( refl = first_refl(this_full_list, &iter);
			  refl != NULL;
			  refl = next_refl(refl, iter) )
		{
			double var;
			int red;

			red = get_redundancy(refl);
			if ( red < min_measurements ) {
				set_redundancy(refl, 0);
				continue;
			}
			var = get_temp2(refl) / get_temp1(refl);
			set_esd_intensity(refl, sqrt(var)/sqrt(red));
		}

		reflist_free( full_list );
		full_list=copy_reflist( this_full_list );
		reflist_free( this_full_list );

		n_iter++;
	}

	// even / odd list
	/* set min_measurements and esd */
	for ( refl = first_refl(even_list, &iter);
		  refl != NULL;
		  refl = next_refl(refl, iter) )
	{
		double var;
		int red;

		red = get_redundancy(refl);
		if ( red < min_measurements ) {
			set_redundancy(refl, 0);
			continue;
		}
		var = get_temp2(refl) / get_temp1(refl);
		set_esd_intensity(refl, sqrt(var)/sqrt(red));
	}

	/* set min_measurements and esd */
	for ( refl = first_refl(odd_list, &iter);
		  refl != NULL;
		  refl = next_refl(refl, iter) )
	{
		double var;
		int red;

		red = get_redundancy(refl);
		if ( red < min_measurements ) {
			set_redundancy(refl, 0);
			continue;
		}
		var = get_temp2(refl) / get_temp1(refl);
		set_esd_intensity(refl, sqrt(var)/sqrt(red));
	}

	return full_list;

}


static void reindex_reflections(FILE *fh, FILE *ofh, int assignment,
                                SymOpList *amb)
{
	/*
		If assignment == 0, amb should be NULL
	*/

	int first = 1;

	do {

		char *rval;
		char line[1024];
		int n;
		signed int h, k, l;
		int r;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) break;

		if ( strcmp(line, REFLECTION_END_MARKER"\n") == 0 ) {
			fputs(line, ofh);
			return;
		}

		if ( first ) {
			fputs(line, ofh);
			first = 0;
			continue;
		}

		r = sscanf(line, "%i %i %i%n", &h, &k, &l, &n);

		/* See scanf() manual page about %n to see why <3 is used */
		if ( (r < 3) && !first ) return;

		if ( assignment > 0 ) {
			get_equiv(amb, NULL, 0, h, k, l, &h, &k, &l);
		}

		fprintf(ofh, "%4i %4i %4i%s", h, k, l, line+n);

	} while ( 1 );
}


/* This is nasty, but means the output includes absolutely everything in the
 * input, even stuff ignored by read_chunk() */
static void write_reindexed_stream(const char *infile, const char *outfile,
                                   int *assignments, int argc, char *argv[])
{
	FILE *fh;
	FILE *ofh;
	int i;
	struct rvec as, bs, cs;
	int have_as = 0;
	int have_bs = 0;
	int have_cs = 0;
	int done = 0;
	SymOpList *amb = NULL;

	fh = fopen(infile, "r");
	if ( fh == NULL ) {
		ERROR("Failed to open '%s'\n", infile);
		return;
	}

	ofh = fopen(outfile, "w");
	if ( ofh == NULL ) {
		ERROR("Failed to open '%s'\n", outfile);
		return;
	}

	/* Copy the header */
	do {

		char line[1024];
		char *rval;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) {
			ERROR("Failed to read stream audit info.\n");
			return;
		}

		if ( strncmp(line, "-----", 5) == 0 ) {

			done = 1;

			/* Add our own header */
			fprintf(ofh, "Re-indexed by EM-detwin "CRYSTFEL_VERSIONSTRING"\n");
			if ( argc > 0 ) {
				for ( i=0; i<argc; i++ ) {
					if ( i > 0 ) fprintf(ofh, " ");
					fprintf(ofh, "%s", argv[i]);
				}
				fprintf(ofh, "\n");
			}

		}

		fputs(line, ofh);

	} while  ( !done );

	i = 0;
	done = 0;
	do {

		char *rval;
		char line[1024];
		int d = 0;
		float u, v, w;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) break;

		// good crystal ?
		if ( assignments[i] < 0 ) {
			if ( strcmp(line, CRYSTAL_START_MARKER"\n") == 0 ) done = 1;
			else if ( strcmp(line, CRYSTAL_END_MARKER"\n") == 0 && done == 1) {
				done = 0;
				i++;
				continue;
			}
		}
		// bad crystal, skip
		if ( done != 0 ) continue;

		if ( strncmp(line, "Cell parameters ", 16) == 0 ) {
			d = 1;
		}

		if ( sscanf(line, "astar = %f %f %f", &u, &v, &w) == 3 ) {
			as.u = u*1e9;  as.v = v*1e9;  as.w = w*1e9;
			have_as = 1;
			d = 1;
		}

		if ( sscanf(line, "bstar = %f %f %f", &u, &v, &w) == 3 ) {
			bs.u = u*1e9;  bs.v = v*1e9;  bs.w = w*1e9;
			have_bs = 1;
			d = 1;
		}

		if ( sscanf(line, "cstar = %f %f %f", &u, &v, &w) == 3 ) {
			cs.u = u*1e9;  cs.v = v*1e9;  cs.w = w*1e9;
			have_cs = 1;
			d = 1;
		}

		if ( have_as && have_bs && have_cs ) {

			UnitCell *cell;
			double asx, asy, asz;
			double bsx, bsy, bsz;
			double csx, csy, csz;
			double a, b, c, al, be, ga;

			cell = cell_new_from_reciprocal_axes(as, bs, cs);
			assert(cell != NULL);

			amb = get_ith_twin_op(assignments[i]);

			if ( amb != NULL ) {

				signed int h, k, l;
				struct rvec na, nb, nc;

				get_equiv(amb, NULL, 0, 1, 0, 0, &h, &k, &l);
				na.u = as.u*h + bs.u*k + cs.u*l;
				na.v = as.v*h + bs.v*k + cs.v*l;
				na.w = as.w*h + bs.w*k + cs.w*l;

				get_equiv(amb, NULL, 0, 0, 1, 0, &h, &k, &l);
				nb.u = as.u*h + bs.u*k + cs.u*l;
				nb.v = as.v*h + bs.v*k + cs.v*l;
				nb.w = as.w*h + bs.w*k + cs.w*l;

				get_equiv(amb, NULL, 0, 0, 0, 1, &h, &k, &l);
				nc.u = as.u*h + bs.u*k + cs.u*l;
				nc.v = as.v*h + bs.v*k + cs.v*l;
				nc.w = as.w*h + bs.w*k + cs.w*l;

				cell_set_reciprocal(cell, na.u, na.v, na.w,
				                          nb.u, nb.v, nb.w,
				                          nc.u, nc.v, nc.w);

			}

			/* The cell parameters might change, so update them.
			 * Unique axis, centering and lattice type can't change,
			 * though. */
			cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
			fprintf(ofh, "Cell parameters %7.5f %7.5f %7.5f nm,"
			        " %7.5f %7.5f %7.5f deg\n",
			        a*1.0e9, b*1.0e9, c*1.0e9,
			        rad2deg(al), rad2deg(be), rad2deg(ga));

			cell_get_reciprocal(cell, &asx, &asy, &asz,
			                   &bsx, &bsy, &bsz,
			                   &csx, &csy, &csz);
			fprintf(ofh, "astar = %+9.7f %+9.7f %+9.7f nm^-1\n",
			        asx/1e9, asy/1e9, asz/1e9);
			fprintf(ofh, "bstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
			        bsx/1e9, bsy/1e9, bsz/1e9);
			fprintf(ofh, "cstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
			        csx/1e9, csy/1e9, csz/1e9);

			cell_free(cell);
			have_as = 0;  have_bs = 0;  have_cs = 0;

		}

		/* Not a bug: REFLECTION_START_MARKER gets passed through */
		if ( !d ) fputs(line, ofh);

		if ( strcmp(line, REFLECTION_START_MARKER"\n") == 0 ) {
			reindex_reflections(fh, ofh, assignments[i++], amb);
		}

	} while ( 1 );

	if ( !feof(fh) ) {
		ERROR("Error reading stream.\n");
	}

	fclose(fh);
	fclose(ofh);
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *output = NULL;
	char output_even[999];
	char output_odd[999];
	char *output_assignments = NULL;
	char *output_stream = NULL;
	char *stat_output = NULL;
	char *sym_str = NULL;
	SymOpList *sym;
	Stream *st;
	RefList *model;
	int tmp = 0;
	int config_scale = 0;
	char *histo = NULL;
	signed int hist_h, hist_k, hist_l;
	signed int hist_nbins=50;
	float hist_min=0.0, hist_max=0.0;
	double *hist_vals = NULL;
	int hist_i;
	char *histo_params = NULL;
	int config_nopolar = 0;
	char *rval;
	int min_measurements = 2;
	int start_after = 0;
	int stop_after = 0;
	double min_snr = -INFINITY;
	double max_adu = +INFINITY;
	double min_res = 0.0;
	double push_res = +INFINITY;
	double min_cc = -INFINITY;
	int twopass = 0;
	char *audit_info;

	int max_n_iter=30;
	int n_crystals_all=0;
	int WINNER_TAKES_ALL = 0;
	int cc_only = 0;
	float highres, lowres;
	double rmin = 0.0;  /* m^-1 */
	double rmax = INFINITY;  /* m^-1 */
	char *new_operaters = NULL;

	struct image_sm **image_array = NULL;
	RefList *full_list = NULL;
	RefList *even_list = NULL;
	RefList *odd_list = NULL;
	int *assignments = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"start-after",        1, NULL,               's'},
		{"stop-after",         1, NULL,               'f'},
		{"max-niter",          1, NULL,               'm'},
		{"winner-takes-all",   0, &WINNER_TAKES_ALL,   1},
		{"scale",              0, &config_scale,       1},
		{"no-polarisation",    0, &config_nopolar,     1},
		{"cc-only",            0, &cc_only,            1},
		{"symmetry",           1, NULL,               'y'},
		{"spacegroupNum",      1, NULL,               'k'},
		{"histogram",          1, NULL,               'g'},
		{"hist-parameters",    1, NULL,               'z'},
		{"min-measurements",   1, NULL,                2},
		{"min-snr",            1, NULL,                3},
		{"max-adu",            1, NULL,                4},
		{"min-res",            1, NULL,                5},
		{"push-res",           1, NULL,                6},
		{"version",            0, NULL,                7},
		{"min-cc",             1, NULL,                8},
		{"stat",               1, NULL,                9},
		{"highres",            1, NULL,               10},
		{"lowres",             1, NULL,               11},
		{"write-assignments",  1, NULL,               12},
		{"add-operators",      1, NULL,               13},
		{"write-stream",       1, NULL,               14},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:o:s:f:m:y:k:g:w:z:",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'i' :
			filename = strdup(optarg);
			break;

			case 'o' :
			output = strdup(optarg);
			break;

			case 's' :
			start_after = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --start-after.\n");
				return 1;
			}
			break;

			case 'f' :
			stop_after = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --stop-after.\n");
				return 1;
			}
			break;

			case 'm' :
			max_n_iter = strtod( optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --max-niter.\n");
				return 1;
			}
			break;

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 'k' :
			space_group_num = strtod( optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --spacegroupNum.\n");
				return 1;
			}
			get_twin_num(space_group_num);
			break;

			case 'g' :
			histo = strdup(optarg);
			break;

			case 'z' :
			histo_params = strdup(optarg);
			break;

			case 2 :
			min_measurements = strtol(optarg, &rval, 10);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-measurements.\n");
				return 1;
			}
			break;

			case 3 :
			min_snr = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-snr.\n");
				return 1;
			}
			ERROR("WARNING: You have used --min-snr.\n");
			ERROR("WARNING: Please read the manual carefully to "
			      "learn about possible detrimental effects of this"
			      " option.\n");
			break;

			case 4 :
			max_adu = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --max-adu.\n");
				return 1;
			}
			break;

			case 5 :
			min_res = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-res.\n");
				return 1;
			}
			min_res = 1e10/min_res;  // m^(-1)
			break;

			case 6 :
			push_res = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --push-res.\n");
				return 1;
			}
			push_res = push_res*1e9;
			break;

			case 7 :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			case 8 :
			min_cc = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-cc.\n");
				return 1;
			}
			twopass = 1;
			break;

			case 9 :
			stat_output = strdup(optarg);
			twopass = 1;
			break;

			case 10 :
			if ( sscanf(optarg, "%e", &highres) != 1 ) {
				ERROR("Invalid value for --highres\n");
				return 1;
			}
			rmax = 1.0 / (highres/1e10);
			break;

			case 11 :
			if ( sscanf(optarg, "%e", &lowres) != 1 ) {
				ERROR("Invalid value for --lowres\n");
				return 1;
			}
			rmin = 1.0 / (lowres/1e10);
			break;

			case 12 :
			output_assignments = strdup(optarg);
			break;

			case 13 :
			new_operaters = strdup(optarg);
			parse_operators(new_operaters);
			break;

			case 14:
			output_stream = strdup(optarg);
			break;

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( filename == NULL ) {
		ERROR("Please specify filename using the -i option\n");
		return 1;
	}

	if ( output == NULL ) {
		output = strdup("processed.hkl");
	}
	strcpy(output_odd, output);
	strcat(output_odd, ".odd");
	strcpy(output_even, output);
	strcat(output_even, ".even");

	if ( sym_str == NULL ) {
		ERROR("Please specify point group using the -y option\n");
		return 1;
	}

	if ( N_TWINS <= 0 ) {
		ERROR("Please specify twinning operators using either -k or --add-operators option\n");
		return 1;
	}
	sym = get_pointgroup(sym_str);
	free(sym_str);

	if( space_group_num <= 0 ){
		STATUS("\nSpace group number is not given, ");
	}
	else{
		STATUS("\nSpace group number is %i, ",space_group_num);
	}
	if( new_operaters != NULL )
		STATUS("with manually added operators '%s', ", new_operaters);
	STATUS("so there are totally %i twinning mode(s).\n", N_TWINS);

	/* only output cc information */
	if(cc_only){
		double compare_cc[N_TWINS];
		RefList* ref1 = read_reflections_2(filename, &sym_str);
		RefList* ref2 = read_reflections_2(output, &sym_str);
		stat_pearson_i_sp2(ref1, ref2, compare_cc, sym, rmin, rmax, NULL);
		STATUS("\nCC between RefList -i and RefList -o for all possible twins : \n")
		for(int tw=0;tw<N_TWINS;tw++)
		{
			printf( "    twin_%d = %f \n", tw, compare_cc[tw] );
		}
		return 0;
	}

	/* Open the data stream */
	STATUS("\nLoading chunks and merging ...\n")
	st = open_stream_for_read(filename);
	if ( st == NULL ) {
		ERROR("Failed to open stream.\n");
		return 1;
	}

	model = reflist_new();

	if ( histo != NULL ) {

		int r;

		r = sscanf(histo, "%i,%i,%i", &hist_h, &hist_k, &hist_l);
		if ( r != 3 ) {
			ERROR("Invalid indices for '--histogram'\n");
			return 1;
		}
		free(histo);

		/* Allocate enough space that hist_vals isn't NULL.
		 * check_hist_size will realloc it straight away */
		hist_vals = malloc(1*sizeof(double));
		STATUS("Histogramming %i %i %i -> ", hist_h, hist_k, hist_l);

		/* Put into the asymmetric cell for the target group */
		get_asymm(sym, hist_h, hist_k, hist_l,
		          &hist_h, &hist_k, &hist_l);
		STATUS("%i %i %i\n", hist_h, hist_k, hist_l);

	}

	if ( histo_params != NULL ) {

		int rr;

		rr = sscanf(histo_params, "%f,%f,%i", &hist_min, &hist_max,
		                                      &hist_nbins);
		if ( rr != 3 ) {
			ERROR("Invalid parameters for '--hist-parameters'\n");
			return 1;
		}
		free(histo_params);
		if ( hist_max <= hist_min ) {
			ERROR("Invalid range for '--hist-parameters'. "
			      "Make sure that 'max' is greater than 'min'.\n");
			return 1;
		}

	}

	/* Need to do a second pass if we are scaling */
	if ( config_scale ) twopass = 1;

	/* only first pass */
	if ( ! twopass ) {

		/* get number of crystals */
		n_crystals_all = get_crystals_num(filename);
		if (n_crystals_all <= 0) return 1;
		else STATUS("There are %d crystals in total.", n_crystals_all);

		/* init image array, assignments */
		image_array = (struct image_sm **)malloc(n_crystals_all*sizeof(struct image_sm*));
		if(image_array == NULL){
			ERROR("Fail to allocate memory for crystals.");
			return 1;
		}
		assignments = (int *)malloc(n_crystals_all*sizeof(int));
		for(tmp=0; tmp<n_crystals_all; tmp++){
			image_array[tmp] = (struct image_sm *)malloc(sizeof(struct image_sm));
			assignments[tmp] = -1;
		}

		hist_i = 0;
		n_crystals_all = add_all(st, model, NULL, sym, image_array, &hist_vals, hist_h, 
			hist_k, hist_l, &hist_i, config_nopolar, min_measurements, min_snr,
			max_adu, start_after, stop_after, min_res, push_res, 
			min_cc, config_scale, stat_output, assignments);

		fprintf(stderr, "\n");
	}

	/* has second pass */
	else {

		hist_i = 0;

		// first pass
		n_crystals_all = add_all(st, model, NULL, sym, NULL, &hist_vals, hist_h, 
			hist_k, hist_l, &hist_i, config_nopolar, min_measurements, min_snr,
			max_adu, start_after, stop_after, min_res, push_res, 
			min_cc, config_scale, stat_output, NULL);

		fprintf(stderr, "\n");

		RefList *reference;

		if ( rewind_stream(st) ) {

			ERROR("Couldn't rewind stream - scaling cannot be "
			      "performed.\n");
			return 1;

		}

		STATUS("\nExtra pass for scaling and/or CCs ...\n");

		reference = model;
		model = reflist_new();

		if ( hist_vals != NULL ) {
			free(hist_vals);
			hist_vals = malloc(1*sizeof(double));
			hist_i = 0;
		}

		// allocate image_array and assignments
		image_array = (struct image_sm **)malloc(n_crystals_all*sizeof(struct image_sm*));
		if(image_array == NULL){
			ERROR("Fail to allocate memory for crystals.");
			return 1;
		}
		assignments = (int *)malloc(n_crystals_all*sizeof(int));
		for(tmp=0; tmp<n_crystals_all; tmp++){
			image_array[tmp] = (struct image_sm *)malloc(sizeof(struct image_sm));
			assignments[tmp] = -1;
		}

		// second pass
		n_crystals_all = add_all(st, model, reference, sym, image_array, &hist_vals, hist_h, 
				hist_k, hist_l, &hist_i, config_nopolar, min_measurements, min_snr, 
				max_adu, start_after, stop_after, min_res, push_res, 
				min_cc, config_scale, stat_output, assignments);

		fprintf(stderr, "\n");
		reflist_free(reference);

	}

	if ( hist_vals != NULL ) {
		STATUS("%i %i %i was seen %i times.\n", hist_h, hist_k, hist_l,
		                                        hist_i);
		plot_histogram(hist_vals, hist_i, hist_min, hist_max, hist_nbins);
		free(hist_vals);
	}

	audit_info = stream_audit_info(st);
	close_stream(st);
	
	// Now start anti-ambiguity
	even_list = reflist_new();
	odd_list = reflist_new();
	full_list = make_reflections_for_uc_from_asymm( model, false, sym, min_measurements );

	printf("\nExpectation Maximization\n");
	full_list = emc( image_array, full_list, even_list, odd_list, 
			sym, n_crystals_all, max_n_iter, WINNER_TAKES_ALL, min_snr, 
			min_measurements, rmin, rmax, assignments );

	printf("\nWriting results ...\n");
	// full list
	reflist_add_command_and_version(full_list, argc, argv);
	reflist_add_notes(full_list, "Audit information from stream:");
	reflist_add_notes(full_list, audit_info);
	write_reflist_2(output, full_list, sym);
	// odd list
	reflist_add_command_and_version(odd_list, argc, argv);
	reflist_add_notes(odd_list, "Audit information from stream:");
	reflist_add_notes(odd_list, audit_info);
	write_reflist_2(output_odd, odd_list, sym);
	// even list
	reflist_add_command_and_version(even_list, argc, argv);
	reflist_add_notes(even_list, "Audit information from stream:");
	reflist_add_notes(even_list, audit_info);
	write_reflist_2(output_even, even_list, sym);
	// wirte assignments
	if ( output_assignments != NULL ){
		FILE *fp = fopen(output_assignments, "w");
		if ( fp != NULL ){
			fprintf(fp, "%-15s%-20s\n", "crystal-rank", "reindexed-manner");
			for(tmp=0; tmp<n_crystals_all; tmp++){
				if(assignments[tmp] < 0) continue;
				fprintf(fp, "%-15d%-20d\n", tmp, assignments[tmp]);
			}
		}
		fclose(fp);
	}
	// write stream file
	if ( output_stream != NULL ){
		write_reindexed_stream(filename, output_stream, assignments, argc, argv);
	}

	// clean
	for(tmp=0; tmp<n_crystals_all; tmp++) {
		if(assignments[tmp] < 0) continue;
		free( image_array[tmp]->serial );
		free( image_array[tmp]->intensity );
		free( image_array[tmp]->sigma );
		cell_free( image_array[tmp]->cell );
		free(image_array[tmp]);
	}
	free(image_array);
	reflist_free(model);
	free_symoplist(sym);
	reflist_free(full_list);
	reflist_free(odd_list);
	reflist_free(even_list);
	free(filename);
	free(output);
	free(output_stream);
	free(output_assignments);
	if(assignments != NULL) free(assignments);

	printf("Done.\n");

	return 0;
	
}
