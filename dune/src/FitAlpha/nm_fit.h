/*
 * nm_fit.h
 * This file is part of ExCeCo
 *
 * Copyright (C) 2011 - Christian Diener, Wolfgang Giese
 *
 * ExCeCo is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * ExCeCo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ExCeCo. If not, see <http://www.gnu.org/licenses/>.
 */



#ifndef NM_FIT_H
#define NM_FIT_H

// GSL includes
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>

// Additional includes
#include <ctime>
#include <fstream>

#define N_P 1
#define TOL 1e-5

double *alpha;

extern double rssq(const gsl_vector *p, void *params);

unsigned int read_config(const char* config_name)
{
	std::ifstream cf(config_name);
	if( !cf.is_open() ) return 0;
	
	cf >> N_ALPHA;
	cf >> N_MATA;
	cf >> ALPHA_SEC;
	cf >> KBAR0;
	cf >> KBARI;
	cf >> BAR1_DIR;
	
	alpha = new double[N_MATA];
	
	for(unsigned int i=0; i<N_MATA; i++) cf >> alpha[i];
	
	cf.close();
	
	return N_MATA;
}

unsigned int run_fit(unsigned int max_iter)
{
	gsl_vector *p = gsl_vector_alloc(N_P);
	gsl_vector *step = gsl_vector_alloc(N_P);
	
	gsl_vector_set(p, 0, ALPHA_SEC);
	
	gsl_vector_set(step, 0, 0.1*ALPHA_SEC);
	
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    
    gsl_multimin_function opt_f;
    
    opt_f.n = N_P;
    opt_f.f = &rssq;
    opt_f.params = NULL;
    
    s = gsl_multimin_fminimizer_alloc (T, opt_f.n);
    gsl_multimin_fminimizer_set (s, &opt_f, p, step);
    
    unsigned int iter = 0;
    int status;
    double size, init = rssq(s->x, NULL);
    
    timespec cpu_timer;
	cpu_timer.tv_sec = 0;
	cpu_timer.tv_nsec = 0;
	clock_settime(CLOCK_PROCESS_CPUTIME_ID, &cpu_timer);
	
    std::cout<<"Start = ";
    for(unsigned int i=0; i<N_P; i++) std::cout<<gsl_vector_get(s->x, i)<<'\t';
    std::cout<<"f() = "<<init<<std::endl;
    
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status) break;
        
        size = gsl_multimin_fminimizer_size (s);
        
        status = gsl_multimin_test_size (size, TOL);
        
        std::cout << "================ ITER #" << iter << "======================" << std::endl;
		for(unsigned int i=0; i<N_P; i++) std::cout<<gsl_vector_get(s->x, i)<<'\t';
		std::cout<<" f() = "<<s->fval<<'\t';
		std::cout<<"size = "<<size<<std::endl<<std::endl;
		  

    }
    while (status == GSL_CONTINUE && iter < max_iter);
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_timer);
	std::cout<<"Done after "<<(double) cpu_timer.tv_sec*1.0e3 + cpu_timer.tv_nsec/1.0e6<<" ms."<<std::endl;
    
	if(status == GSL_SUCCESS) return iter;
	else return 0;
}



#endif
