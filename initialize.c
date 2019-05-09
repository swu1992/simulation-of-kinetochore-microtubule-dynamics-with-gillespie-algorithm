#include <stdio.h>
#include <stdlib.h>
#include "tools.h"


void init_state(int igrp, t_state *state) 
{
    int i, j;
    int nPrt[6][5] = {{50, 50, 50, 50, 50}, //control
                      {50, 2, 50, 50, 50}, //knowdown k67a
                      {50, 50, 2, 50, 50}, //59c
                      {50, 50, 50, 2, 50},  // mast
                      {50, 50, 50, 50, 2}, //eb1
                      {2, 50, 50, 50, 50}}; //msps

    state->totT = 300;
    // state->outT = 0.1;
    state->nProt = 5;
    state->nMT = 13;

    for(i=0; i<state->nProt; i++) {
        state->protein[i] = nPrt[igrp][i];
    }
    // state->MT = (t_MT *) malloc (state->nMT * sizeof(t_MT));
    for(i=0; i<state->nMT; i++) {
        (state->MT[i]).dist_end[0] = 50;
        (state->MT[i]).dist_end[1] = 50;
        (state->MT[i]).end_status[0] = 1;
        (state->MT[i]).end_status[1] = 1;        
        (state->MT[i]).ncomplex = 0;
        (state->MT[i]).position = 0;
        for(j=0; j<state->nProt; j++) {
            (state->MT[i]).complex_prt[j] = 0;
        }
    }
    // state->hc_m = (float *) malloc ((10*state->nMT) * sizeof(float));
    // state->hc_p = (float *) malloc ((18*state->nMT) * sizeof(float));
}

void init_react(t_reaction *reaction) 
{
    int i, j;
    float k_association[5] = {0.04, 0.065, 0.009, 3.3, 0.9};
    float k_dissociation[5] = {0.1, 0.85, 0.1, 0.00235, 0.001};
    float k_sw_m[6] = {0.32, 0.32, 0.32, 0.1, 0.1, 0.1};
    float k_sw_p[6] = {0.0011, 0.56, 0.0002, 0.099, 0.05, 0.061};
    float k_grow_m[3] = {15.75, 0, 0};
    float k_shrink_m[3] = {0, 31.5, -0.63};
    float k_grow_p[3] = {31.5, 0.63, 0.63};
    float k_shrink_p[3] = {15.75, 0, 0};
    int end_status[6][2] = {{0, 2}, {0, 1}, {1, 2}, {1, 0}, {2, 0}, {2, 1}};

    reaction->k_asso_m = 0.04;
    reaction->k_diss_m = 0.1;
    for(i=0; i<5; i++) {
        reaction->k_asso_p[i] = k_association[i];
        reaction->k_diss_p[i] = k_dissociation[i];
    }
    for(i=0; i<6; i++) {
        reaction->k_sw_m[i] = k_sw_m[i];
        reaction->k_sw_p[i] = k_sw_p[i];
        reaction->end_status[i][0] = end_status[i][0];
        reaction->end_status[i][1] = end_status[i][1];        
    }
    for(i=0; i<3; i++) {
        reaction->k_grow_m[i] = k_grow_m[i];
        reaction->k_grow_p[i] = k_grow_p[i];
        reaction->k_shrink_m[i] = k_shrink_m[i];
        reaction->k_shrink_p[i] = k_shrink_p[i];
    }
}
