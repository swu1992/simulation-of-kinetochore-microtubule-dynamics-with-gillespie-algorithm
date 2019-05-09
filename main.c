/* Target: Simulation of microtubule dynamics with Gillespie's algorithm 
 * Reference: Fernandez, Nicolas, et al. "A model for the regulatory network 
 * controlling the dynamics of kinetochore microtubule plus-ends and poleward 
 * flux in metaphase." Proceedings of the National Academy of Sciences 106.19 (2009)
 * Shanshan Wu
 * 05-06-2019
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>                                                                  
#include "tools.h"
#include "initialize.h"

int main(int argc, char *argv[])
{
  t_state state;
  t_reaction reactions;
  int flux, ind_react[3], iRun;
  float t, tau, sum_hc, rand1, rand2;
  char fn_traj[50], fn_summary[50];
  FILE *fp_traj, *fp_summary;
  int i, j, igrp;

  srand(1);
  igrp = strtol(argv[1], NULL, 10);
  sprintf(fn_summary, "./data/summary_%d.dat", igrp);
  fp_summary = fopen(fn_summary, "w");

  for(iRun=0; iRun<15; iRun++) {
    t = 0;
    init_state(igrp, &state);
    init_react(&reactions);
    sprintf(fn_traj, "./data/traj_%d_%d.dat", igrp, iRun+1);
    fp_traj = fopen(fn_traj, "w");
     j=0;
    while(t<state.totT) {
      rand1 = rand() / (RAND_MAX+1.);
      rand2 = rand() / (RAND_MAX+1.);
      sum_hc = get_sum_hc(10, 18, &state, &reactions);
      select_reaction(i, 10, 18, &state, sum_hc*rand1, ind_react);
      update(ind_react, &state, &reactions);
      tau = -log(rand2) / sum_hc;
      t += tau;
      flux = get_flux(&state);  
      if(j%5 == 0) {
        fprintf(fp_traj, "%f %f ", t, 0.008*flux/13);
        for(i=0; i<state.nMT; i++) {
          fprintf(fp_traj, "%d %d %d ", state.MT[i].position, state.MT[i].ncomplex, state.MT[i].complex_prt[0]);
        }
        fprintf(fp_traj, "\n");
      }
    }
    flux = get_flux(&state);  
    fprintf(fp_summary, "%d %f %f \n", iRun+1, t, flux/t);
  }

  fclose(fp_summary);
  return 0;
}
