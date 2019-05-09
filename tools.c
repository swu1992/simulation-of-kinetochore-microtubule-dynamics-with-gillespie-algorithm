#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

//hc_m:10:  0:association, 1; dissociation, 2: grow; 3: shrink; 4-9:switch
//hc_p:18:  0-4:association, 5-9; dissociation, 10:grow; 11:shrink; 12-17: switch
float calc_k(float dist, float *f)
{
  float rate;

  if(dist > 50) {
    rate = f[0];
  }
  else {
    rate = f[1] + dist*f[2];
  }
  return rate;
}

void calc_hc_m(int iStart, t_MT *MT, t_state *state, t_reaction *reactions)
{
  int i, j, indType;

  for(i=0; i<10; i++) {
    state->hc_m[iStart+i] = 0;
  }
  indType = MT->end_status[0];
  state->hc_m[iStart] = (MT->complex_prt[0] > 0 ? 0: reactions->k_asso_m * state->protein[0]);
  state->hc_m[iStart+1] = ((MT->ncomplex == 1 && MT->complex_prt[0] > 0) ? reactions->k_diss_m : 0);
  if(indType == 0 && MT->ncomplex > 0) {
    if(MT->dist_end[0] > 0 || MT->dist_end[1] > 0) {
      state->hc_m[iStart+2] = reactions->k_grow_m[0];
    }
  }
  if(indType == 2 && MT->ncomplex > 0) {
    state->hc_m[iStart+3] = calc_k(MT->dist_end[0], reactions->k_shrink_m);
  }
  state->hc_m[iStart+4+indType*2] = reactions->k_sw_m[indType*2];
  state->hc_m[iStart+4+indType*2+1] = reactions->k_sw_m[indType*2+1];
}

void calc_hc_p(int iStart, t_MT *MT, t_state *state, t_reaction *reactions)
{
  int i, j, indType;
  int *nPrt;
  float *k1, *k2, *ks;

  indType = MT->end_status[1];
  nPrt = MT->complex_prt;
  k1 = reactions->k_asso_p;
  k2 = reactions->k_diss_p;
  ks = reactions->k_sw_p;

  for(i=0; i<18; i++) {
    state->hc_p[iStart+i] = 0;
  }
  state->hc_p[iStart] = (MT->ncomplex > 0 ? 0: k1[0] * state->protein[0]);
  //Association
  if(MT->ncomplex > 0) {
    state->hc_p[iStart+1] = (nPrt[1] == 0 ? k1[1] * state->protein[1] : 0);
    if(nPrt[2] == 0) {
      state->hc_p[iStart+2] = (nPrt[3] > 0 ?
                (0.001*k1[2]*state->protein[2]) : (k1[2]*state->protein[2]));
    }
    if(nPrt[3] == 0) {
      state->hc_p[iStart+3] = (nPrt[2] > 0 ?
                (0.001*k1[3]*state->protein[3]) : (k1[3]*state->protein[3]));
    }
    state->hc_p[iStart+4] = (nPrt[4] == 0 ? k1[4] * state->protein[4] : 0);
  //Dissociation
    state->hc_p[iStart+5] = ((nPrt[0] > 0 && MT->ncomplex==1) ? k2[0] : 0);
    state->hc_p[iStart+6] = (nPrt[1] > 0 ? k2[1] : 0);
    state->hc_p[iStart+7] = (nPrt[2] > 0 ? k2[2] : 0);
    state->hc_p[iStart+8] = (nPrt[3] > 0 ? k2[3] : 0);
    state->hc_p[iStart+9] = (nPrt[4] > 0 ? k2[4] : 0);
  }
  //Growth & Shrink
  if(indType == 0 && MT->ncomplex > 0) {
    if(MT->dist_end[0]>0 || MT->dist_end[1]>0) {
      state->hc_p[iStart+10] = calc_k(MT->dist_end[1], reactions->k_grow_p);
    }
  }
  else if(indType == 2 && MT->ncomplex > 0) {
    state->hc_p[iStart+11] = reactions->k_shrink_p[0];
  }
  //Switch
  switch(indType) {
    case 0:
      state->hc_p[iStart+12] = (nPrt[1] > 0 ? 6.9*ks[0] : ks[0]);
      state->hc_p[iStart+13] = (nPrt[1] > 0 ? 6.9*ks[1] : ks[1]);
      break;
    case 1:
      state->hc_p[iStart+14] = (nPrt[1] > 0 ? 6.9*ks[2] : ks[2]);
      state->hc_p[iStart+15] = (nPrt[2] > 0 ? 0.00119*ks[3] : ks[3]);
      if(nPrt[3] > 0 && nPrt[4] == 0) {
        state->hc_p[iStart+15] *= 1.4;
      }
      else if(nPrt[4] > 0 && nPrt[3] == 0) {
        state->hc_p[iStart+15] *= 2.1;
      }
      else if(nPrt[3] > 0 && nPrt[4] > 0) {
        state->hc_p[iStart+15] *= 2.65;
      }
      break;
    case 2:
      state->hc_p[iStart+16] = (nPrt[2] > 0 ? 0.00119*ks[4] : ks[4]);
      state->hc_p[iStart+17] = (nPrt[2] > 0 ? 0.00119*ks[5] : ks[5]);
      if(nPrt[3] > 0 && nPrt[4] == 0) {
        state->hc_p[iStart+16] *= 1.4;
        state->hc_p[iStart+17] *= 1.4;
      }
      else if(nPrt[4] > 0 && nPrt[3] == 0) {
        state->hc_p[iStart+16] *= 2.1;
        state->hc_p[iStart+17] *= 2.1;
      }
      else if(nPrt[3] > 0 && nPrt[4] > 0) {
        state->hc_p[iStart+16] *= 2.65;
        state->hc_p[iStart+17] *= 2.65;
      }
      break;
  }
}

float get_sum_hc(int n_minus, int n_plus, t_state *state, t_reaction *reactions)
{
  float sum_hc = 0;
  int i, iMT, nMT;
  t_MT *MT;

  nMT = state->nMT;
  MT = state->MT;
  //minus end
  for(iMT=0; iMT<nMT; iMT++) {
    calc_hc_m(iMT*n_minus, &(MT[iMT]), state, reactions);
    calc_hc_p(iMT*n_plus, &(MT[iMT]), state, reactions);
  }
  for(i=0; i<nMT*n_minus; i++) {
    sum_hc += state->hc_m[i];
  }
  for(i=0; i<nMT*n_plus; i++) {
    sum_hc += state->hc_p[i];
  }
  if(sum_hc <= 0) {
    fprintf(stderr, "Wrong calculation of sum_hc %f \n", sum_hc);
    exit(1);
  }
  return sum_hc;
}

void select_reaction(int nstep, int n_minus, int n_plus, t_state *state, float r, int *ind_react)
{
  int i, j, nMT;
  float sp = 0.;

  ind_react[0] = -1;
  ind_react[1] = -1;
  ind_react[2] = -1;
  nMT = state->nMT;
  for(i=0; i<nMT; i++) {
    for(j=0; j<n_minus && ind_react[0] == -1; j++) {
      sp += state->hc_m[i*n_minus+j];
      if(r < sp) {
        // if(i == 0)
        // printf("step %d %d minus %d  %f %f %f %d %d %d \n", nstep, i, j, r, sp,  state->hc_m[i*n_minus+j], state->MT[i].dist_end[0], state->MT[i].dist_end[1], state->MT[i].position);
        ind_react[0] = i;
        ind_react[1] = -1;
        ind_react[2] = j;
        break;
      }
    }
    for(j=0; j<n_plus && ind_react[0] == -1; j++) {
      sp += state->hc_p[i*n_plus+j];
      if(r < sp) {
        // if(i==0)
        // printf("step %d %d plus %d %f %f %f %d %d %d %f %f %d %d %d %f %f %f\n", nstep,  i, j, r, sp, state->hc_p[i*n_plus+j], state->MT[i].dist_end[0], state->MT[i].dist_end[1],
        //     state->MT[i].position, state->hc_p[i*n_plus+16],state->hc_p[i*n_plus+17], state->MT[i].complex_prt[2],
        //     state->MT[i].complex_prt[3], state->MT[i].complex_prt[4], state->hc_p[j*18+2],
        //     state->hc_p[j*18+3],state->hc_p[j*18+4]);
        ind_react[0] = i;
        ind_react[1] = 1;
        ind_react[2] = j;
        break;
      }
    }
  }
  if(ind_react[0]<0) {
    fprintf(stderr, "Mistake in selecting reactions %d %f \n", ind_react[0], r);
    exit(1);
  }
}

void update_dist(int iEnd, int iGrow, t_MT *MT)
{
/*Plus End: grow phase yplus == 0, position+1, yminus-1
          shrink phase yplus == 50, position-1, yminus+1
  Minus End: grow phase yminus==0, position-1, yplus-1
          shrink phase yminus<50 and yplus<50, position+1,yplus+1
          shrink phase yminus<50 and yplus>=50, yminus+1
          shrink phase yminus>=50, yminus+1
  Both end == 0, stop grow or shrink
*/
  if(iEnd == 1 ) { //plus end
      if(iGrow == -1) {//shrink
        if(MT->dist_end[1] >= 50) {
          MT->dist_end[0]++;
          MT->position--;
        }
        else {
          MT->dist_end[1]++;
        }
      }
      else if(iGrow == 1) {//grow
        if(MT->dist_end[1] == 0) {
          MT->dist_end[0]--;
          MT->position++;
        }
        else {
          MT->dist_end[1]--;
        }
      }
  }
  else if(iEnd ==-1) {//minus End
    if(iGrow == -1) {
      if(MT->dist_end[0] < 50 && MT->dist_end[1] < 50) {
        MT->position++;
        MT->dist_end[1]++;
      }
      else {
        MT->dist_end[0]++;
      }
    }
    else if(iGrow == 1){
      if(MT->dist_end[0] == 0) {
        MT->dist_end[1]--;
        MT->position--;
      }
      else {
        MT->dist_end[0]--;
      }
    }
  }
  if(MT->dist_end[0]<0 || MT->dist_end[1]<0) {
    fprintf(stderr, "Mistake in updteing distance\n");
  }
}

void update_association(int iPrt, int iMT, t_state *state)
{
  state->protein[iPrt]--;
  state->MT[iMT].ncomplex++;
  state->MT[iMT].complex_prt[iPrt]++;
}

void update_dissociation(int iPrt, int iMT, t_state *state) {
  state->protein[iPrt]++;
  state->MT[iMT].ncomplex--;
  state->MT[iMT].complex_prt[iPrt]--;
}

void update(int *ind, t_state *state, t_reaction *reaction)
{
  int i, j, iMT, iEnd, iReact, type, nMT;

  nMT = state->nMT;
  iMT = ind[0];
  iEnd = ind[1];
  iReact = ind[2];
  switch(iEnd) {
    case -1:
      switch (iReact) {
        case 0:
          update_association(0, iMT, state);
          break;
        case 1:
          update_dissociation(0, iMT, state);
          break;
        case 2:
          update_dist(-1, 1, &(state->MT[iMT]));
          if(state->MT[iMT].dist_end[0] >200) {
          // printf("%d minus end, grow %d %d %d %f\n", iMT, state->MT[iMT].dist_end[0],
          //   state->MT[iMT].dist_end[1], state->MT[iMT].end_status[0], state->hc_m[iMT*10+iReact]);
          }
          break;
        case 3:
          update_dist(-1, -1, &(state->MT[iMT]));
          break;
        default:
          state->MT[iMT].end_status[0] = reaction->end_status[iReact-4][1];
          break;
      }
      break;
    case 1:
      for(i=0; i<5; i++) {
        if(i == iReact) {
          update_association(i, iMT, state);
        }
      }
      for(i=5; i<10; i++) {
        if(i == iReact) {
          update_dissociation(i-5, iMT, state);
        }
      }
      if(iReact == 10) {
        update_dist(1, 1, &(state->MT[iMT]));
      }
      else if(iReact == 11) {
        update_dist(1, -1, &(state->MT[iMT]));
        if(state->MT[iMT].dist_end[0] >200) {
          // printf("%d plus end, shrink %d %d %d %f\n", iMT, state->MT[iMT].dist_end[0],
          //   state->MT[iMT].dist_end[1], state->MT[iMT].end_status[1], state->hc_p[iMT*18+iReact]);
          }
      }
      else if(iReact > 11) {
        state->MT[iMT].end_status[1] = reaction->end_status[iReact-12][1];
      }
      break;
  }
}

int get_flux(t_state *state)
{
  int i, flux=0;

  for(i=0; i<state->nMT; i++) {
    flux += (state->MT[i]).position;
  }
  return flux;
}