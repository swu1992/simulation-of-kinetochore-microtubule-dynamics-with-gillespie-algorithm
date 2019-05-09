typedef struct
{
  int dist_end[2];
  int end_status[2];  // 0:grow; 1: pause; 2:shrink
  int ncomplex;
  int complex_prt[20]; //0: msps; 1:k67a; 2:k59c; 3:mast; 4eb1;
  int position;
} t_MT;

typedef struct
{
  float totT;
  float outT;
  int nProt;
  int nMT;
  int protein[20]; //number of proteins
  t_MT MT[30];
  float hc_p[1000];
  float hc_m[1000];
} t_state;

typedef struct
{
  float k_asso_m;
  float k_diss_m;
  float k_grow_m[3];
  float k_shrink_m[3];
  float k_sw_m[6];
  float k_asso_p[5];
  float k_diss_p[5];
  float k_grow_p[3];
  float k_shrink_p[3];
  float k_sw_p[6];
  int end_status[6][2];
} t_reaction;

void init_MT(int nProtein, t_state *state);
int init_par(char *fn, t_state *state, t_reaction *(*allReaction));
void select_reaction(int nstep, int n_minus, int n_plus, t_state *state, float r, int *ind_react);
void update(int *ind, t_state *state, t_reaction *reaction);
float get_sum_hc(int nm, int np, t_state *state, t_reaction *allReaction);
void output(FILE *fp, t_state *state);
int get_flux(t_state *state);
