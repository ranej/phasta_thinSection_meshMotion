#include <phastaChef_empty.h>
#include <cstdio>

#ifdef __cplusplus
extern"C"{
#endif
void core_get_pos_on_surf (double dx[], double dy[], double dz[], int numnp,
                           int f[], double px[], double py[], double pz[]) {
  printf("This is a dummy function. Please compile with core at phastaChef level.\n");
}

void core_is_in_closure (int e_dim, int e_tag, int t_dim, int t_tag, int answer) {
  printf("This is a dummy function. Please compile with core at phastaChef level.\n");
}

void core_measure_mesh (double x1[], double x2[], double x3[], int numnp,
                        double minvq, double minfq) {
  printf("This is a dummy function. Please compile with core at phastaChef level.\n");
}

void core_update_rbms (double tx[], double ty[], double tz[],
                         double ax[], double ay[], double az[],
                         double px[], double py[], double pz[],
                         double ag[], double sc[],
                         int tags[],  int numRbm) {
  printf("This is a dummy function. Please compile with core at phastaChef level.\n");
}

void core_get_centroid (int r_tag, double ct[]) {
  printf("This is a dummy function. Please compile with core at phastaChef level.\n");
}


#ifdef __cplusplus
}
#endif

