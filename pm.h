#ifndef PM_H
#define PM_H 1

void pm_init(const int nc_pm, const int pm_factor, Mem* const mem_density, Mem* const mem_force, const float_t boxsize);
void pm_compute_forces(Particles* particles);

#endif
