/* star functions */

void reallocate_memory_maxpartstars(void);
void domain_resize_storage_stars(int count_get_star);

void star_density(void);
void update_SNII(void);
void star_ngb_feedback(void);

double SN_events(void);

integertime get_timestep_star(int p);
void update_star_timesteps(void);
void reconstruct_star_timebins(void);
void update_list_of_active_star_particles(void);