/*black hole functions*/

void reallocate_memory_maxpartbh(void);
void domain_resize_storage_bh(int count_get_bh);
void bh_density(void);
void update_bh_accretion_rate(void);
void bh_ngb_feedback(void);
void update_bh_timesteps(void);
void update_list_of_active_bh_particles(void);
integertime get_timestep_bh(int p, integertime ti_step);

