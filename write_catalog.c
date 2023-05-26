//Designed to be included in make_sf_catalog.c

void stack_weak_lensing(int64_t n, int64_t j);

void write_catalog(int64_t n) {
  char buffer[1024];
  int64_t chunk_num = (catalog_client) ? (catalog_client - 1) : loaded_chunk;
  FILE *output = NULL;
  int64_t i, j, wl_j, output_sfh=0, calc_smhm=-1;
  
  if (PROFILING_MODE) return;
  if (!catalog_client || CATALOG_OUTPUT_MODE) {
    snprintf(buffer, 1024, "%s/sfr_catalog_%f.%"PRId64".bin", OUTBASE, scales[n], chunk_num);
    output = check_fopen(buffer, "w");
    int64_t i;
    for (i=offsets[n]; i<offsets[n+1]; i++) {
      struct catalog_halo h = halos[i];
      if (h.upid > -1) h.upid = halos[h.upid].id;
      if (h.descid > -1) h.descid = halos[h.descid].id;
      check_fwrite(&h, sizeof(struct catalog_halo), 1, output);
    }
    fclose(output);
  }
  
  //printf("Reached catalog %"PRId64"\n", n);
  //return;
  if (!SFH_OUTPUT_MODE && !POSTPROCESSING_MODE) return;

#define NUM_OUTPUTS SMHM_NUM_ZS //10
  int64_t snaps[NUM_OUTPUTS] = {0}, wl_snaps[NUM_WL_ZS] = {0};
  double zs[NUM_OUTPUTS] = SMHM_ZS; 
  double wl_zs[NUM_WL_ZS] = WL_ZS;
  //0.4, 0.616, 0.8756, 1.17, 1.526, 1.92};
  
  for (i=0; i<num_scales; i++) {
    for (j=0; j<NUM_OUTPUTS; j++)
      if (fabs(scales[i]-1.0/(1.0+zs[j])) < fabs(scales[snaps[j]]-1.0/(1.0+zs[j])))
	snaps[j] = i;
    for (j=0; j<NUM_WL_ZS; j++)
      if (fabs(scales[i]-1.0/(1.0+wl_zs[j])) < fabs(scales[wl_snaps[j]]-1.0/(1.0+wl_zs[j])))
	wl_snaps[j] = i;
  }
  
  if (POSTPROCESSING_MODE) {
    for (j=0; j<NUM_OUTPUTS; j++) smhm_scales[j] = scales[snaps[j]]/(double)NUM_BLOCKS;
    for (j=0; j<NUM_WL_ZS; j++) wl_scales[j] = scales[wl_snaps[j]]/(double)NUM_BLOCKS;
  }
    
  for (j=0; j<NUM_OUTPUTS; j++) if (snaps[j] == n) break;
  for (wl_j=0; wl_j<NUM_WL_ZS; wl_j++) if (wl_snaps[wl_j] == n) break;
  if (j<NUM_OUTPUTS) { calc_smhm = j; }
  if (POSTPROCESSING_MODE && (wl_j < NUM_WL_ZS)) stack_weak_lensing(n, wl_j);
  if ((calc_smhm>-1) //&& (scales[n] > 0.32)
      && SFH_OUTPUT_MODE) output_sfh = 1;
  if (!output_sfh && !POSTPROCESSING_MODE) return;
#undef NUM_OUTPUTS

  //Allocate treerootids and sm/icl/m/sfr of main progenitor histories
  int64_t *trids=NULL, *trids2=NULL;
  check_realloc_s(trids, sizeof(int64_t), max_halos);
  check_realloc_s(trids2, sizeof(int64_t), max_halos);
  for (j=offsets[n]; j<offsets[n+1]; j++) trids[j-offsets[n]] = j-offsets[n];

  float *hist = NULL, *sfrh, *iclh, *mh, *smh, *vmph, *dvmaxh,
    *a_infall = NULL, *a_first_infall, *mpeak_cen;
  int64_t num_elems = (n+1)*num_halos[n];
  check_realloc_s(hist, sizeof(float), 6*num_elems);
  check_realloc_s(a_infall, sizeof(float), 3*num_halos[n]);
  a_first_infall = a_infall + num_halos[n];
  mpeak_cen = a_infall + num_halos[n]*2;
  memset(hist, 0, sizeof(float)*6*num_elems);
  for (j=0; j<num_halos[n]; j++) {
    mpeak_cen[j] = 0;
    a_infall[j] = a_first_infall[j] = scales[n];
  }
  sfrh = hist;
  iclh = sfrh+num_elems;
  mh = iclh+num_elems;
  smh = mh+num_elems;
  vmph = smh+num_elems;
  dvmaxh = vmph+num_elems;
  for (i=n; i>=0; i--) {
    if (i<n) {
      for (j=offsets[i]; j<offsets[i+1]; j++) {
	assert(halos[j].descid >= offsets[i+1] && halos[j].descid < offsets[i+2]);
	int64_t trid = trids[halos[j].descid-offsets[i+1]];
	assert(trid >=-1 && trid < (offsets[n+1]-offsets[n]));
	int64_t is_mmp = (halos[j].flags & MMP_FLAG) ? 1 : 0;
	struct catalog_halo *dh = halos+halos[j].descid;
	if (dh->flags & ORPHAN_FLAG)
	  is_mmp = 1; //Orphans do not have MMP flag set, but are MMPs if their desc. is an orphan.
	if (is_mmp && !(halos[i].flags & IGNORE_FLAG))
	  trids2[j-offsets[i]] = trid;
	else
	  trids2[j-offsets[i]] = -1;
      }
      int64_t *tmp = trids;
      trids = trids2;
      trids2 = tmp;
    }
    for (j=offsets[i]; j<offsets[i+1]; j++) {
      if (halos[j].flags & IGNORE_FLAG) continue;
      int64_t trid = trids[j-offsets[i]];
      if (trid < 0) continue;
      int64_t elem = trid*(n+1)+i;
      sfrh[elem] = halos[j].sfr;
      iclh[elem] = halos[j].icl;
      mh[elem] = halos[j].mp;
      smh[elem] = halos[j].sm;
      vmph[elem] = halos[j].vmp;
      dvmaxh[elem] = halos[j].rank1;
      if ((halos[j].upid < 0) && (halos[j].m > mpeak_cen[trid])) mpeak_cen[trid] = halos[j].m;
      if (i<n) {
	assert(halos[j].descid >= offsets[i+1] && halos[j].descid < offsets[i+2]);
	if (a_infall[trid]==scales[i+1]) {
	  if ((halos[j].upid > -1) || (halos[halos[j].descid].upid > -1)) {
	    a_first_infall[trid] = a_infall[trid] = scales[i];
	  }
	}
	if ((halos[j].upid < 0) && (halos[halos[j].descid].upid > -1) &&
	    (halos[j].mp*2.0 > mpeak_cen[trid]) && (i>0))
	  a_first_infall[trid] = scales[i];
      }
    }
  }
  free(trids);
  free(trids2);

  if (POSTPROCESSING_MODE) {
    //SMHM+SSFR Maps
    if (calc_smhm > -1) {
      for (i=offsets[n]; i<offsets[n+1]; i++) {
	if (halos[i].flags & IGNORE_FLAG) continue;	
	//HM Bin:
	double hm = log10(halos[i].mp/h0);
	int64_t hm_b = (hm-SMHM_HM_MIN)*SMHM_HM_BPDEX;
	if ((hm_b>=0) && (hm_b <SMHM_HM_NBINS)) {
	  //SMHM Bin:
	  double shr = halos[i].obs_sm / (halos[i].mp/h0);
	  shr = (shr > 0) ? log10(shr) : SMHM_SMHM_MIN;
	  int64_t shr_b = (shr-SMHM_SMHM_MIN)*SMHM_SMHM_BPDEX;
	  if (shr_b < 0) shr_b = 0;
	  if (shr_b >= SMHM_SMHM_NBINS) shr_b = SMHM_SMHM_NBINS-1;

	  double shr_true = halos[i].sm / (halos[i].mp/h0);
	  shr_true = (shr_true > 0) ? log10(shr_true) : SMHM_SMHM_MIN;
	  int64_t shr_true_b = (shr_true-SMHM_SMHM_MIN)*SMHM_SMHM_BPDEX;
	  if (shr_true_b < 0) shr_true_b = 0;
	  if (shr_true_b >= SMHM_SMHM_NBINS) shr_true_b = SMHM_SMHM_NBINS-1;

	  double shr_icl_true = (halos[i].sm+halos[i].icl) / (halos[i].mp/h0);
	  shr_icl_true = (shr_icl_true > 0) ? log10(shr_icl_true) : SMHM_SMHM_MIN;
	  int64_t shr_icl_true_b = (shr_icl_true-SMHM_SMHM_MIN)*SMHM_SMHM_BPDEX;
	  if (shr_icl_true_b < 0) shr_icl_true_b = 0;
	  if (shr_icl_true_b >= SMHM_SMHM_NBINS) shr_icl_true_b = SMHM_SMHM_NBINS-1;

	  
	  int64_t bin = calc_smhm*SMHM_NBINS + hm_b*SMHM_SMHM_NBINS + shr_b;
	  int64_t bin_true = calc_smhm*SMHM_NBINS + hm_b*SMHM_SMHM_NBINS + shr_true_b;
	  int64_t bin_true_icl = calc_smhm*SMHM_NBINS + hm_b*SMHM_SMHM_NBINS + shr_icl_true_b;

	  smhm[bin]++;
	  if (halos[i].obs_sfr < halos[i].obs_sm*QUENCHED_THRESHOLD) smhm_q[bin]++;
	  else smhm_sf[bin]++;
	  if (halos[i].upid > -1) smhm_sat[bin]++;
	  else {
	    smhm_cen[bin]++;
	    if (halos[i].obs_sfr < halos[i].obs_sm*QUENCHED_THRESHOLD) smhm_cen_q[bin]++;
	    else smhm_cen_sf[bin]++;
	  }

	  smhm_true[bin_true]++;
	  smhm_true_icl[bin_true_icl]++;
	  if (halos[i].sfr < halos[i].sm*QUENCHED_THRESHOLD) smhm_true_q[bin_true]++;
	  else smhm_true_sf[bin_true]++;
	  if (halos[i].upid > -1) smhm_true_sat[bin_true]++;
	  else {
	    smhm_true_cen[bin_true]++;
	    smhm_true_cen_icl[bin_true_icl]++;
	    if (halos[i].sfr < halos[i].sm*QUENCHED_THRESHOLD) smhm_true_cen_q[bin_true]++;
	    else smhm_true_cen_sf[bin_true]++;
	  }
	}

	double sm = (halos[i].obs_sm > 0) ? log10(halos[i].obs_sm) : 0;
	int64_t sm_b = (sm-SSFR_SM_MIN)*SSFR_SM_BPDEX;
	if ((sm_b >= 0) && (sm_b < SSFR_SM_NBINS)) {
	  //SSFR bin:
	  double ssfr = halos[i].obs_sfr / halos[i].obs_sm;
	  ssfr = (ssfr > 0) ? log10(ssfr) : 0;
	  int64_t ssfr_b = (ssfr-SSFR_SSFR_MIN)*SSFR_SSFR_BPDEX;
	  if (ssfr_b < 0) ssfr_b = 0;
	  if (ssfr_b >= SSFR_SSFR_NBINS) ssfr_b = SSFR_SSFR_NBINS-1;

	  int64_t bin = calc_smhm*SSFR_NBINS + sm_b*SSFR_SSFR_NBINS + ssfr_b;
	  ssfr_dist[bin]++;
	  if (halos[i].upid > -1) ssfr_dist_sat[bin]++;
	  else ssfr_dist_cen[bin]++;
	}
      }
    }

    //Other distributions
    if (!ehtags) {
      snprintf(buffer, 1024, "%s/tags.box%"PRId64".dat", INBASE, chunk_num);
      FILE *tags = check_fopen(buffer, "r");
      check_realloc_s(ehtags, sizeof(struct extra_halo_tags), offsets[num_scales]);
      check_fread(ehtags, sizeof(struct extra_halo_tags), offsets[num_scales], tags);
      fclose(tags);
    }

    for (i=offsets[n]; i<offsets[n+1]; i++) {
      if (halos[i].flags & IGNORE_FLAG) continue;
      true_csfr[n] += halos[i].sfr;
      
      int64_t group = 0;
      double hostm = ehtags[i].host_mass;
      if (hostm > 5.62e11*h0 && hostm < 1.78e12*h0) group = 1;
      if (hostm > 5.62e12*h0 && hostm < 1.78e13*h0) group = 2;
      if (hostm > 5.62e13*h0 && hostm < 1.78e14*h0) group = 3;

      double hm = log10(halos[i].mp/h0);
      int64_t hmb = (hm-HM_M_MIN)*HM_M_BPDEX;
      if (hmb < 0) hmb = 0;
      if (hmb >= HM_M_NBINS) hmb = HM_M_NBINS-1;
      hmb += n*HM_M_NBINS;      
      
      hm_dist[hmb]++;
      if (halos[i].upid < 0) hm_cen_dist[hmb]++;
      if (halos[i].sfr < QUENCHED_THRESHOLD*halos[i].sm) {
	hm_q[hmb]++;
	if (halos[i].upid > -1) hm_q_sat[hmb]++;
	else hm_q_cen[hmb]++;
	if (group==1) hm_q_sat_mw[hmb]++;
	if (group==2) hm_q_sat_gr[hmb]++;
	if (group==3) hm_q_sat_cl[hmb]++;
      }

      double hostm_all = (hostm > 0) ? hostm : halos[i].m;
      hostm_all = (hostm_all > 0) ? log10(hostm_all/h0) : 0;
      int64_t hostm_bin = (hostm_all - HM_M_MIN)*HM_M_BPDEX;
      if (hostm_bin < 0) hostm_bin = 0;
      if (hostm_bin >= HM_M_NBINS) hostm_bin = HM_M_NBINS-1;
      hostm_bin += n*HM_M_NBINS;
      if (halos[i].upid < 0) {
	hostm_dist[hostm_bin]++;
	total_sm_central_enclosed[hostm_bin] += halos[i].sm;
	total_sm_icl_central_enclosed[hostm_bin] += halos[i].sm + halos[i].icl;
      }
      total_sm_all_enclosed[hostm_bin] += halos[i].sm + halos[i].icl;
      
      double sm = (halos[i].obs_sm > 0) ? log10(halos[i].obs_sm) : 0;
      int64_t rsmb = (sm-SM_MIN)*SM_RED_BPDEX;
      if (rsmb < 0) rsmb = 0;
      if (rsmb >= SM_RED_NBINS) rsmb = SM_RED_NBINS-1;
      rsmb += n*SM_RED_NBINS;

      sm_red[rsmb]++;
      
      if (halos[i].upid < 0)
	sm_cen[rsmb]++;

      if (group==1) { sm_sat_mw[rsmb]++; hm_sat_mw[hmb]++; }
      if (group==2) { sm_sat_gr[rsmb]++; hm_sat_gr[hmb]++; }
      if (group==3) { sm_sat_cl[rsmb]++; hm_sat_cl[hmb]++; }
      
      if (halos[i].obs_sfr < QUENCHED_THRESHOLD*halos[i].obs_sm) {
	sm_red_q2[rsmb]++;
	if (halos[i].upid > -1) sm_q_sat[rsmb]++;
	else sm_q_cen[rsmb]++;
	if (group==1) sm_q_sat_mw[rsmb]++;
	if (group==2) sm_q_sat_gr[rsmb]++;
	if (group==3) sm_q_sat_cl[rsmb]++;
      }

      total_log_sm[hmb] += sm;
      total_log_smhm_ratio[hmb] += sm-hm;
      double true_sm = (halos[i].sm > 0) ? log10(halos[i].sm) : 0;
      double true_sm_icl = ((halos[i].sm+halos[i].icl) > 0) ? log10(halos[i].sm+halos[i].icl) : 0;
      total_log_true_sm[hmb] += true_sm;
      total_log_true_sm_icl[hmb] += true_sm_icl;
      total_log_true_smhm_ratio[hmb] += true_sm-hm;
      total_log_true_smhm_icl_ratio[hmb] += true_sm_icl-hm;

      double hm_nolog = halos[i].mp/h0;
      double wl_hm = pow(halos[i].mp/h0, 2.0/3.0);
      total_hm[rsmb] += hm_nolog;
      total_log_hm[rsmb] += hm;
      total_wl_hm[rsmb] += wl_hm;

      if (halos[i].upid < 0) {
	total_cen_hm[rsmb] += hm_nolog;
	total_cen_log_hm[rsmb] += hm;
	total_cen_wl_hm[rsmb] += wl_hm;
	if (halos[i].obs_sfr < QUENCHED_THRESHOLD*halos[i].obs_sm) {
	  total_cen_q_hm[rsmb] += hm_nolog;
	  total_cen_q_log_hm[rsmb] += hm;
	  total_cen_q_wl_hm[rsmb] += wl_hm;
	} else {
	  total_cen_sf_hm[rsmb] += hm_nolog;
	  total_cen_sf_log_hm[rsmb] += hm;
	  total_cen_sf_wl_hm[rsmb] += wl_hm;
	}
      } else {
	total_sat_hm[rsmb] += hm_nolog;
	total_sat_log_hm[rsmb] += hm;
	total_sat_wl_hm[rsmb] += wl_hm;
      }

      if (halos[i].obs_sfr < QUENCHED_THRESHOLD*halos[i].obs_sm) {
	total_q_hm[rsmb] += hm_nolog;
	total_q_log_hm[rsmb] += hm;
	total_q_wl_hm[rsmb] += wl_hm;
      } else {
	total_sf_hm[rsmb] += hm_nolog;
	total_sf_log_hm[rsmb] += hm;
	total_sf_wl_hm[rsmb] += wl_hm;
      }
      
      int64_t toff = (i-offsets[n])*(n+1);
      int64_t sat_q_after_infall = 0;
      int64_t scale_index_at_infall = -1;
      double q_tdelay = -1;
      double ssfr_at_infall = -1;
      float *tsfh = sfh_new+(i-offsets[n])*(n+1);
      float tot_in_situ = 0;
      float tot_all = 0;
      
      int64_t sf = 1;
      int64_t rejuv = 0;
      double qt=0, sft=0, qs=0, sfs=0;

      for (j=0; j<=n; j++) {
	if (halos[i].upid > -1 && (scale_index_at_infall > -1) &&
	    (ssfr_at_infall > QUENCHED_THRESHOLD) &&
	    (sfrh[toff+j] < QUENCHED_THRESHOLD*smh[toff+j]) &&
	    (q_tdelay < 0)) {
	  sat_q_after_infall = 1;
	  q_tdelay = (scale_to_years(scales[j]) - scale_to_years(a_first_infall[i-offsets[n]]))/1e9;

	  int64_t tbin = TS_BPGYR*q_tdelay;
	  if (tbin < 0) tbin = 0;
	  if (tbin >= TS_NBINS) tbin = TS_NBINS-1;
	  tbin += rsmb*TS_NBINS;
	  if (group==1) sm_td_sat_mw_infall[tbin]++;
	  if (group==2) sm_td_sat_gr_infall[tbin]++;
	  if (group==3) sm_td_sat_cl_infall[tbin]++;
	}

	  
	if (halos[i].upid > -1 && (fabs(scales[j]-a_first_infall[i-offsets[n]]) < 0.001)) {
	  int64_t bin = rsmb*num_scales+j;
	  if (group==1) sm_ti_sat_mw_infall[bin]++;
	  if (group==2) sm_ti_sat_gr_infall[bin]++;
	  if (group==3) sm_ti_sat_cl_infall[bin]++;
	  scale_index_at_infall = j;
	  ssfr_at_infall = (smh[toff+j] > 0) ? sfrh[toff+j]/smh[toff+j] : -1;
	  if (ssfr_at_infall > 0 && calc_smhm>-1) {
	    double ssfr = log10(ssfr_at_infall);
	    int64_t ssfr_b = (ssfr-SSFR_SSFR_MIN)*SSFR_SSFR_BPDEX;
	    if (ssfr_b < 0) ssfr_b = 0;
	    if (ssfr_b >= SSFR_SSFR_NBINS) ssfr_b = SSFR_SSFR_NBINS-1;
	    int64_t bin = calc_smhm*SSFR_NBINS + (rsmb-n*SM_RED_NBINS)*SSFR_SSFR_NBINS + ssfr_b;
	    if (group==1) sm_ssfr_sat_mw_infall[bin]++;
	    if (group==2) sm_ssfr_sat_gr_infall[bin]++;
	    if (group==3) sm_ssfr_sat_cl_infall[bin]++;
	  }
	}

	int64_t hbin = hmb*num_scales + j;
	int64_t sbin = rsmb*num_scales + j;
	hm_sfh[hbin] += tsfh[j];
	if (halos[i].sfr > halos[i].sm*QUENCHED_THRESHOLD) hm_sfh_sf[hbin] += tsfh[j];
	else hm_sfh_q[hbin] += tsfh[j];
	if (halos[i].upid > -1) hm_sfh_sat[hbin] += tsfh[j];
	else hm_sfh_cen[hbin] += tsfh[j];

	sm_sfh[sbin] += tsfh[j];
	if (halos[i].obs_sfr > halos[i].sm*QUENCHED_THRESHOLD) sm_sfh_sf[sbin] += tsfh[j];
	else sm_sfh_q[sbin] += tsfh[j];
	if (halos[i].upid > -1) sm_sfh_sat[sbin] += tsfh[j];
	else sm_sfh_cen[sbin] += tsfh[j];

	if (n==num_scales-1) {
	  int64_t hm_mpbin =  (hmb - n*HM_M_NBINS)*num_scales+j;
	  int64_t sm_mpbin =  (rsmb - n*SM_RED_NBINS)*num_scales+j;
	  if (halos[i].sfr > halos[i].sm*QUENCHED_THRESHOLD) {
	    mp_sf_counts[hm_mpbin]++;
	    if (sfrh[toff+j] < QUENCHED_THRESHOLD*smh[toff+j]) mp_sf_qf[hm_mpbin]++;
	  }
	  else {
	    mp_q_counts[hm_mpbin]++;
	    if (sfrh[toff+j] < QUENCHED_THRESHOLD*smh[toff+j]) mp_q_qf[hm_mpbin]++;
	  }

	  if (halos[i].obs_sfr > halos[i].obs_sm*QUENCHED_THRESHOLD) {
	    sm_mp_sf_counts[sm_mpbin]++;
	    if (sfrh[toff+j] < QUENCHED_THRESHOLD*smh[toff+j]) sm_mp_sf_qf[sm_mpbin]++;
	  }
	  else {
	    sm_mp_q_counts[sm_mpbin]++;
	    if (sfrh[toff+j] < QUENCHED_THRESHOLD*smh[toff+j]) sm_mp_q_qf[sm_mpbin]++;
	  }
	}

	double dtn = dt[j];
	if (j==n) dtn = dt2[j];
	tot_in_situ += sfrh[toff+j]*rem[n*num_scales+j]*dtn;
	tot_all += tsfh[j]*rem[n*num_scales+j]*dtn/dt[j];

	int64_t elem = (i-offsets[n])*(n+1) + j;
	double sfr = sfrh[elem];
	double sm = smh[elem];
	if (sm) {
	  int64_t cq = (sfr < QUENCHED_THRESHOLD*sm) ? 1 : 0;
	  if (cq) { qs += dt[j]; sfs = 0; }
	  else { sfs += dt[j]; qs = 0; }

	  if (sf && qs > 300e6) {
	    sf = 0;
	    qt = qs;
	    sft -= qs;
	  }
	  else if (!sf && sfs > 300e6) {
	    rejuv++;
	    sf = 1;
	    sft = sfs;
	    qt -= sfs;
	  }

	  if (sf) sft += dt[j];
	  else qt += dt[j];
	}
      }

      if (rejuv) {
	hm_rejuv[hmb]++;
	sm_rejuv[rsmb]++;
      }

      if (fabs(halos[i].sm - tot_all) > 1e-5*halos[i].sm) {
	fprintf(stderr, "[Error] Sm1: %e; Sm2: %e; Halo: %"PRId64"; snap: %"PRId64"\n", halos[i].sm, tot_all, halos[i].id, n);
	//assert(0);
      }
      
      double ex_situ_frac = (tot_in_situ < halos[i].sm && halos[i].sm > 0) ? 1.0-(tot_in_situ / halos[i].sm) : 0;
      hm_ex_situ[hmb] += ex_situ_frac;
      sm_ex_situ[rsmb] += ex_situ_frac;
            
      if ((halos[i].upid > -1) && (halos[i].obs_sfr < QUENCHED_THRESHOLD*halos[i].obs_sm)) {
	if (group==1) {
	  sm_fq_sat_mw_infall[rsmb]+=sat_q_after_infall;
	  hm_fq_sat_mw_infall[hmb]+=sat_q_after_infall;
	}
	if (group==2) {
	  sm_fq_sat_gr_infall[rsmb]+=sat_q_after_infall;
	  hm_fq_sat_gr_infall[hmb]+=sat_q_after_infall;
	}
	if (group==3) {
	  sm_fq_sat_cl_infall[rsmb]+=sat_q_after_infall;
	  hm_fq_sat_cl_infall[hmb]+=sat_q_after_infall;
	}
      }
    }
  }

  if (output_sfh) {
    snprintf(buffer, 1024, "%s/sfh_catalog_%f.%"PRId64".txt", OUTBASE, scales[n], chunk_num);
    output = check_fopen(buffer, "w");
    fprintf(output, "#ID UPID X Y Z VX VY VZ Mpeak Mnow V@Mpeak Vnow Rvir Tidal_Tdyn Rank_DVmax(Z-score) Random_Rank(Z-score) SM ICL SFR Obs_SM Obs_SFR Obs_UV A_first_infall A_last_infall SFH(1..num_scales) ICLH(1..num_scales) SM_main_progenitor(1..num_scales) ICL_main_progenitor(1..num_scales) M_main_progenitor(1..num_scales) SFR_main_progenitor(1..num_scales) V@Mpeak(1..num_scales) Rank_DVmax(Z-score;1..num_scales)\n");
    fprintf(output, "#a = %f\n", scales[n]);
    fprintf(output, "#Units: all masses in Msun (no h).\n");
    fprintf(output, "#Units: all velocities in km/s (physical, not comoving).\n");
    fprintf(output, "#Units: Rvir in kpc/h comoving.\n");
    fprintf(output, "#num_scales: %"PRId64"\n", n+1);
    fprintf(output, "#scale list:");
    for (i=0; i<n+1; i++) fprintf(output, " %f", scales[i]);
    fprintf(output, "\n");
    for (i=offsets[n]; i<offsets[n+1]; i++) {
      if (halos[i].flags & IGNORE_FLAG) continue;
      int64_t upid = halos[i].upid;
      if (upid > -1) upid = halos[upid].id;
      fprintf(output, "%"PRId64" %"PRId64" %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g", halos[i].id, upid, halos[i].pos[0], halos[i].pos[1], halos[i].pos[2], halos[i].pos[3], halos[i].pos[4], halos[i].pos[5], halos[i].mp/h0, halos[i].m/h0, halos[i].vmp, halos[i].v, halos[i].r, halos[i].t_tdyn, halos[i].rank1, halos[i].rarank, halos[i].sm, halos[i].icl, halos[i].sfr, halos[i].obs_sm, halos[i].obs_sfr, halos[i].obs_uv, a_first_infall[i-offsets[n]], a_infall[i-offsets[n]]);
      float *tsfh = sfh_new+(i-offsets[n])*(n+1);
      float *ticl = icl_new+(i-offsets[n])*(n+1);
      for (j=0; j<n+1; j++) fprintf(output, " %g", tsfh[j]/dt[j]);
      for (j=0; j<n+1; j++) fprintf(output, " %g", (CALC_ICL) ? ticl[j]/dt[j] : 0);
      int64_t toff = (i-offsets[n])*(n+1);
      for (j=0; j<n+1; j++) fprintf(output, " %g", smh[toff+j]);
      for (j=0; j<n+1; j++) fprintf(output, " %g", iclh[toff+j]);
      for (j=0; j<n+1; j++) fprintf(output, " %g", mh[toff+j]/h0);
      for (j=0; j<n+1; j++) fprintf(output, " %g", sfrh[toff+j]);
      for (j=0; j<n+1; j++) fprintf(output, " %g", vmph[toff+j]);
      for (j=0; j<n+1; j++) fprintf(output, " %g", dvmaxh[toff+j]);
      fprintf(output, "\n");
    }
    fclose(output);
  }

  if ((SHARED_DATA_OUTPUT_MODE) && (n==(num_scales-1))) {
    snprintf(buffer, 1024, "%s/shared_data.%"PRId64".bin", OUTBASE, chunk_num);
    output = check_fopen(buffer, "w");
    shared_data_to_file(the_model,output);
    fclose(output);
  }
  
  free(a_infall);
  free(hist);
}


void stack_weak_lensing(int64_t n, int64_t j) {
  int64_t chunk_num = (catalog_client) ? (catalog_client - 1) : loaded_chunk;
  int64_t wl_counts[NUM_WL_THRESHES*NUM_WL_TYPES*NUM_WL_BINS] = {0};
  int64_t wl_halo_counts[NUM_WL_THRESHES*NUM_WL_TYPES] = {0};
  int64_t i, k;
  char buffer[1024];
  assert(j<NUM_WL_ZS);
  if (!wl_data) check_calloc_s(wl_data, sizeof(struct wl_data *), NUM_WL_ZS);
  if (!wl_data[j]) {
    check_realloc_s(wl_data[j], sizeof(struct wl_data), num_halos[n]);
    snprintf(buffer, 1024, "%s/pcounts_%f.%"PRId64".bin", INBASE, scales[n], chunk_num);
    FILE *input = check_fopen(buffer, "r");
    check_fread(wl_data[j], sizeof(struct wl_data), num_halos[n], input);
    fclose(input);
  }

  int64_t thresh_off = NUM_WL_TYPES*NUM_WL_BINS;
#define ADD_COUNTS(off) { wl_halo_counts[(off)/NUM_WL_BINS]++; for (k=0; k<NUM_WL_BINS; k++) wl_counts[off+k] += wl_data[j][i-offsets[n]].pcounts[k]; }
  for (i=offsets[n]; i<offsets[n+1]; i++) {
    if (halos[i].flags & IGNORE_FLAG) continue;
    int64_t sf_off = (halos[i].obs_sfr > halos[i].obs_sm*QUENCHED_THRESHOLD) ? NUM_WL_BINS : (2*NUM_WL_BINS);
    if (halos[i].obs_sm > WL_MIN_SM) {
      ADD_COUNTS(0*thresh_off);
      ADD_COUNTS(0*thresh_off+sf_off);
    }
    if (halos[i].obs_sm > WL_MIN_SM2) {
      ADD_COUNTS(1*thresh_off);
      ADD_COUNTS(1*thresh_off+sf_off);
    }
    if (halos[i].obs_sm > WL_MIN_SM3) {
      ADD_COUNTS(2*thresh_off);
      ADD_COUNTS(2*thresh_off+sf_off);
    }
  }
#undef ADD_COUNTS

  float *fwl = wl_pcounts + (thresh_off*NUM_WL_THRESHES)*j;
  for (i=0; i<thresh_off*NUM_WL_THRESHES; i++) fwl[i] = wl_counts[i];
  float *fwh = wl_hcounts + NUM_WL_THRESHES*NUM_WL_TYPES*j;
  for (i=0; i<NUM_WL_THRESHES*NUM_WL_TYPES; i++) fwh[i] = wl_halo_counts[i];
}
