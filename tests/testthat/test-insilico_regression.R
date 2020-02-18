context("test-insilico_regression")

test_that("regression comparison", {
  
  load("insilico_regression.RData")
  
  # run without sub-population
  
  # first test that a warning of missing symptoms exists and is equivalent to expectation
  new_run <- expect_warning(insilico( data, subpop = NULL,
                       Nsim = 400, burnin = 200, thin = 10 , seed = 13,
                       auto.length = FALSE, warning.write = FALSE), regexp="78 symptom missing completely and added to missing list 
List of missing symptoms: 
 child, under5, infant, neonate, magegp1, died_d1, died_d23, died_d36, died_w1, whoop, chest_in, eye_sunk, e_bleed, s_bleed, p_bleed, placent_r, born_early, born_3437, born_38, ab_size, born_small, born_big, twin, comdel, cord, waters, move_lb, cyanosis, baby_br, born_nobr, cried, no_life, mushy, fed_d1, st_suck, ab_posit, conv_d1, conv_d2, arch_b, font_hi, font_lo, unw_d1, unw_d2, cold, umbinf, b_yellow, devel, born_malf, mlf_bk, mlf_lh, mlf_sh, mttv, b_norm, b_assist, b_caes, b_first, b_more4, b_mbpr, b_msmds, b_mcon, b_mbvi, b_mvbl, b_bfac, b_bhome, b_bway, b_bprof, o_trans, fall, smoking, married, vaccin, treat, t_ort, t_iv, blood_tr, t_ngt, antib_i, surgery",fixed=T)
  
  # regression test to see if the results match expectation
  expect_equal(new_run, regress)
  
  
})
