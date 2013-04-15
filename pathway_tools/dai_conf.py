
import dai

config_map = {
"BP":                             (dai.BP,dict(inference="SUMPROD",updates="SEQMAX",logdomain="0",tol="1e-9",maxiter="10000",damping="0.0")),
}

"""
"BP_SEQFIX":                      (dai.BP,dict(inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"BP_SEQRND":                      (dai.BP,dict(inference=SUMPROD,updates=SEQRND,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"BP_SEQMAX":                      (dai.BP,dict(inference=SUMPROD,updates=SEQMAX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"BP_PARALL":                      (dai.BP,dict(inference=SUMPROD,updates=PARALL,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"BP_SEQFIX_LOG":                  (dai.BP,dict(inference=SUMPROD,updates=SEQFIX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"BP_SEQRND_LOG":                  (dai.BP,dict(inference=SUMPROD,updates=SEQRND,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"BP_SEQMAX_LOG":                  (dai.BP,dict(inference=SUMPROD,updates=SEQMAX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"BP_PARALL_LOG":                  (dai.BP,dict(inference=SUMPROD,updates=PARALL,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"MP_SEQFIX":                      (dai.BP,dict(inference=MAXPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"MP_SEQRND":                      (dai.BP,dict(inference=MAXPROD,updates=SEQRND,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"MP_SEQMAX":                      (dai.BP,dict(inference=MAXPROD,updates=SEQMAX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"MP_PARALL":                      (dai.BP,dict(inference=MAXPROD,updates=PARALL,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"MP_SEQFIX_LOG":                  (dai.BP,dict(inference=MAXPROD,updates=SEQFIX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"MP_SEQRND_LOG":                  (dai.BP,dict(inference=MAXPROD,updates=SEQRND,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"MP_SEQMAX_LOG":                  (dai.BP,dict(inference=MAXPROD,updates=SEQMAX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"MP_PARALL_LOG":                  (dai.BP,dict(inference=MAXPROD,updates=PARALL,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),

# --- FBP ---------------------

"FBP":                            (dai.FBP,dict(inference=SUMPROD,updates=SEQMAX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),

"FBP_SEQFIX":                     (dai.FBP,dict(inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"FBP_SEQRND":                     (dai.FBP,dict(inference=SUMPROD,updates=SEQRND,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"FBP_SEQMAX":                     (dai.FBP,dict(inference=SUMPROD,updates=SEQMAX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"FBP_PARALL":                     (dai.FBP,dict(inference=SUMPROD,updates=PARALL,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"FBP_SEQFIX_LOG":                 (dai.FBP,dict(inference=SUMPROD,updates=SEQFIX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"FBP_SEQRND_LOG":                 (dai.FBP,dict(inference=SUMPROD,updates=SEQRND,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"FBP_SEQMAX_LOG":                 (dai.FBP,dict(inference=SUMPROD,updates=SEQMAX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"FBP_PARALL_LOG":                 (dai.FBP,dict(inference=SUMPROD,updates=PARALL,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"FMP_SEQFIX":                     (dai.FBP,dict(inference=MAXPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"FMP_SEQRND":                     (dai.FBP,dict(inference=MAXPROD,updates=SEQRND,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"FMP_SEQMAX":                     (dai.FBP,dict(inference=MAXPROD,updates=SEQMAX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"FMP_PARALL":                     (dai.FBP,dict(inference=MAXPROD,updates=PARALL,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0)),
"FMP_SEQFIX_LOG":                 (dai.FBP,dict(inference=MAXPROD,updates=SEQFIX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"FMP_SEQRND_LOG":                 (dai.FBP,dict(inference=MAXPROD,updates=SEQRND,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"FMP_SEQMAX_LOG":                 (dai.FBP,dict(inference=MAXPROD,updates=SEQMAX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),
"FMP_PARALL_LOG":                 (dai.FBP,dict(inference=MAXPROD,updates=PARALL,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0)),

# --- TRWBP -------------------

"TRWBP":                          (dai.TRWBP,dict(updates=SEQFIX,tol=1e-9,maxiter=10000,logdomain=0,nrtrees=0)),

"TRWBP_SEQFIX":                   (dai.TRWBP,dict(inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWBP_SEQRND":                   (dai.TRWBP,dict(inference=SUMPROD,updates=SEQRND,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWBP_SEQMAX":                   (dai.TRWBP,dict(inference=SUMPROD,updates=SEQMAX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWBP_PARALL":                   (dai.TRWBP,dict(inference=SUMPROD,updates=PARALL,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWBP_SEQFIX_LOG":               (dai.TRWBP,dict(inference=SUMPROD,updates=SEQFIX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWBP_SEQRND_LOG":               (dai.TRWBP,dict(inference=SUMPROD,updates=SEQRND,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWBP_SEQMAX_LOG":               (dai.TRWBP,dict(inference=SUMPROD,updates=SEQMAX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWBP_PARALL_LOG":               (dai.TRWBP,dict(inference=SUMPROD,updates=PARALL,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWMP_SEQFIX":                   (dai.TRWBP,dict(inference=MAXPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWMP_SEQRND":                   (dai.TRWBP,dict(inference=MAXPROD,updates=SEQRND,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWMP_SEQMAX":                   (dai.TRWBP,dict(inference=MAXPROD,updates=SEQMAX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWMP_PARALL":                   (dai.TRWBP,dict(inference=MAXPROD,updates=PARALL,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWMP_SEQFIX_LOG":               (dai.TRWBP,dict(inference=MAXPROD,updates=SEQFIX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWMP_SEQRND_LOG":               (dai.TRWBP,dict(inference=MAXPROD,updates=SEQRND,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWMP_SEQMAX_LOG":               (dai.TRWBP,dict(inference=MAXPROD,updates=SEQMAX,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),
"TRWMP_PARALL_LOG":               (dai.TRWBP,dict(inference=MAXPROD,updates=PARALL,logdomain=1,tol=1e-9,maxiter=10000,damping=0.0,nrtrees=0)),

# --- JTREE -------------------

"JTREE_HUGIN":                    (dai.JTREE,dict(inference=SUMPROD,updates=HUGIN)),
"JTREE_SHSH":                     (dai.JTREE,dict(inference=SUMPROD,updates=SHSH)),

"JTREE_MINFILL_HUGIN":            (dai.JTREE,dict(inference=SUMPROD,heuristic=MINFILL,updates=HUGIN)),
"JTREE_MINFILL_SHSH":             (dai.JTREE,dict(inference=SUMPROD,heuristic=MINFILL,updates=SHSH)),
"JTREE_MINFILL_HUGIN_MAP":        (dai.JTREE,dict(inference=MAXPROD,heuristic=MINFILL,updates=HUGIN)),
"JTREE_MINFILL_SHSH_MAP":         (dai.JTREE,dict(inference=MAXPROD,heuristic=MINFILL,updates=SHSH)),
"JTREE_WEIGHTEDMINFILL_HUGIN":    (dai.JTREE,dict(inference=SUMPROD,heuristic=WEIGHTEDMINFILL,updates=HUGIN)),
"JTREE_WEIGHTEDMINFILL_SHSH":     (dai.JTREE,dict(inference=SUMPROD,heuristic=WEIGHTEDMINFILL,updates=SHSH)),
"JTREE_WEIGHTEDMINFILL_HUGIN_MAP":(dai.JTREE,dict(inference=MAXPROD,heuristic=WEIGHTEDMINFILL,updates=HUGIN)),
"JTREE_WEIGHTEDMINFILL_SHSH_MAP": (dai.JTREE,dict(inference=MAXPROD,heuristic=WEIGHTEDMINFILL,updates=SHSH)),
"JTREE_MINWEIGHT_HUGIN":          (dai.JTREE,dict(inference=SUMPROD,heuristic=MINWEIGHT,updates=HUGIN)),
"JTREE_MINWEIGHT_SHSH":           (dai.JTREE,dict(inference=SUMPROD,heuristic=MINWEIGHT,updates=SHSH)),
"JTREE_MINWEIGHT_HUGIN_MAP":      (dai.JTREE,dict(inference=MAXPROD,heuristic=MINWEIGHT,updates=HUGIN)),
"JTREE_MINWEIGHT_SHSH_MAP":       (dai.JTREE,dict(inference=MAXPROD,heuristic=MINWEIGHT,updates=SHSH)),
"JTREE_MINNEIGHBORS_HUGIN":       (dai.JTREE,dict(inference=SUMPROD,heuristic=MINNEIGHBORS,updates=HUGIN)),
"JTREE_MINNEIGHBORS_SHSH":        (dai.JTREE,dict(inference=SUMPROD,heuristic=MINNEIGHBORS,updates=SHSH)),
"JTREE_MINNEIGHBORS_HUGIN_MAP":   (dai.JTREE,dict(inference=MAXPROD,heuristic=MINNEIGHBORS,updates=HUGIN)),
"JTREE_MINNEIGHBORS_SHSH_MAP":    (dai.JTREE,dict(inference=MAXPROD,heuristic=MINNEIGHBORS,updates=SHSH)),

# --- MF ----------------------

"MF":                             (dai.MF,dict(tol=1e-9,maxiter=10000,damping=0.0,init=UNIFORM,updates=NAIVE)),

"MF_NAIVE_UNI":                   (dai.MF,dict(tol=1e-9,maxiter=10000,damping=0.0,init=UNIFORM,updates=NAIVE)),
"MF_NAIVE_RND":                   (dai.MF,dict(tol=1e-9,maxiter=10000,damping=0.0,init=RANDOM,updates=NAIVE)),
"MF_HARDSPIN_UNI":                (dai.MF,dict(tol=1e-9,maxiter=10000,damping=0.0,init=UNIFORM,updates=HARDSPIN)),
"MF_HARDSPIN_RND":                (dai.MF,dict(tol=1e-9,maxiter=10000,damping=0.0,init=RANDOM,updates=HARDSPIN)),

# --- TREEEP ------------------

"TREEEP":                         (dai.TREEEP,dict(type=ORG,tol=1e-9,maxiter=10000)),
"TREEEPWC":                       (dai.TREEEP,dict(type=ALT,tol=1e-9,maxiter=10000)),

# --- MR ----------------------

"MR_DEFAULT":                     (dai.MR,dict(updates=FULL,inits=RESPPROP,tol=1e-9)),
"MR_RESPPROP_FULL":               (dai.MR,dict(updates=FULL,inits=RESPPROP,tol=1e-9)),
"MR_RESPPROP_LINEAR":             (dai.MR,dict(updates=LINEAR,inits=RESPPROP,tol=1e-9)),
"MR_CLAMPING_FULL":               (dai.MR,dict(updates=FULL,inits=CLAMPING,tol=1e-9)),
"MR_CLAMPING_LINEAR":             (dai.MR,dict(updates=LINEAR,inits=CLAMPING,tol=1e-9)),
"MR_EXACT_FULL":                  (dai.MR,dict(updates=FULL,inits=EXACT,tol=1e-9)),
"MR_EXACT_LINEAR":                (dai.MR,dict(updates=LINEAR,inits=EXACT,tol=1e-9)),

# --- HAK ---------------------

"GBP_MIN":                        (dai.HAK,dict(doubleloop=0,clusters=MIN,init=UNIFORM,tol=1e-9,maxiter=10000)),
"GBP_BETHE":                      (dai.HAK,dict(doubleloop=0,clusters=BETHE,init=UNIFORM,tol=1e-9,maxiter=10000)),
"GBP_DELTA":                      (dai.HAK,dict(doubleloop=0,clusters=DELTA,init=UNIFORM,tol=1e-9,maxiter=10000)),
"GBP_LOOP3":                      (dai.HAK,dict(doubleloop=0,clusters=LOOP,init=UNIFORM,loopdepth=3,tol=1e-9,maxiter=10000)),
"GBP_LOOP4":                      (dai.HAK,dict(doubleloop=0,clusters=LOOP,init=UNIFORM,loopdepth=4,tol=1e-9,maxiter=10000)),
"GBP_LOOP5":                      (dai.HAK,dict(doubleloop=0,clusters=LOOP,init=UNIFORM,loopdepth=5,tol=1e-9,maxiter=10000)),
"GBP_LOOP6":                      (dai.HAK,dict(doubleloop=0,clusters=LOOP,init=UNIFORM,loopdepth=6,tol=1e-9,maxiter=10000)),
"GBP_LOOP7":                      (dai.HAK,dict(doubleloop=0,clusters=LOOP,init=UNIFORM,loopdepth=7,tol=1e-9,maxiter=10000)),
"GBP_LOOP8":                      (dai.HAK,dict(doubleloop=0,clusters=LOOP,init=UNIFORM,loopdepth=8,tol=1e-9,maxiter=10000)),

"HAK_MIN":                        (dai.HAK,dict(doubleloop=1,clusters=MIN,init=UNIFORM,tol=1e-9,maxiter=10000)),
"HAK_BETHE":                      (dai.HAK,dict(doubleloop=1,clusters=BETHE,init=UNIFORM,tol=1e-9,maxiter=10000)),
"HAK_DELTA":                      (dai.HAK,dict(doubleloop=1,clusters=DELTA,init=UNIFORM,tol=1e-9,maxiter=10000)),
"HAK_LOOP3":                      (dai.HAK,dict(doubleloop=1,clusters=LOOP,init=UNIFORM,loopdepth=3,tol=1e-9,maxiter=10000)),
"HAK_LOOP4":                      (dai.HAK,dict(doubleloop=1,clusters=LOOP,init=UNIFORM,loopdepth=4,tol=1e-9,maxiter=10000)),
"HAK_LOOP5":                      (dai.HAK,dict(doubleloop=1,clusters=LOOP,init=UNIFORM,loopdepth=5,tol=1e-9,maxiter=10000)),
"HAK_LOOP6":                      (dai.HAK,dict(doubleloop=1,clusters=LOOP,init=UNIFORM,loopdepth=6,tol=1e-9,maxiter=10000)),
"HAK_LOOP7":                      (dai.HAK,dict(doubleloop=1,clusters=LOOP,init=UNIFORM,loopdepth=7,tol=1e-9,maxiter=10000)),
"HAK_LOOP8":                      (dai.HAK,dict(doubleloop=1,clusters=LOOP,init=UNIFORM,loopdepth=8,tol=1e-9,maxiter=10000)),

# --- LC ----------------------
"LCBP_FULLCAVin_SEQFIX":          (dai.LC,dict(cavity=FULL,reinit=1,updates=SEQFIX,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_FULLCAVin_SEQRND":          (dai.LC,dict(cavity=FULL,reinit=1,updates=SEQRND,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_FULLCAVin_NONE":            (dai.LC,dict(cavity=FULL,reinit=1,updates=SEQFIX,maxiter=0,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_FULLCAV_SEQFIX":            (dai.LC,dict(cavity=FULL,reinit=0,updates=SEQFIX,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_FULLCAV_SEQRND":            (dai.LC,dict(cavity=FULL,reinit=0,updates=SEQRND,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_FULLCAV_NONE":              (dai.LC,dict(cavity=FULL,reinit=0,updates=SEQFIX,maxiter=0,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIRCAVin_SEQFIX":          (dai.LC,dict(cavity=PAIR,reinit=1,updates=SEQFIX,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIRCAVin_SEQRND":          (dai.LC,dict(cavity=PAIR,reinit=1,updates=SEQRND,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIRCAVin_NONE":            (dai.LC,dict(cavity=PAIR,reinit=1,updates=SEQFIX,maxiter=0,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIRCAV_SEQFIX":            (dai.LC,dict(cavity=PAIR,reinit=0,updates=SEQFIX,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIRCAV_SEQRND":            (dai.LC,dict(cavity=PAIR,reinit=0,updates=SEQRND,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIRCAV_NONE":              (dai.LC,dict(cavity=PAIR,reinit=0,updates=SEQFIX,maxiter=0,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIR2CAVin_SEQFIX":         (dai.LC,dict(cavity=PAIR2,reinit=1,updates=SEQFIX,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIR2CAVin_SEQRND":         (dai.LC,dict(cavity=PAIR2,reinit=1,updates=SEQRND,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIR2CAVin_NONE":           (dai.LC,dict(cavity=PAIR2,reinit=1,updates=SEQFIX,maxiter=0,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIR2CAV_SEQFIX":           (dai.LC,dict(cavity=PAIR2,reinit=0,updates=SEQFIX,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIR2CAV_SEQRND":           (dai.LC,dict(cavity=PAIR2,reinit=0,updates=SEQRND,maxiter=10000,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_PAIR2CAV_NONE":             (dai.LC,dict(cavity=PAIR2,reinit=0,updates=SEQFIX,maxiter=0,cavainame=BP,cavaiopts=dai.dict(updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0)),,tol=1e-9)),
"LCBP_UNICAV_SEQFIX":             (dai.LC,dict(cavity=UNIFORM,updates=SEQFIX,maxiter=10000,tol=1e-9,cavaiopts=dai.dict()),,cavainame=NONE)),
"LCBP_UNICAV_SEQRND":             (dai.LC,dict(cavity=UNIFORM,updates=SEQRND,maxiter=10000,tol=1e-9,cavaiopts=dai.dict()),,cavainame=NONE)),

"LCTREEEP":                       (dai.LC,dict(cavity=FULL,reinit=1,updates=SEQFIX,maxiter=10000,cavainame=TREEEP,cavaiopts="[type=ORG,tol=1e-9,maxiter=10000]",tol=1e-9)),
"LCMF":                           (dai.LC,dict(cavity=FULL,reinit=1,updates=SEQFIX,maxiter=10000,cavainame=MF,cavaiopts=dai.dict(tol=1e-9,maxiter=10000)),,tol=1e-9)),
"LCBP":                           LCBP_FULLCAVin_SEQRND

# --- GIBBS -------------------

"GIBBS":                          (dai.GIBBS,dict(maxiter=10000,burnin=100,restart=10000)),

# --- CBP ---------------------

"CBP":                            (dai.CBP,dict(max_levels=12,updates=SEQMAX,tol=1e-9,rec_tol=1e-9,maxiter=500,choose=CHOOSE_RANDOM,recursion=REC_FIXED,clamp=CLAMP_VAR,min_max_adj=1.0e-9,bbp_cfn=CFN_FACTOR_ENT,rand_seed=0,bbp_props="[tol=1.0e-9,maxiter=10000,damping=0,updates=SEQ_BP_REV]", clamp_outfile="")),
"BBP":                            (dai.CBP,dict(choose=CHOOSE_BBP)),

# --- DECMAP ------------------

"DECMAP":				(dai.DECMAP,dict(ianame=BP,iaopts="[inference=MAXPROD,updates=SEQRND,logdomain=1,tol=1e-9,maxiter=10000,damping=0.1,verbose=0]", reinit=1,verbose=0)),

#-------GLC(+)--------------
"GLC_UNICAV_SINGLE":              (dai.GLC, dict=(rgntype=SINGLE,cavity=UNIFORM,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame="EXACT",inaiopts="[]",tol=1e-9,verbose=0)),
"GLC_PAIRCAV_SINGLE":             (dai.GLC, dict=(rgntype=SINGLE,cavity=PAIR,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame="EXACT",inaiopts=[],tol=1e-9)),
"GLC_PAIR2CAV_SINGLE":            (dai.GLC, dict=(rgntype=SINGLE,cavity=PAIR2,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame="EXACT",inaiopts=[],tol=1e-9)),
"GLC_FULLCAV_SINGLE":             (dai.GLC, dict=(rgntype=SINGLE,cavity=FULL,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame=EXACT,inaiopts=[],tol=1e-9)),
"GLC+_UNICAV_FACTOR":             (dai.GLC, dict=(rgntype=OVFACTOR,cavity=UNIFORM,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame=JTREE,inaiopts=[inference=SUMPROD,updates=HUGIN],tol=1e-9)),
"GLC+_FULLCAV_FACTOR":            (dai.GLC, dict=(rgntype=OVFACTOR,cavity=FULL,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame=EXACT,inaiopts=[],tol=1e-9)),
"GLC+_UNICAV_LOOP3":              (dai.GLC, dict=(rgntype=OVLOOP,loopdepth=3,cavity=UNIFORM,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame=JTREE,inaiopts=[inference=SUMPROD,updates=HUGIN],tol=1e-9)),
"GLC+_FULLCAV_LOOP3":             (dai.GLC, dict=(rgntype=OVLOOP,loopdepth=3,cavity=FULL,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame=EXACT,inaiopts=[],tol=1e-9)),
"GLC_UNICAV_LOOP4":               (dai.GLC, dict=(rgntype=LOOP,loopdepth=4,cavity=UNIFORM,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame=JTREE,inaiopts=[inference=SUMPROD,updates=HUGIN],tol=1e-9)),
"GLC_FULLCAV_LOOP4":              (dai.GLC, dict=(rgntype=LOOP,loopdepth=4,cavity=UNIFORM,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame=EXACT,inaiopts=[],tol=1e-9])),
"GLC+_UNICAV_LOOP4":              (dai.GLC, dict=(rgntype=OVLOOP,loopdepth=4,cavity=UNIFORM,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame=JTREE,inaiopts=[inference=SUMPROD,updates=HUGIN],tol=1e-9)),
"GLC+_FULLCAV_LOOP4":             (dai.GLC, dict=(rgntype=OVLOOP,loopdepth=4,cavity=FULL,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame=EXACT,inaiopts=[],tol=1e-9)),
"GLC+_UNICAV_LOOP5":              (dai.GLC, dict=(rgntype=OVLOOP,loopdepth=5,cavity=UNIFORM,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame=JTREE,inaiopts=[inference=SUMPROD,updates=HUGIN],tol=1e-9)),
"GLC+_FULLCAV_LOOP5":             (dai.GLC, dict=(rgntype=OVLOOP,loopdepth=5,cavity=FULL,updates=SEQRND,maxiter=10000,tol=1e-9,cavainame=BP,cavaiopts="[updates=SEQMAX,tol=1e-9,maxiter=10000,logdomain=0]",inainame=EXACT,inaiopts=[],tol=1e-9))
}

"""
