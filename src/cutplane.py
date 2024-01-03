###############################################################################
#
# This code was written and is being maintained by Matias Villagra,
# PhD Student in Operations Research @ Columbia, supervised by 
# Daniel Bienstock.
#
# Please report any bugs or issues (for sure there will be) to
#                         mjv2153@columbia.edu
#
# Oct 2023
###############################################################################

import sys
import math
from log import danoLogger
from gurobipy import *
import numpy as np
import bisect
from myutils import breakexit
import reader
import time
import math
from cuts import *
from tools import *
import random

      
def gocutplane(log, all_data):

  formulation_start = time.time()
  themodel          = Model("Cutplane")
  buses             = all_data['buses']
  numbuses          = all_data['numbuses']
  branches          = all_data['branches']
  numbranches       = all_data['numbranches']
  gens              = all_data['gens']
  IDtoCountmap      = all_data['IDtoCountmap']
  FeasibilityTol    = all_data['FeasibilityTol']
  threshold         = all_data['threshold']

  ############################ LOAD SOLUTION ##################################

  if all_data['matpower_sol']:
    getsol_matpower(log,all_data)

  if all_data['knitro_sol']:
    #getsol_knitro(log,all_data)
    getsol_knitro2(log,all_data)

  ################################ VARIABLES ##################################

  cvar    = {}
  svar    = {}
  Pvar_f  = {}
  Qvar_f  = {}
  Pvar_t  = {}
  Qvar_t  = {}
  Pinjvar = {}
  Qinjvar = {}
  GenPvar = {}
  GenQvar = {}
  GenTvar = {}

  log.joint(' creating variables...\n')

  varcount = 0

  #cbus, bus-injection, GenP, GenQ, and GenT variables
  for bus in buses.values():
    maxprod = bus.Vmax*bus.Vmax
    minprod = bus.Vmin*bus.Vmin
          
    ubound = maxprod
    lbound = minprod

    Pubound, Plbound, Qubound, Qlbound = computebalbounds(log, all_data, bus)  
          
    cvar[bus] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                name = "c_" + str(bus.nodeID) + "_" 
                                + str(bus.nodeID))
    Pinjvar[bus] = themodel.addVar(obj = 0.0, lb = Plbound, ub = Pubound, 
                                   name = "IP_"+str(bus.nodeID))
    Qinjvar[bus] = themodel.addVar(obj = 0.0, lb = Qlbound, ub = Qubound, 
                                   name = "IQ_"+str(bus.nodeID))
    
    varcount += 3

    for genid in bus.genidsbycount:
      gen = gens[genid]

      lower = gen.Pmin*gen.status
      upper = gen.Pmax*gen.status

      GenPvar[gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, 
                                     name = "GP_" + str(gen.count) + "_" 
                                     + str(gen.nodeID))
      lower = gen.Qmin*gen.status
      upper = gen.Qmax*gen.status

      if bus.nodetype == 3:
        upper = GRB.INFINITY
        lower = -GRB.INFINITY

      GenQvar[gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, 
                                     name = "GQ_" + str(gen.count) + "_" 
                                     + str(gen.nodeID))

      varcount += 2

      #if model has quadratic objective and we linearize
      if( gen.costdegree == 2 and gen.costvector[0] != 0 and 
          (all_data['linear_objective'] or all_data['hybrid']) ):

        GenTvar[gen] = themodel.addVar(obj = 0.0, lb = 0.0, ub = GRB.INFINITY,
                                       name = 't_g_' + str(gen.count) + '_' 
                                       + str(gen.nodeID))         
        varcount += 1



  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]
    maxprod     = buses[count_of_f].Vmax*buses[count_of_t].Vmax
    minprod     = buses[count_of_f].Vmin*buses[count_of_t].Vmin

    #c variables
    ubound =  maxprod
    lbound =  minprod*math.cos(branch.maxangle_rad)

    if branch.upperanglenone == 1: 
      ubound = maxprod
      lbound = 0

    if branch.loweranglenone == 1:
      ubound = maxprod
      lbound = 0                                                                         

    cvar[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                   name = "c_" + str(branchcount) + "_" 
                                   + str(f) + "_" + str(t))

    #s variables
    if branch.maxangle_rad > 0:      
      ubound = maxprod*math.sin(branch.maxangle_rad)
    else:
      ubound = minprod*math.sin(branch.maxangle_rad)        
    if branch.minangle_rad <= 0:
      lbound = maxprod*math.sin(branch.minangle_rad)
    else:
      lbound = minprod*math.sin(branch.minangle_rad)

    if branch.upperanglenone == 1:
      ubound = maxprod
    if branch.loweranglenone == 1:
      lbound = -maxprod

    svar[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                   name = "s_" + str(branchcount) + "_" 
                                   + str(f) + "_" + str(t))

    varcount +=2
    

  
  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    ubound = branch.limit 
    lbound = -branch.limit

    Pvar_f[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                     name = "P_" + str(branch.count) + "_" 
                                     + str(f) + "_" + str(t))
    Pvar_t[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                     name = "P_" + str(branch.count) + "_" 
                                     + str(t) + "_" + str(f))
    Qvar_f[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                     name = "Q_" + str(branch.count) + "_" 
                                     + str(f) + "_" + str(t))
    Qvar_t[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                     name = "Q_" + str(branch.count) + "_" 
                                     + str(t) + "_" + str(f))

    varcount +=4 

  #i2 variables
  if all_data['i2']:
    i2var_f = {}
    #i2var_t = {}

    for branch in branches.values():
      branchcount = branch.count
      f           = branch.f
      t           = branch.t
      count_of_f  = IDtoCountmap[f]
      count_of_t  = IDtoCountmap[t]
      bus_f       = buses[count_of_f]
      #bus_t       = buses[count_of_t]

      upperbound_f = branch.limit**2 / (bus_f.Vmin * bus_f.Vmin)
      #upperbound_t = branch.limit**2 / (bus_t.Vmin * bus_t.Vmin)

      i2var_f[branch] = themodel.addVar(obj = 0.0, lb = 0, ub = upperbound_f ,
                                        name = "i2_" + str(branch.count) + "_"
                                        + str(f) + "_" + str(t))

      #i2var_t[branch] = themodel.addVar(obj = 0.0, lb = 0, ub = upperbound_t,
                                        #name = "i2_" + str(branch.count) + "_"
                                        #+ str(t) + "_" + str(f))
    
      varcount += 1
      
  themodel.update()

  log.joint('   %d variables added\n' %varcount)

  all_data['themodel']    = themodel
  all_data['cvar']        = cvar
  all_data['svar']        = svar
  all_data['GenPvar']     = GenPvar
  all_data['GenTvar']     = GenTvar
  all_data['Pvar_f']      = Pvar_f
  all_data['Pvar_t']      = Pvar_t
  all_data['Qvar_f']      = Qvar_f
  all_data['Qvar_t']      = Qvar_t

  if all_data['i2']:
    all_data['i2var_f']   = i2var_f
    #all_data['i2var_t']  = i2var_f

  ############################## OBJECTIVE ####################################

  log.joint(' creating objective...\n')
  #varobjcount = 0

  #constant term
  constobjval = 0  
  for gen in gens.values():  
    if gen.status > 0:
      constobjval += gen.costvector[gen.costdegree]
  
  constvar = themodel.addVar(obj = constobjval, lb = 1.0, ub = 1.0, 
                             name = "constant")
  
  #varobjcount +=1

  #quad terms
  objvar   = themodel.addVar(obj = 1.0, lb = -GRB.INFINITY, ub = GRB.INFINITY,
                             name = "objvar")
  qcostvar = themodel.addVar(obj = 1.0, lb = 0, ub = GRB.INFINITY, 
                             name = "qcostvar") 
  
  #varobjcount +=2

  if all_data['linear_objective'] == 0 or all_data['hybrid']:  
    objvar.setAttr("Obj",0)
    qcostexpr = QuadExpr()
    for gen in gens.values():
      if gen.costdegree == 2 and gen.costvector[0] != 0:
        qcostexpr += gen.costvector[0]*GenPvar[gen]*GenPvar[gen]
    qcost = themodel.addConstr(qcostexpr <= qcostvar, name = "qcost")
    all_data['qcost'] = qcost
    
  #linear terms
  if all_data['linear_objective'] or all_data['hybrid']:
    lincostvar = themodel.addVar(obj = 0.0, lb = -GRB.INFINITY, 
                                 ub = GRB.INFINITY, name = "lincostvar")
    lincost    = themodel.addConstr(lincostvar <= objvar, name= "lincost")
    all_data['lincost'] = lincost
  else:
    lincostvar = themodel.addVar(obj = 1.0, lb = -GRB.INFINITY, 
                                 ub = GRB.INFINITY, name = "lincostvar")

  #varobjcount += 1

  coeff       = [gen.costvector[gen.costdegree-1] for gen in gens.values()]
  variables   = [GenPvar[gen] for gen in gens.values()]
  lincostexpr = LinExpr(coeff, variables)
  lincostdef  = themodel.addConstr(lincostexpr == lincostvar, 
                                   name= "lincostdef")

  #sumTvars if using linear_objective
  if all_data['linear_objective'] or all_data['hybrid']:
    qcostvar.setAttr("Obj",0.0)
    sumTvars = LinExpr()
    for gen in gens.values():
      if gen.costdegree == 2 and gen.costvector[0] != 0:
        sumTvars += GenTvar[gen]
    sumTconstr = themodel.addConstr(sumTvars <= objvar, name='obj_var_quad')
    all_data['sumTconstr'] = sumTconstr

  if all_data['hybrid']: #we start with full objective
    objvar.setAttr("Obj",0)
    lincostvar.setAttr("Obj",1)
    qcostvar.setAttr("Obj",1)
    themodel.remove(sumTconstr)
    themodel.remove(lincost)
    
  if all_data['linear_objective']:
    log.joint('  linear objective added\n')
  elif all_data['hybrid']:
    log.joint('  hybrid algorithm, linear and quadratic objectives created\n')
  else:
    log.joint('  quadratic objective added\n')

  themodel.update()
  
  all_data['objvar']     = objvar
  all_data['lincostvar'] = lincostvar
  all_data['qcostvar']   = qcostvar

  ############################# CONSTRAINTS ###################################

  log.joint(' creating constraints...\n')

  constrcount = 0
  count       = 0

  #definition flow variables
  log.joint('  active power flow variables definition\n')

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    if branch.status == 0:
      log.joint(' branch ' + str(branch.count) + ' f ' + str(f) + ' t ' 
                + str(t) + ' is OFF\n')
      breakexit('check, reader does not include branches with status = 0')

    #  Gff cff + Gft cft + Bft sft
    constrname = "Pdef_"+str(branch.count)+"_"+str(f)+"_"+str(t)
    expr = LinExpr()
    expr += branch.Gff*cvar[buses[count_of_f]]
    expr += branch.Gft*cvar[branch]
    expr += branch.Bft*svar[branch]
    
    themodel.addConstr(expr == Pvar_f[branch], name = constrname)

    #  Gtt ctt + Gtf cft + Btf stf = Gtt ctt + Gtf cft - Btf sft
    constrname = "Pdef_"+str(branch.count)+"_"+str(t)+"_"+str(f)
    expr = LinExpr()
    expr += branch.Gtt*cvar[buses[count_of_t]]
    expr += branch.Gtf*cvar[branch]
    expr += -branch.Btf*svar[branch]

    themodel.addConstr(expr == Pvar_t[branch], name = constrname)
    
    constrcount += 2
    count       += 2

  log.joint('   %d active power flow definition constraints added\n'%count)

  log.joint('  reactive power flow variables definition\n')
  count = 0

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    constrname = "Qdef_"+str(branch.count)+"_"+str(f)+"_"+str(t)
    
    # -Bff cff - Bft cft + Gft sft
    expr = LinExpr()
    expr += -branch.Bff*cvar[buses[count_of_f]]
    expr += -branch.Bft*cvar[branch]
    expr += +branch.Gft*svar[branch]

    themodel.addConstr(expr == Qvar_f[branch], name = constrname)

    # -Btt ctt - Btf cft + Gtf stf = -Btt ctt - Btf cft - Gtf sft 
    constrname = "Qdef_"+str(branch.count)+"_"+str(t)+"_"+str(f)
    expr = LinExpr()
    expr += -branch.Btt*cvar[buses[count_of_t]]
    expr += -branch.Btf*cvar[branch]
    expr += -branch.Gtf*svar[branch]

    themodel.addConstr(expr == Qvar_t[branch], name = constrname)

    constrcount += 2
    count       += 2
  log.joint('   %d reactive power flow definition constraints added\n'%count)

  #balance constraints
  log.joint('  active power injection constraints\n')
  count = 0

  for bus in buses.values():
    constrname = "PBaldef"+str(bus.nodeID)
    expr = LinExpr()
    for branchid in bus.frombranchids.values():
      expr += Pvar_f[ branches[branchid] ]
    for branchid in bus.tobranchids.values():
      expr += Pvar_t[ branches[branchid] ]

    if ( (bus.Gs != 0) and ( ( len(bus.frombranchids) != 0 ) 
                             or ( len(bus.tobranchids) != 0 ) ) ):
      expr += bus.Gs*cvar[bus]

    themodel.addConstr(expr == Pinjvar[bus], name = constrname)

    constrcount += 1
    count       += 1
  
  log.joint('   %d active power injection constraints added\n'%count)
  log.joint('  reactive power injection constraints\n')
  count = 0

  for bus in buses.values():
    constrname = "QBaldef"+str(bus.nodeID)
    expr = LinExpr()
    for branchid in bus.frombranchids.values():
      expr += Qvar_f[ branches[branchid] ]
    for branchid in bus.tobranchids.values():
      expr += Qvar_t[ branches[branchid] ]
 
    if ( (bus.Bs != 0) and ( ( len(bus.frombranchids) != 0 ) 
                             or ( len(bus.tobranchids) != 0 ) ) ):
      expr += (-bus.Bs)*cvar[bus]

    themodel.addConstr(expr == Qinjvar[bus], name = constrname)
    
    constrcount += 1
    count       += 1
  log.joint('   %d reactive power injection constraints added\n'%count)

  #definition bus-injection variables

  log.joint('  adding injection definition constraints...\n')

  count = 0

  for bus in buses.values():
    constrname = "Bus_PInj_"+str(bus.nodeID)
    expr = LinExpr()

    if len(bus.genidsbycount) > 0:
      for genid in bus.genidsbycount:
        gen = gens[genid]
        expr += GenPvar[gen]

    themodel.addConstr(Pinjvar[bus] == expr - bus.Pd, name = constrname)

    constrname = "Bus_QInj_"+str(bus.nodeID)
    expr = LinExpr()

    if len(bus.genidsbycount) > 0:
      for genid in bus.genidsbycount:
        gen = gens[genid]
        expr += GenQvar[gen]

    themodel.addConstr(Qinjvar[bus] == expr - bus.Qd, name = constrname)

    constrcount += 2
    count       += 2

  log.joint('   %d power injection definitions added\n'%count)
  
  #definition i2 variables
  if all_data['i2']:
    constrcount += i2_def(log,all_data)

  #active power loss-inequalities
  if all_data['loss_inequalities']:    
    constrcount += loss_inequalities(log,all_data)

  #reactive power loss-inequalities
  if all_data['qloss_inequalities']:    
    constrcount += qloss_inequalities(log,all_data)

  #jabr inequalities
  if all_data['jabr_inequalities']:
    constrcount += jabr_inequalities(log,all_data)

  #i2 inequalities
  if all_data['i2_inequalities']:
    constrcount += i2_inequalities(log,all_data)

  #limit constraints
  if all_data['limit_inequalities']:
    constrcount += limit_inequalities(log,all_data)
  
  log.joint('  %d constraints added\n'%constrcount)
    
  themodel.update()

  formulation_end = time.time()

  all_data['formulation_time'] = formulation_end - formulation_start

  log.joint(' formulation time: %g\n' % all_data['formulation_time'])

  # Write model to .lp file
  if all_data['writelps']:
    log.joint(' writing to lpfile ' + all_data['lpfilename'] + '\n')  
    themodel.write(all_data['lpfilename'])

  ###################### INIT DATA STRUCTURES FOR CUTS ########################

  if all_data['jabrcuts']:
    jabr_cuts_info = all_data['jabr_cuts_info']

    for branch in branches.values():
      jabr_cuts_info[branch] = {}

  if all_data['i2cuts']:
    i2_cuts_info = all_data['i2_cuts_info']

    for branch in branches.values():
      i2_cuts_info[branch] = {}

  if all_data['limitcuts']:
    limit_cuts_info = all_data['limit_cuts_info']

    for branch in branches.values():
      limit_cuts_info[branch] = {}

  if all_data['linear_objective'] or all_data['hybrid']:
    dicGenPvalues = {}
    for gen in gens.values():
      dicGenPvalues[gen] = []
    all_data['dicGenPvalues'] = dicGenPvalues

  ######################## FIXING/WRITING AN AC SOLUTION ######################

  if all_data['fixflows']:
    fixflows(log,all_data)
    if all_data['fixcs'] == 0:
      return None

  if all_data['fixcs']:
    fixcs(log,all_data)
    return None

  if all_data['writeACsol']:
    writeACsol(log,all_data)
    return None

  ########################### SOLVER PARAMETERS ###############################

  themodel.Params.Method    = all_data['solver_method']
  themodel.Params.Crossover = all_data['crossover'] 


  if all_data['solver_method'] == 2:
    themodel.Params.BarHomogeneous = 1
    themodel.Params.BarConvTol     = 1e-6
    themodel.Params.FeasibilityTol = 1e-6
    themodel.Params.OptimalityTol  = 1e-6

  themodel.Params.NumericFocus = 1
  themodel.Params.OutPutFlag = 1
  
  #themodel.Params.Logfile = 'newlogs/' + all_data['casename'] + '_gurobi.log'
  themodel.Params.Logfile = all_data['casename'] + '_gurobi.log'

  ######################### READING AND LOADING CUTS ##########################

  if all_data['addcuts']:

    t0_cuts = time.time()

    if all_data['max_rounds'] > 1:
      add_cuts_ws(log,all_data)
    else:
      add_cuts(log,all_data)

    themodel.update()

    t1_cuts = time.time()

    all_data['addcuts_time'] = t1_cuts - t0_cuts

    log.joint(' pre-computed cuts added and model updated\n')

    log.joint(' reading and loading cuts time = '
              + str(all_data['addcuts_time']) + '\n')

    if all_data['writelps']:
      themodel.write(all_data['casename']+'_precomputed_cuts.lp')
      log.joint(' model with precomputed written to .lp file\n')
      breakexit('c?')

  ########################## CUTPLANE MAIN LOOP ###############################

  all_data['round']                  = 1
  all_data['runtime']                = time.time() - all_data['T0']
  all_data['cumulative_solver_time'] = 0
  all_data['ftol_counter']           = 0
  oldobj                             = 1
  gap                                = 1e20

  #  all_data['losstest'] 
  all_data['losstest_dic']        = {}
  all_data['losstest_loss']       = {}
  all_data['losstest_branchloss'] = {}
  all_data['losstest_obj']        = {}

  while ((all_data['round'] <= all_data['max_rounds']) and 
         (all_data['runtime'] <= all_data['max_time']) and 
         (all_data['ftol_counter'] <= all_data['ftol_iterates'])):
    

    #new tolerances
    if ((all_data['ftol_counter'] == 4) or (all_data['max_rounds'] == 1) 
        or (all_data['runtime'] > 150)):

      themodel.Params.BarHomogeneous = 1
      themodel.Params.NumericFocus   = 1 #off and then on doesnt help
      themodel.Params.BarConvTol     = 1e-6
      themodel.Params.FeasibilityTol = 1e-6
      themodel.Params.OptimalityTol  = 1e-6


    ########################### HYBRID ALGORITHM ##############################

    if all_data['hybrid']:
      hybrid(log,all_data)
  
    ############################ SOLVING MODEL ################################

    cutplane_optimize(log,all_data)

    ########################### STORING SOLUTION ##############################

    if all_data['writesol'] == 0 and all_data['writesol_wtol'] == 0:

      log.joint(' storing current solution ...\n')
      
      all_data['Pfvalues'] = themodel.getAttr("X",Pvar_f)
      all_data['Qfvalues'] = themodel.getAttr("X",Qvar_f)
      all_data['Ptvalues'] = themodel.getAttr("X",Pvar_t)
      all_data['Qtvalues'] = themodel.getAttr("X",Qvar_t)
      all_data['cvalues']  = themodel.getAttr("X",cvar)
      all_data['svalues']  = themodel.getAttr("X",svar)

      # numpy
      # all_data['cvalues_array']  = np.zeros(numbranches, dtype = 'float')
      # all_data['svalues_array']  = np.zeros(numbranches, dtype = 'float')
      # all_data['cfvalues_array'] = np.zeros(numbranches, dtype = 'float')
      # all_data['ctvalues_array']  = np.zeros(numbranches, dtype = 'float')

      # count = 0 
      # for branch in branches.values():
      #   all_data['cvalues_array'][count] = all_data['cvalues'][branch]
      #   all_data['svalues_array'][count] = all_data['svalues'][branch]
        
      #   busf = buses[IDtoCountmap[branch.f]]
      #   bust = buses[IDtoCountmap[branch.t]]
        
      #   all_data['cfvalues_array'][count] = all_data['cvalues'][busf]
      #   all_data['ctvalues_array'][count] = all_data['cvalues'][bust]        
        
      #   count += 1

      # all_data['Pfvalues_array'] = np.array(all_data['Pfvalues'].values(), dtype = 'float')
      # all_data['Qfvalues_array'] = np.array(all_data['Qfvalues'].values(), dtype = 'float')
      # all_data['Ptvalues_array'] = np.array(all_data['Ptvalues'].values(), dtype = 'float')
      # all_data['Qtvalues_array'] = np.array(all_data['Qtvalues'].values(), dtype = 'float')
      ##########


      #this is only **needed** when using objective cuts...
      #all_data['GenPvalues'] = themodel.getAttr("X",GenPvar)
      #all_data['GenQvalues'] = themodel.getAttr("X",GenQvar)
      #linear objective and objective cuts and hybrid disable if writesol = 0

      if all_data['i2']:
        all_data['i2fvalues'] = themodel.getAttr("X",i2var_f)
      
      log.joint(' done storing values\n')
      
    ####################### WRITING SOL TO A FILE #############################
    
    if all_data['writesol'] or all_data['writesol_wtol']:

      eps = all_data['writesol_wtol']

      namesolfile = 'sol_ws_' + all_data['casename'] + '.txt'
      #namesolfile = 'sols_warmstarted/sol_ws_' + all_data['casename'] + '.txt'
      #namesolfile = 'cutplane_sols/cutplane_sol_' + all_data['casename'] + '.txt'

      if all_data['newsol']:
        namesolfile = 'cutplanesol_'+ all_data['casename'] +'.txt'


      solfile     = open(namesolfile,'a+')
      solfile.write('round' + str(all_data['round']) + '\n' )
      solfile.write('obj ' + str(all_data['objval']) + '\n')

      #values
      log.joint(' storing current solution ...\n')

      cvalues       = {}
      svalues       = {}
      GenPvalues    = {}
      GenQvalues    = {}
      plossvalues   = {}
      qlossvalues   = {}
      qgains        = {}
      qlossvalues_3 = {}
      Pfvalues      = {}
      Qfvalues      = {}
      Ptvalues      = {}
      Qtvalues      = {}
      i2fvalues     = {}
      #i2tvalues     = {}    

      for bus in buses.values():
        cvalues[bus] = cvar[bus].x
        varname      = cvar[bus].varname
        if all_data['writesol']:
          lines        = [varname,' = ',str(cvalues[bus]),'\n']
        elif all_data['writesol_wtol']:
          lines        = [str((cvalues[bus]-eps)),' <= ' + varname 
                          + ' <= ',str((cvalues[bus]+eps)),'\n']
        solfile.writelines(lines)

      ##### LOSS
      LOSS = False

      for branch in branches.values():
        f                   = branch.f
        t                   = branch.t
        count_of_f          = IDtoCountmap[f]
        count_of_t          = IDtoCountmap[t]
        bus_f               = buses[count_of_f]
        bus_t               = buses[count_of_t]
        bc                  = branch.bc
        ratio               = branch.ratio

        cvalues[branch]       = cvar[branch].x
        svalues[branch]       = svar[branch].x
        plossvalues[branch]   = Pvar_f[branch].x + Pvar_t[branch].x
        qlossvalues[branch]   = Qvar_f[branch].x + Qvar_t[branch].x
        qgains[branch]        = - (bc/2) * ( (cvalues[bus_f] / (ratio**2) ) + cvalues[bus_t] )
        Pfvalues[branch]      = Pvar_f[branch].x
        Ptvalues[branch]      = Pvar_t[branch].x
        Qfvalues[branch]      = Qvar_f[branch].x
        Qtvalues[branch]      = Qvar_t[branch].x

        cvarname = cvar[branch].varname 
        svarname = svar[branch].varname 
        
        if all_data['writesol'] and all_data['newsol'] == 0:
          cslines  = [cvarname,' = ',str(cvalues[branch]),'\n',svarname,' = ',str(svalues[branch]),'\n']
          Plines   = [Pvar_f[branch].varname + ' = ',str(Pvar_f[branch].x),'\n',Pvar_t[branch].varname + ' = ',str(Pvar_t[branch].x),'\n']
          Qlines   = [Qvar_f[branch].varname + ' = ',str(Qvar_f[branch].x),'\n',Qvar_t[branch].varname + ' = ',str(Qvar_t[branch].x),'\n']

          solfile.writelines(cslines)
          solfile.writelines(Plines)
          solfile.writelines(Qlines)

          if (Pvar_f[branch].x + Pvar_t[branch].x) < 0:
            LOSS = True


        elif all_data['writesol_wtol'] and all_data['newsol'] == 0:
          cslines  = [str((cvalues[branch]-eps)), ' <= ' + cvarname + ' <= ', str((cvalues[branch]+eps)),'\n',str((svalues[branch]-eps)), ' <= ' + svarname + ' <= ', str((svalues[branch]+eps)),'\n']
          #Plines   = [str((Pvar_f[branch].x - eps)) ,' <= ',Pvar_f[branch].varname,' <= ',str((Pvar_f[branch].x + eps)),'\n',str((Pvar_t[branch].x-eps)),' <= ' + Pvar_t[branch].varname + ' <= ',str((Pvar_t[branch].x+eps)),'\n']
          #Qlines   = [str((Qvar_f[branch].x - eps)) ,' <= ',Qvar_f[branch].varname,' <= ',str((Qvar_f[branch].x + eps)),'\n',str((Qvar_t[branch].x-eps)),' <= ' + Qvar_t[branch].varname + ' <= ',str((Qvar_t[branch].x+eps)),'\n']
    
          solfile.writelines(cslines)

        #solfile.writelines(Plines)
        #solfile.writelines(Qlines)

        

        if all_data['i2']:
          i2fvalues[branch] = i2var_f[branch].x
          i2line            = [i2var_f[branch].varname + ' = ',str(i2fvalues[branch]),'\n']   
          #solfile.writelines(i2line)
          #i2tvalues[branch] = i2var_t[branch].x 

      if LOSS:
        log.joint('free generation\n') 

      if all_data['newsol']:
        all_data['cvalues'] = cvalues
        all_data['svalues'] = svalues
        getangles(log,all_data)
        angles = all_data['angles']
        for busid in angles.keys():
          varname = 'theta_' + str(busid)
          angle   = angles[busid]
          lines = [varname,' = ',str(angle),'\n']
          solfile.writelines(lines)

      for gen in gens.values():
        GenPvalues[gen] = GenPvar[gen].x
        GenQvalues[gen] = GenQvar[gen].x
        GenPvarname     = GenPvar[gen].varname + ' = '
        GenQvarname     = GenQvar[gen].varname + ' = '
        GenPlines       = [GenPvarname,str(GenPvalues[gen]),'\n']
        GenQlines       = [GenQvarname,str(GenQvalues[gen]),'\n']
      
        solfile.writelines(GenPlines)
        solfile.writelines(GenQlines)

        if all_data['linear_objective'] or all_data['hybrid']:
          if gen.costdegree == 2 and gen.costvector[0] != 0:
            gen_values = dicGenPvalues[gen]
            bisect.insort(gen_values,GenPvar[gen].x)

      solfile.close()
        
      all_data['cvalues']                 = cvalues
      all_data['svalues']                 = svalues
      all_data['GenPvalues']              = GenPvalues
      all_data['GenQvalues']              = GenQvalues
      all_data['plossvalues']             = plossvalues
      all_data['qlossvalues']             = qlossvalues
      all_data['total_active_gen']        = sum(GenPvalues.values())
      all_data['total_active_losses']     = sum(plossvalues.values())
      all_data['total_reactive_gen']      = sum(GenQvalues.values())
      all_data['total_reactive_losses']   = sum(qlossvalues.values())
      all_data['total_reactive_gains']    = sum(qgains.values())
      all_data['Pfvalues']                = Pfvalues
      all_data['Ptvalues']                = Ptvalues
      all_data['Qfvalues']                = Qfvalues
      all_data['Qtvalues']                = Qtvalues

      if all_data['i2']:
        all_data['i2fvalues'] = i2fvalues
        #all_data['i2tvalues'] = i2tvalues

      if all_data['linear_objective'] or all_data['hybrid']:
        all_data['dicGenPvalues'] = dicGenPvalues



      log.joint(' done storing and writing down values\n')


    ########################### ROUND STATISTICS ##############################

    cutplane_stats(log,all_data)

    ########################### Flow Decomposition ############################

    # log.joint(' ploss ' + str(sum(plossvalues.values())) + '\n')

    # if all_data['writesol']:
    #   all_data['sol_Pfvalues']   = all_data['Pfvalues']
    #   all_data['sol_Ptvalues']   = all_data['Ptvalues']
    #   all_data['sol_GenPvalues'] = {}

    #   for gen in gens.values():
    #     all_data['sol_GenPvalues'][gen.count] = all_data['GenPvalues'][gen]
    
    #   log.joint(' sol_GenPvalues ' + str(all_data['sol_GenPvalues']) + '\n') 

    #   flowdecomp(log,all_data)

    ######################### SUMMARY EXPERIMENTS #############################

    all_data['runtime'] = time.time() - all_data['T0']

    log.joint("\n writing casename, opt stauts, obj and runtime to summary_ws.log\n")

    summary_ws = open("summary_ws.log","a+") #later add feasibility errors, etc                            

    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(all_data['optstatus']) + ' obj ' 
                     + str(all_data['objval']) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + '\n')

    summary_ws.close()

    #if all_data['max_rounds'] == 1:
    #  log.joint(' bye!\n')
    #  sys.exit(0)

    ############################### CUTS ######################################

    # Cut computations and management
    cutplane_cuts(log,all_data)

    # Cut statistics
    cutplane_cutstats(log,all_data)
    
    themodel.update()

    log.joint(' model updated\n')
    log.joint('\n')

    ############################### WRITE CUTS ################################

    if all_data['writecuts']:
      write_cuts(log,all_data)

    ############################### WRITE LPS #################################
    
    if all_data['losstest']:
      themodel.write('perturbedjabr.lp')
      constr = themodel.getQConstrs()[all_data['rbranch']]
      constr.setAttr("QCRHS",0)
  
      themodel.update()
 
    ###########

    if all_data['writelps'] and ( all_data['round'] > 0 ):
      name = 'post_cuts' + '_' + str(all_data['round']) + '.lp'
      themodel.write(name)
      log.joint(' model with new cuts written to .lp file\n')

    # valloss = 0
    # for branch in branches.values():
    #   loss = all_data['Pfvalues'][branch] + all_data['Ptvalues'][branch]
    #   log.joint(' loss at branch ' + str(branch.count) + ' = ' 
    #             + str(loss) + '\n')
    #   valloss += loss

    # log.joint('ploss ' + str(sum(plossvalues.values())) + '\n')
    # log.joint('value loss ' + str(valloss) + '\n')

    ################# losses test

    if all_data['losstest']:
      rbranch = all_data['branches'][all_data['rbranch']]

      rbranchloss = all_data['Pfvalues'][rbranch] + all_data['Ptvalues'][rbranch]
      all_data['losstest_dic'][all_data['round']] = (all_data['rbranch'],
                                                     all_data['objval'],
                                                     sum(plossvalues.values()),
                                                     rbranchloss)

      all_data['losstest_loss'][all_data['round']] = sum(plossvalues.values())
      all_data['losstest_branchloss'][all_data['round']] = rbranchloss
      all_data['losstest_obj'][all_data['round']] = all_data['objval']

      avg        = sum(all_data['losstest_loss'].values()) / all_data['round']
      avg_branch = sum(all_data['losstest_branchloss'].values()) / all_data['round']

      avg_obj = sum(all_data['losstest_obj'].values()) / all_data['round'] 
      maxloss = max(all_data['losstest_loss'].values())
      maxloss_branch = max(all_data['losstest_branchloss'].values())

      minloss = min(all_data['losstest_loss'].values())
      minloss_branch = min(all_data['losstest_branchloss'].values())

      log.joint(' round ' + str(all_data['round']) + '\n')
      log.joint(' random branch ' + str(all_data['rbranch']) + '\n')
      log.joint(' current obj avg ' + str(avg_obj) + '\n' )
      log.joint(' current losstest avg ' + str(avg) + '\n' )
      log.joint(' current losstest avg (branch)' + str(avg_branch) + '\n' )
      log.joint(' current max loss ' + str(maxloss) + '\n')
      log.joint(' current max loss (at the branch)' + str(maxloss_branch) + '\n')
      log.joint(' current min loss ' + str(minloss) + '\n')
      log.joint(' current min loss (at the branch)' + str(minloss_branch) + '\n')
      #log.joint(' losstest ' + str(all_data['losstest_loss']) + '\n')
      #log.joint(' losstest (branch) ' + str(all_data['losstest_branchloss']) + '\n')
      #breakexit('c')
      if all_data['round'] == 50:
        return None


    ########################## CHECK OBJ IMPROVEMENT ##########################

    if ((all_data['objval'] - oldobj)/oldobj) < all_data['ftol']:
      all_data['ftol_counter'] += 1
    else:
      all_data['ftol_counter'] = 0

    oldobj              = all_data['objval']
    all_data['runtime'] = time.time() - all_data['T0']


    all_data['round'] += 1
    
    ###########################################################################



def fixflows(log,all_data):

  if (all_data['matpower_sol'] == 0) and (all_data['knitro_sol'] == 0):
    log.joint(' cannot fix flows since no solution has been loaded!\n')
    return None

  themodel    = all_data['themodel']
  branches    = all_data['branches']
  tolerance   = all_data['tol_fix']
  Pvar_f      = all_data['Pvar_f']
  Pvar_t      = all_data['Pvar_t']
  Qvar_f      = all_data['Qvar_f']
  Qvar_t      = all_data['Qvar_t']

  sol_Pfvalues = all_data['sol_Pfvalues']
  sol_Ptvalues = all_data['sol_Ptvalues']
  sol_Qfvalues = all_data['sol_Qfvalues']
  sol_Qtvalues = all_data['sol_Qtvalues']

  for branch in branches.values():
        
    sol_Pf = sol_Pfvalues[branch]
    sol_Pt = sol_Ptvalues[branch]
    sol_Qf = sol_Qfvalues[branch]
    sol_Qt = sol_Qtvalues[branch]

    #Pf
    ubound_Pf = sol_Pf + tolerance
    lbound_Pf = sol_Pf - tolerance

    Pvar_f[branch].setAttr("ub",ubound_Pf)
    Pvar_f[branch].setAttr("lb",lbound_Pf)

    #Pt
    ubound_Pt = sol_Pt + tolerance
    lbound_Pt = sol_Pt - tolerance

    Pvar_t[branch].setAttr("ub",ubound_Pt)
    Pvar_t[branch].setAttr("lb",lbound_Pt)

    #Qf
    ubound_Qf = sol_Qf + tolerance
    lbound_Qf = sol_Qf - tolerance

    Qvar_f[branch].setAttr("ub",ubound_Qf)
    Qvar_f[branch].setAttr("lb",lbound_Qf)

    #Qt
    ubound_Qt = sol_Qt + tolerance
    lbound_Qt = sol_Qt - tolerance

    Qvar_t[branch].setAttr("ub",ubound_Qt)
    Qvar_t[branch].setAttr("lb",lbound_Qt)

  themodel.update()
  themodel.write('fixflows.lp')
  log.joint('check fixflows.lp\n')  


def fixcs(log,all_data):

  if (all_data['matpower_sol'] == 0) and (all_data['knitro_sol'] == 0):
    log.joint(' cannot fix flows since no solution has been loaded!\n')
    return None

  themodel   = all_data['themodel']
  tolerance  = all_data['tol_fix']
  buses      = all_data['buses']
  branches   = all_data['branches']
  cvar       = all_data['cvar']
  svar       = all_data['svar']
  sol_cvalues = all_data['sol_cvalues']
  sol_svalues = all_data['sol_svalues']

  for bus in buses.values():

    sol_v2 = all_data['sol_cvalues'][bus]

    ubound_v2 = sol_v2 + tolerance
    lbound_v2 = sol_v2 - tolerance

    cvar[bus].setAttr("ub",ubound_v2)
    cvar[bus].setAttr("lb",lbound_v2)

  for branch in branches.values():

    sol_c = all_data['sol_cvalues'][branch]
    sol_s = all_data['sol_svalues'][branch]

    #c
    ubound_c = sol_c + tolerance
    lbound_c = sol_c - tolerance

    cvar[branch].setAttr("ub",ubound_c)
    cvar[branch].setAttr("lb",lbound_c)

    #s
    ubound_s = sol_s + tolerance
    lbound_s = sol_s - tolerance

    svar[branch].setAttr("ub",ubound_s)
    svar[branch].setAttr("lb",lbound_s)

  themodel.update()
  themodel.write('fixCS.lp')
  log.joint('check fixCS.lp\n')


def computebalbounds(log, all_data, bus):

  #first let's get max/min generations

  loud = 0

  baseMVA = all_data['baseMVA']
  gens = all_data['gens']

  Pubound = Plbound = 0
  Qubound = Qlbound = 0


  for gencounter in bus.genidsbycount:
    if gens[gencounter].status:
      Pubound += gens[gencounter].Pmax
      Plbound += gens[gencounter].Pmin
      Qubound += gens[gencounter].Qmax
      Qlbound += gens[gencounter].Qmin

    if loud:
     #log.joint(" Pubound for " + str(bus.nodeID) + " " + str(Pubound) + " genc " + str(gencounter) + "\n")
     #log.joint(" Plbound for " + str(bus.nodeID) + " " + str(Plbound) + " genc " + str(gencounter) + "\n")
     log.joint(" Qubound for " + str(bus.nodeID) + " " + str(Qubound) + " genc " + str(gencounter) + "\n")
     log.joint(" Qlbound for " + str(bus.nodeID) + " " + str(Qlbound) + " genc " + str(gencounter) + "\n")

  Pubound -= bus.Pd
  Plbound -= bus.Pd
  Qubound -= bus.Qd
  Qlbound -= bus.Qd

  if bus.nodetype == 4:
    Pubound = Plbound = Qubound = Qlbound = 0

  if loud:
     #log.joint(" Pubound for " + str(bus.nodeID) + " final " + str(Pubound) + "\n")
     #log.joint(" (Pd was %g)\n" %bus.Pd)
     #log.joint(" Plbound for " + str(bus.nodeID) + " final " + str(Plbound) + "\n")
     log.joint(" Qubound for " + str(bus.nodeID) + " final " + str(Qubound) + "\n")
     log.joint(" (Qd was %g)\n" %bus.Qd)
     log.joint(" Qlbound for " + str(bus.nodeID) + " final " + str(Qlbound) + "\n")
     breakexit(" ")
  
  return Pubound, Plbound, Qubound, Qlbound


def computeangles(log,all_data):
  
  cvalues  = all_data['cvalues']
  svalues  = all_data['svalues']
  branches = all_data['branches']

  for branch in branches.values():
    s = svalues[branch]
    c = cvalues[branch]

    log.joint(' branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' angle ' + str(math.atan(s/c)) + '\n')

def getsol_knitro(log,all_data):

  casename    = all_data['casename']
  filename    = 'knitro_sols/ksol_'+ casename +'.txt'

  try:
    thefile   = open(filename, "r")
    lines     = thefile.readlines()
    lenlines  = len(lines)
    thefile.close()
  except:
    log.stateandquit("cannot open file " + datafilename)
    sys.exit("failure")

  branches     = all_data['branches']
  buses        = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap'] 

  sol_obj      = 0
  sol_vm       = {}
  sol_angle    = {}
  sol_cvalues  = {}
  sol_svalues  = {}
  sol_Pfvalues = {}
  sol_Ptvalues = {}
  sol_Qfvalues = {}
  sol_Qtvalues = {}

  linenum = 0
  log.joint(' reading file\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    if thisline[0] == 'value':
      sol_obj              = float(lines[0].split()[1])
    elif thisline[0] == 'bus':
      buscount             = int(thisline[1])
      bus                  = buses[buscount]
      sol_vm[bus]          = float(thisline[3])
      sol_angle[bus]       = float(thisline[5]) 
      sol_cvalues[bus]     = sol_vm[bus]**2
    elif thisline[0] == 'branch':
      branchcount          = int(thisline[1])
      branch               = branches[branchcount]
      sol_Pfvalues[branch] = float(thisline[7])
      sol_Ptvalues[branch] = float(thisline[9])
      sol_Qfvalues[branch] = float(thisline[11])
      sol_Qtvalues[branch] = float(thisline[13])

    linenum += 1

  for branch in branches.values():
    f          = branch.f
    t          = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    bus_f      = buses[count_of_f]
    bus_t      = buses[count_of_t]
    vm_f       = sol_vm[bus_f]
    vm_t       = sol_vm[bus_t]

    sol_cvalues[branch] = vm_f * vm_t * math.cos(sol_angle[bus_f] - sol_angle[bus_t])
    sol_svalues[branch] = vm_f * vm_t * math.sin(sol_angle[bus_f] - sol_angle[bus_t])

  all_data['sol_vm']       = sol_vm
  all_data['sol_angle']    = sol_angle
  all_data['sol_cvalues']  = sol_cvalues
  all_data['sol_svalues']  = sol_svalues
  all_data['sol_Pfvalues'] = sol_Pfvalues
  all_data['sol_Ptvalues'] = sol_Ptvalues
  all_data['sol_Qfvalues'] = sol_Qfvalues
  all_data['sol_Qtvalues'] = sol_Qtvalues

  log.joint(' knitro solution loaded\n')
  

def getsol_matpower(log,all_data):

  casename = all_data['casename']
  filename = 'mp_sols/solution_va_'+ casename +'.txt'

  log.joint(" reading file with matpower solution volt magnitudes and angles " + filename + "\n")

  try:
      thefile  = open(filename, "r")
      lines    = thefile.readlines()
      lenlines = len(lines)
      thefile.close()
  except:
      log.stateandquit("cannot open file", filename)
      sys.exit("failure")

  sol_vm       = {}
  sol_angle    = {}
  sol_cvalues  = {}
  sol_svalues  = {}

  buses        = all_data['buses']
  branches     = all_data['branches']
  IDtoCountmap = all_data['IDtoCountmap']


  linenum  = 1
  while linenum < lenlines:
    thisline = lines[linenum].split(',')
    buscount = int(thisline[0]) #bug, not bus_id -> matpower uses buscount
    bus      = buses[buscount]

    sol_vm[bus]        = float(thisline[1])
    sol_angle[bus]     = float(thisline[2]) * math.pi / 180
    sol_cvalues[bus]   = sol_vm[bus]**2

    linenum += 1

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]
    bus_f       = buses[count_of_f]
    bus_t       = buses[count_of_t]


    vm_f        = sol_vm[bus_f]
    vm_t        = sol_vm[bus_t]

    sol_c      = vm_f * vm_t * math.cos(sol_angle[bus_f] - sol_angle[bus_t])
    sol_s      = vm_f * vm_t * math.sin(sol_angle[bus_f] - sol_angle[bus_t])

    sol_cvalues[branch] = sol_c
    sol_svalues[branch] = sol_s

  all_data['sol_vm']      = sol_vm
  all_data['sol_angle']   = sol_angle
  all_data['sol_cvalues'] = sol_cvalues
  all_data['sol_svalues'] = sol_svalues

  log.joint(' done loading volts, angles, and cs values\n')

  #########

  filename = 'mp_sols/solution_'+ casename +'.txt'
  log.joint(" reading file with matpower power flows " + filename + "\n")

  try:
      thefile  = open(filename, "r")
      lines    = thefile.readlines()
      lenlines = len(lines)
      thefile.close()
  except:
      log.stateandquit("cannot open file", filename)
      sys.exit("failure")

  sol_Pfvalues    = {}
  sol_Ptvalues    = {}
  sol_Qfvalues    = {}
  sol_Qtvalues    = {}

  linenum = 1

  while linenum < lenlines:
    thisline            = lines[linenum].split(',')
    branchcount         = int(thisline[0])
    branch              = branches[branchcount]
    f                   = branch.f
    t                   = branch.t

    if int(thisline[1]) != f or int(thisline[2]) != t:
      breakexit('check')

    sol_Pfvalues[branch]      = float(thisline[3])/all_data['baseMVA']
    sol_Ptvalues[branch]      = float(thisline[4])/all_data['baseMVA']
    sol_Qfvalues[branch]      = float(thisline[5])/all_data['baseMVA']
    sol_Qtvalues[branch]      = float(thisline[6])/all_data['baseMVA']

    linenum   += 1
  
  all_data['sol_Pfvalues'] = sol_Pfvalues
  all_data['sol_Ptvalues'] = sol_Ptvalues
  all_data['sol_Qfvalues'] = sol_Qfvalues
  all_data['sol_Qtvalues'] = sol_Qtvalues

  log.joint(' done loading power flows\n')
  log.joint(' matpower solution loaded\n')

def writeACsol(log,all_data):

  if (all_data['matpower_sol'] == 0) and (all_data['knitro_sol'] == 0):
    log.joint(' cannot write AC solution since no solution has been loaded!\n')
    return None

  branches     = all_data['branches']
  buses        = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap']
  sol_vm       = all_data['sol_vm']
  sol_angle    = all_data['sol_angle']
  sol_cvalues  = all_data['sol_cvalues']
  sol_svalues  = all_data['sol_svalues']
  sol_Pfvalues = all_data['sol_Pfvalues']
  sol_Ptvalues = all_data['sol_Ptvalues']
  sol_Qfvalues = all_data['sol_Qfvalues']
  sol_Qtvalues = all_data['sol_Qtvalues']
  tolerance    = all_data['tol_fix']

  filename   = 'ACsol_' + all_data['casename'] + '.lp'
  thefile    = open(filename,"w+")

  log.joint(' writing voltages to ' + filename + ' ...\n')

  for bus in buses.values():
    f        = bus.nodeID
    sol_v2   = sol_cvalues[bus]
    linelp_v = str(sol_v2 - tolerance) + ' <= c_' + str(f) + '_' + str(f) + ' <= ' + str(sol_v2 + tolerance) + '\n'
    thefile.write(linelp_v)

  log.joint(' writing cs values to ' + filename + ' ...\n')

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]
    bus_f       = buses[count_of_f]
    bus_t       = buses[count_of_t]

    vm_f        = sol_vm[bus_f]
    vm_t        = sol_vm[bus_t]

    sol_c      = vm_f * vm_t * math.cos(sol_angle[bus_f] - sol_angle[bus_t])
    sol_s      = vm_f * vm_t * math.sin(sol_angle[bus_f] - sol_angle[bus_t])

    linelp_c = str(sol_c - tolerance) + ' <= c_' + str(f) + '_' + str(t) + ' <= ' + str(sol_c + tolerance) + '\n'
    thefile.write(linelp_c)

    linelp_s = str(sol_s - tolerance) + ' <= s_' + str(f) + '_' + str(t) + ' <= ' + str(sol_s + tolerance) + '\n'
    thefile.write(linelp_s)

  log.joint(' writing flows to ' + filename + ' ...\n')

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    Pfval       = sol_Pfvalues[branch]
    Ptval       = sol_Ptvalues[branch]
    Qfval       = sol_Qfvalues[branch]
    Qtval       = sol_Qtvalues[branch]
  
    linelp_Pf = str(Pfval-tolerance) + ' <= P_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' <= ' + str(Pfval+tolerance) + '\n'
    thefile.write(linelp_Pf)
    
    linelp_Pt = str(Ptval-tolerance) + ' <= P_' + str(branchcount) + '_' + str(t) + '_' + str(f) + ' <= ' + str(Ptval+tolerance) + '\n'
    thefile.write(linelp_Pt)
    
    linelp_Qf = str(Qfval-tolerance) + ' <= Q_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' <= ' + str(Qfval+tolerance) + '\n'
    thefile.write(linelp_Qf)

    linelp_Qt = str(Qtval-tolerance) + ' <= Q_' + str(branchcount) + '_' + str(t) + '_' + str(f) + ' <= ' + str(Qtval+tolerance) + '\n'
    thefile.write(linelp_Qt)

  thefile.close()
  log.joint(' done writing AC solution to .lp file\n')
  

def computei2value(log,all_data,branch,mp_c,mp_s,mp_cbusf,mp_cbust):

  ratio  = branch.ratio
  y      = branch.y
  g      = y.real
  b      = y.imag
  bshunt = branch.bc
  angle  = branch.angle_rad

  #first f                                                                                                                                                                                           
  i2f  = 0
  i2f += (g*g + b*b)/(ratio*ratio) * ( (mp_cbusf/(ratio*ratio)) + mp_cbust - (2/ratio) * ( mp_c * math.cos(angle) + mp_s * math.sin(angle) ) )
  i2f += b*bshunt/(ratio**3) * ( (mp_cbusf/ratio) - (mp_c * math.cos(angle) + mp_s * math.sin(angle) ))
  i2f += g*bshunt/(ratio**3) * ( mp_s * math.cos(angle) - mp_c * math.sin(angle) )
  i2f += (bshunt*bshunt*mp_cbusf/(4*(ratio**4)) )

  return i2f
  
def loss_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  Pvar_f         = all_data['Pvar_f']
  Pvar_t         = all_data['Pvar_t']
  FeasibilityTol = all_data['FeasibilityTol']

  if ( all_data['matpower_sol'] or all_data['knitro_sol'] ) and all_data['loss_validity']:
    FeasibilityTol = all_data['FeasibilityTol']
    log.joint('  adding and checking validity of loss inequalities wrt a solution\n')
    log.joint('   feasibility tolerance ' + str(FeasibilityTol) + '\n')
  else:
    log.joint('  loss inequalities\n')

  all_data['loss_cuts'] = {}
  counter_loss = 0

  for branch in branches.values():
    if branch.r < 0: ##!!!
      continue

    branchcount = branch.count
    f           = branch.f
    t           = branch.t

    if all_data['loss_validity']:
     sol_Pf     = all_data['sol_Pfvalues'][branch]
     sol_Pt     = all_data['sol_Ptvalues'][branch]
     violation  = - (sol_Pf + sol_Pt)

     if violation > FeasibilityTol:
       log.joint('   WARNING, the loss inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
       log.joint('   violation ' + str(violation) + '\n')
       log.joint('   values (AC solution) ' + ' Pf ' + str(sol_Pf) + ' Pt ' + str(sol_Pt) + '\n')
       breakexit('check!')
     else:
       log.joint('   AC solution satisfies loss inequality at branch ' + str(branchcount) + ' with slack ' + str(violation) + '\n')

    counter_loss                 += 1
    all_data['loss_cuts'][branch] = (0,0,FeasibilityTol)

    lossexp    = LinExpr()
    constrname = "loss_ineq_"+str(branch.count)+"_"+str(f)+"_"+str(t)
    lossexp   += Pvar_f[branch] + Pvar_t[branch]
    themodel.addConstr(lossexp >= 0, name = constrname)

  all_data['num_loss_cuts']         = counter_loss
  all_data['num_loss_cuts_dropped'] = 0
  all_data['dropped_loss']          = []

  log.joint('   %d loss inequalities added\n'%counter_loss) 
  
  return counter_loss
  
def qloss_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  Pvar_f         = all_data['Pvar_f']
  Pvar_t         = all_data['Pvar_t']
  Qvar_f         = all_data['Qvar_f']
  Qvar_t         = all_data['Qvar_t']
  FeasibilityTol = all_data['FeasibilityTol']

  if ( all_data['matpower_sol'] or all_data['knitro_sol'] ) and all_data['qloss_validity']:
    FeasibilityTol = all_data['FeasibilityTol']
    log.joint('  adding and checking validity of loss inequalities wrt a solution\n')
    log.joint('   feasibility tolerance ' + str(FeasibilityTol) + '\n')
  else:
    log.joint('  qloss inequalities\n')

  all_data['qloss_cuts'] = {}
  qcounter_loss = 0

  for branch in branches.values():
    if branch.r <= 0 or (branch.bc != 0): ##!!!
      continue

    branchcount = branch.count
    f           = branch.f
    t           = branch.t

    if all_data['qloss_validity']:
     sol_Pf     = all_data['sol_Qfvalues'][branch]
     sol_Pt     = all_data['sol_Qtvalues'][branch]
     violation  = - (sol_Qf + sol_Qt)

     if violation > FeasibilityTol:
       log.joint('   WARNING, the qloss inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
       log.joint('   violation ' + str(violation) + '\n')
       log.joint('   values (AC solution) ' + ' Qf ' + str(sol_Qf) + ' Pt ' + str(sol_Qt) + '\n')
       breakexit('check!')
     else:
       log.joint('   AC solution satisfies loss inequality at branch ' + str(branchcount) + ' with slack ' + str(violation) + '\n')

    qcounter_loss                 += 1
    all_data['qloss_cuts'][branch] = (0,0,FeasibilityTol)

    qlossexp    = LinExpr()
    constrname = "qloss_ineq_"+str(branch.count)+"_"+str(f)+"_"+str(t)
    qlossexp   += (Qvar_f[branch] + Qvar_t[branch]) - (branch.x/branch.r)*(Pvar_f[branch] + Pvar_t[branch])
    themodel.addConstr(qlossexp == 0, name = constrname)

  all_data['num_qloss_cuts']         = qcounter_loss
  all_data['num_qloss_cuts_dropped'] = 0
  all_data['dropped_qloss']          = []

  log.joint('   %d qloss inequalities added\n'%qcounter_loss) 
  
  return qcounter_loss

def jabr_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  buses          = all_data['buses']
  cvar           = all_data['cvar']
  svar           = all_data['svar']
  IDtoCountmap   = all_data['IDtoCountmap']
  FeasibilityTol = all_data['FeasibilityTol']

  if ( all_data['matpower_sol'] or all_data['knitro_sol'] ) and all_data['jabr_validity']:
    log.joint('  adding and checking validity of Jabr inequalities wrt a solution\n')
    log.joint('   feasibility tolerance ' + str(FeasibilityTol) + '\n')
  else:
    log.joint('  Jabr inequalities\n')

  maxviolation = 0
  violated     = 0
  maxbranch    = -1
  maxbusf      = -1
  maxbust      = -1
  counter_jabr = 0

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]

    if ( all_data['matpower_sol'] or all_data['knitro_sol'] ) and all_data['jabr_validity']:
      sol_c         = all_data['sol_cvalues'][branch]
      sol_s         = all_data['sol_svalues'][branch]
      sol_cbusf     = all_data['sol_cvalues'][buses[count_of_f]]
      sol_cbust     = all_data['sol_cvalues'][buses[count_of_t]]
      relviolation = violation = sol_c * sol_c + sol_s * sol_s - sol_cbusf * sol_cbust
      #relviolation = violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 )
      if relviolation > maxviolation:
        maxviolation = relviolation
        maxbranch    = branch.count
        maxbusf      = f
        maxbust      = t

      if relviolation > FeasibilityTol:
        violated += 1
        log.joint('   WARNING, the Jabr-inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
        log.joint('   violation ' + str(violation) + '\n')
        log.joint('   relative violation ' + str(relviolation) + '\n')
        log.joint('   values (AC solution) ' + ' cft ' + str(sol_c) + ' sft ' + str(sol_s) + ' cff ' + str(sol_cbusf) + ' ctt ' + str(sol_cbust) + '\n' )
        #breakexit('check!')
      else:
        log.joint('   AC solution satisfies loss inequality at branch ' + str(branchcount) + ' with slack ' + str(relviolation) + '\n')

    counter_jabr += 1
    trigexp       = QuadExpr()
    constrname    = "jabr_"+str(branchcount)+"_"+str(f)+"_"+str(t)
    trigexp      += cvar[branch]*cvar[branch] + svar[branch]*svar[branch] - cvar[buses[count_of_f]]*cvar[buses[count_of_t]]

    themodel.addConstr(trigexp <= 0, name = constrname)

  if ( all_data['matpower_sol'] or all_data['knitro_sol'] ) and all_data['jabr_validity']:
    log.joint('  max violation of Jabr-inequalities by AC solution ' + str(maxviolation) + ' at branch ' + str(maxbranch) + ' f ' + str(maxbusf) + ' t ' + str(maxbust) + '\n')
    log.joint('  number of violated Jabr-inequalities ' + str(violated) + '\n')
    breakexit('  check Jabr violation')

  log.joint('   %d Jabr inequalities added\n'%counter_jabr) 

  return counter_jabr



def i2_def(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  buses          = all_data['buses']
  i2var_f        = all_data['i2var_f']
  cvar           = all_data['cvar']
  svar           = all_data['svar']
  IDtoCountmap   = all_data['IDtoCountmap']
  
  counter_i2def  = 0

  log.joint('  i2 variables definition\n')

  for branch in branches.values():
    expr_f = LinExpr()
    #expr_t = LinExpr()

    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]
    ratio       = branch.ratio
    y           = branch.y
    g           = y.real
    b           = y.imag
    bshunt      = branch.bc
    angle       = branch.angle_rad

    # if f == 549 and t == 5002:
    #   log.joint(' g ' + str(g) + ' b ' + str(b) + ' bshunt ' + str(bshunt) + ' angle ' + str(angle) + ' ratio ' + str(ratio) + '\n')
    #   breakexit('c')

    #first f
    expr_f += (g*g + b*b)/(ratio*ratio) * ( (cvar[buses[count_of_f]]/(ratio*ratio)) + cvar[buses[count_of_t]] - (2/ratio) * ( cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) ) )
    expr_f += b*bshunt/(ratio**3) * ( (cvar[buses[count_of_f]]/ratio) - (cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) )) #since gsh is 0
    expr_f += g*bshunt/(ratio**3) * ( svar[branch] * math.cos(angle) - cvar[branch] * math.sin(angle) )
    expr_f += (bshunt*bshunt*cvar[buses[count_of_f]]/(4*(ratio**4)) )

    # #now t
    # expr_t += (g*g + b*b) * ( cvar[buses[count_of_t]] + cvar[buses[count_of_f]]/(ratio*ratio) - (2/ratio) * ( cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) ) )
    # expr_t += b*bshunt * ( cvar[buses[count_of_t]] - (1/ratio) * (cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) ))
    # expr_t += (g*bshunt/ratio) * ( svar[branch] * math.cos(angle) - cvar[branch] * math.sin(angle) )
    # expr_t += ( bshunt*bshunt*cvar[buses[count_of_t]]/4 )
    
    counter_i2def  += 1

    constrname_f    = 'i2def_'+str(branch.count)+"_"+str(f) + "_" + str(t)
    themodel.addConstr(expr_f == i2var_f[branch],name = constrname_f) 

    #constrname_t = 'i2def_'+str(branch.count)+"_"+str(t) + "_" + str(f)
    #themodel.addConstr(expr_t == i2var_t[branch],name = constrname_t) 

  log.joint('   %d i2 definition constraints added\n'%counter_i2def) 

  return counter_i2def


def i2_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  buses          = all_data['buses']
  cvar           = all_data['cvar']
  i2var_f        = all_data['i2var_f']
  Pvar_f         = all_data['Pvar_f']
  Qvar_f         = all_data['Qvar_f']
  IDtoCountmap   = all_data['IDtoCountmap']
  FeasibilityTol = all_data['FeasibilityTol']

  if ( all_data['matpower_sol'] or all_data['knitro_sol'] ) and all_data['i2_validity']:
    log.joint('  adding and checking validity of i2 inequalities wrt a solution\n')
    log.joint('   feasibility tolerance ' + str(FeasibilityTol) + '\n')
  else:
    log.joint('  i2 inequalities\n')

  maxviolation = 0
  violated     = 0
  maxbranch    = -1
  maxbusf      = -1
  maxbust      = -1
  maxi2f       = -1
  maxcff       = -1
  maxPf        = -1
  maxQf        = -1

  counter_i2   = 0

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]

    if ( all_data['matpower_sol'] or all_data['knitro_sol'] ) and all_data['i2_validity']:
      sol_Pf        = all_data['sol_Pfvalues'][branch]
      sol_Qf        = all_data['sol_Qfvalues'][branch]
      sol_c         = all_data['sol_cvalues'][branch]
      sol_s         = all_data['sol_svalues'][branch]
      sol_cbusf     = all_data['sol_cvalues'][buses[count_of_f]]
      sol_cbust     = all_data['sol_cvalues'][buses[count_of_t]]
      sol_i2f       = computei2value(log,all_data,branch,sol_c,sol_s,sol_cbusf,sol_cbust)
      relviolation = violation = sol_Pf * sol_Pf + sol_Qf * sol_Qf - sol_cbusf * sol_i2f
      #relviolation = violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 )
      if relviolation > maxviolation:
        maxviolation = relviolation
        maxbranch    = branch.count
        maxbusf      = f
        maxbust      = t
        maxi2f       = sol_i2f
        maxcff       = sol_cbusf
        maxPf        = sol_Pf
        maxQf        = sol_Qf

      if relviolation > FeasibilityTol:
        violated += 1
        log.joint('   WARNING, the i2 inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
        log.joint('   violation ' + str(violation) + '\n')
        log.joint('   relative violation ' + str(relviolation) + '\n')
        log.joint('   values (AC solution) ' + ' Pft ' + str(sol_Pf) + ' Qft ' + str(sol_Qf) + ' cff ' + str(sol_cbusf) + ' i2ft ' + str(sol_i2f) + '\n' )
        #breakexit('check!')
      else:
        log.joint('   AC solution satisfies i2 inequality at branch ' + str(branchcount) + ' with slack ' + str(relviolation) + '\n')

    counter_i2 += 1
    trigexp     = QuadExpr()
    constrname  = "i2_"+str(branchcount)+"_"+str(f)+"_"+str(t)
    trigexp    += Pvar_f[branch]**2 + Qvar_f[branch]**2 - cvar[buses[count_of_f]] * i2var_f[branch]
    themodel.addConstr(trigexp <= 0, name = constrname)

  if ( all_data['matpower_sol'] or all_data['knitro_sol'] ) and all_data['i2_validity']:
    log.joint('  max violation of i2 inequalities by AC solution ' + str(maxviolation) + ' at branch ' + str(maxbranch) + ' f ' + str(maxbusf) + ' t ' + str(maxbust) + '\n')
    log.joint('  values (AC solution) ' + ' Pft ' + str(maxPf) + ' Qft ' + str(maxQf) + ' cff ' + str(maxcff) + ' i2ft ' + str(maxi2f) + '\n' )
    log.joint('  number of violated i2 inequalities ' + str(violated) + '\n')
    breakexit('  check i2 violation')

  log.joint('   %d i2 inequalities added\n'%counter_i2) 

  return counter_i2  


def limit_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  Pvar_f         = all_data['Pvar_f']
  Pvar_t         = all_data['Pvar_t']
  Qvar_f         = all_data['Qvar_f']
  Qvar_t         = all_data['Qvar_t']

  counter_limit  = 0

  log.joint('  limit inequalities\n')

  for branch in branches.values():
    if branch.constrainedflow:
      branchcount = branch.count
      f           = branch.f
      t           = branch.t

      constrname = "limit_f_"+str(branchcount)+"_"+str(f)+"_"+str(t)
      limexp = QuadExpr()
      limexp += Pvar_f[branch]*Pvar_f[branch] + Qvar_f[branch]*Qvar_f[branch]
      themodel.addConstr(limexp <= branch.limit**2, name = constrname)

      constrname = "limit_t_"+str(branchcount)+"_"+str(t)+"_"+str(f)
      limexp = QuadExpr()
      limexp += Pvar_t[branch]*Pvar_t[branch] + Qvar_t[branch]*Qvar_t[branch]
      themodel.addConstr(limexp <= branch.limit**2, name = constrname)

      counter_limit += 2

  log.joint('   %d limit inequalities added\n'%counter_limit) 

  return counter_limit


def hybrid(log,all_data):

  themodel   = all_data['themodel']
  gens       = all_data['gens']
  GenPvar    = all_data['GenPvar']
  objvar     = all_data['objvar']
  lincostvar = all_data['lincostvar']
  qcostvar   = all_data['qcostvar']
  qcost      = all_data['qcost']
  sumTconstr = all_data['sumTconstr']
  lincost    = all_data['lincost']

  QP                  = 50
  QP_to_LP            = 3
  no_objective_cuts   = 11

  log.joint(' running iteration ' + str(all_data['round']) + ' of hybrid algorithm\n')

  #QP
  if all_data['round'] < QP_to_LP or all_data['round'] >= QP:
    if all_data['linear_objective']:  #switching from linear to convex QP
      themodel.params.method       = 2
      themodel.Params.BarHomogeneous = 1
      themodel.Params.BarConvTol = 1e-06
      all_data['linear_objective'] = 0
      all_data['objective_cuts']   = 0
      objvar.setAttr("Obj",0)
      lincostvar.setAttr("Obj",1)
      qcostvar.setAttr("Obj",1)
      themodel.remove(lincost)
      themodel.remove(sumTconstr)
      all_data['qcost'] = qcost = themodel.addConstr(qcostexpr <= qcostvar, name = "qcost")

  #LP
  if QP_to_LP <= all_data['round'] and all_data['round'] < QP:
    if all_data['linear_objective'] == 0: #switching from convex QP to linear
      all_data['linear_objective'] = 1
      all_data['objective_cuts']   = 1
      themodel.params.method       = 1
      objvar.setAttr("Obj",1)
      lincostvar.setAttr("Obj",0)
      qcostvar.setAttr("Obj",0)
      themodel.remove(qcost)
      all_data['sumTconstr'] = sumTconstr = themodel.addConstr(sumTvars <= objvar, name= 'objvar_quad')
      all_data['lincost']    = lincost = themodel.addConstr(lincostvar <= objvar, name= 'objvar_linear')

    if all_data['round'] > no_objective_cuts:
      all_data['objective_cuts'] = 0


  themodel.update()
  hybridlpname = 'hybrid_' + all_data['casename'] + "_" + str(all_data['round']) + '.lp'
  log.joint(' writing down .lp file (hybrid)...\n')
  themodel.write(hybridlpname)

  #breakexit('check.lp')


def check13659(log,all_data):

  try:
    f = open('mp_sols/case13659pegase.out', "r")
    lines = f.readlines()
    f.close()
  except:
    log.stateandquit("cannot open file " + datafilename)
    sys.exit("failure")

  lenlines  = len(lines)
  afile     = open('mp_sols/SOL.txt',"w+")
  afile_kit = open('mp_sols/SOL_kit.txt',"w+")

  branches     = all_data['branches']
  buses        = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap'] 

  mp_vm      = {}
  mp_angle   = {}
  mp_cvalues = {}
  mp_svalues = {}

  mp_Pfvalues = {}
  mp_Ptvalues = {}
  mp_Qfvalues = {}
  mp_Qtvalues = {}

  linenum = 3
  log.joint(' reading file\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    log.joint(' the line ' + str(thisline) + '\n')
    if thisline[0] == 'bus':
      bus_id          = int(thisline[1])
      buscount        = IDtoCountmap[bus_id]
      bus             = buses[buscount]
      mp_vm[bus]      = float(thisline[3])
      mp_angle[bus]   = float(thisline[5]) * math.pi / 180
      mp_cvalues[bus] = mp_vm[bus]**2
    elif thisline[0] == 'branch':
      branchcount     = int(thisline[1])
      branch          = branches[branchcount]

      mp_Pfvalues[branch] = float(thisline[7])/all_data['baseMVA']
      mp_Ptvalues[branch] = float(thisline[9])/all_data['baseMVA']
      mp_Qfvalues[branch] = float(thisline[11])/all_data['baseMVA']
      mp_Qtvalues[branch] = float(thisline[13])/all_data['baseMVA']

    linenum += 1

  
  log.joint(' writing on ../mp_sols/SOL.txt\n')
  for bus in buses.values():
    f = bus.nodeID
    c = mp_cvalues[bus]
    log.joint('c_' + str(f) + '_' + str(f) + ' = ' + str(c) + '\n')
    afile.write('c_' + str(f) + '_' + str(f) + ' = ' + str(c) + '\n')


  for branch in branches.values():

    branchcount = branch.count
    f          = branch.f
    t          = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    bus_f  = buses[count_of_f]
    bus_t  = buses[count_of_t]

    vm_f   = mp_vm[bus_f]
    vm_t   = mp_vm[bus_t]

    mp_cvalues[branch] = c = vm_f * vm_t * math.cos(mp_angle[bus_f] - mp_angle[bus_t])
    mp_svalues[branch] = s = vm_f * vm_t * math.sin(mp_angle[bus_f] - mp_angle[bus_t])

    log.joint('c_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(c) + '\n')
    afile.write('c_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(c) + '\n')

    log.joint('s_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(s) + '\n')
    afile.write('s_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(s) + '\n')


  afile.close()


  all_data['mp_Pfvalues'] = mp_Pfvalues
  all_data['mp_Ptvalues'] = mp_Ptvalues
  all_data['mp_Qfvalues'] = mp_Qfvalues
  all_data['mp_Qtvalues'] = mp_Qtvalues

  all_data['mp_vm']      = mp_vm
  all_data['mp_angle']   = mp_angle
  all_data['mp_cvalues'] = mp_cvalues
  all_data['mp_svalues'] = mp_svalues


def cutplane_stats(log,all_data):

  themodel = all_data['themodel']
  
  log.joint('\n ******************** round statistics **********************\n') 
    
  log.joint(' casename = %s\n' % all_data['casename'] )
  log.joint(' round = %g\n' % all_data['round'] )
  log.joint(' objective = %g\n' % all_data['objval'])
  log.joint(' solver status = ' + str(all_data['optstatus']) 
            + ' solver method = ' + str(themodel.params.method) + '\n')
  if themodel.params.method == 2:
    log.joint(' crossover ' + str(themodel.params.crossover) + '\n')
  log.joint(' BarConTol = ' + str(themodel.params.BarConvTol) 
            + ' FeasTol = ' + str(themodel.params.FeasibilityTol) 
            + ' OptTol = ' + str(themodel.params.OptimalityTol) + '\n') 

  #log.joint(' -- active power --\n')
  #log.joint(' total active power generation = %g\n' % all_data['total_active_gen'] )
  #log.joint(' active power losses = %g\n' % all_data['total_active_losses'] )

  #log.joint(' -- reactive power --\n')
  #log.joint(' total reactive power generation = %g\n' % all_data['total_reactive_gen'] )
  #log.joint(' reactive power net losses = %g\n' % all_data['total_reactive_losses'] )
  #log.joint(' reactive power gains = %g\n' % - all_data['total_reactive_gains'] )

  if all_data['addcuts'] and all_data['round'] == 1:
    log.joint(' -- precomputed cuts --\n')
    log.joint(' Jabr-envelope cuts = %d\n' 
               % all_data['addcuts_numjabrcuts'])
    log.joint(' i2-envelope cuts = %d\n' 
               % all_data['addcuts_numi2cuts'])
    log.joint(' Limit-envelope cuts = %d\n' 
               % all_data['addcuts_numlimitcuts'])
  else:
    log.joint(' -- cut parameters --\n')
    log.joint(' cut age limit = ' 
              + str(all_data['cut_age_limit']) + '\n')
    log.joint(' parallel-cuts threshold = ' 
              + str(all_data['threshold_dotprod']) + '\n')
    log.joint(' initial threshold = ' 
              + str(all_data['initial_threshold']) + '\n') #check what to do with this...
    log.joint(' threshold = ' 
              + str(all_data['tolerance']) + '\n') #check what to do with this...

    if all_data['jabrcuts']:
      log.joint(' -- Jabr-envelope cuts --\n')
      log.joint(' cuts = %d\n' 
                % all_data['num_jabr_cuts'])
      log.joint(' top percent of most violated cuts added = %g\n' 
                % (100*all_data['most_violated_fraction_jabr']) )
      log.joint(' max error (abs) = %g\n' 
                % all_data['max_error_jabr'])
      log.joint(' cut threshold = %g\n' 
                % all_data['threshold']) #check this, one threshold for all?

    if all_data['i2cuts']:
      log.joint(' -- i2-envelope cuts --\n')
      log.joint(' cuts = %d\n' 
                % all_data['num_i2_cuts'])
      log.joint(' top percent of most violated cuts added = %g\n' 
                % (100*all_data['most_violated_fraction_i2']))
      log.joint(' max error (abs) = %g\n'
                % all_data['max_error_i2'])
      log.joint(' cut threshold = %g\n' 
                % all_data['threshold_i2']) #check this as well

    if all_data['limitcuts']:
      log.joint(' -- Limit-envelope cuts --\n')
      log.joint(' cuts = %d\n' 
                % all_data['num_limit_cuts'])
      log.joint(' top percent of most violated added = %g\n' 
                % (100*all_data['most_violated_fraction_limit']))
      log.joint(' max error (abs) = %g\n' 
                % all_data['max_error_limit'])
      log.joint(' cut threshold = %g\n' 
                % all_data['threshold_limit'])

    if all_data['loss_inequalities']:
      log.joint(' -- loss inequalities --\n')
      log.joint(' inequalities = %d\n'
                % all_data['num_loss_cuts'])
      log.joint(' top percent of most violated cuts added = %g\n'
                % (100*all_data['most_violated_fraction_loss']))
      log.joint(' max error (abs) = %g\n'
                % all_data['max_error_loss'])

  if all_data['linear_objective'] or all_data['hybrid']:
    log.joint(' -- objective cuts --\n')
    log.joint(' objective-cuts = %d\n'
              % all_data['num_objective_cuts'])
    log.joint(' threshold = %g\n'
              % all_data['threshold_objcuts'])
    
  log.joint(' -- runtimes --\n')
  log.joint(' solver runtime (current round) = %g\n'
            % all_data['solvertime'])
  log.joint(' cumulative solver time = %g\n' 
            % all_data['cumulative_solver_time'])
  
  timenow = time.time()

  log.joint(' time so far (overall - formulation time) = %g\n'
            % (timenow - all_data['T0'] - all_data['formulation_time']))
  log.joint(' time so far (overall) = %g\n'
            % (timenow - all_data['T0']))
  log.joint(' max runtime = %g\n'
            % all_data['max_time'])
  if all_data['ftol_counter']:
    log.joint(' minimum obj improvement (ftol) = %g\n'
              % all_data['ftol'])
    log.joint(' consecutive rounds with obj improvement < ftol = %d\n' 
              % all_data['ftol_counter'])

  log.joint(' ************************************************************\n') 


def cutplane_cutstats(log,all_data):

  themodel = all_data['themodel']

  log.joint('\n *********************** cuts statistics ********************\n') 
  if all_data['jabrcuts']:
    log.joint(' -- Jabr-envelope cuts --\n')
    log.joint(' cuts = %d\n'
              % all_data['num_jabr_cuts'])
    if all_data['addcuts']:
      log.joint('  precomputed cuts = %d\n' 
              % all_data['addcuts_numjabrcuts'])
    log.joint('  added in current round = %d\n' 
              % all_data['num_jabr_cuts_added'])
    if all_data['dropjabr']:
      log.joint('  dropped in current round = %d\n' 
                % all_data['num_jabr_cuts_dropped'])
    log.joint(' added (overall) = %d\n'
              % all_data['ID_jabr_cuts'])
    if all_data['dropjabr']:
      log.joint(' dropped (overall) = %d\n'
                % all_data['total_jabr_dropped'])

  if all_data['i2cuts']:
    log.joint(' -- i2-envelope cuts --\n')
    log.joint(' cuts = %d\n'
              % all_data['num_i2_cuts'])
    if all_data['addcuts']:
      log.joint('  precomputed cuts = %d\n' 
                % all_data['addcuts_numi2cuts'])
    log.joint('  added in current round = %d\n' 
              % all_data['num_i2_cuts_added'])
    if all_data['dropi2']:
      log.joint('  dropped in current round = %d\n'
                % all_data['num_i2_cuts_dropped'] )
  
    log.joint(' added (overall) = %d\n'
              % all_data['ID_i2_cuts'])

    if all_data['dropi2']:
      log.joint(' dropped (overall) = %d\n'
                % all_data['total_i2_dropped'])

  if all_data['limitcuts']:
    log.joint(' -- Limit-envelope cuts --\n')
    log.joint(' cuts = %d\n'
              % all_data['num_limit_cuts'])
    if all_data['addcuts']:
      log.joint('  precomputed cuts = %d\n' 
                % all_data['addcuts_numlimitcuts'])
    log.joint('  added in current round = %d\n'
              % all_data['num_limit_cuts_added'])
    if all_data['droplimit']:
      log.joint('  dropped in current round = %d\n'
                % all_data['num_limit_cuts_dropped'] )
    log.joint(' added (overall) = %d\n'
               % all_data['ID_limit_cuts'])
    if all_data['droplimit']:
      log.joint(' dropped (overall) = %d\n'
                % all_data['total_limit_dropped'])

  if all_data['loss_inequalities']:
    log.joint(' -- loss inequalities --\n')
    log.joint(' Loss-inequalities = %d\n'
              % all_data['num_loss_cuts'])
    log.joint('  added in current round = %d\n'
              % all_data['num_loss_cuts_added'])
    if all_data['droploss']:
      log.joint('  dropped in current round = %d\n'
                % all_data['num_loss_cuts_dropped'])

  if all_data['linear_objective'] or all_data['hybrid']:
    log.joint(' -- objective cuts --\n')
    log.joint(' objective-cuts = %g\n' % all_data['num_objective_cuts'])
    log.joint(' objective-cuts threshold = %g\n' % all_data['threshold_objcuts'])    

  log.joint(' ************************************************************\n\n') 


def cutplane_cuts(log,all_data):

  log.joint('\n starting cut procedure ...\n')

  t0_cuts = time.time()

  t0_jabr = time.time()

  if all_data['jabrcuts']:
    jabr_cuts(log,all_data)

    if all_data['NO_jabrs_violated']:
      if all_data['threshold'] > all_data['tolerance']:
        all_data['threshold'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold']) + '\n' )
        #continue
      # else:
      #   log.joint(' threshold below ' + str(all_data['tolerance']) + ', we are done\n' )
      #   if all_data['hybrid'] == 0:
      #     log.joint(' writing to lpfile ' + all_data['lpfilename_cuts'] + '\n')  
      #     themodel.write(all_data['lpfilename_cuts'])
      #     log.joint(' total runtime = ' + str(time.time() - all_data['T0']) )
      #     log.joint(' bye.\n')
      #     return None

  t1_jabr = time.time()

  log.joint(' time spent on Jabr-cuts ' + str(t1_jabr - t0_jabr) + '\n')
  #breakexit(' jabrs ')

  t0_i2 = time.time()

  if all_data['i2cuts']:
    i2_cuts(log,all_data)

    if all_data['NO_i2_cuts_violated']:
      if all_data['threshold_i2'] > all_data['tolerance']:
        all_data['threshold_i2'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold_i2']) + '\n' )
        #continue
      # else:
      #   log.joint(' threshold below ' + str(all_data['tolerance']) + ', we are done\n' )
      #   if all_data['hybrid'] == 0:
      #     log.joint(' writing to lpfile ' + all_data['lpfilename_cuts'] + '\n')  
      #     themodel.write(all_data['lpfilename_cuts'])
      #     log.joint(' total runtime = ' + str(time.time() - all_data['T0']) )
      #     log.joint(' bye.\n')
      #     return None

  t1_i2 = time.time()
  log.joint(' time spent on i2-cuts ' + str(t1_i2 - t0_i2) + '\n')


  t0_lim = time.time()

  if all_data['limitcuts']:
    limit_cuts(log,all_data)
    if all_data['NO_limit_cuts_violated']:
      if all_data['threshold_limit'] > all_data['tolerance']:
        all_data['threshold_limit'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold_limit']) + '\n' )
        #continue
      else:
        log.joint(' threshold below ' + str(all_data['tolerance']) + ', we continue\n' )

  t1_lim = time.time()
  log.joint(' time spent on lim-cuts ' + str(t1_lim - t0_lim) + '\n')


  if all_data['losscuts']:
    loss_cuts(log,all_data)
    if all_data['NO_loss_violated']:
      if all_data['threshold_limit'] > all_data['tolerance']:
        all_data['threshold_limit'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold_limit']) + '\n' )
        #continue
      else:
        log.joint(' threshold below ' + str(all_data['tolerance']) + ', we continue\n' )


  if all_data['linear_objective']:
    if all_data['objective_cuts']:
      objective_cuts(log,all_data)

  t1_cuts = time.time()

  log.joint('\n time spent on cuts ' + str(t1_cuts - t0_cuts) + '\n')

  #breakexit('time adding cuts')


def cutplane_cutmanagement(log,all_data):

  t0_cutmanag = time.time()

  if all_data['jabrcuts'] and all_data['dropjabr'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_jabr(log,all_data)

  if all_data['i2cuts'] and all_data['dropi2'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_i2(log,all_data)

  if all_data['limitcuts'] and all_data['droplimit'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_limit(log,all_data)

  if all_data['losscuts'] and all_data['droploss'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_loss(log,all_data)

  if all_data['loss_inequalities'] and all_data['droploss'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_loss(log,all_data)

  if all_data['cut_analysis']:
    cut_analysis(log,all_data)

  log.joint('\n')

  t1_cutmanag = time.time()

  #log.joint(' time spent on cut management ' + str(t1_cutmanag - t0_cutmanag) + '\n')


def cutplane_optimize(log,all_data):

  themodel = all_data['themodel']

  log.joint(' solving model with method ' + str(themodel.params.method) + '\n')
  log.joint(' crossover ' + str(themodel.params.crossover) + '\n')

  if all_data['losstest']:
    all_data['rbranch'] = randombranch(log,all_data)
  
  t0_solve = time.time()
  themodel.optimize()
  t1_solve = time.time()

  if themodel.status == GRB.status.INF_OR_UNBD:
    log.joint(' -> LP infeasible or unbounded\n')
    log.joint(' turning presolve off and reoptimizing\n')
    themodel.params.presolve = 0

    t0_solve = time.time()
    themodel.optimize()
    t1_solve = time.time()

    breakexit('for now')

    themodel.computeIIS()
    themodel.write('model.ilp')

    all_data['runtime'] = t1_solve - all_data['T0']

    log.joint(' writing casename, opt status, and runtime to summary_ws.log\n')

    summary_ws = open("summary_ws.log","a+") 
    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(themodel.status) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + '\n')

    summary_ws.close()

    log.joint(' optimization status ' + str(themodel.status) + '\n')
    log.joint(' solver runtime current round = %g\n' % (t1_solve - t0_solve) )
    log.joint(' overall time = %g\n' % all_data['runtime'] )
    log.joint(' bye.\n')
    exit(0)

  elif themodel.status == GRB.status.INFEASIBLE:
    log.joint(' -> LP infeasible\n')

    all_data['runtime'] = time.time() - all_data['T0']

    breakexit('for now')

    themodel.computeIIS()
    model.write('model.ilp')

    log.joint(' writing casename, opt status, and runtime to summary_ws.log\n')

    summary_ws = open("summary_ws.log","a+") #later add feasibility errors, etc                            
    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(themodel.status) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + '\n')

    summary_ws.close()

    log.joint(' solver runtime current round = %g\n' % (t1_solve - t0_solve) )
    log.joint(' overall time = %g\n' % all_data['runtime'] )
    log.joint(' bye.\n')
    exit(0)

    
  elif ( themodel.status != GRB.status.OPTIMAL and 
         themodel.status != GRB.status.SUBOPTIMAL ):

    log.joint(' -> solver terminated with status ' + str(themodel.status) + '\n')

    all_data['runtime'] = time.time() - all_data['T0']

    log.joint(' writing casename, opt status and runtime to summary_ws.log\n')

    summary_ws = open("summary_ws.log","a+") #later add feasibility errors, etc                            
    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(themodel.status) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + '\n')

    summary_ws.close()

    log.joint(' solver runtime current round = %g\n' % (t1_solve - t0_solve) )
    log.joint(' overall time = %g\n' % all_data['runtime'] )
    log.joint(' bye.\n')
    exit(0)

  all_data['objval']                  = themodel.ObjVal
  all_data['optstatus']               = themodel.status
  all_data['solvertime']              = t1_solve - t0_solve
  all_data['cumulative_solver_time'] += (t1_solve - t0_solve)  


def getsol_knitro2(log,all_data):

  casename    = all_data['casename']
  filename    = 'knitro_sols2/ksol_'+ casename +'.txt'

  try:
    thefile   = open(filename, "r")
    lines     = thefile.readlines()
    lenlines  = len(lines)
    thefile.close()
  except:
    log.stateandquit("cannot open file " + datafilename)
    sys.exit("failure")

  branches     = all_data['branches']
  buses        = all_data['buses']
  gens         = all_data['gens']
  IDtoCountmap = all_data['IDtoCountmap'] 

  sol_obj        = 0
  sol_vm         = {}
  sol_angle      = {}
  sol_cvalues    = {}
  sol_svalues    = {}
  sol_evalues    = {}
  sol_fvalues    = {}
  sol_Pfvalues   = {}
  sol_Ptvalues   = {}
  sol_Qfvalues   = {}
  sol_Qtvalues   = {}
  sol_GenPvalues = {}
  sol_GenQvalues = {}
  sol_IPvalues   = {}
  sol_IQvalues   = {}
  sol_allvars    = {}


  linenum = 0
  log.joint(' reading file\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    if thisline[0] == 'value':
      sol_obj              = float(lines[0].split()[1])
    elif thisline[0] == 'bus':
      buscount             = int(thisline[1])
      bus                  = buses[buscount]
      busid                = bus.nodeID
      sol_vm[bus]          = float(thisline[3])
      sol_angle[bus]       = float(thisline[5])
      v2value              = sol_vm[bus]**2
      evalue               = sol_vm[bus] * math.cos(sol_angle[bus])
      fvalue               = sol_vm[bus] * math.sin(sol_angle[bus])
      
      sol_cvalues[bus]     = v2value
      sol_evalues[bus]     = evalue
      sol_fvalues[bus]     = fvalue

      cname = 'c_' + str(busid) + '_' + str(busid)
      ename = 'e_' + str(busid)
      fname = 'f_' + str(busid)

      sol_allvars[cname] = v2value
      sol_allvars[ename] = evalue
      sol_allvars[fname] = fvalue

    elif thisline[0] == 'branch':
      branchid             = int(thisline[1])
      branch               = branches[branchid]
      f                    = branch.f
      t                    = branch.t
      sol_Pfvalues[branch] = float(thisline[7])
      sol_Ptvalues[branch] = float(thisline[9])
      sol_Qfvalues[branch] = float(thisline[11])
      sol_Qtvalues[branch] = float(thisline[13])

      Pfname  = 'P_' + str(branchid) + '_' + str(f) + '_' + str(t)
      Ptname  = 'P_' + str(branchid) + '_' + str(t) + '_' + str(f)
      Qfname  = 'Q_' + str(branchid) + '_' + str(f) + '_' + str(t)
      Qtname  = 'Q_' + str(branchid) + '_' + str(t) + '_' + str(f)

      sol_allvars[Pfname]  = float(thisline[7])
      sol_allvars[Ptname]  = float(thisline[9])
      sol_allvars[Qfname]  = float(thisline[11])
      sol_allvars[Qtname]  = float(thisline[13])

    elif thisline[0] == 'genid':
      genid      = int(thisline[1])
      gen_nodeID = int(thisline[3])
      sol_GenPvalues[genid] = float(thisline[5])
      sol_GenQvalues[genid] = float(thisline[7])

      GPname = "GP_" + str(genid) + "_" + str(gen_nodeID)
      GQname = "GQ_" + str(genid) + "_" + str(gen_nodeID)

      sol_allvars[GPname] = float(thisline[5])
      sol_allvars[GQname] = float(thisline[7])

    linenum += 1

  for branch in branches.values(): #this could be added to ksol.txt (1)
    branchid   = branch.count
    f          = branch.f
    t          = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    bus_f      = buses[count_of_f]
    bus_t      = buses[count_of_t]
    vm_f       = sol_vm[bus_f]
    vm_t       = sol_vm[bus_t]

    sol_cvalues[branch] = vm_f * vm_t * math.cos(sol_angle[bus_f] - sol_angle[bus_t])
    sol_svalues[branch] = vm_f * vm_t * math.sin(sol_angle[bus_f] - sol_angle[bus_t])
    
    cftname = 'c_' + str(branchid) + '_' + str(f) + '_' + str(t)
    sftname = 's_' + str(branchid) + '_' + str(f) + '_' + str(t)

    sol_allvars[cftname] = sol_cvalues[branch]
    sol_allvars[sftname] = sol_svalues[branch]
  
  for bus in buses.values(): #this could be added to ksol.txt (2)
    IPvalue = - bus.Pd
    IQvalue = - bus.Qd
    IPname  = 'IP_' + str(bus.nodeID)
    IQname  = 'IQ_' + str(bus.nodeID)

    for gencounter in bus.genidsbycount:
      if gens[gencounter].status:
        IPvalue += sol_GenPvalues[gencounter]
        IQvalue += sol_GenQvalues[gencounter]

    sol_IPvalues[bus] = IPvalue
    sol_IQvalues[bus] = IQvalue

    sol_allvars[IPname] = IPvalue
    sol_allvars[IQname] = IQvalue


  all_data['sol_vm']       = sol_vm
  all_data['sol_angle']    = sol_angle
  all_data['sol_cvalues']  = sol_cvalues
  all_data['sol_svalues']  = sol_svalues
  all_data['sol_evalues']  = sol_evalues
  all_data['sol_fvalues']  = sol_fvalues
  all_data['sol_Pfvalues']   = sol_Pfvalues
  all_data['sol_Ptvalues']   = sol_Ptvalues
  all_data['sol_Qfvalues']   = sol_Qfvalues
  all_data['sol_Qtvalues']   = sol_Qtvalues
  all_data['sol_GenPvalues'] = sol_GenPvalues
  all_data['sol_GenQvalues'] = sol_GenQvalues
  all_data['sol_IPvalues']   = sol_IPvalues
  all_data['sol_IQvalues']   = sol_IQvalues
  all_data['sol_allvars']    = sol_allvars

  log.joint(' knitro solution loaded\n')

  #log.joint(' sol_GenPvalues ' + str(sol_GenPvalues) + '\n')
  #flowdecomp(log,all_data)

  strictchecker(log,all_data)
  exit(0)
  

def flowdecomp(log,all_data):
  
  buses          = all_data['buses']
  branches       = all_data['branches']
  gens           = all_data['gens']
  IDtoCountmap   = all_data['IDtoCountmap']
  all_data['fd_threshold'] = 1e-8

  all_data['fd_flux'] = {}
  paths                 = {}
  Pd                    = {}
  losses                = {}

  for bus in buses.values():
    buscount     = bus.count
    Pd[buscount] = bus.Pd

  all_data['fd_Pd']     = Pd
  all_data['fd_loss']   = False
  all_data['fd_cycle']  = False

  sol_Pfvalues          = all_data['sol_Pfvalues']
  sol_Ptvalues          = all_data['sol_Ptvalues']
  sol_GenPvalues        = all_data['sol_GenPvalues']
  
  #log.joint(' GenPvalues before ' + str(sol_GenPvalues) + '\n')
  #log.joint(' Pd before ' + str(Pd) + '\n')

  totalgen_initial      = sum(sol_GenPvalues.values())
  totaldemand_initial   = sum(Pd.values())

  #log.joint(' initial total (net) generation ' + str(totalgen_initial) + '\n')
  #log.joint(' initial total (net) demand ' + str(totaldemand_initial) + '\n')

  ##########NET GENERATION

  for genid in sol_GenPvalues.keys():
    
    gen      = gens[genid]
    busid    = gen.nodeID 
    buscount = IDtoCountmap[busid]
    
    if sol_GenPvalues[genid] >= Pd[buscount]:
      sol_GenPvalues[genid] -= Pd[buscount]
      Pd[buscount]           = 0
    else:
      Pd[buscount]          -= sol_GenPvalues[genid]
      sol_GenPvalues[genid]  = 0

  #log.joint(' GenPvalues after ' + str(sol_GenPvalues) + '\n')
  #log.joint(' Pd after ' + str(Pd) + '\n')

  ##########################

  totalgen_initial      = sum(sol_GenPvalues.values())
  totaldemand_initial   = sum(Pd.values())
  
  log.joint(' checking losses ...\n')

  for branch in branches.values():

    losses[branch] = sol_Pfvalues[branch] + sol_Ptvalues[branch]
    log.joint(' loss at branch ' + str(branch.count) + ' : ' 
              + str((sol_Pfvalues[branch] + sol_Ptvalues[branch])) + '\n')
  
  all_data['fd_losses'] = losses  
  totallosses_initial = sum(losses.values())

  log.joint(' initial total (net) generation ' + str(totalgen_initial) + '\n')
  log.joint(' initial total (net) demand ' + str(totaldemand_initial) + '\n')
  log.joint(' initial losses ' + str(totallosses_initial) + '\n')

  #breakexit('check')

  ################

  log.joint(' initializing flow decomposition\n')

  GP_values = np.fromiter(sol_GenPvalues.values(), dtype=float) 
  GP_genids = np.fromiter(sol_GenPvalues.keys(), dtype=int)

  #log.joint(' GP values ' + str(GP_values) + '\n')
  #log.joint(' GP genids ' + str(GP_genids) + '\n')

  sorting_array    = np.argsort(GP_values)[::-1]
  sorted_GP_values = GP_values[sorting_array]
  sorted_GP_genids = GP_genids[sorting_array]

  all_data['fd_iteration'] = 0
  threshold                  = 0.1
  all_data['fd_numcycles'] = 0

  while ((sum(all_data['fd_losses'].values()) 
          + sum(Pd.values()) ) > threshold): # before loop sum(sorted_GP_values)
    
    all_data['fd_iteration'] += 1
    
    maxgen_genid = sorted_GP_genids[0]
    maxgen_value = sorted_GP_values[0]
    gen          = gens[maxgen_genid]
    gen_nodeID   = gen.nodeID
    buscount     = IDtoCountmap[gen_nodeID]
    bus          = buses[buscount]


    all_data['fd_maxgen_value'] = maxgen_value #new 
    all_data['fd_done']         = False
    all_data['fd_cycle']        = False
    prepath                     = {}
    pathlen                     = 0
    nocycling                   = []
    all_data['fd_push']         = 1e20

    nocycling.append(gen_nodeID)

    log.joint(' iter ' + str(all_data['fd_iteration']) + ' maxg: ' 
              + str(maxgen_value) + ' at bus ' + str(gen_nodeID) + ' genid ' 
              + str(maxgen_genid) + '\n')
    
    # Computing largest flow
    while (all_data['fd_done'] == False):
      pathlen += 1
      nextbus_id, prepath, paths, nocycling = heavybranch(log,all_data,paths,bus,prepath,pathlen,threshold,nocycling)
      if nextbus_id != 'cycle':
        bus = buses[nextbus_id]

    # Tracing paths/cycles and update demands and gen
    if all_data['fd_cycle']:
      all_data['fd_numcycles'] += 1
      tracecycle(log,all_data)
    else:
      tracepath(log,all_data,Pd,prepath,gen)

    # Update GP_values
    if all_data['fd_cycle'] == False:
      sorted_GP_genids, sorted_GP_values = updategeneration(log,all_data,gen,sorted_GP_values,sorted_GP_genids)

    #breakexit('done w iteration')
    log.joint(' -------------\n')

  #log.joint('  sorted GP values ' + str(sorted_GP_values) + '\n')
  #log.joint('  sorted GP genids ' + str(sorted_GP_genids) + '\n')

  log.joint('  total gen ' + str(totalgen_initial) + '\n')
  log.joint('  total processed gen (flowdecomp) ' + str(sum(sorted_GP_values)) + '\n')
  log.joint('\n  total demand ' + str(totaldemand_initial) + '\n')
  log.joint('  total processed demand (flowdecomp) ' + str(sum(all_data['fd_Pd'].values())) + '\n')
  log.joint('\n  total losses ' + str(totallosses_initial) + '\n')
  log.joint('  total processed losses (flowdecomp) ' + str(sum(all_data['fd_losses'].values())) + '\n')
  #log.joint('  paths (flux) ' + str(all_data['fd_flux']) + '\n')

  log.joint('\n')
  flux    = sum(all_data['fd_flux'].values())
  avgflux = flux/all_data['fd_iteration']

  log.joint('  number of paths ' + str(all_data['fd_iteration']) + ' avg flux ' 
            + str(avgflux) + '\n')

  log.joint('  number of cycles ' + str(all_data['fd_numcycles']) + '\n')

  #breakexit('c')

def tracepath(log,all_data,Pd,prepath,gen):

  IDtoCountmap          = all_data['IDtoCountmap']
  branches              = all_data['branches']
  sol_Pfvalues          = all_data['sol_Pfvalues']
  sol_Ptvalues          = all_data['sol_Ptvalues']
  prepathlen            = len(prepath)

  log.joint('  now tracing path of length ' + str(prepathlen) + '\n')
  #log.joint('  the prepath is ' + str(prepath) + '\n')

  for item in prepath.values():
    branchid = item[0]
    kind     = item[2]
    loss     = item[3]
    branch   = branches[branchid]
    
    if kind == 'f' and loss == 0:
      newPfvalue = sol_Pfvalues[branch] - all_data['fd_push']
      newPtvalue = sol_Ptvalues[branch] + all_data['fd_push']
      initialbus = branch.f
      finalbus   = branch.t
    elif kind == 't' and loss == 0:
      newPfvalue = sol_Pfvalues[branch] + all_data['fd_push']
      newPtvalue = sol_Ptvalues[branch] - all_data['fd_push']
      initialbus = branch.t
      finalbus = branch.f
    elif kind == 'f' and loss:
      newPfvalue = sol_Pfvalues[branch] - all_data['fd_push']
      newPtvalue = sol_Ptvalues[branch]
      initialbus = branch.f
      finalbus = branch.t
    elif kind == 't' and loss:
      newPtvalue = sol_Ptvalues[branch] - all_data['fd_push']
      newPfvalue = sol_Pfvalues[branch]
      initialbus = branch.t
      finalbus = branch.f
      
    log.joint('   bus ' + str(initialbus) + ' branch ' + str(branchid) + ' ' + str(branch.f)
              + ' -> ' + str(branch.t) + ', kind: ' + kind + '\n')
    #log.joint(' oldflows: ' + str(sol_Pfvalues[branch]) + ' ' + str(sol_Ptvalues[branch]) + '\n')
    newloss = newPfvalue + newPtvalue

    log.joint('   newflows: ' + str(newPfvalue) + ' ' + str(newPtvalue) 
              + ' newloss ' + str(newloss) + '\n')

    all_data['sol_Pfvalues'][branch] = newPfvalue
    all_data['sol_Ptvalues'][branch] = newPtvalue
    all_data['fd_losses'][branch]    = newloss

  all_data['finalbus'] = finalbus

  log.joint('  final bus: ' + str(finalbus) + '\n')


def tracecycle(log,all_data):

  IDtoCountmap          = all_data['IDtoCountmap']
  branches              = all_data['branches']
  sol_Pfvalues          = all_data['sol_Pfvalues']
  sol_Ptvalues          = all_data['sol_Ptvalues']
  cycle                 = all_data['thecycle']

  log.joint('  now tracing cycle of length ' + str(len(cycle)) + '\n')
  log.joint('  the cycle is ' + str(cycle) + '\n')

  for item in cycle.values():
    branchid = item[0]
    kind     = item[2]
    branch   = branches[branchid]
    
    if kind == 'f':
      newPfvalue = sol_Pfvalues[branch] - all_data['cycle_flow']
      newPtvalue = sol_Ptvalues[branch] + all_data['cycle_flow']
      initialbus = branch.f
      finalbus   = branch.t
    elif kind == 't':
      newPfvalue = sol_Pfvalues[branch] + all_data['cycle_flow']
      newPtvalue = sol_Ptvalues[branch] - all_data['cycle_flow'] #HOMOGENEIZE
      initialbus = branch.t
      finalbus = branch.f
      
    log.joint('   bus ' + str(initialbus) + ' branch ' + str(branchid) + ' ' + str(branch.f)
              + ' -> ' + str(branch.t) + ', kind: ' + kind + '\n')
    #log.joint(' oldflows: ' + str(sol_Pfvalues[branch]) + ' ' + str(sol_Ptvalues[branch]) + '\n')
    newloss = newPfvalue + newPtvalue

    log.joint('   newflows: ' + str(newPfvalue) + ' ' + str(newPtvalue) + '\n')

    all_data['sol_Pfvalues'][branch] = newPfvalue
    all_data['sol_Ptvalues'][branch] = newPtvalue
    all_data['fd_push'] = all_data['cycle_flow']

  breakexit('c')

def updategeneration(log,all_data,gen,sorted_GP_values,sorted_GP_genids):

  IDtoCountmap          = all_data['IDtoCountmap']
  branches              = all_data['branches']
  Pd                    = all_data['fd_Pd']
  finalbus              = all_data['finalbus']
  finalbus_id           = IDtoCountmap[finalbus]
  
  maxgen_value = sorted_GP_values[0]
  maxgen_genid = gen.count
  
  newmaxgen_value  = maxgen_value - all_data['fd_push']
  if newmaxgen_value < 0:
    log.joint(' negative maxgen\n')
    breakexit(' check!!')

  if all_data['fd_loss'] == False: #check this!!!!
    proxy_newdemand = Pd[finalbus_id] - all_data['fd_push']
    if proxy_newdemand < 0:
      Pd[finalbus_id] = 0
    else:
      Pd[finalbus_id] = proxy_newdemand
  else:
    all_data['fd_loss'] = False

  log.joint('  now generation at ' + str(gen.nodeID) + ' : ' 
            + str(newmaxgen_value) + ', demand at ' 
            + str(all_data['finalbus']) + ' : ' + str(Pd[finalbus_id]) 
            + '\n')

  sorted_GP_values[0] = newmaxgen_value 
  sorting_array       = np.argsort(sorted_GP_values)[::-1]
  sorted_GP_values    = sorted_GP_values[sorting_array]
  sorted_GP_genids    = sorted_GP_genids[sorting_array]


  return sorted_GP_genids, sorted_GP_values


def catchcycle(log,all_data,bus,prepath,pathlen,threshold,nocycling):

  IDtoCountmap          = all_data['IDtoCountmap']
  branches              = all_data['branches']
  branches_from         = bus.frombranchids.values()
  branches_to           = bus.tobranchids.values()
  sol_Pfvalues          = all_data['sol_Pfvalues']
  sol_Ptvalues          = all_data['sol_Ptvalues']
  Pd                    = all_data['fd_Pd']

  largestflow_value    = 0
  largestflow_branchid = -1
  largestflow_kind     = 'foo'
  largestflow_nextbus  = 'foo2'

  log.joint(' pass to catch cycle!\n')

  for branchid in branches_from:
    branch = branches[branchid]
    log.joint(' branch ' + str(branchid) + ' f ' + str(branch.f) 
              + ' t ' + str(branch.t) + ' Pf ' + str(sol_Pfvalues[branch]) 
              + ' Pt ' + str(sol_Ptvalues[branch]) + '\n')

    if (sol_Pfvalues[branch] > largestflow_value  
        and branch.t in nocycling):

      largestflow_value    = sol_Pfvalues[branch]
      largestflow_branchid = branchid
      largestflow_kind     = 'f'
      largestflow_nextbus  = branch.t
      largestflow_nextbus_id = IDtoCountmap[branch.t]

  for branchid in branches_to:
    branch = branches[branchid]

    log.joint(' branch ' + str(branchid) + ' f ' + str(branch.f) 
              + ' t ' + str(branch.t) + ' Pt ' + str(sol_Ptvalues[branch]) 
              + ' Pf ' + str(sol_Pfvalues[branch]) + '\n')

    if (sol_Ptvalues[branch] > largestflow_value 
        and branch.f in nocycling):

      largestflow_value     = sol_Ptvalues[branch]
      largestflow_branchid  = branchid
      largestflow_kind      = 't'
      largestflow_nextbus   = branch.f
      largestflow_nextbus_id = IDtoCountmap[branch.f] 
  

  cycle    = {}  
  cycle[0] = (largestflow_branchid,largestflow_value,largestflow_kind,0)
  all_data['cycle_flow'] = largestflow_value

  log.joint(' cycle is given by:\n')  
  log.joint(' pathlen ' + str(pathlen) + '\n')
  log.joint(' first bus in cycle, bus id ' + str(largestflow_nextbus) + '\n') 
  
  branch = branches[largestflow_branchid]
  
  log.joint(' arc0 branchid ' + str(largestflow_branchid) + ' Pf ' + str(sol_Pfvalues[branch])
            + ' Pt ' + str(sol_Ptvalues[branch]) + '\n')

  i           = pathlen - 1
  cycle_count = 1

  while i <= pathlen:
    
    busid = nocycling[i]

    if busid != largestflow_nextbus:
      arc    = prepath[i]
      branch = branches[arc[0]]
      flow   = arc[1]
      log.joint(' arc ' + str(arc) + ' branchid ' + str(arc[0]) 
                + ' f ' + str(branch.f) + ' t ' + str(branch.t) + '\n')
      
      cycle[cycle_count] = arc
      if all_data['cycle_flow'] > flow:
        all_data['cycle_flow'] = flow

      log.joint(' ' + str(pathlen - i + 1) + 'th bus in cycle, bus id ' + str(busid) + '\n') 
      i           -= 1
      cycle_count += 1
      
    else:
      break

  log.joint(' cycle ' + str(cycle) + '\n')

  breakexit('c')

  if largestflow_value < 0:
    log.joint(' check, something off\n')
    breakexit('something off')

  log.joint('  push will be ' + str(all_data['cycle_flow']) + ' through cycle\n')

  all_data['fd_done']   = True
  all_data['thecycle']  = cycle

def heavybranch(log,all_data,paths,bus,prepath,pathlen,threshold,nocycling):

  IDtoCountmap          = all_data['IDtoCountmap']
  branches              = all_data['branches']
  branches_from         = bus.frombranchids.values()
  branches_to           = bus.tobranchids.values()
  sol_Pfvalues          = all_data['sol_Pfvalues']
  sol_Ptvalues          = all_data['sol_Ptvalues']
  Pd                    = all_data['fd_Pd']

  largestflow_value    = 0
  largestflow_branchid = -1
  largestflow_kind     = 'foo'
  largestflow_nextbus  = 'foo2'

  # log.joint(' branches from ' + str(branches_from) +  ' at bus ' 
  #           + str(bus.nodeID) + '\n')

  # log.joint(' branches to ' + str(branches_to) +  ' at bus ' 
  #           + str(bus.nodeID) + '\n')

  for branchid in branches_from:
    branch = branches[branchid]
    log.joint(' branch ' + str(branchid) + ' f ' + str(branch.f) 
              + ' t ' + str(branch.t) + ' Pf ' + str(sol_Pfvalues[branch]) + '\n')

    if (sol_Pfvalues[branch] > largestflow_value  
        and branch.t not in nocycling):

      largestflow_value    = math.fabs(sol_Pfvalues[branch])
      largestflow_branchid = branchid
      largestflow_kind     = 'f'
      largestflow_nextbus  = branch.t
      largestflow_nextbus_id = IDtoCountmap[branch.t]

  for branchid in branches_to:
    branch = branches[branchid]

    log.joint(' branch ' + str(branchid) + ' f ' + str(branch.f) 
              + ' t ' + str(branch.t) + ' Pt ' + str(sol_Ptvalues[branch]) + '\n')

    if (sol_Ptvalues[branch] > largestflow_value 
        and branch.f not in nocycling):

      largestflow_value     = math.fabs(sol_Ptvalues[branch])
      largestflow_branchid  = branchid
      largestflow_kind      = 't'
      largestflow_nextbus   = branch.f
      largestflow_nextbus_id = IDtoCountmap[branch.f] 
      
  if largestflow_branchid == -1:
    log.joint(' deal with the cycle!\n')
    all_data['fd_cycle'] = True
    log.joint(' prepath ' + str(prepath) + '\n')
    log.joint(' no cycling list ' + str(nocycling) + '\n')
      
    breakexit(' there is a cycle!' )
    catchcycle(log,all_data,bus,prepath,pathlen,threshold,nocycling)
    
    breakexit('check again')
    return 'cycle', prepath, paths, nocycling

    
  log.joint(' largestflow branchid ' + str(largestflow_branchid) + '\n')
  largestflow_branch = branches[largestflow_branchid]
  largestflow_loss   = sol_Pfvalues[largestflow_branch] + sol_Ptvalues[largestflow_branch]

  nocycling.append(largestflow_nextbus)

  if largestflow_value < 0:
    log.joint(' check, something off\n')
    breakexit('something off')

  if largestflow_value < all_data['fd_push']:
    all_data['fd_push'] = largestflow_value

  if all_data['fd_maxgen_value'] < all_data['fd_push']: # new
    all_data['fd_push'] = all_data['fd_maxgen_value']

  log.joint('  largest flow ' + str(largestflow_value) + ' at branch '
            + str(largestflow_branchid) + ' kind ' + str(largestflow_kind)
            + '\n')

  log.joint('   next bus: ' + str(largestflow_nextbus) + ' d: ' 
            + str(Pd[largestflow_nextbus_id]) 
            + ' loss : ' + str(largestflow_loss) + '\n')

  if largestflow_loss > all_data['fd_threshold']:
    log.joint('  hit loss at branch ' + str(largestflow_branchid) 
              + ' currloss: ' + str(largestflow_loss) + '\n')

    if all_data['fd_push'] - largestflow_loss > 0:
      all_data['fd_push'] = largestflow_loss

    log.joint('  push will be ' + str(all_data['fd_push']) + '\n')

    all_data['fd_flux'][all_data['fd_iteration']] = pathlen * all_data['fd_push']
    prepath[pathlen] = (largestflow_branchid,largestflow_value,largestflow_kind,1)

    all_data['fd_done'] = True
    all_data['fd_loss'] = True

    return largestflow_nextbus_id, prepath, paths, nocycling

  elif Pd[largestflow_nextbus_id] > 0: 
    log.joint('  hit load at bus ' + str(largestflow_nextbus) 
              + ' currdemand: ' + str(Pd[largestflow_nextbus_id]) + '\n')

    if all_data['fd_push'] - Pd[largestflow_nextbus_id] > 0:
      all_data['fd_push'] = Pd[largestflow_nextbus_id]

    log.joint('  push will be ' + str(all_data['fd_push']) + '\n')

    all_data['fd_flux'][all_data['fd_iteration']] = pathlen * all_data['fd_push']
   
    prepath[pathlen] = (largestflow_branchid,largestflow_value,largestflow_kind,0)
                       
    all_data['fd_done'] = True

    return largestflow_nextbus_id, prepath, paths, nocycling
  else:
    log.joint('   moved to bus ' + str(largestflow_nextbus) + '\n')

    prepath[pathlen] = (largestflow_branchid,largestflow_value,largestflow_kind,0)
                         
    return largestflow_nextbus_id, prepath, paths, nocycling
    

def randombranch(log,all_data):

  IDtoCountmap = all_data['IDtoCountmap']
  branches     = all_data['branches']
  numbranches  = all_data['numbranches'] 
  buses        = all_data['buses']
  themodel     = all_data['themodel']

  #pick a random branch
  randombranch = random.randrange(1,numbranches+1)
  branch       = branches[randombranch]

  while branch.status == 0:
      randombranch = random.randrange(1,numbranches+1)
      branch       = branches[randombranch]

  f          = branch.f
  t          = branch.t
  constrname = "jabr_"+str(randombranch)+"_"+str(f)+"_"+str(t)
  busf       = buses[IDtoCountmap[f]]
  bust       = buses[IDtoCountmap[t]]
  Vmaxf      = busf.Vmax
  Vmaxt      = bust.Vmax
  jabrbound  = 2*( (Vmaxf**2) * (Vmaxt**2) )

  log.joint(' changing RHS of Jabr ' + str(randombranch) + ' to ' 
            + str(jabrbound) + '\n')

  themodel.write("check.lp")

  constrs = themodel.getQConstrs()

  log.joint(' constrs ' + str(constrs) + '\n')
  
  constr = constrs[randombranch]

  log.joint(' the constr ' + str(constr) + '\n')
  #constr = themodel.getqconstrbyname(constrname)

  constr.setAttr("QCRHS",jabrbound)
  #constr.setAttr("QCRHS",0)

  themodel.update()

  themodel.write("check.lp")

  log.joint(' resolving SOCP ...\n')

  return randombranch



def strictchecker(log,all_data):
  
  log.joint(' initializing strict checker ...\n')

  buses          = all_data['buses']
  branches       = all_data['branches']
  gens           = all_data['gens']
  IDtoCountmap   = all_data['IDtoCountmap']
  sol_Pfvalues   = all_data['sol_Pfvalues']
  sol_Ptvalues   = all_data['sol_Ptvalues']
  sol_Qfvalues   = all_data['sol_Qfvalues']
  sol_Qtvalues   = all_data['sol_Qtvalues']
  sol_GenPvalues = all_data['sol_GenPvalues']
  sol_GenQvalues = all_data['sol_GenQvalues']
  sol_vm         = all_data['sol_vm']    
  sol_angle      = all_data['sol_angle']
  sol_cvalues    = all_data['sol_cvalues']
  sol_svalues    = all_data['sol_svalues']
  sol_evalues    = all_data['sol_evalues'] 
  sol_fvalues    = all_data['sol_fvalues'] 
  sol_allvars    = all_data['sol_allvars'] 

  log.joint(' now reading QCQP ...\n')

  qcqpfilename = 'qcqps/qcqp_' + all_data['casename'] + '.lp'
  qcqpmodel    = all_data['qcqpmodel'] = read(qcqpfilename)

  qcqp_LContrs  = qcqpmodel.getConstrs()
  qcqp_QConstrs = qcqpmodel.getQConstrs()

  max_violation_string = 'none'
  max_violation_value  = 0

  all_data['violation'] = {}

  vmagviol        = all_data['violation']['vmagviol']    = {}
  GPviol          = all_data['violation']['GPviol']      = {}
  GQviol          = all_data['violation']['GQviol']      = {}
  branchlimitviol = all_data['violation']['branchlimit'] = {}
  Pdefviol        = all_data['violation']['Pdef']        = {}
  Qdefviol        = all_data['violation']['Qdef']        = {}
  PBalviol        = all_data['violation']['PBalviol']    = {}
  QBalviol        = all_data['violation']['QBalviol']    = {}

  log.joint(' checking violation of variable bounds ...\n')

  #log.joint('  voltages\n')
  
  max_vmagviol_string = 'none'
  max_vmagviol_value  = 0

  for bus in buses.values():
    
    name     = 'bus_' + str(bus.nodeID)
    value    = sol_vm[bus]
    viol_max = 0
    viol_min = 0

    if value > bus.Vmax:
      viol_max = value - bus.Vmax

    if value < bus.Vmin:
      viol_min = value - bus.Vmin

    if viol_max > viol_min:
      viol = viol_max
    else:
      viol = viol_min

    vmagviol[bus] = viol
    
    if viol > max_vmagviol_value:
      max_vmagviol_value  = viol
      max_vmagviol_string = name

  log.joint('  max violation of voltage magnitude ' 
            + str(max_vmagviol_value) 
            + ' at ' + max_vmagviol_string + '\n')

  if max_vmagviol_value > max_violation_value:
    max_violation_value  = max_vmagviol_value
    max_violation_string = max_vmagviol_string

  #log.joint('  GP and GQ\n')

  max_GPviol_string = 'none'
  max_GPviol_value  = 0

  max_GQviol_string = 'none'
  max_GQviol_value  = 0

  for gen in gens.values():
    genid    = gen.count
    busid    = gen.nodeID
    bus      = buses[IDtoCountmap[busid]]

    if bus.nodetype == 3 or gen.status == 0: # slack bus, we skip it
      continue

    gen_name = 'gen_' + str(genid) + '_' + str(busid)
    GPvalue  = sol_GenPvalues[genid]
    GQvalue  = sol_GenQvalues[genid]

    GPviol_max = max(GPvalue - gen.Pmax,0)
    GPviol_min = max(gen.Pmin - GPvalue,0)

    if GPviol_max > GPviol_min:
      GPviol_value = GPviol_max
    else:
      GPviol_value = GPviol_min

    GPviol[gen] = GPviol_value

    if GPviol_value > max_GPviol_value:
      max_GPviol_value  = GPviol_value
      max_GPviol_string = gen_name 

    GQviol_max = max(GQvalue - gen.Qmax,0)
    GQviol_min = max(gen.Qmin - GQvalue,0)

    if GQviol_max > GQviol_min:
      GQviol_value = GQviol_max
    else:
      GQviol_value = GQviol_min

    GQviol[gen] = GQviol_value

    if GQviol_value > max_GQviol_value:
      max_GQviol_value  = GQviol_value
      max_GQviol_string = gen_name 
    
    if gen.status == 0 and (GPvalue != 0 or GQvalue != 0):
      log.joint(' gen status is 0 but gen has output, BUG\n')
      breakexit('check')
      
  log.joint('  max violation active power gen ' 
            + str(max_GPviol_value) 
            + ' at ' + max_GPviol_string + '\n')
  log.joint('  max violation reactive power gen ' 
            + str(max_GQviol_value) 
            + ' at ' + max_GQviol_string + '\n')

  if max_GPviol_value > max_violation_value:
    max_violation_value  = max_GPviol_value
    max_violation_string = max_GPviol_string

  if max_GQviol_value > max_violation_value:
    max_violation_value  = max_GQviol_value
    max_violation_string = max_GQviol_string

  # for gen in gens.values():
  #   log.joint(' gen ' + str(gen.count) + ' at bus ' + str(gen.nodeID)
  #             + ' GP ' + str(sol_GenPvalues[gen.count])
  #             + ' GQ ' + str(sol_GenQvalues[gen.count]) + '\n')

  #   log.joint(' gen ' + str(gen.count) + ' at bus ' + str(gen.nodeID) 
  #             + ' Pmax ' + str(gen.Pmax) + ' Pmin ' + str(gen.Pmin) 
  #             + ' Qmax ' + str(gen.Qmax) + ' Qmin ' + str(gen.Qmin) + '\n')
  # log.joint(' Pdic ' + str(sol_GenPvalues) + '\n')
  # log.joint(' Qdic ' + str(sol_GenQvalues) + '\n')
              

  log.joint(' checking violations of constraints ...\n')

  #log.joint('  power flow defintions\n')

  max_Pdefviol_string = 'none'
  max_Pdefviol_value  = 0

  max_Qdefviol_string = 'none'
  max_Qdefviol_value  = 0

  for branch in branches.values():
    branchid       = branch.count
    branch_f       = branch.f
    branch_t       = branch.t
    count_of_f     = IDtoCountmap[branch_f]
    bus_f          = buses[count_of_f]

    Pdefconstrname_f = 'Pdef_' + str(branchid) + '_' + str(branch_f) + '_' + str(branch_t)
    Pdefconstrname_t = 'Pdef_' + str(branchid) + '_' + str(branch_t) + '_' + str(branch_f)
    Qdefconstrname_f = 'Qdef_' + str(branchid) + '_' + str(branch_f) + '_' + str(branch_t)
    Qdefconstrname_t = 'Qdef_' + str(branchid) + '_' + str(branch_t) + '_' + str(branch_f)
    
    Pdef_f = qcqpmodel.getConstrByName(Pdefconstrname_f)    
    Pdefviol_f, Pdefviol_f_string = checkviol_linear(log,all_data,Pdef_f,Pdefconstrname_f)

    Pdef_t = qcqpmodel.getConstrByName(Pdefconstrname_t)
    Pdefviol_t, Pdefviol_t_string = checkviol_linear(log,all_data,Pdef_t,Pdefconstrname_t)
    
    if Pdefviol_f > Pdefviol_t:
      Pdefviol_value   = Pdefviol_f
      Pdefviol_string  = Pdefviol_f_string
      Pdefviol[branch] = Pdefviol_value
    else:
      Pdefviol_value   = Pdefviol_t
      Pdefviol_string  = Pdefviol_t_string
      Pdefviol[branch] = Pdefviol_value

    if Pdefviol_value > max_Pdefviol_value:
      max_Pdefviol_string  = Pdefviol_string
      max_Pdefviol_value   = Pdefviol_value

    Qdef_f = qcqpmodel.getConstrByName(Qdefconstrname_f)
    Qdefviol_f, Qdefviol_f_string = checkviol_linear(log,all_data,Qdef_f,Qdefconstrname_f)

    Qdef_t = qcqpmodel.getConstrByName(Qdefconstrname_t)
    Qdefviol_t, Qdefviol_t_string = checkviol_linear(log,all_data,Qdef_t,Qdefconstrname_t)

    if Qdefviol_f > Qdefviol_t:
      Qdefviol_value   = Qdefviol_f
      Qdefviol_string  = Qdefviol_f_string
      Qdefviol[branch] = Qdefviol_value
    else:
      Qdefviol_value   = Qdefviol_t
      Qdefviol_string  = Qdefviol_t_string
      Qdefviol[branch] = Qdefviol_value

    if Qdefviol_value > max_Qdefviol_value:
      max_Qdefviol_string  = Qdefviol_string
      max_Qdefviol_value   = Qdefviol_value

  log.joint('  max violation of active power flow definitions ' 
            + str(max_Pdefviol_value) + ' at ' + max_Pdefviol_string + '\n')
  log.joint('  max violation of reactive power flow definitions ' 
            + str(max_Qdefviol_value) + ' at ' + max_Qdefviol_string + '\n')

  if max_Pdefviol_value > max_violation_value:
    max_violation_value  = max_Pdefviol_value
    max_violation_string = max_Pdefviol_string

  if max_Qdefviol_value > max_violation_value:
    max_violation_value  = max_Qdefviol_value
    max_violation_string = max_Qdefviol_string
  
  #log.joint('  power balance\n')

  max_PBalviol_string  = 'none'
  max_PBalviol_value = 0

  max_QBalviol_string  = 'none'
  max_QBalviol_value = 0

  for bus in buses.values():
    busid       = bus.nodeID
    buscount    = bus.count

    PBalconstrname = 'PBaldef' + str(buscount) + '_' + str(busid) 
    QBalconstrname = 'QBaldef' + str(buscount) + '_' + str(busid) 

    #log.joint(' constrname ' + PBalconstrname + '\n')
    PBalconstr = qcqpmodel.getConstrByName(PBalconstrname)

    PBalviol_value, PBalviol_string = checkviol_linear(log,all_data,PBalconstr,PBalconstrname)

    PBalviol[bus] = PBalviol_value

    if PBalviol_value > max_PBalviol_value:
      max_PBalviol_value = PBalviol_value
      max_PBalviol_string = PBalviol_string

    QBalconstr = qcqpmodel.getConstrByName(QBalconstrname)
    QBalviol_value, QBalviol_string = checkviol_linear(log,all_data,QBalconstr,QBalconstrname)

    QBalviol[bus] = QBalviol_value

    if QBalviol_value > max_QBalviol_value:
      max_QBalviol_value = QBalviol_value
      max_QBalviol_string = QBalviol_string

  log.joint('  max violation of active power balance constraints ' 
            + str(max_PBalviol_value) + ' at ' + max_PBalviol_string + '\n')
  log.joint('  max violation of reactive power balance constraints ' 
            + str(max_QBalviol_value) + ' at ' + max_QBalviol_string + '\n')

  if max_PBalviol_value > max_violation_value:
    max_violation_value  = max_PBalviol_value
    max_violation_string = max_PBalviol_string

  if max_QBalviol_value > max_violation_value:
    max_violation_value  = max_QBalviol_value
    max_violation_string = max_QBalviol_string

  #log.joint('  branch limits\n')

  max_branchviol_string = 'none'
  max_branchviol_value  = 0

  for branch in branches.values():
    
    fromvalue = math.sqrt(sol_Pfvalues[branch]*sol_Pfvalues[branch] + 
                          sol_Qfvalues[branch]*sol_Qfvalues[branch])
    fromviol  = max(fromvalue - branch.limit,0)

    tovalue   = math.sqrt(sol_Ptvalues[branch]*sol_Ptvalues[branch] + 
                          sol_Qtvalues[branch]*sol_Qtvalues[branch])
    toviol    = max(tovalue - branch.limit,0)

    if fromviol > toviol:
      thisviol   = fromviol
      branchname = 'branch_'+str(branch.count)+'_from'
      branchlimitviol[branch] = thisviol
    else:
      thisviol = toviol
      branchname = 'branch_'+str(branch.count)+'_to'
      branchlimitviol[branch] = thisviol

    if thisviol > max_branchviol_value:
      max_branchviol_value  = thisviol
      max_branchviol_string = branchname
  
  log.joint('  max violation branchlimits ' 
            + str(max_branchviol_value) 
            + ' at ' + max_branchviol_string + '\n')

  if max_branchviol_value > max_violation_value:
    max_violation_value  = max_branchviol_value
    max_violation_string = max_branchviol_string

  log.joint(' max violation ' + str(max_violation_value) + ' at ' 
            + max_violation_string + '\n')


  breakexit('end checker')

  # log.joint(' checking volt2 variables ...\n')
  
  # max_v2viol      = 0
  # max_v2viol_name = ''

  # for bus in buses.values():
  #   expr_val = math.fabs(sol_cvalues[bus] - sol_evalues[bus] * sol_evalues[bus] - sol_fvalues[bus] * sol_fvalues[bus])

  #   if expr_val > max_v2viol:
  #     max_v2viol      = expr_val
  #     max_v2viol_name = 'cbusdef_' + str(bus.nodeID) + '_' + str(bus.nodeID)
  
  # log.joint('  max violation of volt2 ' + str(max_v2viol)
  #           + ' at ' + max_v2viol_name + '\n')

  # log.joint(' checking e,f variables ...\n')
  
  # max_efviol      = 0
  # max_efviol_name = ''

  # for branch in branches.values():
    
  #   branch_f   = branch.f
  #   branch_t   = branch.t
  #   count_of_f = IDtoCountmap[branch_f]
  #   count_of_t = IDtoCountmap[branch_t]
  #   bus_f      = buses[count_of_f]
  #   bus_t      = buses[count_of_t]
    
  #   expr_val_c = math.fabs(sol_cvalues[branch] - sol_evalues[bus_f] * sol_evalues[bus_t] - sol_fvalues[bus_f] * sol_fvalues[bus_t])
  #   expr_val_s = math.fabs(sol_svalues[branch] + sol_evalues[bus_f] * sol_fvalues[bus_t] - sol_evalues[bus_t] * sol_fvalues[bus_f])

  #   if expr_val_c > max_efviol:
  #     max_efviol      = expr_val_c
  #     max_efviol_name = 'cdef_' + str(branch.count) + '_' + str(bus.nodeID) + '_' + str(bus.nodeID)
  
  #   if expr_val_s > max_efviol:
  #     max_efviol      = expr_val_s
  #     max_efviol_name = 'sdef_' + str(branch.count) + '_' + str(bus.nodeID) + '_' + str(bus.nodeID)

  # log.joint('  max violation of cs variables (ef) ' + str(max_efviol)
  #           + ' at ' + max_efviol_name + '\n')
  

def checkviol_linear(log,all_data,constr,constr_string):

  sol_allvars = all_data['sol_allvars']
  qcqpmodel   = all_data['qcqpmodel']
  lhs         = qcqpmodel.getRow(constr)
  rhs         = constr.RHS
  expr_size   = lhs.size()
  expr_val    = 0 
    
  for k in range(expr_size):
    var       = lhs.getVar(k)
    varname   = var.varname
    coeff     = lhs.getCoeff(k)
    if varname not in sol_allvars.keys():
      log.joint(' the status of ' + varname + ' is 0, we skip it\n')
      continue
    expr_val += coeff * sol_allvars[varname]
  
  if constr.Sense == '=':
    viol      = math.fabs(expr_val - rhs)
  elif constr.Sense == '>':
    viol      = max(0,rhs - expr_val)
  elif constr.Sense == '<':
    viol      = max(0,expr_val - rhs)

  viol_string = constr_string

  return viol, constr_string




