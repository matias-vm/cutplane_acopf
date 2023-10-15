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
    getsol_knitro(log,all_data)

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
    #themodel.Params.FeasibilityTol = 1e-4
    #themodel.Params.OptimalityTol  = 1e-4

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

  ########################## CUTPLANE MAIN LOOP ###############################

  all_data['round']                  = 1
  all_data['runtime']                = time.time() - all_data['T0']
  all_data['cumulative_solver_time'] = 0
  all_data['ftol_counter']           = 0
  oldobj                             = 1
  gap                                = 1e20


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

    if all_data['writesol'] == 0:

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
    
    if all_data['writesol']:

      namesolfile = 'sol_ws_' + all_data['casename'] + '.txt'
      #namesolfile = 'sols_warmstarted/sol_ws_' + all_data['casename'] + '.txt'
      #namesolfile = 'cutplane_sols/cutplane_sol_' + all_data['casename'] + '.txt'
      solfile     = open(namesolfile,'a')
      solfile.write('round' + str(all_data['round']) + '\n' )
      solfile.write('obj ' + str(all_data['objval']) + '\n')

      #values
      log.joint(' storing current solution...\n')

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
        varname      = cvar[bus].varname + ' = '
        lines        = [varname,str(cvalues[bus]),'\n']
        solfile.writelines(lines)
      
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

        cvarname = cvar[branch].varname + ' = '
        svarname = svar[branch].varname + ' = '
        cslines  = [cvarname,str(cvalues[branch]),'\n',svarname,str(svalues[branch]),'\n']
        Plines   = [Pvar_f[branch].varname + ' = ',str(Pvar_f[branch].x),'\n',Pvar_t[branch].varname + ' = ',str(Pvar_t[branch].x),'\n']
        Qlines   = [Qvar_f[branch].varname + ' = ',str(Qvar_f[branch].x),'\n',Qvar_t[branch].varname + ' = ',str(Qvar_t[branch].x),'\n']

        solfile.writelines(cslines)
        solfile.writelines(Plines)
        solfile.writelines(Qlines)

        if all_data['i2']:
          i2fvalues[branch] = i2var_f[branch].x
          i2line            = [i2var_f[branch].varname + ' = ',str(i2fvalues[branch]),'\n']   
          solfile.writelines(i2line)
          #i2tvalues[branch] = i2var_t[branch].x 


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

    if all_data['max_rounds'] == 1:
      log.joint(' bye!\n')
      sys.exit(0)

    ############################### CUTS ######################################

    # Cuts' computations and management
    cutplane_cuts(log,all_data)

    # Cuts' statistics
    cutplane_cutstats(log,all_data)
    
    themodel.update()

    log.joint(' model updated\n')
    log.joint('\n')

    ############################### WRITE CUTS ################################

    if all_data['writecuts']:
      write_cuts(log,all_data)

    ############################### WRITE LPS #################################

    if all_data['writelps'] and ( all_data['round'] > 0 ):
      name = 'post_cuts' + '_' + str(all_data['round']) + '.lp'
      themodel.write(name)
      log.joint(' model with new cuts written to .lp file\n')

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
    bus_id   = int(thisline[0])
    buscount = IDtoCountmap[bus_id]
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


