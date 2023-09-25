import sys
import math
from log import danoLogger
from gurobipy import *
import numpy as np
import bisect
from myutils import breakexit
import time
import math
from cuts import *
from mincut import *
from grbgraphical import *

      
def gocutplane(log, all_data):
  log.joint("Initializing AC-OPF Cutting-plane algorithm\n")

  starttime = time.time()

  buses         = all_data['buses']
  numbuses      = all_data['numbuses']
  branches      = all_data['branches']
  numbranches   = all_data['numbranches']
  gens          = all_data['gens']
  IDtoCountmap  = all_data['IDtoCountmap']
  themodel      = Model("csmodel")

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
  
  FeasibilityTol  = all_data['FeasibilityTol']
  threshold       = all_data['threshold']


  log.joint(' creating variables...\n')

  varcount = 0
  for bus in buses.values():
    maxprod = bus.Vmax*bus.Vmax
    minprod = bus.Vmin*bus.Vmin
          
    ubound = maxprod
    lbound = minprod

    Pubound, Plbound, Qubound, Qlbound = computebalbounds(log, all_data, bus)    
          

    cvar[bus]    = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "c_"+str(bus.nodeID)+","+str(bus.nodeID))
    Pinjvar[bus] = themodel.addVar(obj = 0.0, lb = Plbound, ub = Pubound, name = "IP_"+str(bus.nodeID))
    Qinjvar[bus] = themodel.addVar(obj = 0.0, lb = Qlbound, ub = Qubound, name = "IQ_"+str(bus.nodeID))
    
    varcount += 3

    for genid in bus.genidsbycount:
      gen = gens[genid]

      lower = gen.Pmin*gen.status
      upper = gen.Pmax*gen.status

      GenPvar[gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, name = "GP_"+str(gen.count)+"_"+str(gen.nodeID))

      lower = gen.Qmin*gen.status
      upper = gen.Qmax*gen.status

      if bus.nodetype == 3:
        upper = GRB.INFINITY
        lower = -GRB.INFINITY

      GenQvar[gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, name = "GQ_"+str(gen.count)+"_"+str(gen.nodeID))

      varcount += 2

      if gen.costdegree == 2 and gen.costvector[0] != 0 and (all_data['linear_objective'] or all_data['hybrid']):
        GenTvar[gen] = themodel.addVar(obj = 0.0, lb = 0.0, ub = GRB.INFINITY, name = 't_g_' + str(gen.count) + '_' + str(gen.nodeID)) 
        
        varcount += 1

  for branch in branches.values():
      f = branch.f
      t = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]
      maxprod = buses[count_of_f].Vmax*buses[count_of_t].Vmax
      minprod = buses[count_of_f].Vmin*buses[count_of_t].Vmin

      #c #check kit for bounds, might be tighter
      ubound =  maxprod
      lbound =  minprod*math.cos(branch.maxangle_rad)

      if branch.upperanglenone == 1: 
        ubound = maxprod
        lbound = 0

      if branch.loweranglenone == 1:
        ubound = maxprod
        lbound = 0                                                                         
                                          
      cvar[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "c_"+str(branch.count)+","+str(f) + "," + str(t))

      #s
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
                                  
      svar[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "s_"+str(branch.count)+","+str(f) + "," + str(t))
    
      varcount +=2
      
  themodel.update()

  for branch in branches.values():
      f = branch.f
      t = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]

      #P and Q both ways  

      ubound = branch.limit  #if constrainedflow? check kit
      lbound = -branch.limit

      Pvar_f[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "P_"+str(branch.count)+","+str(f) + "," + str(t))
      Pvar_t[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "P_"+str(branch.count)+","+str(t) + "," + str(f))
      Qvar_f[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "Q_"+str(branch.count)+","+str(f) + "," + str(t))
      Qvar_t[branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, name = "Q_"+str(branch.count)+","+str(t) + "," + str(f))

      varcount +=4 

  themodel.update()

  log.joint('   %d variables added\n' %varcount)

  #objective

  log.joint(' creating objective...\n')
  #varobjcount = 0

  #constant term
  constobjval = 0  
  for gen in gens.values():  
    if gen.status > 0:
      constobjval += gen.costvector[gen.costdegree]
  
  constvar = themodel.addVar(obj = constobjval, lb = 1.0, ub = 1.0, name = "constant")
  
  #varobjcount +=1

  #quad terms
  objvar   = themodel.addVar(obj = 1.0, lb = -GRB.INFINITY, ub = GRB.INFINITY, name = "objvar")
  qcostvar = themodel.addVar(obj = 1.0, lb = 0, ub = GRB.INFINITY, name = "qcostvar") 
  
  #varobjcount +=2

  if all_data['linear_objective'] == 0 or all_data['hybrid']:  
    objvar.setAttr("Obj",0)
    qcostexpr = QuadExpr()
    for gen in gens.values():
      if gen.costdegree == 2 and gen.costvector[0] != 0:
        qcostexpr += gen.costvector[0]*GenPvar[gen]*GenPvar[gen]
    qcost = themodel.addConstr(qcostexpr <= qcostvar, name = "qcost")
    
  #linear terms
  if all_data['linear_objective']:
    lincostvar = themodel.addVar(obj = 0.0, lb = -GRB.INFINITY, ub = GRB.INFINITY, name = "lincostvar")
    lincost    = themodel.addConstr(lincostvar <= objvar, name= "lincost")
  else:
    lincostvar = themodel.addVar(obj = 1.0, lb = -GRB.INFINITY, ub = GRB.INFINITY, name = "lincostvar")

  #varobjcount += 1

  coeff       = [gen.costvector[gen.costdegree-1] for gen in gens.values()]
  variables   = [GenPvar[gen] for gen in gens.values()]
  lincostexpr = LinExpr(coeff, variables)
  lincostdef  = themodel.addConstr(lincostexpr == lincostvar, name= "lincostdef")


  #sumTvars if using linear_objective
  if all_data['linear_objective'] or all_data['hybrid']:
    qcostvar.setAttr("Obj",0.0)
    sumTvars = LinExpr()
    for gen in gens.values():
      if gen.costdegree == 2 and gen.costvector[0] != 0:
        sumTvars += GenTvar[gen]
    sumTconstr = themodel.addConstr(sumTvars <= objvar, name='obj_var_quad')

  
  if all_data['hybrid']: #we start with full objective
    objvar.setAttr("Obj",0)
    lincostvar.setAttr("Obj",1)
    qcostvar.setAttr("Obj",1)
    themodel.remove(sumTconstr)

  themodel.update()
  
  #constraints

  log.joint(' creating constraints...\n')

  constrcount = 0
  count       = 0
  #define flow variables
  log.joint('  active power flow variables definition\n')

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    if branch.status == 0:
      log.joint(' branch ' + str(branch.count) + ' f ' + str(f) + ' t ' + str(t) + ' is OFF\n')
      breakexit('check')

    #  Gff cff + Gft cft + Bft sft
    constrname = "Pdef_"+str(branch.count)+","+str(f)+","+str(t)
    expr = LinExpr()
    expr += branch.Gff*cvar[buses[count_of_f]]
    expr += branch.Gft*cvar[branch]
    expr += branch.Bft*svar[branch]
    
    themodel.addConstr(expr == Pvar_f[branch], name = constrname)

    #  Gtt ctt + Gtf cft + Btf stf = Gtt ctt + Gtf cft - Btf sft
    constrname = "Pdef_"+str(branch.count)+","+str(t)+","+str(f)
    expr = LinExpr()
    expr += branch.Gtt*cvar[buses[count_of_t]]
    expr += branch.Gtf*cvar[branch]
    expr += -branch.Btf*svar[branch] #minus because svarft = -svartf

    themodel.addConstr(expr == Pvar_t[branch], name = constrname)
    
    constrcount += 2
    count       += 2

  log.joint('   %d active power flow variables defined\n'%count)

  log.joint('  reactive power flow variables definition\n')
  count = 0

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    constrname = "Qdef_"+str(branch.count)+","+str(f)+","+str(t)
    
    # -Bff cff - Bft cft + Gft sft
    expr = LinExpr()
    expr += -branch.Bff*cvar[buses[count_of_f]]
    expr += -branch.Bft*cvar[branch]
    expr += +branch.Gft*svar[branch]

    themodel.addConstr(expr == Qvar_f[branch], name = constrname)

    # -Btt ctt - Btf cft + Gtf stf = -Btt ctt - Btf cft - Gtf sft 
    constrname = "Qdef_"+str(branch.count)+","+str(t)+","+str(f)
    expr = LinExpr()
    expr += -branch.Btt*cvar[buses[count_of_t]]
    expr += -branch.Btf*cvar[branch]
    expr += -branch.Gtf*svar[branch] #again, same minus

    themodel.addConstr(expr == Qvar_t[branch], name = constrname)

    constrcount += 2
    count        = 2
  log.joint('   %d reactive power flow variables defined\n'%count)

  #balance constraints
  log.joint('  adding bus injection constraints...\n')
  log.joint('  active power injection constraints\n')
  count = 0

  for bus in buses.values():
    constrname = "PBaldef"+str(bus.nodeID)
    expr = LinExpr()
    for branchid in bus.frombranchids.values():
      expr += Pvar_f[ branches[branchid] ]
    for branchid in bus.tobranchids.values():
      expr += Pvar_t[ branches[branchid] ]

    if bus.Gs != 0 and ( ( len(bus.frombranchids) != 0 ) or ( len(bus.tobranchids) != 0 ) ):
      expr += bus.Gs*cvar[bus]

    themodel.addConstr(expr == Pinjvar[bus], name = constrname)

    constrcount += 2
    count       += 2
  
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
 
    if bus.Bs != 0 and ( ( len(bus.frombranchids) != 0 ) or ( len(bus.tobranchids) != 0 ) ):
      expr += (-bus.Bs)*cvar[bus]

    themodel.addConstr(expr == Qinjvar[bus], name = constrname)
    
    constrcount += 2
    count       += 2
  log.joint('   %d reactive power injection constraints added\n'%count)

  #bus-injection definitions

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
  
  log.joint('  %d constraints added\n'%constrcount)
  ################## matpower solution ##############
  
  if all_data['mp_txt']:
    getpower(log,all_data)

  if all_data['getsol']:
    if all_data['casename'] == 'case_ACTIVSg70k':
      check70(log,all_data)
    elif all_data['casename'] == 'case_ACTIVSg10k':
      check10(log,all_data)
    elif all_data['casename'] == 'case13659pegase':
      check13659(log,all_data)
    else:
      get_flows(log,all_data)
      get_va(log,all_data)

  ####################################################

  #loss inequalities

  #active loss

  if all_data['loss_ineqs']:
    log.joint(' adding and checking validity of loss inequalities (wrt a matpower solution)\n')
    log.joint(' feasibility tolerance ' + str(FeasibilityTol) + '\n')
    all_data['loss_cuts'] = {}
    counter_loss = 0
    for branch in branches.values():

      if branch.r < 0: ##!!!
        continue

      if all_data['loss_validity']:
       mp_Pf     = all_data['mp_Pfvalues'][branch]
       mp_Pt     = all_data['mp_Ptvalues'][branch]
       violation = - (mp_Pf + mp_Pt)

       if violation > FeasibilityTol:
         log.joint(' WARNING, the loss inequality associated to branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
         log.joint(' violation ' + str(violation) + '\n')
         log.joint(' values (AC solution) ' + ' Pf ' + str(mp_Pf) + ' Pt ' + str(mp_Pt) + '\n')
         breakexit('check!')
       else:
         log.joint(' AC solution satisfies loss inequality at branch ' + str(branch.count) + ' with slack ' + str(violation) + '\n')
 
      counter_loss += 1
      all_data['loss_cuts'][branch] = (0,0,threshold)

      f = branch.f
      t = branch.t
      lossexp = LinExpr()
      constrname = "loss_ineq_"+str(branch.count)+"_"+str(f)+","+str(t)
      lossexp += Pvar_f[branch] + Pvar_t[branch]
      themodel.addConstr(lossexp >= 0, name = constrname)

    all_data['num_loss_cuts']         = counter_loss
    all_data['num_loss_cuts_dropped'] = 0
    all_data['dropped_loss']          = []

  #qloss
  if all_data['qloss_cuts']:
    qloss_cuts(log,all_data)

  #jabr rotated cone inequalities
  if all_data['jabrs']:
    maxviolation = 0
    violated     = 0
    maxbranch    = -1
    maxbusf      = -1
    maxbust      = -1
    for branch in branches.values():
      f = branch.f
      t = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]

      if all_data['jabr_validity']:
        mp_c         = all_data['mp_cvalues'][branch]
        mp_s         = all_data['mp_svalues'][branch]
        mp_cbusf     = all_data['mp_cvalues'][buses[count_of_f]]
        mp_cbust     = all_data['mp_cvalues'][buses[count_of_t]]
        relviolation = violation    =  mp_c * mp_c + mp_s * mp_s - mp_cbusf * mp_cbust
        #relviolation = violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 )
        if relviolation > maxviolation:
          maxviolation = relviolation
          maxbranch    = branch.count
          maxbusf      = f
          maxbust      = t

        if relviolation > FeasibilityTol:
          violated += 1
          log.joint(' WARNING, the Jabr-inequality associated to branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
          log.joint(' violation ' + str(violation) + '\n')
          log.joint(' relative violation ' + str(relviolation) + '\n')
          log.joint(' values (AC solution) ' + ' cft ' + str(mp_c) + ' sft ' + str(mp_s) + ' cff ' + str(mp_cbusf) + ' ctt ' + str(mp_cbust) + '\n' )
          #breakexit('check!')
        else:
          log.joint(' AC solution satisfies loss inequality at branch ' + str(branch.count) + ' with slack ' + str(relviolation) + '\n')


      trigexp = QuadExpr()
      constrname = "jabr_"+str(branch.count)+"_"+str(f)+","+str(t)
      trigexp += cvar[branch]*cvar[branch] + svar[branch]*svar[branch] - cvar[buses[count_of_f]]*cvar[buses[count_of_t]]
      themodel.addConstr(trigexp <= 0, name = constrname)

    if all_data['jabr_validity']:
      log.joint(' max violation of Jabr-inequalities by AC solution ' + str(maxviolation) + ' at branch ' + str(maxbranch) + ' f ' + str(maxbusf) + ' t ' + str(maxbust) + '\n')
      log.joint(' number of violated Jabr-inequalities ' + str(violated) + '\n')
      breakexit(' check Jabr violation')

  #i2-envelope cuts 
  if all_data['sqd_current_cuts']:
    i2var_f = {}
    #i2var_t = {}

    for branch in branches.values():
      f          = branch.f
      t          = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]

      #upperbound_f = upperbound_t = GRB.INFINITY
      upperbound_f = branch.limit**2 / (buses[count_of_f].Vmin*buses[count_of_f].Vmin)
      #upperbound_t = branch.limit / (buses[count_of_t].Vmin*buses[count_of_t].Vmin)

      i2var_f[branch] = themodel.addVar(obj = 0.0, lb = 0, ub = upperbound_f , name = "i2_"+str(branch.count)+","+str(f) + "," + str(t) )
      expr_f = LinExpr()

      #i2var_t[branch] = themodel.addVar(obj = 0.0, lb = 0, ub = upperbound_t, name = "i2_"+str(branch.count)+","+str(t) + "," + str(f) )
      #expr_t = LinExpr()

      ratio  = branch.ratio
      y      = branch.y
      g      = y.real
      b      = y.imag
      bshunt = branch.bc
      angle = branch.angle_rad
      
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

      # if (ratio != 1 and ratio != 0):
      #   angle = branch.angle_rad
      #   #first f
      #   expr_f += (g*g + b*b)/(ratio*ratio) * ( (cvar[buses[count_of_f]]/(ratio*ratio)) + cvar[buses[count_of_t]] - (2/ratio) * ( cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) ) )
      #   expr_f += b*bshunt/(ratio**3) * ( (cvar[buses[count_of_f]]/ratio) - (cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) )) #since gsh is 0
      #   expr_f += g*bshunt/(ratio**3) * ( svar[branch] * math.cos(angle) - cvar[branch] * math.sin(angle) )
      #   expr_f += (bshunt*bshunt*cvar[buses[count_of_f]]/(4*(ratio**4)) )

      #   #now t
      #   #expr_t += (g*g + b*b) * ( cvar[buses[count_of_t]] + cvar[buses[count_of_f]]/(ratio*ratio) - (2/ratio) * ( cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) ) )
      #   #expr_t += b*bshunt * ( cvar[buses[count_of_t]] - (1/ratio) * (cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) ))
      #   #expr_t += (g*bshunt/ratio) * ( svar[branch] * math.cos(angle) - cvar[branch] * math.sin(angle) )
      #   #expr_t += ( bshunt*bshunt*cvar[buses[count_of_t]]/4 )
      # else:
      #   expr_f += (g*g + b*b) * ( cvar[buses[count_of_f]] + cvar[buses[count_of_t]] - 2 * cvar[branch] ) + bshunt * ( cvar[buses[count_of_f]] * ( (bshunt/4) + b ) - b * cvar[branch] + g * svar[branch])
      #   #now t
      #   #expr_t += (g*g + b*b) * ( cvar[buses[count_of_f]] + cvar[buses[count_of_t]] - 2 * cvar[branch] ) + bshunt * ( cvar[buses[count_of_t]] * ( (bshunt/4) + b ) - b * cvar[branch] - g * svar[branch]) #Matias, changed - g * svar[branch] to + g * svar[branch]

      constrname_f = 'i2def_'+str(branch.count)+","+str(f) + "," + str(t)
      themodel.addConstr(expr_f == i2var_f[branch],name = constrname_f) 

      constrname_t = 'i2def_'+str(branch.count)+","+str(t) + "," + str(f)
      #themodel.addConstr(expr_t == i2var_t[branch],name = constrname_t) 

  #i2 rotated cone inequalities | notice that we are not adding i2var_tf ... 
  if all_data['i2']:
    maxviolation = 0
    violated     = 0
    maxbranch    = -1
    maxbusf      = -1
    maxbust      = -1
    maxi2f       = -1
    maxcff       = -1
    maxPf        = -1
    maxQf        = -1

    for branch in branches.values():
      f = branch.f
      t = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]

      #
      if all_data['i2ineq_validity']:
        mp_Pf        = all_data['mp_Pfvalues'][branch]
        mp_Qf        = all_data['mp_Qfvalues'][branch]
        mp_c         = all_data['mp_cvalues'][branch]
        mp_s         = all_data['mp_svalues'][branch]
        mp_cbusf     = all_data['mp_cvalues'][buses[count_of_f]]
        mp_cbust     = all_data['mp_cvalues'][buses[count_of_t]]
        mp_i2f       = computei2value(log,all_data,branch,mp_c,mp_s,mp_cbusf,mp_cbust)
        relviolation = violation    = mp_Pf * mp_Pf + mp_Qf * mp_Qf - mp_cbusf * mp_i2f
        #relviolation = violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 )
        if relviolation > maxviolation:
          maxviolation = relviolation
          maxbranch    = branch.count
          maxbusf      = f
          maxbust      = t
          maxi2f       = mp_i2f
          maxcff       = mp_cbusf
          maxPf        = mp_Pf
          maxQf        = mp_Qf


        if relviolation > FeasibilityTol:
          violated += 1
          log.joint(' WARNING, the i2 inequality associated to branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
          log.joint(' violation ' + str(violation) + '\n')
          log.joint(' relative violation ' + str(relviolation) + '\n')
          log.joint(' values (AC solution) ' + ' Pft ' + str(mp_Pf) + ' Qft ' + str(mp_Qf) + ' cff ' + str(mp_cbusf) + ' i2ft ' + str(mp_i2f) + '\n' )
          #breakexit('check!')
        else:
          log.joint(' AC solution satisfies loss inequality at branch ' + str(branch.count) + ' with slack ' + str(relviolation) + '\n')

        #if branch.count == 14080: #knitro sol off
        #  breakexit("check line")
      trigexp = QuadExpr()
      constrname = "i2_"+str(branch.count)+"_"+str(f)+","+str(t)
      trigexp += Pvar_f[branch]**2 + Qvar_f[branch]**2 - cvar[buses[count_of_f]] * i2var_f[branch]
      themodel.addConstr(trigexp <= 0, name = constrname)

    log.joint(' max violation of i2-inequalities by AC solution ' + str(maxviolation) + ' at branch ' + str(maxbranch) + ' f ' + str(maxbusf) + ' t ' + str(maxbust) + '\n')
    log.joint(' values (AC solution) ' + ' Pft ' + str(maxPf) + ' Qft ' + str(maxQf) + ' cff ' + str(maxcff) + ' i2ft ' + str(maxi2f) + '\n' )
    if maxbranch == 14080:
      branch    = branches[maxbranch]
      y         = branch.y
      g         = y.real
      b         = y.imag
      busf      = buses[IDtoCountmap[branch.f]]
      bust      = buses[IDtoCountmap[branch.t]]
      againi2ft = (g**2 + b**2) * ( all_data['mp_cvalues'][busf] + all_data['mp_cvalues'][bust] - 2 * all_data['mp_cvalues'][branch] )
    log.joint(' again i2f ' + str(againi2ft) + '\n')
    log.joint(' number of violated i2-inequalities ' + str(violated) + '\n')
    breakexit(' check i2 violation')
  #limit constraints
  if all_data['limits']:
    for branch in branches.values():
      if branch.constrainedflow:
        f = branch.f
        t = branch.t
        constrname = "limit_f_"+str(branch.count)+","+str(f)+","+str(t)
        limexp = QuadExpr()
        limexp += Pvar_f[branch]*Pvar_f[branch] + Qvar_f[branch]*Qvar_f[branch]
        themodel.addConstr(limexp <= branch.limit**2, name = constrname)

        constrname = "limit_t_"+str(branch.count)+","+str(t)+","+str(f)
        limexp = QuadExpr()
        limexp += Pvar_t[branch]*Pvar_t[branch] + Qvar_t[branch]*Qvar_t[branch]
        themodel.addConstr(limexp <= branch.limit**2, name = constrname)


  #reactive loss inequalities -> implied by pkm + pmk >= 0 | easy fact
      
  themodel.update()

  endtime = time.time()

  log.joint(' formulation time: %g\n' %(endtime - starttime))

  log.joint(' we are not writing down the base LP model\n')
  #log.joint(" writing to lpfile " + all_data['lpfilename'] + "\n")  
  #themodel.write(all_data['lpfilename'])

  #data
  all_data['themodel']    = themodel
  all_data['cvar']        = cvar
  all_data['svar']        = svar
  all_data['GenPvar']     = GenPvar
  all_data['GenTvar']     = GenTvar
  all_data['Pvar_t']      = Pvar_t
  all_data['Qvar_f']      = Qvar_f
  all_data['Qvar_t']      = Qvar_t

  if all_data['sqd_current_cuts']:
    all_data['Pvar_f']  = Pvar_f
    all_data['Qvar_f']  = Qvar_f
    all_data['i2var_f'] = i2var_f
    all_data['i2var_t'] = i2var_f

  dicGenPvalues = {}  #check if we are using these (objective cuts)
  for gen in gens.values():
    dicGenPvalues[gen] = []
  all_data['dicGenPvalues'] = dicGenPvalues

  jabr_cuts_info         = all_data['jabr_cuts_info']  #these were created in main.py
  jabr_cuts_info_updated = all_data['jabr_cuts_info_updated']

  for branch in branches.values():
    jabr_cuts_info[branch]         = {}
    jabr_cuts_info_updated[branch] = {}

  if all_data['sqd_current_cuts']:
    i2_cuts_info         = all_data['i2_cuts_info']
    i2_cuts_info_updated = all_data['i2_cuts_info_updated']

    for branch in branches.values():
      i2_cuts_info[branch]         = {}
      i2_cuts_info_updated[branch] = {}

  if all_data['lim_cuts']:
    limit_cuts_info         = all_data['limit_cuts_info']
    limit_cuts_info_updated = all_data['limit_cuts_info_updated']

    for branch in branches.values():
      limit_cuts_info[branch]         = {}
      limit_cuts_info_updated[branch] = {}

  if all_data['mincut']:
    all_data['size_supernodes'] = math.floor(0.1 * numbuses) 
    all_data['new_mincut']      = 0
    all_data['mincut_cuts']     = {}
    all_data['mincut_allcuts']  = []
    all_data['num_mincut_cuts'] = 0

  ################## fixing a solution ##############

  if all_data['fixflows']:
    fixflows(log,all_data)
    if all_data['fixCS'] == 0:
      return None

  if all_data['fixCS']:
    fixCS(log,all_data)
    return None

  #if all_data['fixCS_knitro']:
  if all_data['knitro_sol']:
     getsol_knitro(log,all_data)
  
  #   fixCS_knitro(log,all_data)
  #   fixflows_knitro(log,all_data)
  #   return None

  ################## param setting before loop ##############

  themodel.Params.Method = all_data['solver_method']
  themodel.Params.Crossover = all_data['crossover'] 
  if all_data['solver_method'] == 2:
    themodel.Params.BarHomogeneous = 1
    themodel.Params.BarConvTol = 1e-06

  themodel.Params.NumericFocus = 1
  themodel.Params.OutPutFlag = 1
  themodel.Params.Logfile = 'newlogs/' + all_data['casename'] + '_gurobi.log'


  ################ writing solutions to a file  ##############
    
  namesolfile = 'cutplane_sols/cutplane_sol_' + all_data['casename'] + '.txt'
  solfile     = open(namesolfile,'a')

  ####################### main loop ##########################

  all_data['round']                  = 0
  all_data['runtime']                = time.time() - all_data['T0']
  all_data['T0_cutplane']            = time.time()
  all_data['cumulative_solver_time'] = 0
  gap                                = 1e20
  all_data['ftol_counter']           = 0
  oldobj                             = 1

  mc = 0

  if all_data['addcuts']:
    if add_cuts(log,all_data):
      themodel.update()

    themodel.write(all_data['casename']+'_precomputed_cuts.lp')
    log.joint(' pre-computed cuts added!\n') #fix renaming of new cuts 


  while (all_data['round'] <= all_data['max_rounds']) and (all_data['runtime'] <= all_data['max_time']) and (all_data['ftol_counter'] <= all_data['ftol_iterates']):
    
    if all_data['knitro_sol'] and (all_data['round'] >= 30):
      getsol_knitro(log,all_data)
      fixCS_knitro(log,all_data)
      fixflows_knitro(log,all_data)
      return None

    all_data['round'] += 1

    ########## HYBRID ############
    if all_data['hybrid']:
      QP                  = 50
      QP_to_LP            = 3
      no_objective_cuts   = 11

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
          qcost = themodel.addConstr(qcostexpr <= qcostvar, name = "qcost")

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
          sumTconstr = themodel.addConstr(sumTvars <= objvar, name= 'objvar_quad')
          lincost    = themodel.addConstr(lincostvar <= objvar, name= 'objvar_linear')
        if all_data['round'] > no_objective_cuts:
          all_data['objective_cuts'] = 0


      themodel.update()
      hybridlpname = 'hybrid' + all_data['casename'] + str(all_data['round']) + '.lp'
      log.joint(' writing down .lp file for hybrid ...\n')
      themodel.write(hybridlpname)
      
      #breakexit('check.lp')
    
    #################### solving model #########################

    log.joint(' solving model with method ' + str(themodel.params.method) + '\n')
    t0_solve = time.time()
    themodel.optimize()
    t1_solve = time.time()

    if themodel.status == GRB.status.INF_OR_UNBD:
      log.joint('->LP infeasible or unbounded\n')
      log.joint(' solver runtime = %g\n' % (t1_solve - t0_solve) )
      return themodel.status, 0

    if themodel.status == GRB.status.INFEASIBLE:
      log.joint('->LP infeasible\n')
      log.joint(' solver runtime = %g\n' % (t1_solve - t0_solve) )
      return themodel.status, 0

    if themodel.status == GRB.status.UNBOUNDED:
      log.joint('->LP unbounded\n')
      log.joint(' solver runtime = %g\n' % (t1_solve - t0_solve) )
      return themodel.status, 0


    all_data['cumulative_solver_time'] += (t1_solve - t0_solve)  

    ############################################################

    solfile     = open(namesolfile,'a')
    solfile.write('round' + str(all_data['round']) + '\n' )
    solfile.write('obj ' + str(themodel.ObjVal) + '\n')

    #### values ####
    log.joint(' storing current solution values ...\n')

    newobj        = themodel.ObjVal
    
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
    #Pinjvalues    = {}
    #Qinjvalues    = {}

    dicGenPvalues = all_data['dicGenPvalues']

    for bus in buses.values():
      cvalues[bus] = cvar[bus].x
      varname      = cvar[bus].varname + ' = '
      lines        = [varname,str(cvalues[bus]),'\n']
      solfile.writelines(lines)
      #Pinjvalues[bus] = Pinjvar[bus].x
      #Qinjvalues[bus] = Qinjvar[bus].x
      
    for branch in branches.values():
      f                   = branch.f
      t                   = branch.t
      count_of_f          = IDtoCountmap[f]
      count_of_t          = IDtoCountmap[t]
      bus_f               = buses[count_of_f]
      bus_t               = buses[count_of_t]
      bc                  = branch.bc

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

      if all_data['sqd_current_cuts']:
        i2fvalues[branch] = i2var_f[branch].x
        i2line            = [i2var_f[branch].varname + ' = ',str(i2fvalues[branch]),'\n']   
        solfile.writelines(i2line)
        #i2tvalues[branch] = i2var_t[branch].x #are we actually using this variable?


    for gen in gens.values():
      GenPvalues[gen] = GenPvar[gen].x
      GenQvalues[gen] = GenQvar[gen].x
      GenPvarname     = GenPvar[gen].varname + ' = '
      GenQvarname     = GenQvar[gen].varname + ' = '
      GenPlines       = [GenPvarname,str(GenPvalues[gen]),'\n']
      GenQlines       = [GenQvarname,str(GenQvalues[gen]),'\n']
      
      solfile.writelines(GenPlines)
      solfile.writelines(GenQlines)

      if gen.costdegree == 2 and gen.costvector[0] != 0:
        gen_values = dicGenPvalues[gen]
        bisect.insort(gen_values,GenPvar[gen].x)

    solfile.close()
        
    all_data['cvalues']                 = cvalues
    all_data['svalues']                 = svalues
    all_data['GenPvalues']              = GenPvalues
    all_data['GenQvalues']              = GenQvalues
    all_data['dicGenPvalues']           = dicGenPvalues
    all_data['plossvalues']             = plossvalues
    all_data['qlossvalues']             = qlossvalues
    all_data['total_active_gen']        = sum(GenPvalues.values())
    all_data['total_active_losses']     = sum(plossvalues.values())
    all_data['total_reactive_gen']      = sum(GenQvalues.values())
    all_data['total_reactive_losses']   = sum(qlossvalues.values())
    all_data['total_reactive_gains']    = sum(qgains.values())
    #all_data['Pinjvalues']              = Pinjvalues
    #all_data['Qinjvalues']              = Qinjvalues
    all_data['Pfvalues']                = Pfvalues
    all_data['Ptvalues']                = Ptvalues
    all_data['Qfvalues']                = Qfvalues
    all_data['Qtvalues']                = Qtvalues

    if all_data['sqd_current_cuts']:
      all_data['i2fvalues'] = i2fvalues
      #all_data['i2tvalues'] = i2tvalues

    log.joint(' done storing values\n')

    ### adding cuts if there are ... 
    # if all_data['addcuts']:
    #   if add_cut(log,all_data):
    #     themodel.update()
  
    #   all_data['addcuts'] = 0
    #   themodel.write('new_lp.lp')
    #   #breakexit('check new lp')
    #   log.joint(' we skip this iteration and resolve with new cuts') #fix renaming of new cuts 
    #   continue


    #### CUTS ####
    
    log.joint(' adding cuts ...\n')

    jabr_cuts(log,all_data)
    
    if all_data['NO_jabrs_violated']:
      if all_data['threshold'] > all_data['tolerance']:
        all_data['threshold'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold']) + '\n' )
        #continue
      else:
        log.joint(' threshold below ' + str(all_data['tolerance']) + ', we are done\n' )
        if all_data['hybrid'] == 0:
          log.joint(' writing to lpfile ' + all_data['lpfilename_cuts'] + '\n')  
          themodel.write(all_data['lpfilename_cuts'])
          log.joint(' total runtime = ' + str(time.time() - all_data['T0']) )
          log.joint(' bye.\n')
          return None

    if all_data['sqd_current_cuts']:
      i2_cuts(log,all_data)

      if all_data['NO_i2_cuts_violated']:
        if all_data['threshold_i2'] > all_data['tolerance']:
          all_data['threshold_i2'] *= 1e-01
          log.joint(' threshold updated to ' + str(all_data['threshold_i2']) + '\n' )
          #continue
        else:
          log.joint(' threshold below ' + str(all_data['tolerance']) + ', we are done\n' )
          if all_data['hybrid'] == 0:
            log.joint(' writing to lpfile ' + all_data['lpfilename_cuts'] + '\n')  
            themodel.write(all_data['lpfilename_cuts'])
            log.joint(' total runtime = ' + str(time.time() - all_data['T0']) )
            log.joint(' bye.\n')
            return None


    if all_data['lim_cuts']:
      limit_cuts(log,all_data)
      if all_data['NO_limit_cuts_violated']:
        if all_data['threshold_limit'] > all_data['tolerance']:
          all_data['threshold_limit'] *= 1e-01
          log.joint(' threshold updated to ' + str(all_data['threshold_limit']) + '\n' )
          #continue
        else:
          log.joint(' threshold below ' + str(all_data['tolerance']) + ', we continue\n' )

    if all_data['loss_inequalities']:
      loss_cuts(log,all_data)
      
    if all_data['linear_objective']:
      if all_data['objective_cuts']:
        objective_cuts(log,all_data)

    if all_data['sqd_current_cuts'] and all_data['dropi2'] and (all_data['round'] >= all_data['cut_age_limit']):
      drop_i2(log,all_data)

    if all_data['lim_cuts'] and all_data['droplimit'] and (all_data['round'] >= all_data['cut_age_limit']):
      drop_limit(log,all_data)
      
    #if all_data['loss_ineqs'] and all_data['droploss'] and (all_data['round'] >= all_data['cut_age_limit']):
    #  drop_loss(log,all_data)

    if all_data['loss_inequalities'] and all_data['droploss'] and (all_data['round'] >= all_data['cut_age_limit']):
      drop_loss(log,all_data)
      
    if all_data['dropjabrs'] and (all_data['round'] >= all_data['cut_age_limit']):
      drop_jabr(log,all_data)

    if all_data['cut_analysis']:
      cut_analysis(log,all_data)
    
    log.joint('\n')

    log.joint(' ******************************************\n')
    
    log.joint(' casename = %s\n' % all_data['casename'] )
    log.joint(' round = %g\n' % all_data['round'] )
    log.joint(' objective = %g\n' % themodel.ObjVal)
    if type(all_data['primal_bound']) == float:
      primal_bound = all_data['primal_bound']
      gap = 100 * ((primal_bound - themodel.ObjVal) / primal_bound)
      log.joint(' primal_bound = %g\n' % primal_bound)
      log.joint(' gap (percent) = %g\n' % gap)
    log.joint(' solver status = ' + str(themodel.status) + ' solver method = ' + str(themodel.params.method) + '\n')
    log.joint(' threshold relative improvement obj = %g\n' % all_data['ftol']) 
    log.joint(' -- active power --\n')
    log.joint(' total active power generation = %g\n' % all_data['total_active_gen'] )
    log.joint(' active power losses = %g\n' % all_data['total_active_losses'] )
    if all_data['mp_txt']:
      mp_ploss = all_data['mp_ploss']
      mp_sum_plosses = (1/100)*sum(mp_ploss.values())
      log.joint(' primal bound active power losses = %g\n' % mp_sum_plosses )
      ploss_gap = 100*(mp_sum_plosses - all_data['total_active_losses'])/mp_sum_plosses
      log.joint(' active loss gap (percent) = %g\n' % ploss_gap)
    log.joint(' -- reactive power --\n')
    log.joint(' total reactive power generation = %g\n' % all_data['total_reactive_gen'] )
    log.joint(' reactive power net losses = %g\n' % all_data['total_reactive_losses'] )
    log.joint(' reactive power gains = %g\n' % - all_data['total_reactive_gains'] )
    if all_data['mp_txt']:
      mp_qloss          = all_data['mp_qloss']
      mp_qfchg          = all_data['mp_qfchg']
      mp_qtchg          = all_data['mp_qtchg']
      mp_sum_qlosses    = (1/100)*sum(mp_qloss.values())
      mp_sum_qgains     = (1/100)*(sum(mp_qfchg.values()) + sum(mp_qtchg.values()))
      mp_sum_qnetlosses = mp_sum_qlosses - mp_sum_qgains
      log.joint(' primal bound reactive power net losses = %g\n' % mp_sum_qnetlosses )
      #log.joint(' primal bound reactive power losses = %g\n' % mp_sum_qlosses )
      log.joint(' primal bound reactive power gains = %g\n' % mp_sum_qgains )
      qloss_gap = 100*(mp_sum_qnetlosses - all_data['total_reactive_losses'])/mp_sum_qnetlosses
      log.joint(' net reactive loss gap (percent) = %g\n' % qloss_gap)

    log.joint(' -- Jabr-envelope cuts --\n')
    log.joint(' current number of Jabr-envelope cuts = %g\n' % all_data['num_jabr_cuts'])
    log.joint(' number of Jabr-envelope cuts added in current round = %g\n' % all_data['num_jabr_cuts_added'] )
    log.joint(' top percent of most violated Jabr-envelope cuts added = %g\n' % (100*all_data['most_violated_fraction']) )
    log.joint(' max error (Jabrs) in current round = %g\n' % all_data['max_error'] )
    log.joint(' initial Jabr-envelope threshold = %g\n' % all_data['initial_threshold'] )
    log.joint(' current Jabr-envelope cuts threshold = %g\n' % all_data['threshold'] )
    log.joint(' Jabr-envelope cuts tolerance = %g\n' % all_data['tolerance'] )
    log.joint(' parallel-cuts threshold = %g\n' % all_data['threshold_dotprod'] )
    if all_data['dropjabrs']:
      log.joint(' number of Jabr-envelope cuts dropped in current round = %g\n' % all_data['num_jabr_cuts_dropped'] )
      log.joint(' cut age limit = %g\n' % all_data['cut_age_limit'])

    if all_data['sqd_current_cuts']:
      log.joint(' -- i2-envelope cuts --\n')
      log.joint(' current number of i2-envelope cuts = %g\n' % all_data['num_i2_cuts'])
      log.joint(' number of i2-envelope cuts added in current round = %g\n' % all_data['num_i2_cuts_added'] )
      log.joint(' top percent of most violated i2-envelope cuts added = %g\n' % (100*all_data['most_violated_fraction_i2']) )
      log.joint(' max error (i2) in current round = %g\n' % all_data['max_error_i2'] )
      log.joint(' initial i2-envelope threshold = %g\n' % all_data['threshold_i2'] )
      if all_data['dropi2']:
        log.joint(' number of i2-envelope cuts dropped in current round = %g\n' % all_data['num_i2_cuts_dropped'] )

    if all_data['lim_cuts']:
      log.joint(' -- Limit-envelope cuts --\n')
      log.joint(' current number of limit-envelope cuts = %g\n' % all_data['num_limit_cuts'])
      log.joint(' number of limit-envelope cuts added in current round = %g\n' % all_data['num_limit_cuts_added'] )
      log.joint(' top percent of most violated limit-envelope cuts added = %g\n' % (100*all_data['most_violated_fraction_limit']) )
      log.joint(' max error (limit) in current round = %g\n' % all_data['max_error_limit'] )
      log.joint(' initial limit-envelope threshold = %g\n' % all_data['threshold_limit'] )
      if all_data['droplimit']:
        log.joint(' number of limit-envelope cuts dropped in current round = %g\n' % all_data['num_limit_cuts_dropped'] )

    if all_data['loss_inequalities']:
      log.joint(' -- loss inequalities --\n')
      log.joint(' current number of loss inequalities = %g\n' % all_data['num_loss_cuts'])
      log.joint(' number of loss inequalities added in current round = %g\n' % all_data['num_loss_cuts_added'] )
      log.joint(' top percent of most violated loss inequalities added = %g\n' % (100*all_data['most_violated_fraction_loss'])
 )
      log.joint(' max error (loss) in current round = %g\n' % all_data['max_error_loss'] )
    if all_data['droploss']:
      log.joint(' number of loss ineqs dropped in current round = %g\n' % all_data['num_loss_cuts_dropped'] )
    if all_data['mincut']:
      log.joint(' current N = ' + str(all_data['size_supernodes']) + '\n')
    if all_data['linear_objective'] or all_data['hybrid']:
      log.joint(' -- objective cuts --\n')
      log.joint(' current number of objective-cuts = %g\n' % all_data['num_objective_cuts'])
      log.joint(' objective-cuts threshold = %g\n' % all_data['threshold_objcuts'])    
    log.joint(' -- runtimes --\n')
    log.joint(' solver runtime current round = ' + str((t1_solve - t0_solve)) + '\n')
    log.joint(' cumulative solver time = %g\n' % all_data['cumulative_solver_time']) 
    log.joint(' time so far (cutplane) = %g\n' % (time.time() - all_data['T0_cutplane']) )
    log.joint(' time so far (overall) = %g\n' % (time.time() - all_data['T0']) )
    log.joint(' max running time = %g\n' % all_data['max_time'] )
    if all_data['ftol_counter']:
      log.joint(' relative improvement in obj less than %g\n' % all_data['ftol']) 
    log.joint(' ******************************************\n')
    
    themodel.update()

    log.joint('\n')
    log.joint(' model updated with new Jabr,i2, and limit-cuts and dropped slack-and-old cuts (if any)\n' )
    log.joint('\n')

    ### writing down cuts

    if all_data['writecuts']:
      write_cuts(log,all_data)

    #### max-flow #####
    IT = 1

    if all_data['mincut'] and (all_data['round'] >= IT): #>=
      mc = 1
      supernodes(log,all_data)
      mincut(log,all_data)

      cut = all_data['mincut_cuts'][1]

      log.joint(' adding loss ineqs associated to mincut in the network\n')
      if all_data['loss_ineqs'] == 0:
        breakexit('ch')
        count_loss = 0
        #all_data['loss_cuts'] = {}
        for branchid in cut:
          count_loss += 1
          branch = branches[branchid]
          #all_data['loss_cuts'][count_loss] = branch
          f = branch.f
          t = branch.t
          lossexp = LinExpr()
          constrname = "loss_ineq_"+str(branchid)+"_"+str(f)+","+str(t)
          lossexp += Pvar_f[branch] + Pvar_t[branch]
          themodel.addConstr(lossexp >= 0, name = constrname)
        all_data['num_loss_cuts']         = count_loss
        all_data['num_loss_cuts_dropped'] = 0
        all_data['dropped_loss']          = []
        #all_data['loss_cuts'] = []
      
        themodel.update()

      if all_data['dographics']:
        textlist = []
        #plottype = all_data['plottype'] = 'active_losses'
        plottype = all_data['plottype'] = 'losses'
        grbgraphical(log,all_data, plottype, textlist)

    #### obbt #####
    if all_data['obbt'] and mc and all_data['round'] >= IT:
        obbt(log,all_data)

        #mincutjabrs(log,all_data)
        #mincuti2(log,all_data)

        expr = LinExpr()
        expr += constvar
        expr += objvar
        if all_data['objective_cuts'] == 0:
          expr += qcost
        themodel.setObjective(expr)
        themodel.update()
        themodel.write('newmodel.lp')

    if all_data['loss_analysis']:
      loss_analysis(log,all_data)

    if all_data['writelps'] and all_data['round'] > 0:
      name = 'post_cuts'+'_'+str(all_data['round'])+'.lp'
      themodel.write(name)


    if ((newobj - oldobj)/oldobj) < all_data['ftol']:
      all_data['ftol_counter'] += 1
    else:
      all_data['ftol_counter'] = 0

    oldobj = newobj
    all_data['runtime'] = time.time() - all_data['T0']

    #volts2file(log,all_data)
    #log.joint(' writing down voltages to a file')


def volts2file(log,all_data):

  buses   = all_data['buses'] 
  cvalues = all_data['cvalues']

  casefilename = all_data['casefilename']
  casename     = casefilename.split('/')[2].split('.')[0]
  filename     = 'voltsol_'+ casename +'.txt'

  log.joint(" reading file with volts sol " + filename + "\n")
  
  volts   = open(filename,"w+")

  volts.write('buscount vm\n')
  for bus in buses.values():
    volts.write(str(bus.count) + ' ' + str(cvalues[bus]**0.5) + '\n')

  volts.close()

def fixflows(log,all_data):
  tolerance   = all_data['tol_fix']
  themodel    = all_data['themodel']
  branches    = all_data['branches']
  Pvar_f      = all_data['Pvar_f']
  Pvar_t      = all_data['Pvar_t']
  Qvar_f      = all_data['Qvar_f']
  Qvar_t      = all_data['Qvar_t']

  mp_Pfvalues = all_data['mp_Pfvalues']
  mp_Ptvalues = all_data['mp_Ptvalues']
  mp_Qfvalues = all_data['mp_Qfvalues']
  mp_Qtvalues = all_data['mp_Qtvalues']

  for branch in branches.values():
        
    mp_Pf = mp_Pfvalues[branch]
    mp_Pt = mp_Ptvalues[branch]
    mp_Qf = mp_Qfvalues[branch]
    mp_Qt = mp_Qtvalues[branch]

    #Pf
    ubound_Pf = mp_Pf + tolerance
    lbound_Pf = mp_Pf - tolerance

    Pvar_f[branch].setAttr("ub",ubound_Pf)
    Pvar_f[branch].setAttr("lb",lbound_Pf)

    #Pt
    ubound_Pt = mp_Pt + tolerance
    lbound_Pt = mp_Pt - tolerance

    Pvar_t[branch].setAttr("ub",ubound_Pt)
    Pvar_t[branch].setAttr("lb",lbound_Pt)

    #Qf
    ubound_Qf = mp_Qf + tolerance
    lbound_Qf = mp_Qf - tolerance

    Qvar_f[branch].setAttr("ub",ubound_Qf)
    Qvar_f[branch].setAttr("lb",lbound_Qf)

    #Qt
    ubound_Qt = mp_Qt + tolerance
    lbound_Qt = mp_Qt - tolerance

    Qvar_t[branch].setAttr("ub",ubound_Qt)
    Qvar_t[branch].setAttr("lb",lbound_Qt)

  themodel.update()
  themodel.write('fixflows.lp')
  log.joint('check fixflows.lp\n')  


def fixCS(log,all_data):
  tolerance  = all_data['tol_fix']
  themodel   = all_data['themodel']
  buses      = all_data['buses']
  branches   = all_data['branches']
  cvar       = all_data['cvar']
  svar       = all_data['svar']
  mp_cvalues = all_data['mp_cvalues']
  mp_svalues = all_data['mp_svalues']

  for bus in buses.values():

    mp_v2 = all_data['mp_cvalues'][bus]

    ubound_v2 = mp_v2 + tolerance
    lbound_v2 = mp_v2 - tolerance

    cvar[bus].setAttr("ub",ubound_v2)
    cvar[bus].setAttr("lb",lbound_v2)

  for branch in branches.values():

    mp_c = all_data['mp_cvalues'][branch]
    mp_s = all_data['mp_svalues'][branch]

    #c
    ubound_c = mp_c + tolerance
    lbound_c = mp_c - tolerance

    cvar[branch].setAttr("ub",ubound_c)
    cvar[branch].setAttr("lb",lbound_c)

    #s
    ubound_s = mp_s + tolerance
    lbound_s = mp_s - tolerance

    svar[branch].setAttr("ub",ubound_s)
    svar[branch].setAttr("lb",lbound_s)

  themodel.update()
  themodel.write('fixCS.lp')
  log.joint('check fixCS.lp\n')

def fixflows_knitro(log,all_data):
  tolerance   = all_data['tol_fix']
  themodel    = all_data['themodel']
  branches    = all_data['branches']
  Pvar_f      = all_data['Pvar_f']
  Pvar_t      = all_data['Pvar_t']
  Qvar_f      = all_data['Qvar_f']
  Qvar_t      = all_data['Qvar_t']

  knitro_Pfvalues = all_data['knitro_Pfvalues']
  knitro_Ptvalues = all_data['knitro_Ptvalues']
  knitro_Qfvalues = all_data['knitro_Qfvalues']
  knitro_Qtvalues = all_data['knitro_Qtvalues']

  for branch in branches.values():
        
    knitro_Pf = knitro_Pfvalues[branch]
    knitro_Pt = knitro_Ptvalues[branch]
    knitro_Qf = knitro_Qfvalues[branch]
    knitro_Qt = knitro_Qtvalues[branch]

    #Pf
    ubound_Pf = knitro_Pf + tolerance
    lbound_Pf = knitro_Pf - tolerance

    Pvar_f[branch].setAttr("ub",ubound_Pf)
    Pvar_f[branch].setAttr("lb",lbound_Pf)

    #Pt
    ubound_Pt = knitro_Pt + tolerance
    lbound_Pt = knitro_Pt - tolerance

    Pvar_t[branch].setAttr("ub",ubound_Pt)
    Pvar_t[branch].setAttr("lb",lbound_Pt)

    #Qf
    ubound_Qf = knitro_Qf + tolerance
    lbound_Qf = knitro_Qf - tolerance

    Qvar_f[branch].setAttr("ub",ubound_Qf)
    Qvar_f[branch].setAttr("lb",lbound_Qf)

    #Qt
    ubound_Qt = knitro_Qt + tolerance
    lbound_Qt = knitro_Qt - tolerance

    Qvar_t[branch].setAttr("ub",ubound_Qt)
    Qvar_t[branch].setAttr("lb",lbound_Qt)

  themodel.update()
  fixmodelname = 'fixflows_knitro.lp'
  themodel.write(fixmodelname)
  log.joint('check fixed model\n')  


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
  
  cvalues = all_data['cvalues']
  svalues = all_data['svalues']
  branches = all_data['branches']

  for branch in branches.values():
    s = svalues[branch]
    c = cvalues[branch]

    log.joint(' branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' angle ' + str(math.atan(s/c)) + '\n')

def check70(log,all_data):

  try:
    f = open('mp_sols/case_ACTIVSg70k.out', "r")
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
    #log.joint(' the line ' + str(thisline) + '\n')
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

  
  log.joint(' writing solution to ../mp_sols/SOL.txt\n')
  for bus in buses.values():
    f = bus.nodeID
    c = mp_cvalues[bus]
    #log.joint('c_' + str(f) + ',' + str(f) + ' = ' + str(c) + '\n')
    afile.write('c_' + str(f) + ',' + str(f) + ' = ' + str(c) + '\n')
    afile_kit.write('c_' + str(f) + '_' + str(f) + ' = ' + str(c) + '\n')

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

    #log.joint('c_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(c) + '\n')
    afile.write('c_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(c) + '\n')
    afile_kit.write('c_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(c) + '\n')

    #log.joint('s_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(s) + '\n')
    afile.write('s_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(s) + '\n')
    afile_kit.write('s_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(s) + '\n')


  afile.close()
  afile_kit.close()


  all_data['mp_Pfvalues'] = mp_Pfvalues
  all_data['mp_Ptvalues'] = mp_Ptvalues
  all_data['mp_Qfvalues'] = mp_Qfvalues
  all_data['mp_Qtvalues'] = mp_Qtvalues

  all_data['mp_vm']      = mp_vm
  all_data['mp_angle']   = mp_angle
  all_data['mp_cvalues'] = mp_cvalues
  all_data['mp_svalues'] = mp_svalues

def check10(log,all_data):

  try:
    f = open('mp_sols/case_ACTIVSg10k.out', "r")
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
    log.joint('c_' + str(f) + ',' + str(f) + ' = ' + str(c) + '\n')
    afile.write('c_' + str(f) + ',' + str(f) + ' = ' + str(c) + '\n')
    afile_kit.write('c_' + str(f) + '_' + str(f) + ' = ' + str(c) + '\n')

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

    log.joint('c_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(c) + '\n')
    afile.write('c_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(c) + '\n')
    afile_kit.write('c_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(c) + '\n')

    log.joint('s_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(s) + '\n')
    afile.write('s_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(s) + '\n')
    afile_kit.write('s_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(s) + '\n')


  afile.close()
  afile_kit.close()


  all_data['mp_Pfvalues'] = mp_Pfvalues
  all_data['mp_Ptvalues'] = mp_Ptvalues
  all_data['mp_Qfvalues'] = mp_Qfvalues
  all_data['mp_Qtvalues'] = mp_Qtvalues

  all_data['mp_vm']      = mp_vm
  all_data['mp_angle']   = mp_angle
  all_data['mp_cvalues'] = mp_cvalues
  all_data['mp_svalues'] = mp_svalues


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
    log.joint('c_' + str(f) + ',' + str(f) + ' = ' + str(c) + '\n')
    afile.write('c_' + str(f) + ',' + str(f) + ' = ' + str(c) + '\n')
    afile_kit.write('c_' + str(f) + '_' + str(f) + ' = ' + str(c) + '\n')

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

    log.joint('c_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(c) + '\n')
    afile.write('c_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(c) + '\n')
    afile_kit.write('c_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(c) + '\n')

    log.joint('s_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(s) + '\n')
    afile.write('s_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(s) + '\n')
    afile_kit.write('s_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(s) + '\n')


  afile.close()
  afile_kit.close()


  all_data['mp_Pfvalues'] = mp_Pfvalues
  all_data['mp_Ptvalues'] = mp_Ptvalues
  all_data['mp_Qfvalues'] = mp_Qfvalues
  all_data['mp_Qtvalues'] = mp_Qtvalues

  all_data['mp_vm']      = mp_vm
  all_data['mp_angle']   = mp_angle
  all_data['mp_cvalues'] = mp_cvalues
  all_data['mp_svalues'] = mp_svalues


def get_flows(log,all_data):

  casefilename = all_data['casefilename']
  casename     = casefilename.split('/')[2].split('.')[0]
  filename     = 'mp_sols/solution_'+ casename +'.txt'

  print(casename)
  print(filename)
  log.joint(" reading file with matpower power flows " + filename + "\n")

  try:
      f = open(filename, "r")
      lines = f.readlines()
      f.close()
  except:
      log.stateandquit("cannot open file", filename)
      sys.exit("failure")

  branches       = all_data['branches']    
  mp_Pfvalues    = {}
  mp_Ptvalues    = {}
  mp_Qfvalues    = {}
  mp_Qtvalues    = {}
  lenlines       = len(lines)

  solfilename     = 'mp_sols/' + casename + '_SOL.txt'
  solfilename_kit = 'mp_sols/' + casename + '_SOL_kit.txt'
  afile           = open(solfilename,"w+")
  afile_kit       = open(solfilename_kit,"w+")

  log.joint(' writing on ' + solfilename + ' ...\n')

  linenum = 1
  while linenum < lenlines:
    thisline            = lines[linenum].split(',')
    branchcount         = int(thisline[0])
    branch              = branches[branchcount]
    log.joint(' branch ' + str(branchcount) + '\n')
    f                   = branch.f
    t                   = branch.t
    if int(thisline[1]) != f or int(thisline[2]) != t:
      breakexit('check')

    mp_Pfvalues[branch]      = float(thisline[3])/all_data['baseMVA']
    mp_Ptvalues[branch]      = float(thisline[4])/all_data['baseMVA']
    mp_Qfvalues[branch]      = float(thisline[5])/all_data['baseMVA']
    mp_Qtvalues[branch]      = float(thisline[6])/all_data['baseMVA']
    
    log.joint('P_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(mp_Pfvalues[branch]) + '\n')
    afile.write('P_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(mp_Pfvalues[branch]) + '\n')
    afile_kit.write('P_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(mp_Pfvalues[branch]) + '\n') #kit
    
    log.joint('P_' + str(branchcount) + ',' + str(t) + ',' + str(f) + ' = ' + str(mp_Ptvalues[branch]) + '\n')
    afile.write('P_' + str(branchcount) + ',' + str(t) + ',' + str(f) + ' = ' + str(mp_Ptvalues[branch]) + '\n')
    afile_kit.write('P_' + str(branchcount) + '_' + str(t) + '_' + str(f) + ' = ' + str(mp_Ptvalues[branch]) + '\n') #kit

    log.joint('Q_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(mp_Qfvalues[branch]) + '\n')
    afile.write('Q_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(mp_Qfvalues[branch]) + '\n')
    afile_kit.write('Q_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(mp_Qfvalues[branch]) + '\n') #kit
    
    log.joint('Q_' + str(branchcount) + ',' + str(t) + ',' + str(f) + ' = ' + str(mp_Qtvalues[branch]) + '\n')
    afile.write('Q_' + str(branchcount) + ',' + str(t) + ',' + str(f) + ' = ' + str(mp_Qtvalues[branch]) + '\n')
    afile_kit.write('Q_' + str(branchcount) + '_' + str(t) + '_' + str(f) + ' = ' + str(mp_Qtvalues[branch]) + '\n') #kit

    linenum   += 1

  afile.close()
  afile_kit.close()
  
  all_data['mp_Pfvalues'] = mp_Pfvalues
  all_data['mp_Ptvalues'] = mp_Ptvalues
  all_data['mp_Qfvalues'] = mp_Qfvalues
  all_data['mp_Qtvalues'] = mp_Qtvalues
  

def get_va(log,all_data):
  casefilename = all_data['casefilename']
  casename = casefilename.split('/')[2].split('.')[0]
  filename = 'mp_sols/solution_va_'+ casename +'.txt'

  log.joint(" reading file matpower solution volt magnitudes and angles " + filename + "\n")

  try:
      f = open(filename, "r")
      lines = f.readlines()
      f.close()
  except:
      log.stateandquit("cannot open file", filename)
      sys.exit("failure")

  lenlines    = len(lines)

  mp_vm       = {}
  mp_angle    = {}
  mp_cvalues  = {}
  mp_svalues  = {}

  buses                  = all_data['buses']
  branches               = all_data['branches']
  IDtoCountmap           = all_data['IDtoCountmap']

  solfilename     = 'mp_sols/' + casename + '_SOL_va.txt'
  solfilename_kit = 'mp_sols/' + casename + '_SOL_va_kit.txt'

  afile     = open(solfilename,"w+")
  afile_kit = open(solfilename_kit,"w+")

  log.joint(' writing on ' + solfilename + ' squared-voltages and angles ...\n')

  linenum  = 1
  while linenum < lenlines:
    thisline = lines[linenum].split(',')
    bus_id   = int(thisline[0])
    buscount = IDtoCountmap[bus_id]
    bus      = buses[buscount]

    mp_vm[bus]        = float(thisline[1])
    mp_angle[bus]     = float(thisline[2]) * math.pi / 180
    mp_cvalues[bus]   = mp_vm[bus]**2

    log.joint('c_' + str(bus_id) + ',' + str(bus_id) + ' = ' + str(mp_cvalues[bus]) + '\n')
    afile.write('c_' + str(bus_id) + ',' + str(bus_id) + ' = ' + str(mp_cvalues[bus]) + '\n')
    afile_kit.write('c_' + str(bus_id) + '_' + str(bus_id) + ' = ' + str(mp_cvalues[bus]) + '\n')

    linenum += 1

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]
    bus_f       = buses[count_of_f]
    bus_t       = buses[count_of_t]


    vm_f        = mp_vm[bus_f]
    vm_t        = mp_vm[bus_t]

    c      = vm_f * vm_t * math.cos(mp_angle[bus_f] - mp_angle[bus_t])
    s      = vm_f * vm_t * math.sin(mp_angle[bus_f] - mp_angle[bus_t])

    mp_cvalues[branch] = c
    mp_svalues[branch] = s

    log.joint('c_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(c) + '\n')
    afile.write('c_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(c) + '\n')
    afile_kit.write('c_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(c) + '\n')
    
    log.joint('s_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(s) + '\n')
    afile.write('s_' + str(branchcount) + ',' + str(f) + ',' + str(t) + ' = ' + str(s) + '\n')
    afile_kit.write('s_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' = ' + str(s) + '\n')

  afile.close()
  afile_kit.close()

  all_data['mp_vm']      = mp_vm
  all_data['mp_angle']   = mp_angle
  all_data['mp_cvalues'] = mp_cvalues
  all_data['mp_svalues'] = mp_svalues

def getsol_knitro(log,all_data):

  casefilename = all_data['casefilename']
  casename     = casefilename.split('/')[2].split('.')[0]
  filename     = 'knitro_sols/ksol_'+ casename +'.txt'

  log.joint(" reading file with knitro sol " + filename + "\n")

  try:
      f = open(filename, "r")
      lines = f.readlines()
      f.close()
  except:
      log.stateandquit("cannot open file", filename)
      sys.exit("failure")

  branches           = all_data['branches']    
  buses              = all_data['buses']
  IDtoCountmap       = all_data['IDtoCountmap']
  knitro_Pfvalues    = {}
  knitro_Ptvalues    = {}
  knitro_Qfvalues    = {}
  knitro_Qtvalues    = {}
  knitro_vm          = {}
  knitro_angle       = {}
  knitro_cvalues     = {}
  knitro_svalues     = {}

  numlines                = len(lines)
  linenum                 = 1
  lookingforendofbranches = 1
  
  while linenum < numlines:
    while lookingforendofbranches and linenum < numlines:
      thisline = lines[linenum].split()

      if thisline[0] == 'nodeID':
        lookingforendofbranches = 0
        linenum += 1
        break

      branchcount         = int(thisline[0])
      branch              = branches[branchcount]
      if branch.f != int(thisline[1]) or branch.t != int(thisline[2]):
        breakexit('buggy')

      knitro_Pfvalues[branch]      = float(thisline[3])
      knitro_Ptvalues[branch]      = float(thisline[4])
      knitro_Qfvalues[branch]      = float(thisline[5])
      knitro_Qtvalues[branch]      = float(thisline[6])
      knitro_angle[branch]         = float(thisline[7])
    
      linenum += 1

    while linenum < numlines:
      thisline = lines[linenum].split()

      buscount  = int(thisline[0])
      bus       = buses[buscount]
                  
      knitro_vm[bus] = float(thisline[1])
  
      linenum += 1

  for bus in buses.values():
    knitro_cvalues[bus] = knitro_vm[bus] * knitro_vm[bus]
  
  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    bus_f      = buses[count_of_f]
    bus_t      = buses[count_of_t]

    knitro_cvalues[branch]    = knitro_vm[bus_f] * knitro_vm[bus_t] * math.cos(knitro_angle[branch])
    knitro_svalues[branch]    = knitro_vm[bus_f] * knitro_vm[bus_t] * math.sin(knitro_angle[branch])

  all_data['knitro_Pfvalues'] = knitro_Pfvalues
  all_data['knitro_Ptvalues'] = knitro_Ptvalues
  all_data['knitro_Qfvalues'] = knitro_Qfvalues
  all_data['knitro_Qtvalues'] = knitro_Qtvalues
  all_data['knitro_vm']       = knitro_vm
  all_data['knitro_angle']    = knitro_angle
  all_data['knitro_cvalues']  = knitro_cvalues
  all_data['knitro_svalues']  = knitro_svalues


def fixCS_knitro(log,all_data):
  tolerance      = all_data['tol_fix']
  themodel       = all_data['themodel']
  buses          = all_data['buses']
  branches       = all_data['branches']
  cvar           = all_data['cvar']
  svar           = all_data['svar']
  knitro_cvalues = all_data['knitro_cvalues']
  knitro_svalues = all_data['knitro_svalues']

  for bus in buses.values():

    v2 = all_data['knitro_cvalues'][bus]

    ubound_v2 = v2 + tolerance
    lbound_v2 = v2 - tolerance

    cvar[bus].setAttr("ub",ubound_v2)
    cvar[bus].setAttr("lb",lbound_v2)

  for branch in branches.values():

    c = all_data['knitro_cvalues'][branch]
    s = all_data['knitro_svalues'][branch]

    #c
    ubound_c = c + tolerance
    lbound_c = c - tolerance

    cvar[branch].setAttr("ub",ubound_c)
    cvar[branch].setAttr("lb",lbound_c)

    #s
    ubound_s = s + tolerance
    lbound_s = s - tolerance

    svar[branch].setAttr("ub",ubound_s)
    svar[branch].setAttr("lb",lbound_s)

  themodel.update()
  themodel.write('fixCS_knitro.lp')
  log.joint('check fixCS_knitro.lp\n')


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
  
